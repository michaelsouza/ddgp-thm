#!/usr/bin/env python3
"""Generate DDGP instances with exactly one explicit pruning edge.

The generated dataset is JSONL. Vertices are written as 1-based indices to make
the files easier to inspect by hand. The first K+1 vertices are treated as a
fixed rigid seed and do not carry binary decisions.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from pathlib import Path


Point = tuple[float, ...]
Edge = tuple[int, int]

PRUNING_KINDS = ("last-seed", "last-nonseed", "nonlast-nonseed")


def edge(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


def dist(a: Point, b: Point) -> float:
    return math.dist(a, b)


def triangle_area2(a: Point, b: Point, c: Point) -> float:
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def tetra_volume6(a: Point, b: Point, c: Point, d: Point) -> float:
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    ad = (d[0] - a[0], d[1] - a[1], d[2] - a[2])
    cross = (
        ac[1] * ad[2] - ac[2] * ad[1],
        ac[2] * ad[0] - ac[0] * ad[2],
        ac[0] * ad[1] - ac[1] * ad[0],
    )
    return abs(ab[0] * cross[0] + ab[1] * cross[1] + ab[2] * cross[2])


def simplex_measure(points: tuple[Point, ...]) -> float:
    if len(points) == 3:
        return triangle_area2(points[0], points[1], points[2])
    if len(points) == 4:
        return tetra_volume6(points[0], points[1], points[2], points[3])
    raise ValueError(f"Unsupported simplex size: {len(points)}.")


def generate_points(n: int, k: int, rng: random.Random) -> tuple[Point, ...]:
    return tuple(
        tuple(rng.uniform(-1.0, 1.0) for _ in range(k))
        for _ in range(n)
    )


def is_globally_nondegenerate(points: tuple[Point, ...], k: int, min_measure: float) -> bool:
    if min_measure <= 0.0:
        return True
    return all(
        simplex_measure(tuple(points[index] for index in simplex)) > min_measure
        for simplex in itertools.combinations(range(len(points)), k + 1)
    )


def choose_predecessor_sets(
    points: tuple[Point, ...],
    *,
    k: int,
    mode: str,
    rng: random.Random,
    min_area2: float,
) -> dict[int, tuple[int, ...]] | None:
    if k not in (2, 3):
        raise ValueError("This generator currently supports only K=2 and K=3.")

    predecessor_sets: dict[int, tuple[int, ...]] = {}
    for i in range(k, len(points)):
        if mode == "dmdgp":
            candidates = [tuple(range(i - k, i))]
        else:
            candidates = list(itertools.combinations(range(i), k))
            rng.shuffle(candidates)

        chosen = None
        for pred in candidates:
            simplex = tuple(points[index] for index in (*pred, i))
            if simplex_measure(simplex) > min_area2:
                chosen = tuple(sorted(pred))
                break
        if chosen is None:
            return None
        predecessor_sets[i] = chosen

    return predecessor_sets


def discretization_edges(predecessor_sets: dict[int, tuple[int, ...]], k: int) -> set[Edge]:
    edges = {edge(u, i) for i, predecessors in predecessor_sets.items() for u in predecessors}
    edges |= {edge(u, v) for u, v in itertools.combinations(range(k + 1), 2)}
    return edges


def base_edges(predecessor_sets: dict[int, tuple[int, ...]], k: int) -> set[Edge]:
    edges = discretization_edges(predecessor_sets, k)
    for predecessors in predecessor_sets.values():
        for u, v in itertools.combinations(predecessors, 2):
            edges.add(edge(u, v))
    return edges


def candidate_pruning_edges(n: int, k: int, base: set[Edge], kind: str) -> list[Edge]:
    if kind == "last-seed":
        candidates = [edge(u, n - 1) for u in range(k + 1)]
    elif kind == "last-nonseed":
        candidates = [edge(u, n - 1) for u in range(k, n - 1)]
    elif kind == "nonlast-nonseed":
        candidates = [edge(u, v) for v in range(k + 1, n - 1) for u in range(k, v)]
    else:
        raise ValueError(f"Unknown pruning kind: {kind}")

    return [e for e in candidates if e not in base]


def encode_edge(e: Edge) -> list[int]:
    return [e[0] + 1, e[1] + 1]


def encode_distance_key(e: Edge) -> str:
    return f"{e[0] + 1}-{e[1] + 1}"


def build_record(
    *,
    instance_id: str,
    n: int,
    k: int,
    mode: str,
    pruning_kind: str,
    rng: random.Random,
    min_area2: float,
    global_min_area2: float,
) -> dict:
    for _ in range(2000):
        points = generate_points(n, k, rng)
        if not is_globally_nondegenerate(points, k, global_min_area2):
            continue
        predecessor_sets = choose_predecessor_sets(
            points,
            k=k,
            mode=mode,
            rng=rng,
            min_area2=min_area2,
        )
        if predecessor_sets is None:
            continue

        d_edges = discretization_edges(predecessor_sets, k)
        b_edges = base_edges(predecessor_sets, k)
        candidates = candidate_pruning_edges(n, k, b_edges, pruning_kind)
        if not candidates:
            continue

        pruning_edge = rng.choice(candidates)
        all_edges = b_edges | {pruning_edge}
        distances = {
            encode_distance_key(e): dist(points[e[0]], points[e[1]])
            for e in sorted(all_edges)
        }

        return {
            "id": instance_id,
            "mode": mode,
            "k": k,
            "n": n,
            "pruning_kind": pruning_kind,
            "points": [list(point) for point in points],
            "predecessor_sets": {
                str(i + 1): [u + 1 for u in predecessors]
                for i, predecessors in predecessor_sets.items()
            },
            "discretization_edges": [encode_edge(e) for e in sorted(d_edges)],
            "base_edges": [encode_edge(e) for e in sorted(b_edges)],
            "pruning_edge": encode_edge(pruning_edge),
            "distances": distances,
        }

    raise RuntimeError(f"Could not generate instance {instance_id}.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("experiments/single-pruning/instances.jsonl"))
    parser.add_argument("--samples-per-kind", type=int, default=20)
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--mode", choices=["ddgp", "dmdgp"], default="ddgp")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-area2", type=float, default=1e-3)
    parser.add_argument(
        "--global-min-area2",
        type=float,
        default=0.0,
        help="Reject instances with any nearly degenerate K-simplex; disabled by default.",
    )
    parser.add_argument("--kinds", nargs="+", choices=PRUNING_KINDS, default=list(PRUNING_KINDS))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with args.out.open("w", encoding="utf-8") as fh:
        for kind in args.kinds:
            for sample in range(1, args.samples_per_kind + 1):
                count += 1
                instance_id = f"{args.mode}-{kind}-{sample:04d}"
                record = build_record(
                    instance_id=instance_id,
                    n=args.n,
                    k=args.k,
                    mode=args.mode,
                    pruning_kind=kind,
                    rng=rng,
                    min_area2=args.min_area2,
                    global_min_area2=args.global_min_area2,
                )
                fh.write(json.dumps(record, sort_keys=True) + "\n")

    print(f"wrote {count} instance(s) to {args.out}")


if __name__ == "__main__":
    main()
