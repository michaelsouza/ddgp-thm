#!/usr/bin/env python3
"""Generate DDGP instances with multiple explicit pruning edges."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from ddgp_generate_single_pruning import (
    Edge,
    base_edges,
    choose_predecessor_sets,
    discretization_edges,
    dist,
    edge,
    encode_distance_key,
    encode_edge,
    generate_points,
    is_globally_nondegenerate,
)


def candidate_pruning_edges(n: int, k: int, base: set[Edge]) -> list[Edge]:
    return [
        edge(u, v)
        for v in range(k + 1, n)
        for u in range(v)
        if edge(u, v) not in base
    ]


def build_record(
    *,
    instance_id: str,
    n: int,
    k: int,
    mode: str,
    pruning_count: int,
    rng: random.Random,
    min_area2: float,
    global_min_area2: float,
) -> dict:
    for _ in range(5000):
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
        candidates = candidate_pruning_edges(n, k, b_edges)
        if len(candidates) < pruning_count:
            continue

        pruning_edges = sorted(rng.sample(candidates, pruning_count))
        all_edges = b_edges | set(pruning_edges)
        distances = {
            encode_distance_key(e): dist(points[e[0]], points[e[1]])
            for e in sorted(all_edges)
        }

        return {
            "id": instance_id,
            "mode": mode,
            "k": k,
            "n": n,
            "pruning_count": pruning_count,
            "points": [list(point) for point in points],
            "predecessor_sets": {
                str(i + 1): [u + 1 for u in predecessors]
                for i, predecessors in predecessor_sets.items()
            },
            "discretization_edges": [encode_edge(e) for e in sorted(d_edges)],
            "base_edges": [encode_edge(e) for e in sorted(b_edges)],
            "pruning_edges": [encode_edge(e) for e in pruning_edges],
            "distances": distances,
        }

    raise RuntimeError(f"Could not generate instance {instance_id}.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("research/multi-pruning-instances.jsonl"))
    parser.add_argument("--samples-per-count", type=int, default=50)
    parser.add_argument("--pruning-counts", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--n", type=int, default=9)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--mode", choices=["ddgp", "dmdgp"], default="ddgp")
    parser.add_argument("--seed", type=int, default=100)
    parser.add_argument("--min-area2", type=float, default=1e-3)
    parser.add_argument(
        "--global-min-area2",
        type=float,
        default=1e-4,
        help="Reject instances with any nearly degenerate K-simplex.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with args.out.open("w", encoding="utf-8") as fh:
        for pruning_count in args.pruning_counts:
            for sample in range(1, args.samples_per_count + 1):
                count += 1
                instance_id = f"{args.mode}-p{pruning_count}-n{args.n}-{sample:04d}"
                record = build_record(
                    instance_id=instance_id,
                    n=args.n,
                    k=args.k,
                    mode=args.mode,
                    pruning_count=pruning_count,
                    rng=rng,
                    min_area2=args.min_area2,
                    global_min_area2=args.global_min_area2,
                )
                fh.write(json.dumps(record, sort_keys=True) + "\n")

    print(f"wrote {count} instance(s) to {args.out}")


if __name__ == "__main__":
    main()
