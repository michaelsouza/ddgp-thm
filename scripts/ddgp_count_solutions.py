#!/usr/bin/env python3
"""Count solutions for DDGP JSONL instances.

The input format is produced by scripts/ddgp_generate_single_pruning.py.
The output is a CSV with solution counts and dependency-DAG metrics. Counts use
the convention that the first K+1 vertices are fixed, removing the trivial
global reflection from the binary decision space.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


Point = tuple[float, ...]
Edge = tuple[int, int]

EPS = 1e-7


@dataclass(frozen=True)
class Instance:
    instance_id: str
    mode: str
    pruning_kind: str
    k: int
    n: int
    points: tuple[Point, ...]
    predecessor_sets: dict[int, tuple[int, ...]]
    discretization_edges: frozenset[Edge]
    base_edges: frozenset[Edge]
    pruning_edge: Edge
    distances: dict[Edge, float]


def edge(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


def parse_edge(values: list[int]) -> Edge:
    return edge(values[0] - 1, values[1] - 1)


def parse_distance_key(key: str) -> Edge:
    u, v = key.split("-", 1)
    return edge(int(u) - 1, int(v) - 1)


def load_instance(record: dict) -> Instance:
    return Instance(
        instance_id=record["id"],
        mode=record["mode"],
        pruning_kind=record["pruning_kind"],
        k=record["k"],
        n=record["n"],
        points=tuple(tuple(float(value) for value in point) for point in record["points"]),
        predecessor_sets={
            int(i) - 1: tuple(u - 1 for u in predecessors)
            for i, predecessors in record["predecessor_sets"].items()
        },
        discretization_edges=frozenset(parse_edge(e) for e in record["discretization_edges"]),
        base_edges=frozenset(parse_edge(e) for e in record["base_edges"]),
        pruning_edge=parse_edge(record["pruning_edge"]),
        distances={parse_distance_key(k): float(v) for k, v in record["distances"].items()},
    )


def dist(a: Point, b: Point) -> float:
    return math.dist(a, b)


def dot(a: Point, b: Point) -> float:
    return sum(x * y for x, y in zip(a, b))


def add(a: Point, b: Point) -> Point:
    return tuple(x + y for x, y in zip(a, b))


def sub(a: Point, b: Point) -> Point:
    return tuple(x - y for x, y in zip(a, b))


def scale(value: float, vector: Point) -> Point:
    return tuple(value * x for x in vector)


def norm(vector: Point) -> float:
    return math.sqrt(dot(vector, vector))


def cross3(a: Point, b: Point) -> Point:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def circle_intersections(c1: Point, r1: float, c2: Point, r2: float) -> list[Point]:
    dx = c2[0] - c1[0]
    dy = c2[1] - c1[1]
    d = math.hypot(dx, dy)
    if d < EPS:
        return []

    a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d)
    h2 = r1 * r1 - a * a
    if h2 < -EPS:
        return []

    h = math.sqrt(max(0.0, h2))
    mx = c1[0] + a * dx / d
    my = c1[1] + a * dy / d
    px = -dy / d
    py = dx / d

    p1 = (mx + h * px, my + h * py)
    p2 = (mx - h * px, my - h * py)
    if dist(p1, p2) < EPS:
        return [p1]
    return [p1, p2]


def sphere_intersections_3d(
    c0: Point,
    r0: float,
    c1: Point,
    r1: float,
    c2: Point,
    r2: float,
) -> list[Point]:
    e1_raw = sub(c1, c0)
    d = norm(e1_raw)
    if d < EPS:
        return []

    e1 = scale(1.0 / d, e1_raw)
    c2_raw = sub(c2, c0)
    i = dot(e1, c2_raw)
    e2_raw = sub(c2_raw, scale(i, e1))
    j = norm(e2_raw)
    if j < EPS:
        return []

    e2 = scale(1.0 / j, e2_raw)
    e3 = cross3(e1, e2)
    x = (r0 * r0 - r1 * r1 + d * d) / (2.0 * d)
    y = (r0 * r0 - r2 * r2 + dot(c2_raw, c2_raw) - 2.0 * i * x) / (2.0 * j)
    z2 = r0 * r0 - x * x - y * y
    if z2 < -EPS:
        return []

    z = math.sqrt(max(0.0, z2))
    base = add(c0, add(scale(x, e1), scale(y, e2)))
    p1 = add(base, scale(z, e3))
    p2 = add(base, scale(-z, e3))
    if dist(p1, p2) < EPS:
        return [p1]
    return [p1, p2]


def intersection_candidates(centers: tuple[Point, ...], radii: tuple[float, ...]) -> list[Point]:
    if len(centers) == 2:
        return circle_intersections(centers[0], radii[0], centers[1], radii[1])
    if len(centers) == 3:
        return sphere_intersections_3d(
            centers[0],
            radii[0],
            centers[1],
            radii[1],
            centers[2],
            radii[2],
        )
    raise ValueError(f"Unsupported number of predecessors: {len(centers)}.")


def constraints_by_target(edges: Iterable[Edge]) -> dict[int, list[Edge]]:
    by_target: dict[int, list[Edge]] = {}
    for u, v in edges:
        by_target.setdefault(v, []).append((u, v))
    return by_target


def enumerate_realizations(instance: Instance, active_edges: Iterable[Edge], tol: float) -> int:
    return enumerate_realizations_on_vertices(instance, active_edges, set(range(instance.n)), tol)


def enumerate_realizations_on_vertices(
    instance: Instance,
    active_edges: Iterable[Edge],
    vertices: set[int],
    tol: float,
) -> int:
    active_constraints = set(active_edges) - set(instance.discretization_edges)
    active_constraints = {e for e in active_constraints if e[0] in vertices and e[1] in vertices}
    by_target = constraints_by_target(active_constraints)
    search_vertices = [i for i in range(instance.k + 1, instance.n) if i in vertices]

    coords: list[Point | None] = [None] * instance.n
    for i in range(instance.k + 1):
        coords[i] = instance.points[i]

    def search(position: int) -> int:
        if position == len(search_vertices):
            return 1

        i = search_vertices[position]
        predecessors = instance.predecessor_sets[i]
        if any(predecessor not in vertices for predecessor in predecessors):
            raise ValueError(f"Vertex {i + 1} has a predecessor outside the induced vertex set.")
        centers = tuple(coords[predecessor] for predecessor in predecessors)
        assert all(center is not None for center in centers)
        candidates = intersection_candidates(
            centers,  # type: ignore[arg-type]
            tuple(instance.distances[edge(predecessor, i)] for predecessor in predecessors),
        )

        total = 0
        for candidate in candidates:
            coords[i] = candidate
            feasible = True
            for u, v in by_target.get(i, []):
                assert v == i
                assert coords[u] is not None
                if abs(dist(coords[u], candidate) - instance.distances[(u, v)]) > tol:
                    feasible = False
                    break
            if feasible:
                total += search(position + 1)
            coords[i] = None
        return total

    return search(0)


def dependency_sets(instance: Instance) -> list[set[int]]:
    fix_sets = [set([i]) for i in range(instance.n)]
    for i in range(instance.k + 1, instance.n):
        fixed = {i}
        for u in instance.predecessor_sets[i]:
            fixed |= fix_sets[u]
        fix_sets[i] = fixed
    return fix_sets


def format_vertices(vertices: Iterable[int]) -> str:
    return "{" + ",".join(str(v + 1) for v in sorted(vertices)) + "}"


def count_record(instance: Instance, tol: float) -> dict:
    raw_count = 2 ** (instance.n - instance.k - 1)
    base_count = enumerate_realizations(instance, instance.base_edges, tol)
    valid_count = enumerate_realizations(instance, set(instance.base_edges) | {instance.pruning_edge}, tol)

    fix_sets = dependency_sets(instance)
    u, v = instance.pruning_edge
    fu = fix_sets[u]
    fv = fix_sets[v]
    union = fu | fv
    intersection = fu & fv
    symdiff = fu ^ fv
    noninitial = set(range(instance.k + 1, instance.n))

    union_choices = len(union & noninitial)
    intersection_choices = len(intersection & noninitial)
    symdiff_choices = len(symdiff & noninitial)
    outside_choices = (instance.n - instance.k - 1) - union_choices
    outside_factor = 2 ** outside_choices
    local_base_edges = {e for e in instance.base_edges if e[0] in union and e[1] in union}
    local_base_count = enumerate_realizations_on_vertices(instance, local_base_edges, union, tol)
    local_valid_count = enumerate_realizations_on_vertices(
        instance,
        local_base_edges | {instance.pruning_edge},
        union,
        tol,
    )

    return {
        "id": instance.instance_id,
        "mode": instance.mode,
        "pruning_kind": instance.pruning_kind,
        "k": instance.k,
        "n": instance.n,
        "pruning_edge": f"{u + 1}-{v + 1}",
        "raw_count": raw_count,
        "base_count": base_count,
        "valid_count": valid_count,
        "removed_by_pruning": base_count - valid_count,
        "valid_fraction_of_base": valid_count / base_count if base_count else 0.0,
        "outside_choices": outside_choices,
        "outside_factor": outside_factor,
        "normalized_base_count": base_count / outside_factor,
        "normalized_valid_count": valid_count / outside_factor,
        "local_raw_count": 2**union_choices,
        "local_base_count": local_base_count,
        "local_valid_count": local_valid_count,
        "local_removed_by_pruning": local_base_count - local_valid_count,
        "local_valid_fraction_of_base": local_valid_count / local_base_count if local_base_count else 0.0,
        "fix_u_size": len(fu),
        "fix_v_size": len(fv),
        "fix_union_size": len(union),
        "fix_intersection_size": len(intersection),
        "fix_symdiff_size": len(symdiff),
        "union_choices": union_choices,
        "intersection_choices": intersection_choices,
        "symdiff_choices": symdiff_choices,
        "fix_u": format_vertices(fu),
        "fix_v": format_vertices(fv),
        "fix_union": format_vertices(union),
        "fix_intersection": format_vertices(intersection),
        "fix_symdiff": format_vertices(symdiff),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, default=Path("experiments/single-pruning/instances.jsonl"))
    parser.add_argument("--out", type=Path, default=Path("experiments/single-pruning/counts.csv"))
    parser.add_argument("--tol", type=float, default=1e-5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with args.instances.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            rows.append(count_record(load_instance(json.loads(line)), args.tol))

    if not rows:
        raise RuntimeError(f"No instances found in {args.instances}.")

    with args.out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} row(s) to {args.out}")


if __name__ == "__main__":
    main()
