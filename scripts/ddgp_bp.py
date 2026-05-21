#!/usr/bin/env python3
"""Branch-and-Prune (BP) Algorithm for Discretizable Distance Geometry Problems (DDGP).

This script implements the classic Branch-and-Prune (BP) depth-first search
algorithm to reconstruct full geometric realizations of DDGP instances.
It uses K-lateration to branch at each decision vertex, and prunes infeasible
branches immediately using base and pruning constraints.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

# Reuse the core geometric utilities and Instance definition from ddgp_count_solutions
from ddgp_count_solutions import (
    Edge,
    Instance,
    Point,
    dist,
    edge,
    intersection_candidates,
    load_instance,
    parse_edge,
    parse_distance_key,
)


def load_record_with_pruning(record: dict[str, Any]) -> tuple[Instance, frozenset[Edge]]:
    """Loads a DDGP instance, robust to both single and multiple pruning edges."""
    if "pruning_edges" not in record:
        instance = load_instance(record)
        return instance, frozenset({instance.pruning_edge})

    pruning_edges = frozenset(parse_edge(edge_values) for edge_values in record["pruning_edges"])
    first_pruning_edge = sorted(pruning_edges)[0] if pruning_edges else (0, 1)
    
    instance = Instance(
        instance_id=record["id"],
        mode=record["mode"],
        pruning_kind=f"multi-{record.get('pruning_count', len(pruning_edges))}",
        k=record["k"],
        n=record["n"],
        points=tuple(tuple(float(value) for value in point) for point in record["points"]),
        predecessor_sets={
            int(i) - 1: tuple(u - 1 for u in predecessors)
            for i, predecessors in record["predecessor_sets"].items()
        },
        discretization_edges=frozenset(parse_edge(e) for e in record["discretization_edges"]),
        base_edges=frozenset(parse_edge(e) for e in record["base_edges"]),
        pruning_edge=first_pruning_edge,
        distances={parse_distance_key(k): float(v) for k, v in record["distances"].items()},
    )
    return instance, pruning_edges


def constraints_by_target(edges: Iterable[Edge]) -> dict[int, list[Edge]]:
    """Groups edges by their maximum vertex index (the target for check evaluation)."""
    by_target: dict[int, list[Edge]] = {}
    for u, v in edges:
        # Guarantee that we index by the larger vertex number (the one positioned last)
        target = max(u, v)
        source = min(u, v)
        by_target.setdefault(target, []).append((source, target))
    return by_target


def branch_and_prune(
    instance: Instance,
    pruning_edges: frozenset[Edge],
    tol: float = 1e-5,
    verbose: bool = False,
) -> list[tuple[Point, ...]]:
    """Executes the classic Branch-and-Prune DFS algorithm on a DDGP instance.
    
    Returns a list of all valid full realizations (tuples of coordinates).
    """
    n = instance.n
    k = instance.k

    # Identify all active constraints that need to be evaluated during search.
    # Discretization edges are always satisfied by construction (sphere intersections).
    # Thus, constraints to check are: (base_edges union pruning_edges) minus discretization_edges.
    active_constraints = (instance.base_edges | pruning_edges) - instance.discretization_edges
    by_target = constraints_by_target(active_constraints)

    # Output details about search space
    if verbose:
        print(f"\n[BP Search Init] Instance: {instance.instance_id}")
        print(f"  - Dimension K: {k}, Total Vertices N: {n}")
        print(f"  - Decision vertices: {list(range(k + 1, n))}")
        print(f"  - Active check constraints: {len(active_constraints)} edges")

    # Keep track of coordinates
    coords: list[Point | None] = [None] * n
    for i in range(k + 1):
        coords[i] = instance.points[i]

    realizations: list[tuple[Point, ...]] = []
    nodes_visited = 0
    nodes_pruned = 0

    def dfs(i: int) -> None:
        nonlocal nodes_visited, nodes_pruned
        nodes_visited += 1

        # Base case: successfully positioned all vertices
        if i == n:
            realization = tuple(coords[j] for j in range(n))  # type: ignore
            realizations.append(realization)
            if verbose:
                print(f"  --> [SOLUTION FOUND] Realization #{len(realizations)} completed!")
            return

        predecessors = instance.predecessor_sets[i]
        centers = tuple(coords[p] for p in predecessors)
        
        # Verify that all predecessors have already been placed
        assert all(c is not None for c in centers), f"Predecessors for vertex {i+1} not positioned!"

        # Branching step: K-lateration generates circle or sphere intersection candidates
        distances = tuple(instance.distances[edge(p, i)] for p in predecessors)
        candidates = intersection_candidates(centers, distances)  # type: ignore[arg-type]

        if not candidates:
            # Pruned due to empty intersection (incompatible discretization distances)
            nodes_pruned += 1
            if verbose:
                print(f"  Level {i+1}: Pruned (no intersection candidates)")
            return

        for idx, candidate in enumerate(candidates):
            coords[i] = candidate
            
            # Pruning step: evaluate constraints ending at vertex i
            feasible = True
            for u, v in by_target.get(i, []):
                assert v == i, "Target vertex mismatch in constraints"
                assert coords[u] is not None, f"Source vertex {u+1} of constraint not positioned"
                
                dist_calc = dist(coords[u], candidate)
                dist_real = instance.distances[(u, v)]
                if abs(dist_calc - dist_real) > tol:
                    feasible = False
                    nodes_pruned += 1
                    if verbose:
                        print(f"  Level {i+1} (branch {idx}): Pruned by edge constraint {u+1}-{v+1} (err: {abs(dist_calc - dist_real):.2e} > {tol:.2e})")
                    break

            if feasible:
                if verbose:
                    print(f"  Level {i+1} (branch {idx}): Feasible. Descending...")
                dfs(i + 1)
            
            coords[i] = None

    # Start recursive DFS search at first decision vertex (k + 1)
    dfs(k + 1)

    if verbose:
        print(f"[BP Search End] Found {len(realizations)} realization(s).")
        print(f"  - Total nodes visited: {nodes_visited}")
        print(f"  - Branches pruned: {nodes_pruned}")

    return realizations


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, required=True, help="Path to JSONL instances file")
    parser.add_argument("--instance-id", type=str, default=None, help="ID of a specific instance to run (runs all if omitted)")
    parser.add_argument("--out", type=Path, default=None, help="Optional output path to write coordinates of found solutions (JSON format)")
    parser.add_argument("--tol", type=float, default=1e-5, help="Pruning tolerance (default: 1e-5)")
    parser.add_argument("--verbose", action="store_true", help="Print detailed DFS search steps and pruning messages")
    args = parser.parse_args()

    if not args.instances.exists():
        print(f"Error: Instances file not found: {args.instances}", file=sys.stderr)
        sys.exit(1)

    records_processed = 0
    results_output = {}

    with args.instances.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            
            record = json.loads(line)
            inst_id = record.get("id")
            
            # If a specific instance-id is requested, skip others
            if args.instance_id and inst_id != args.instance_id:
                continue

            instance, pruning_edges = load_record_with_pruning(record)
            
            realizations = branch_and_prune(
                instance=instance,
                pruning_edges=pruning_edges,
                tol=args.tol,
                verbose=args.verbose or (args.instance_id is not None),
            )
            
            print(f"Instance {instance.instance_id}: Found {len(realizations)} valid realization(s).")
            
            # Format and collect coordinates for output if requested
            results_output[instance.instance_id] = {
                "instance_id": instance.instance_id,
                "solution_count": len(realizations),
                "realizations": [
                    [[round(coord, 6) for coord in point] for point in realization]
                    for realization in realizations
                ]
            }
            records_processed += 1

    if records_processed == 0:
        if args.instance_id:
            print(f"Error: Instance ID '{args.instance_id}' not found in {args.instances}", file=sys.stderr)
        else:
            print("Error: No instances were processed.", file=sys.stderr)
        sys.exit(1)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", encoding="utf-8") as f:
            json.dump(results_output, f, indent=2)
        print(f"Wrote coordinates of solutions for {records_processed} instance(s) to {args.out}")


if __name__ == "__main__":
    main()
