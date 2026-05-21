#!/usr/bin/env python3
r"""Verify Rank-Count predictions against Branch-and-Prune (BP).

This script compares:
1. The local Rank-Count prediction against restricted local BP search.
2. The global Rank-Count prediction against full BP search.

The global prediction applies the same labelled-violation rank formula with
L = V, B = V \ V0, and F = base_edges union pruning_edges, so it validates the
solution count estimated by the theory directly against scripts/ddgp_bp.py.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ddgp_count_solutions import Instance
from ddgp_bp import load_record_with_pruning, branch_and_prune
from ddgp_hypothesis_metrics import (
    dependency_sets,
    dependency_children,
    graph_generators,
    labelled_violation_kernel_rank,
    enumerate_bit_realizations,
)


def rank_count_for_vertices(
    instance: Instance,
    active_edges: frozenset[tuple[int, int]] | set[tuple[int, int]],
    vertices: set[int],
) -> int:
    decisions = tuple(i for i in range(instance.k + 1, instance.n) if i in vertices)
    children = dependency_children(instance)
    generator_records = graph_generators(instance, children, vertices, decisions)
    rank = labelled_violation_kernel_rank(
        generator_records,
        active_edges,
        decisions,
    )
    return 2**rank


def run_verification(instances_path: Path, max_instances: int = 10, tol: float = 1e-5) -> None:
    print("=" * 90)
    print("VERIFYING COMPATIBILITY: Rank-Count Prediction vs Branch-and-Prune (BP)")
    print(f"Source instances: {instances_path}")
    print("=" * 90)
    print(
        f"{'Instance ID':<25} | "
        f"{'Local theory':<12} | {'Local BP':<8} | {'Local':<6} | "
        f"{'Global theory':<13} | {'Global BP':<9} | {'Global':<6}"
    )
    print("-" * 100)

    count = 0
    local_matches = 0
    global_matches = 0
    with instances_path.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            
            record = json.loads(line)
            instance, pruning_edges = load_record_with_pruning(record)

            # --- 1. Compute Rank-Count Prediction ---
            fix_sets = dependency_sets(instance)
            
            # Local vertices L_P are the union of the fixing sets of the pruning endpoints
            local_vertices = set()
            for u, v in pruning_edges:
                local_vertices |= fix_sets[u]
                local_vertices |= fix_sets[v]
                
            local_base_edges = {
                edge for edge in instance.base_edges
                if edge[0] in local_vertices and edge[1] in local_vertices
            }
            local_valid_edges = local_base_edges | pruning_edges
            predicted_local_count = rank_count_for_vertices(
                instance,
                local_valid_edges,
                local_vertices,
            )

            # --- 1b. Compute global Rank-Count prediction ---
            global_vertices = set(range(instance.n))
            global_valid_edges = instance.base_edges | pruning_edges
            predicted_global_count = rank_count_for_vertices(
                instance,
                global_valid_edges,
                global_vertices,
            )

            # --- 2. Compute Local BP Count (Restricted Search) ---
            # Using enumerate_bit_realizations which runs a BP restricted to the local decision vertices
            local_bits = enumerate_bit_realizations(instance, local_valid_edges, local_vertices, tol)
            local_bp_count = len(local_bits)

            # --- 3. Compute Global BP Count (Full Search) ---
            global_realizations = branch_and_prune(instance, pruning_edges, tol=tol, verbose=False)
            global_bp_count = len(global_realizations)

            # --- 4. Compare and Report ---
            match_local = (predicted_local_count == local_bp_count)
            match_global = (predicted_global_count == global_bp_count)
            local_status = "PASS" if match_local else "FAIL"
            global_status = "PASS" if match_global else "FAIL"
            local_matches += int(match_local)
            global_matches += int(match_global)
            
            print(
                f"{instance.instance_id:<25} | "
                f"{predicted_local_count:<12} | {local_bp_count:<8} | {local_status:<6} | "
                f"{predicted_global_count:<13} | {global_bp_count:<9} | {global_status:<6}"
            )
            
            count += 1
            if count >= max_instances:
                break

    print("=" * 90)
    print(f"Local matches:  {local_matches}/{count}")
    print(f"Global matches: {global_matches}/{count}")
    print()
    print("Local theory uses L_P from pruning-edge dependency sets.")
    print("Global theory uses all vertices and all active edges: base_edges union pruning_edges.")
    print("=" * 90)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, required=True, help="Path to JSONL instances file")
    parser.add_argument("--count", type=int, default=10, help="Maximum number of instances to verify (default: 10)")
    parser.add_argument("--tol", type=float, default=1e-5, help="Tolerância de poda (default: 1e-5)")
    args = parser.parse_args()

    if not args.instances.exists():
        print(f"Error: Instances file not found: {args.instances}", file=sys.stderr)
        sys.exit(1)

    run_verification(args.instances, max_instances=args.count, tol=args.tol)


if __name__ == "__main__":
    main()
