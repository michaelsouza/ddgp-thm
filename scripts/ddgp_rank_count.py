#!/usr/bin/env python3
"""Predict local DDGP solution counts by a labelled-violation rank formula.

The prediction is graph-combinatorial. It does not enumerate branch strings.
Use ``--verify`` to add exact local counts by enumeration for small instances.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from ddgp_count_solutions import Instance, parse_distance_key, parse_edge, load_instance
from ddgp_hypothesis_metrics import (
    Edge,
    dependency_children,
    dependency_sets,
    enumerate_bit_realizations,
    format_generator,
    format_vertices,
    gf2_rank,
    graph_generators,
    labelled_violation_kernel_rank,
    labelled_violation_vectors,
    log2_int,
    mirror_rule_preserves_edges,
)


def load_record(record: dict) -> tuple[Instance, frozenset[Edge]]:
    if "pruning_edges" not in record:
        instance = load_instance(record)
        return instance, frozenset({instance.pruning_edge})

    pruning_edges = frozenset(parse_edge(edge_values) for edge_values in record["pruning_edges"])
    first_pruning_edge = sorted(pruning_edges)[0]
    instance = Instance(
        instance_id=record["id"],
        mode=record["mode"],
        pruning_kind=f"multi-{record['pruning_count']}",
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


def format_edges(edges: set[Edge] | frozenset[Edge]) -> str:
    return ";".join(f"{u + 1}-{v + 1}" for u, v in sorted(edges))


def local_vertices_for_edges(instance: Instance, pruning_edges: frozenset[Edge]) -> set[int]:
    fix_sets = dependency_sets(instance)
    vertices: set[int] = set()
    for u, v in pruning_edges:
        vertices |= fix_sets[u]
        vertices |= fix_sets[v]
    return vertices


def count_row(
    instance: Instance,
    pruning_edges: frozenset[Edge],
    *,
    verify: bool,
    tol: float,
) -> dict[str, object]:
    local_vertices = local_vertices_for_edges(instance, pruning_edges)
    local_base_edges = {
        edge for edge in instance.base_edges
        if edge[0] in local_vertices and edge[1] in local_vertices
    }
    local_valid_edges = local_base_edges | set(pruning_edges)
    decisions = tuple(i for i in range(instance.k + 1, instance.n) if i in local_vertices)

    children = dependency_children(instance)
    generator_records = graph_generators(instance, children, local_vertices, decisions)
    base_generator_records = [
        generator
        for generator in generator_records
        if mirror_rule_preserves_edges(generator, local_base_edges, decisions)
    ]
    base_generator_masks = sorted({generator[0] for generator in base_generator_records})
    individually_valid_records = [
        generator
        for generator in base_generator_records
        if mirror_rule_preserves_edges(generator, pruning_edges, decisions)
    ]
    individually_valid_masks = sorted({generator[0] for generator in individually_valid_records})
    violation_vectors, labels = labelled_violation_vectors(
        base_generator_records,
        pruning_edges,
        decisions,
    )
    combined_vectors = [
        generator[0] + violation
        for generator, violation in zip(base_generator_records, violation_vectors)
    ]

    base_rank = gf2_rank(base_generator_masks)
    individual_rank = gf2_rank(individually_valid_masks)
    labelled_rank = labelled_violation_kernel_rank(
        base_generator_records,
        pruning_edges,
        decisions,
    )
    violation_rank = gf2_rank(violation_vectors)
    combined_rank = gf2_rank(combined_vectors)

    row: dict[str, object] = {
        "id": instance.instance_id,
        "mode": instance.mode,
        "k": instance.k,
        "n": instance.n,
        "pruning_count": len(pruning_edges),
        "pruning_edges": format_edges(pruning_edges),
        "local_vertices": format_vertices(local_vertices),
        "local_decisions": format_vertices(decisions),
        "local_decision_count": len(decisions),
        "base_generator_count": len(base_generator_records),
        "base_generator_rank": base_rank,
        "predicted_local_base_count": 2**base_rank,
        "individually_valid_generator_count": len(individually_valid_records),
        "individually_valid_generator_rank": individual_rank,
        "predicted_by_individual_compatibility": 2**individual_rank,
        "labelled_violation_count": len(labels),
        "labelled_violation_rank": violation_rank,
        "combined_mask_violation_rank": combined_rank,
        "labelled_kernel_rank": labelled_rank,
        "predicted_local_valid_count": 2**labelled_rank,
        "rank_identity_holds": labelled_rank == combined_rank - violation_rank,
        "base_generators": ";".join(
            format_generator(generator, decisions) for generator in base_generator_records
        ),
        "individually_valid_generators": ";".join(
            format_generator(generator, decisions) for generator in individually_valid_records
        ),
    }

    if verify:
        base_bits = enumerate_bit_realizations(instance, local_base_edges, local_vertices, tol)
        valid_bits = enumerate_bit_realizations(instance, local_valid_edges, local_vertices, tol)
        row.update(
            {
                "exact_local_base_count": len(base_bits),
                "exact_local_valid_count": len(valid_bits),
                "exact_local_base_dim": "" if log2_int(len(base_bits)) is None else log2_int(len(base_bits)),
                "exact_local_valid_dim": "" if log2_int(len(valid_bits)) is None else log2_int(len(valid_bits)),
                "base_prediction_matches": 2**base_rank == len(base_bits),
                "valid_prediction_matches": 2**labelled_rank == len(valid_bits),
                "individual_prediction_matches": 2**individual_rank == len(valid_bits),
            }
        )

    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("research/rank-counts.csv"))
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--tol", type=float, default=1e-5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = []
    with args.instances.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            instance, pruning_edges = load_record(json.loads(line))
            rows.append(count_row(instance, pruning_edges, verify=args.verify, tol=args.tol))

    if not rows:
        raise RuntimeError(f"No instances found in {args.instances}.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    if args.verify:
        matches = sum(row["valid_prediction_matches"] for row in rows)
        print(f"wrote {len(rows)} row(s) to {args.out}; valid matches {matches}/{len(rows)}")
    else:
        print(f"wrote {len(rows)} row(s) to {args.out}")


if __name__ == "__main__":
    main()
