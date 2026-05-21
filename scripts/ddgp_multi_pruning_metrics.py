#!/usr/bin/env python3
"""Test local solution-count hypotheses on DDGP instances with multiple pruning edges."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from ddgp_count_solutions import Instance, parse_distance_key, parse_edge
from ddgp_hypothesis_metrics import (
    Edge,
    affine_hull_rank,
    base_clique_metrics,
    dependency_children,
    dependency_sets,
    enumerate_bit_realizations,
    format_generator,
    format_mask,
    format_vertices,
    gf2_rank,
    graph_generators,
    is_affine_subspace,
    labelled_violation_kernel_masks,
    labelled_violation_kernel_rank,
    log2_int,
    mask_preserves,
    masks_in_span,
    mirror_rule_preserves_edges,
    span_elements,
    translation_stabilizer,
    violation_kernel_masks,
    violation_kernel_rank,
)


def load_multi_instance(record: dict) -> tuple[Instance, frozenset[Edge]]:
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


def local_vertices_for_pruning_edges(instance: Instance, pruning_edges: frozenset[Edge]) -> set[int]:
    fix_sets = dependency_sets(instance)
    local_vertices: set[int] = set()
    for u, v in pruning_edges:
        local_vertices |= fix_sets[u]
        local_vertices |= fix_sets[v]
    return local_vertices


def multi_pruning_row(instance: Instance, pruning_edges: frozenset[Edge], tol: float) -> dict[str, object]:
    local_vertices = local_vertices_for_pruning_edges(instance, pruning_edges)
    local_base_edges = {
        e for e in instance.base_edges if e[0] in local_vertices and e[1] in local_vertices
    }
    local_valid_edges = local_base_edges | set(pruning_edges)
    decisions = tuple(i for i in range(instance.k + 1, instance.n) if i in local_vertices)

    base_bits = enumerate_bit_realizations(instance, local_base_edges, local_vertices, tol)
    valid_bits = enumerate_bit_realizations(instance, local_valid_edges, local_vertices, tol)

    full_base_stabilizer = translation_stabilizer(base_bits, len(decisions))
    full_valid_stabilizer = translation_stabilizer(valid_bits, len(decisions))
    full_base_stabilizer_rank = gf2_rank(full_base_stabilizer)
    full_valid_stabilizer_rank = gf2_rank(full_valid_stabilizer)

    children = dependency_children(instance)
    generator_records = graph_generators(instance, children, local_vertices, decisions)
    edge_rule_base_records = [
        generator
        for generator in generator_records
        if mirror_rule_preserves_edges(generator, local_base_edges, decisions)
    ]
    edge_rule_base_masks = sorted({generator[0] for generator in edge_rule_base_records})
    edge_rule_valid_records = [
        generator
        for generator in edge_rule_base_records
        if mirror_rule_preserves_edges(generator, pruning_edges, decisions)
    ]
    edge_rule_valid_masks = sorted({generator[0] for generator in edge_rule_valid_records})
    edge_rule_valid_span = span_elements(edge_rule_valid_masks, width=len(decisions))
    edge_rule_valid_rank = gf2_rank(edge_rule_valid_masks)
    predicted_valid = 2**edge_rule_valid_rank
    violation_kernel_valid_masks = violation_kernel_masks(
        edge_rule_base_records,
        pruning_edges,
        decisions,
    )
    violation_kernel_valid_rank = gf2_rank(violation_kernel_valid_masks)
    violation_kernel_valid_formula_rank = violation_kernel_rank(
        edge_rule_base_records,
        pruning_edges,
        decisions,
    )
    predicted_valid_by_violation_kernel = 2**violation_kernel_valid_rank
    labelled_violation_kernel_valid_masks = labelled_violation_kernel_masks(
        edge_rule_base_records,
        pruning_edges,
        decisions,
    )
    labelled_violation_kernel_valid_rank = gf2_rank(labelled_violation_kernel_valid_masks)
    labelled_violation_kernel_valid_formula_rank = labelled_violation_kernel_rank(
        edge_rule_base_records,
        pruning_edges,
        decisions,
    )
    predicted_valid_by_labelled_violation_kernel = 2**labelled_violation_kernel_valid_rank

    graph_base_span = span_elements(edge_rule_base_masks, width=len(decisions))
    graph_valid_span_stabilizer = [
        mask for mask in graph_base_span if mask_preserves(valid_bits, mask)
    ]
    graph_valid_span_rank = gf2_rank(graph_valid_span_stabilizer)

    local_base_count = len(base_bits)
    local_valid_count = len(valid_bits)
    local_base_dim = log2_int(local_base_count)
    local_valid_dim = log2_int(local_valid_count)

    row: dict[str, object] = {
        "id": instance.instance_id,
        "mode": instance.mode,
        "k": instance.k,
        "n": instance.n,
        "pruning_count": len(pruning_edges),
        "pruning_edges": format_edges(pruning_edges),
        "local_decisions": len(decisions),
        "local_base_count": local_base_count,
        "local_valid_count": local_valid_count,
        "local_base_dim": "" if local_base_dim is None else local_base_dim,
        "local_valid_dim": "" if local_valid_dim is None else local_valid_dim,
        "base_is_affine_subspace": is_affine_subspace(base_bits),
        "valid_is_affine_subspace": is_affine_subspace(valid_bits),
        "base_affine_hull_rank": affine_hull_rank(base_bits),
        "valid_affine_hull_rank": affine_hull_rank(valid_bits),
        "full_base_stabilizer_rank": full_base_stabilizer_rank,
        "full_valid_stabilizer_rank": full_valid_stabilizer_rank,
        "full_valid_rank_matches": 2**full_valid_stabilizer_rank == local_valid_count,
        "edge_rule_base_generator_rank": gf2_rank(edge_rule_base_masks),
        "edge_rule_valid_generator_rank": edge_rule_valid_rank,
        "predicted_valid_by_edge_rule": predicted_valid,
        "edge_rule_rank_matches": predicted_valid == local_valid_count,
        "edge_rule_rank_matches_full_valid": edge_rule_valid_rank == full_valid_stabilizer_rank,
        "full_valid_stabilizer_in_edge_rule_span": masks_in_span(
            full_valid_stabilizer,
            edge_rule_valid_masks,
        ),
        "violation_kernel_valid_rank": violation_kernel_valid_rank,
        "violation_kernel_valid_formula_rank": violation_kernel_valid_formula_rank,
        "violation_kernel_formula_rank_matches_masks": violation_kernel_valid_formula_rank
        == violation_kernel_valid_rank,
        "predicted_valid_by_violation_kernel": predicted_valid_by_violation_kernel,
        "violation_kernel_rank_matches": predicted_valid_by_violation_kernel == local_valid_count,
        "violation_kernel_rank_matches_full_valid": violation_kernel_valid_rank
        == full_valid_stabilizer_rank,
        "full_valid_stabilizer_in_violation_kernel_span": masks_in_span(
            full_valid_stabilizer,
            violation_kernel_valid_masks,
        ),
        "labelled_violation_kernel_valid_rank": labelled_violation_kernel_valid_rank,
        "labelled_violation_kernel_valid_formula_rank": labelled_violation_kernel_valid_formula_rank,
        "labelled_violation_kernel_formula_rank_matches_masks": (
            labelled_violation_kernel_valid_formula_rank == labelled_violation_kernel_valid_rank
        ),
        "predicted_valid_by_labelled_violation_kernel": predicted_valid_by_labelled_violation_kernel,
        "labelled_violation_kernel_rank_matches": predicted_valid_by_labelled_violation_kernel
        == local_valid_count,
        "labelled_violation_kernel_rank_matches_full_valid": labelled_violation_kernel_valid_rank
        == full_valid_stabilizer_rank,
        "full_valid_stabilizer_in_labelled_violation_kernel_span": masks_in_span(
            full_valid_stabilizer,
            labelled_violation_kernel_valid_masks,
        ),
        "graph_valid_span_rank": graph_valid_span_rank,
        "graph_valid_span_rank_matches": 2**graph_valid_span_rank == local_valid_count,
        "graph_valid_rank_matches_full_valid": graph_valid_span_rank == full_valid_stabilizer_rank,
        "decisions": format_vertices(decisions),
        "fix_union": format_vertices(local_vertices),
        "edge_rule_base_generators": ";".join(
            format_generator(generator, decisions) for generator in edge_rule_base_records
        ),
        "edge_rule_valid_generators": ";".join(
            format_generator(generator, decisions) for generator in edge_rule_valid_records
        ),
        "edge_rule_valid_span_masks": ";".join(
            format_mask(mask, decisions) for mask in sorted(edge_rule_valid_span)
        ),
        "violation_kernel_valid_masks": ";".join(
            format_mask(mask, decisions) for mask in violation_kernel_valid_masks
        ),
        "labelled_violation_kernel_valid_masks": ";".join(
            format_mask(mask, decisions) for mask in labelled_violation_kernel_valid_masks
        ),
        "full_valid_stabilizer_masks": ";".join(
            format_mask(mask, decisions) for mask in full_valid_stabilizer
        ),
    }
    row.update(base_clique_metrics(instance, local_vertices))
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, default=Path("research/multi-pruning-instances.jsonl"))
    parser.add_argument("--out", type=Path, default=Path("research/multi-pruning-metrics.csv"))
    parser.add_argument("--tol", type=float, default=1e-5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = []
    with args.instances.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            instance, pruning_edges = load_multi_instance(json.loads(line))
            rows.append(multi_pruning_row(instance, pruning_edges, args.tol))

    if not rows:
        raise RuntimeError(f"No instances found in {args.instances}.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} row(s) to {args.out}")


if __name__ == "__main__":
    main()
