#!/usr/bin/env python3
"""Plot adjacency matrices and dependency DAGs for DDGP JSONL instances."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap


Edge = tuple[int, int]


@dataclass(frozen=True)
class Instance:
    instance_id: str
    mode: str
    pruning_kind: str
    k: int
    n: int
    predecessor_sets: dict[int, tuple[int, ...]]
    discretization_edges: frozenset[Edge]
    base_edges: frozenset[Edge]
    pruning_edge: Edge


def edge(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


def parse_edge(values: list[int]) -> Edge:
    return edge(values[0] - 1, values[1] - 1)


def load_instance(record: dict) -> Instance:
    return Instance(
        instance_id=record["id"],
        mode=record["mode"],
        pruning_kind=record["pruning_kind"],
        k=record["k"],
        n=record["n"],
        predecessor_sets={
            int(i) - 1: tuple(u - 1 for u in predecessors)
            for i, predecessors in record["predecessor_sets"].items()
        },
        discretization_edges=frozenset(parse_edge(e) for e in record["discretization_edges"]),
        base_edges=frozenset(parse_edge(e) for e in record["base_edges"]),
        pruning_edge=parse_edge(record["pruning_edge"]),
    )


def dependency_sets(instance: Instance) -> list[set[int]]:
    fix_sets = [set([i]) for i in range(instance.n)]
    for i in range(instance.k + 1, instance.n):
        fixed = {i}
        for u in instance.predecessor_sets[i]:
            fixed |= fix_sets[u]
        fix_sets[i] = fixed
    return fix_sets


def dependency_levels(instance: Instance) -> list[int]:
    levels = [0] * instance.n
    for i in range(instance.k + 1, instance.n):
        levels[i] = 1 + max(levels[u] for u in instance.predecessor_sets[i])
    return levels


def base_clique_nodes(instance: Instance) -> dict[tuple[int, ...], list[int]]:
    nodes: dict[tuple[int, ...], list[int]] = {}
    for child in range(instance.k + 1, instance.n):
        clique = tuple(sorted(instance.predecessor_sets[child]))
        nodes.setdefault(clique, []).append(child)
    return nodes


def matrix_codes(instance: Instance, ordered: list[int]) -> list[list[int]]:
    data = []
    clique_edges = set(instance.base_edges) - set(instance.discretization_edges)
    for u in ordered:
        row = []
        for v in ordered:
            e = edge(u, v)
            if u == v:
                row.append(4)
            elif e == instance.pruning_edge:
                row.append(3)
            elif e in instance.discretization_edges:
                row.append(1)
            elif e in clique_edges:
                row.append(2)
            else:
                row.append(0)
        data.append(row)
    return data


def plot_matrix(
    instance: Instance,
    path: Path,
    *,
    title: str,
    subtitle: str,
    vertices: Iterable[int] | None = None,
    fix_u: set[int] | None = None,
    fix_v: set[int] | None = None,
) -> None:
    ordered = list(range(instance.n) if vertices is None else sorted(vertices))
    labels = [str(v + 1) for v in ordered]
    data = matrix_codes(instance, ordered)
    size = len(ordered)
    fig_size = max(5.2, 0.64 * size + 2.2)

    cmap = ListedColormap(["#f8fafc", "#2563eb", "#f59e0b", "#dc2626", "#111827"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)

    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    sns.heatmap(
        data,
        ax=ax,
        cmap=cmap,
        norm=norm,
        cbar=False,
        square=True,
        linewidths=0.7,
        linecolor="#cbd5e1",
        xticklabels=labels,
        yticklabels=labels,
    )

    ax.set_title(f"{title}\n{subtitle}", loc="left", fontsize=10.5, fontweight="bold", pad=13)
    ax.tick_params(axis="x", rotation=0, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=9)
    ax.set_xlabel("vertex")
    ax.set_ylabel("vertex")

    def membership_color(v: int) -> str:
        in_u = fix_u is not None and v in fix_u
        in_v = fix_v is not None and v in fix_v
        if in_u and in_v:
            return "#0f766e"
        if in_u:
            return "#16a34a"
        if in_v:
            return "#f97316"
        return "#334155"

    for tick, v in zip(ax.get_xticklabels(), ordered):
        tick.set_color(membership_color(v))
        tick.set_fontweight("bold")
    for tick, v in zip(ax.get_yticklabels(), ordered):
        tick.set_color(membership_color(v))
        tick.set_fontweight("bold")

    positions = {v: idx for idx, v in enumerate(ordered)}
    u, v = instance.pruning_edge
    if u in positions and v in positions:
        for row, col in [(positions[u], positions[v]), (positions[v], positions[u])]:
            ax.add_patch(
                mpatches.Rectangle(
                    (col, row),
                    1,
                    1,
                    fill=False,
                    edgecolor="#9333ea",
                    linewidth=2.8,
                )
            )

    legend_items = [
        mpatches.Patch(color="#2563eb", label="discretization"),
        mpatches.Patch(color="#f59e0b", label="required clique"),
        mpatches.Patch(color="#dc2626", label="single pruning"),
        mpatches.Patch(color="#f8fafc", label="absent"),
        mpatches.Patch(color="#111827", label="diagonal"),
    ]
    if fix_u is not None and fix_v is not None:
        legend_items.extend(
            [
                mpatches.Patch(color="#16a34a", label="Fix(u)"),
                mpatches.Patch(color="#f97316", label="Fix(v)"),
                mpatches.Patch(color="#0f766e", label="intersection"),
            ]
        )
    ax.legend(
        handles=legend_items,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.11),
        ncol=3,
        frameon=False,
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_dependency_dag(
    instance: Instance,
    path: Path,
    *,
    title: str,
    fix_u: set[int],
    fix_v: set[int],
) -> None:
    levels = dependency_levels(instance)
    by_level: dict[int, list[int]] = {}
    for vertex, level in enumerate(levels):
        by_level.setdefault(level, []).append(vertex)

    max_level = max(by_level)
    max_width = max(len(vertices) for vertices in by_level.values())
    fig_width = max(7.0, 1.55 * (max_level + 1))
    fig_height = max(4.8, 0.85 * max_width + 2.4)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    positions: dict[int, tuple[float, float]] = {}
    for level, vertices in by_level.items():
        vertices = sorted(vertices)
        offset = (len(vertices) - 1) / 2.0
        for idx, vertex in enumerate(vertices):
            positions[vertex] = (level, offset - idx)

    def node_color(vertex: int) -> str:
        in_u = vertex in fix_u
        in_v = vertex in fix_v
        if in_u and in_v:
            return "#0f766e"
        if in_u:
            return "#16a34a"
        if in_v:
            return "#f97316"
        if vertex <= instance.k:
            return "#64748b"
        return "#e2e8f0"

    for child, predecessors in instance.predecessor_sets.items():
        if child <= instance.k:
            continue
        x1, y1 = positions[child]
        for parent in predecessors:
            x0, y0 = positions[parent]
            ax.annotate(
                "",
                xy=(x1 - 0.13, y1),
                xytext=(x0 + 0.13, y0),
                arrowprops={
                    "arrowstyle": "->",
                    "color": "#64748b",
                    "linewidth": 1.2,
                    "shrinkA": 8,
                    "shrinkB": 8,
                },
            )

    u, v = instance.pruning_edge
    if u in positions and v in positions:
        x0, y0 = positions[u]
        x1, y1 = positions[v]
        ax.annotate(
            "",
            xy=(x1, y1 + 0.12),
            xytext=(x0, y0 + 0.12),
            arrowprops={
                "arrowstyle": "-",
                "color": "#dc2626",
                "linewidth": 2.4,
                "linestyle": "--",
                "connectionstyle": "arc3,rad=0.22",
            },
        )

    for vertex, (x, y) in positions.items():
        circle = mpatches.Circle(
            (x, y),
            radius=0.18,
            facecolor=node_color(vertex),
            edgecolor="#0f172a",
            linewidth=1.0,
            zorder=3,
        )
        ax.add_patch(circle)
        ax.text(
            x,
            y,
            str(vertex + 1),
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color="white" if node_color(vertex) in {"#0f766e", "#16a34a", "#f97316", "#64748b"} else "#0f172a",
            zorder=4,
        )

    for level in range(max_level + 1):
        ax.axvline(level, color="#e2e8f0", linewidth=0.8, zorder=0)
        ax.text(
            level,
            -max_width / 2 - 0.55,
            f"level {level}",
            ha="center",
            va="top",
            fontsize=9,
            color="#475569",
        )

    ax.set_title(f"{title}\ndiscretization dependency DAG by level", loc="left", fontsize=10.5, fontweight="bold", pad=12)
    ax.set_xlim(-0.6, max_level + 0.6)
    ax.set_ylim(-max_width / 2 - 0.9, max_width / 2 + 0.7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    legend_items = [
        mpatches.Patch(color="#64748b", label="V0 / level 0"),
        mpatches.Patch(color="#16a34a", label="Fix(u)"),
        mpatches.Patch(color="#f97316", label="Fix(v)"),
        mpatches.Patch(color="#0f766e", label="intersection"),
        mpatches.Patch(color="#e2e8f0", label="outside both"),
        mpatches.Patch(color="#dc2626", label="pruning edge"),
    ]
    ax.legend(
        handles=legend_items,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.04),
        ncol=3,
        frameon=False,
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_base_clique_graph(
    instance: Instance,
    path: Path,
    *,
    title: str,
    fix_u: set[int],
    fix_v: set[int],
) -> None:
    vertex_levels = dependency_levels(instance)
    nodes = base_clique_nodes(instance)
    if not nodes:
        return

    node_order = sorted(nodes, key=lambda clique: (max(vertex_levels[v] for v in clique), clique))
    clique_level = {clique: max(vertex_levels[v] for v in clique) for clique in node_order}
    by_level: dict[int, list[tuple[int, ...]]] = {}
    for clique in node_order:
        by_level.setdefault(clique_level[clique], []).append(clique)

    max_level = max(by_level)
    max_width = max(len(cliques) for cliques in by_level.values())
    fig_width = max(8.2, 1.8 * (max_level + 1))
    fig_height = max(5.0, 1.05 * max_width + 2.8)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    positions: dict[tuple[int, ...], tuple[float, float]] = {}
    for level, cliques in by_level.items():
        cliques = sorted(cliques)
        offset = (len(cliques) - 1) / 2.0
        for idx, clique in enumerate(cliques):
            positions[clique] = (level, offset - idx)

    child_to_base: dict[int, tuple[int, ...]] = {}
    for clique, children in nodes.items():
        for child in children:
            child_to_base[child] = clique

    transition_edges: dict[tuple[tuple[int, ...], tuple[int, ...]], set[int]] = {}
    for target in node_order:
        for vertex in target:
            source = child_to_base.get(vertex)
            if source is not None and source != target:
                transition_edges.setdefault((source, target), set()).add(vertex)

    overlap_edges: list[tuple[tuple[int, ...], tuple[int, ...], int]] = []
    transition_pairs = {edge_pair for edge_pair in transition_edges}
    transition_pairs |= {(target, source) for source, target in transition_edges}
    for idx, left in enumerate(node_order):
        for right in node_order[idx + 1 :]:
            shared = len(set(left) & set(right))
            if shared >= instance.k - 1 and shared > 0 and (left, right) not in transition_pairs:
                overlap_edges.append((left, right, shared))

    def clique_color(clique: tuple[int, ...]) -> str:
        vertices = set(clique)
        in_u = vertices <= fix_u
        in_v = vertices <= fix_v
        if in_u and in_v:
            return "#0f766e"
        if in_u:
            return "#16a34a"
        if in_v:
            return "#f97316"
        if vertices <= set(range(instance.k + 1)):
            return "#64748b"
        return "#e2e8f0"

    for left, right, shared in overlap_edges:
        x0, y0 = positions[left]
        x1, y1 = positions[right]
        ax.plot(
            [x0, x1],
            [y0, y1],
            color="#94a3b8",
            linewidth=0.8 + 0.35 * shared,
            linestyle=":",
            zorder=1,
        )

    for (source, target), carrier_vertices in transition_edges.items():
        x0, y0 = positions[source]
        x1, y1 = positions[target]
        ax.annotate(
            "",
            xy=(x1 - 0.16, y1),
            xytext=(x0 + 0.16, y0),
            arrowprops={
                "arrowstyle": "->",
                "color": "#334155",
                "linewidth": 1.1 + 0.35 * (len(carrier_vertices) - 1),
                "shrinkA": 9,
                "shrinkB": 9,
                "connectionstyle": "arc3,rad=0.06",
            },
            zorder=2,
        )
        if len(carrier_vertices) > 1:
            ax.text(
                (x0 + x1) / 2,
                (y0 + y1) / 2 + 0.12,
                ",".join(str(v + 1) for v in sorted(carrier_vertices)),
                ha="center",
                va="center",
                fontsize=7,
                color="#334155",
                bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
                zorder=5,
            )

    pruning_vertices = set(instance.pruning_edge)
    for clique, (x, y) in positions.items():
        color = clique_color(clique)
        contains_pruning_endpoint = bool(set(clique) & pruning_vertices)
        rect = mpatches.FancyBboxPatch(
            (x - 0.36, y - 0.21),
            0.72,
            0.42,
            boxstyle="round,pad=0.04,rounding_size=0.05",
            facecolor=color,
            edgecolor="#dc2626" if contains_pruning_endpoint else "#0f172a",
            linewidth=2.0 if contains_pruning_endpoint else 1.0,
            zorder=3,
        )
        ax.add_patch(rect)
        label = "{" + ",".join(str(v + 1) for v in clique) + "}"
        children = ",".join(str(child + 1) for child in nodes[clique])
        ax.text(
            x,
            y + 0.045,
            label,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            color="white" if color in {"#0f766e", "#16a34a", "#f97316", "#64748b"} else "#0f172a",
            zorder=4,
        )
        ax.text(
            x,
            y - 0.105,
            f"-> {children}",
            ha="center",
            va="center",
            fontsize=6.8,
            color="white" if color in {"#0f766e", "#16a34a", "#f97316", "#64748b"} else "#475569",
            zorder=4,
        )

    for level in range(max_level + 1):
        ax.axvline(level, color="#e2e8f0", linewidth=0.8, zorder=0)
        ax.text(
            level,
            -max_width / 2 - 0.7,
            f"clique level {level}",
            ha="center",
            va="top",
            fontsize=9,
            color="#475569",
        )

    ax.set_title(
        f"{title}\nbase-clique transition graph",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlim(-0.7, max_level + 0.7)
    ax.set_ylim(-max_width / 2 - 1.05, max_width / 2 + 0.8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    legend_items = [
        mpatches.Patch(color="#64748b", label="seed-only base"),
        mpatches.Patch(color="#16a34a", label="base in Fix(u)"),
        mpatches.Patch(color="#f97316", label="base in Fix(v)"),
        mpatches.Patch(color="#0f766e", label="base in both"),
        mpatches.Patch(color="#e2e8f0", label="outside both"),
        mpatches.Patch(facecolor="white", edgecolor="#dc2626", linewidth=2.0, label="uses pruning endpoint"),
        mpatches.Patch(facecolor="white", edgecolor="#334155", label="arrow: generated vertex reused"),
        mpatches.Patch(facecolor="white", edgecolor="#94a3b8", label="dotted: large overlap"),
    ]
    ax.legend(
        handles=legend_items,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.05),
        ncol=3,
        frameon=False,
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def load_counts(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    with path.open(newline="", encoding="utf-8") as fh:
        return {row["id"]: row for row in csv.DictReader(fh)}


def plot_instance(instance: Instance, counts: dict[str, str], out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    fix_sets = dependency_sets(instance)
    u, v = instance.pruning_edge
    union = fix_sets[u] | fix_sets[v]

    count_text = ""
    if counts:
        count_text = (
            f" | solutions={counts['valid_count']} base={counts['base_count']} "
            f"local={counts.get('local_valid_count', '?')} fraction={counts['valid_fraction_of_base']}"
        )

    title = (
        f"{instance.instance_id}: {instance.pruning_kind}, "
        f"pruning edge ({u + 1},{v + 1}){count_text}"
    )

    paths = []
    full_path = out_dir / f"{instance.instance_id}-full.png"
    plot_matrix(
        instance,
        full_path,
        title=title,
        subtitle="full adjacency matrix",
    )
    paths.append(full_path)

    fix_path = out_dir / f"{instance.instance_id}-fix-union.png"
    plot_matrix(
        instance,
        fix_path,
        title=title,
        subtitle=f"induced matrix on Fix({u + 1}) union Fix({v + 1})",
        vertices=union,
        fix_u=fix_sets[u],
        fix_v=fix_sets[v],
    )
    paths.append(fix_path)

    dag_path = out_dir / f"{instance.instance_id}-dependency-dag.png"
    plot_dependency_dag(
        instance,
        dag_path,
        title=title,
        fix_u=fix_sets[u],
        fix_v=fix_sets[v],
    )
    paths.append(dag_path)

    base_clique_path = out_dir / f"{instance.instance_id}-base-clique-graph.png"
    plot_base_clique_graph(
        instance,
        base_clique_path,
        title=title,
        fix_u=fix_sets[u],
        fix_v=fix_sets[v],
    )
    paths.append(base_clique_path)
    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, required=True)
    parser.add_argument("--counts", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts_by_id = load_counts(args.counts)
    paths = []
    with args.instances.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            instance = load_instance(json.loads(line))
            paths.extend(plot_instance(instance, counts_by_id.get(instance.instance_id, {}), args.out_dir))

    print(f"wrote {len(paths)} image(s) to {args.out_dir}")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
