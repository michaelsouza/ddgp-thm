---
term: BP Sign Cover
tags: [geometry, algorithms, group-theory]
sources: []
---

## Intuition

The BP sign cover is the binary space of Branch-and-Prune choices before
pruning constraints are imposed.

Each non-seed vertex has two generic lateration positions once its predecessors
are fixed. Choosing one of the two positions is encoded by a sign or bit. The
rank-count theory studies transformations of these sign strings rather than
only transformations of final Euclidean coordinates.

## Formal definition

Let $B_P$ be the [local binary decision set](local-binary-decision-set.md). The
local BP sign cover is the vector space

$$
\mathbb F_2^{B_P}.
$$

A sign string $s\in\mathbb F_2^{B_P}$ determines a local BP realization whenever
all required lateration steps are feasible. For a branch mask
$h\in\mathbb F_2^{B_P}$, the BP sign transformation is

$$
\tau_h(s)=s\oplus h.
$$

The [local solution code](local-solution-code.md) $\Xi_F$ is the subset of the
sign cover whose realizations satisfy every edge in the
[active edge set](active-edge-set.md) $F$.
