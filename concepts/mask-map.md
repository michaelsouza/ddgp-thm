---
term: Mask Map
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The mask map sends a generator presentation to the total branch mask it applies
to the BP sign cover.

It forgets the labeled violation data and keeps only the combined binary
effect on local branch decisions.

## Formal definition

Let $\mathcal G$ be the graph-generator family and let
$m_g\in\mathbb F_2^{B_P}$ be the branch mask of generator $g$. The mask map is

$$
M:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{B_P}
$$

with

$$
M\alpha
=
\bigoplus_{g\in\operatorname{supp}(\alpha)}m_g.
$$

The matrix of this map is the [generator mask matrix](generator-mask-matrix.md).
