---
term: Violation Map
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The violation map sends a generator presentation to its total labelled
obstruction vector.

Its kernel consists of generator combinations whose labelled violations cancel
over $\mathbb F_2$. After projection through the mask map, this kernel gives the
rank-count stabilizer space.

## Formal definition

Let $\mathcal G$ be the graph-generator family and let
$\nu_F(g)\in\mathbb F_2^{\mathcal L_F}$ be the labelled violation vector of
generator $g$. The violation map is

$$
V_F:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{\mathcal L_F}
$$

with

$$
V_F\alpha
=
\bigoplus_{g\in\operatorname{supp}(\alpha)}\nu_F(g).
$$

The matrix of this map is the [labelled violation
matrix](labelled-violation-matrix.md).
