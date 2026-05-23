---
term: Labelled Violation Set
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The labelled violation set is the coordinate system used to record active-edge
obstructions.

Instead of labelling a violation only by the edge that is crossed, the
rank-count theory labels it by both the edge and the mirror clique. This keeps
track of the geometric reason for the obstruction.

## Formal definition

Let $F$ be an [active edge set](active-edge-set.md), and let
$g=(m_g,C_g)$ be a DDGP graph generator with support
$S_g=\operatorname{supp}(m_g)$. The generator $g$ violates an edge
$e=\{a,b\}\in F$ when

$$
|e\cap S_g|=1
\quad\text{and}\quad
e\setminus S_g\not\subseteq C_g.
$$

The labelled violation set is

$$
\mathcal L_F
=
\{(e,C_g)\mid g\in\mathcal G,\ e\in F,\ g\text{ violates }e\}.
$$

It indexes the rows of the [labelled violation matrix](labelled-violation-matrix.md).
