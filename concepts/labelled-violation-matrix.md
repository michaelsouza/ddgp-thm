---
term: Labelled Violation Matrix
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The labelled violation matrix records which active edges are generically broken
by each graph generator, with each violation labelled by the mirror clique that
caused it.

This label is essential. Two generator presentations may create the same edge
crossing but through different mirror cliques, and those should not be treated
as the same geometric obstruction.

## Formal definition

Let $\mathcal G$ be the graph-generator family and let $F$ be the
[active edge set](active-edge-set.md). The [labelled violation
set](labelled-violation-set.md) is

$$
\mathcal L_F=\{(e,C_g)\mid g\in\mathcal G,\ e\in F,\ g\text{ violates }e\}.
$$

For each generator $g$, define its labelled violation vector
$\nu_F(g)\in\mathbb F_2^{\mathcal L_F}$ by

$$
\nu_F(g)_{(e,C)}=1
\quad\Longleftrightarrow\quad
C=C_g \text{ and } g\text{ violates }e.
$$

The labelled violation matrix $V_F$ has columns $\nu_F(g)$:

$$
V_F:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{\mathcal L_F}.
$$

Its kernel consists of generator combinations with zero total labelled
obstruction.
