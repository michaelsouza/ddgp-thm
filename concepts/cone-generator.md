---
term: Cone Generator
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

A cone generator records the direct effect of toggling one local branch decision
and propagating that change through the DDGP dependency DAG.

It generalizes the suffix partial reflections used in the DMDGP. Instead of a
consecutive suffix, the support is the local part of the vertex's dependency
cone.

## Formal definition

For a decision vertex $q\in B_P$, the cone generator is

$$
g_q=(m_q,U_q),
$$

where $U_q$ is the predecessor clique of $q$ and

$$
\operatorname{supp}(m_q)
=
B_P\cap\operatorname{Cone}_U(q).
$$

Equivalently,

$$
(m_q)_z=1
\quad\Longleftrightarrow\quad
z\in B_P \text{ and } q\leadsto z
\text{ in the dependency DAG}.
$$

When the vertices in $B_P$ are ordered by the DDGP order, the cone masks form an
upper-triangular matrix with ones on the diagonal. Hence they span the full
local branch space $\mathbb F_2^{B_P}$.

_Avoid_: dependency-cone generator. Use **dependency cone** for the graph set
and **cone generator** for the generator built from it.
