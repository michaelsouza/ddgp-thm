---
term: Dependency Cone
tags: [geometry, combinatorial-optimization]
sources: []
---

## Intuition

The dependency cone of a vertex is the forward set of vertices whose coordinates
depend on that vertex's branch choice.

In the DMDGP, changing a branch choice affects a suffix of the vertex order. In
the general DDGP, predecessor sets can be nonconsecutive, so the affected set is
instead a descendant cone in the discretization dependency DAG.

## Formal definition

Let $D_U$ be the discretization dependency DAG with an arc $u\to i$ whenever
$u\in U_i$. For a vertex $q$, the dependency cone is

$$
\operatorname{Cone}_U(q)
=
\{q\}\cup\{z\mid q\leadsto z \text{ in } D_U\}.
$$

Equivalently, it satisfies the recursion

$$
\operatorname{Cone}_U(q)
=
\{q\}\cup
\bigcup_{w:q\in U_w}\operatorname{Cone}_U(w).
$$

Intersecting this set with the local decision set $B_P$ gives the support of a
[cone generator](cone-generator.md).
