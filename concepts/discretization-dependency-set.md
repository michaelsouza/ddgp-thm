---
term: Discretization Dependency Set
tags: [geometry, combinatorial-optimization]
sources: [goncalves2021new, liberti2011number]
---

## Intuition

The discretization dependency set of a vertex is the set of vertices whose positions must already be fixed in order to fix that vertex through the chosen discretization construction. It is also useful to think of it as the fixing set of a vertex, relative to a selected family of predecessor sets $U_i$.

This concept is especially important in the [DDGP](discretizable-distance-geometry-problem.md). In these notes, the first $K+1$ vertices are fixed as an initial rigid simplex and have no associated binary decisions. In the DDGP, the $K$ predecessors used to place a later vertex $v_i$ may be nonconsecutive, so many earlier non-seed vertices may be irrelevant to the position of $v_i$.

For algorithms inspired by [SBBU](symmetry-based-build-up.md), this replaces linear intervals of vertices by dependency sets. Instead of assuming that every vertex before $v_i$ matters, one can track only the ancestors of $v_i$ in the discretization dependency graph.

A complementary derived object is the [base-clique transition graph](base-clique-transition-graph.md), whose vertices are the predecessor cliques $U_i$ and whose directed edges record when a vertex generated from one base clique is reused by another.

The forward descendants of a vertex in this DAG form a dependency cone. For a pruning edge, cones containing exactly one endpoint are [endpoint-separating cones](endpoint-separating-cone.md).

## Formal definition

Let $G = (V, E, d)$ be a DDGP instance in $\mathbb{R}^K$ with discretization order $v_1, \dots, v_n$. The initial rigid simplex is

$$
V_0 = \{v_1, \dots, v_{K+1}\}.
$$

For each $i > K+1$, let

$$
U_i \subseteq \{v_1, \dots, v_{i-1}\}, \qquad |U_i| = K,
$$

be the predecessor set whose [discretization edges](discretization-edge.md) are used to place $v_i$.

Define the discretization dependency digraph

$$
D_U = (V, A_U),
$$

where

$$
(u, v_i) \in A_U
\quad\Longleftrightarrow\quad
i > K+1 \text{ and } u \in U_i.
$$

Because every arc points from an earlier vertex to a later vertex, $D_U$ is acyclic.

The ancestor set of $v_i$ is

$$
\operatorname{Anc}_U(v_i)
=
\{v_h \in V \mid \text{there is a directed path } v_h \leadsto v_i \text{ in } D_U\}.
$$

The discretization dependency set, or fixing set, of $v_i$ is

$$
\operatorname{Fix}_U(v_i)
=
\operatorname{Anc}_U(v_i) \cup \{v_i\}.
$$

Equivalently, it is defined recursively by

$$
\operatorname{Fix}_U(v_i)
=
\{v_i\}
\cup
\bigcup_{u \in U_i} \operatorname{Fix}_U(u),
\qquad i > K+1,
$$

with base case

$$
\operatorname{Fix}_U(v_i) = \{v_i\},
\qquad i \le K+1,
$$

assuming the initial rigid simplex has already been fixed as the coordinate seed.

For a [pruning edge](pruning-edge.md) $\{v_i, v_j\}$, the vertices needed to evaluate its distance constraint are

$$
\operatorname{Fix}_U(v_i, v_j)
=
\operatorname{Fix}_U(v_i) \cup \operatorname{Fix}_U(v_j).
$$

Only vertices in $\operatorname{Fix}_U(v_i) \setminus V_0$ carry binary decisions relevant to $v_i$; the vertices in $V_0$ are fixed.
