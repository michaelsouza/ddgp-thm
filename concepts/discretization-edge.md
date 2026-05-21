---
term: Discretization Edge
tags: [geometry, combinatorial-optimization]
sources: [goncalves2021new, liberti2011number]
---

## Intuition

A discretization edge is an edge whose distance is used to place a new vertex during the sequential construction of a realization. In discretizable distance geometry, these edges are what turn a continuous search over coordinates into a finite branching process: once the first $K+1$ vertices are fixed as an initial rigid simplex, each later vertex is located from distances to $K$ already placed vertices in $\mathbb{R}^K$.

In the [DMDGP](discretizable-molecular-distance-geometry-problem.md), the discretization edges have a particularly simple form. For each vertex $v_i$ with $i > K$, they connect $v_i$ to its $K$ immediate predecessors $v_{i-1}, \dots, v_{i-K}$. In the adjacency matrix ordered by the DMDGP order, these are the first $K$ diagonals away from the main diagonal.

In the [DDGP](discretizable-distance-geometry-problem.md), the same idea is more general: the $K$ predecessors used to place $v_i$ do not need to be immediate predecessors. A discretization edge is therefore relative to the chosen predecessor set used in the construction, which must satisfy the required simplex nondegeneracy conditions. These choices define the [discretization dependency set](discretization-dependency-set.md) of each vertex.

## Formal definition

Let $G = (V, E, d)$ be a weighted graph embedded in $\mathbb{R}^K$, and let $v_1, \dots, v_n$ be a discretization order. The first $K+1$ vertices form the fixed initial simplex. For each $i > K+1$, choose a set

$$
U_i \subseteq \{v_1, \dots, v_{i-1}\}, \qquad |U_i| = K,
$$

of predecessors used to locate $v_i$. The discretization edges associated with $v_i$ are

$$
E_D(i) = \{\{u, v_i\} \in E \mid u \in U_i\}.
$$

For a DDGP discretization, $U_i$ is chosen so that it satisfies the [Strict Simplex Inequalities](strict-simplex-inequalities.md), ensuring that the sphere intersection used to place $v_i$ is nondegenerate.

The full set of discretization edges is

$$
E_D = \bigcup_{i=K+2}^{n} E_D(i),
$$

together with the clique edges on the fixed initial simplex $\{v_1,\dots,v_{K+1}\}$.

In the DMDGP special case,

$$
U_i = \{v_{i-1}, \dots, v_{i-K}\},
$$

so

$$
E_D(i) = \{\{v_{i-j}, v_i\} \mid j = 1, \dots, K\}.
$$

The initial clique on $v_1, \dots, v_{K+1}$ supplies the seed distances needed to fix the coordinate system and remove the trivial global reflection; depending on convention, those seed edges may be included together with $E_D$.
