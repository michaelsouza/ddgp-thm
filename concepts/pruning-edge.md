---
term: Pruning Edge
tags: [geometry, combinatorial-optimization]
sources: [goncalves2021new, liberti2011number]
---

## Intuition

A pruning edge is an edge whose distance is not used to generate a candidate position, but is used to test whether that candidate is valid. In algorithms such as [Branch-and-Prune](branch-and-prune.md), discretization edges create the binary branching structure, while pruning edges remove branches whose coordinates violate additional distance constraints.

In the [DMDGP](discretizable-molecular-distance-geometry-problem.md), after $v_i$ is placed from its $K$ immediate predecessors, any edge from $v_i$ to an earlier vertex outside that immediate predecessor window acts as a pruning edge. In the [DDGP](discretizable-distance-geometry-problem.md), the same concept is relative to the chosen set of $K$ predecessors: an earlier edge not selected as a [discretization edge](discretization-edge.md) becomes a pruning constraint.

Pruning edges are essential for reducing the number of candidate realizations. Without them, a DMDGP instance generically yields a full binary tree with up to $2^{n-K}$ leaves after the initial $K$ vertices are fixed.

## Formal definition

Let $G = (V, E, d)$ be a weighted graph with discretization order $v_1, \dots, v_n$. For each $i > K$, let $U_i$ be the set of $K$ predecessors used to place $v_i$, and let

$$
E_D(i) = \{\{u, v_i\} \in E \mid u \in U_i\}
$$

be the discretization edges ending at $v_i$. A pruning edge ending at $v_i$ is an edge

$$
\{u, v_i\} \in E
$$

such that

$$
u \in \{v_1, \dots, v_{i-1}\}
\quad\text{and}\quad
u \notin U_i.
$$

The set of pruning edges is therefore

$$
E_P =
\{\{u, v_i\} \in E \mid i > K,\ u \in \{v_1, \dots, v_{i-1}\},\ u \notin U_i\}.
$$

In the DMDGP special case, where $U_i = \{v_{i-1}, \dots, v_{i-K}\}$, pruning edges are exactly the edges

$$
\{v_h, v_i\} \in E
\quad\text{with}\quad
h < i-K.
$$

Given a candidate position $x_i$, such an edge is checked by verifying

$$
\left|\|x_u - x_i\| - d_{ui}\right| < \varepsilon,
$$

for a numerical tolerance $\varepsilon$.
