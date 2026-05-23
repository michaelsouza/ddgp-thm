---
term: Local Binary Decision Set
tags: [geometry, combinatorial-optimization]
sources: []
---

## Intuition

The local binary decision set is the set of non-seed vertices whose BP branch
choices affect a local DDGP subproblem.

It removes the fixed initial simplex from the local vertex set. Each remaining
vertex contributes one binary coordinate to the local sign cover.

## Formal definition

Let $L_P$ be the [local vertex set](local-vertex-set.md) and let

$$
V_0=\{v_1,\dots,v_{K+1}\}
$$

be the [initial rigid simplex](initial-rigid-simplex.md). The local binary
decision set is

$$
B_P=L_P\setminus V_0.
$$

The local BP sign cover is $\mathbb F_2^{B_P}$, and a local solution code is a
subset

$$
\Xi_F\subseteq\mathbb F_2^{B_P}.
$$
