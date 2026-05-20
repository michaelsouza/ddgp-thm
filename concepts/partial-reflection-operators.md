---
term: Partial Reflection Operators
tags: [geometry, algebra]
sources: [goncalves2021new]
---

## Intuition

In a discretizable distance geometry problem (DMDGP), every new vertex position $x_i$ (for $i > K$) is found by intersecting spheres centered at its $K$ predecessors, which gives two possible candidate positions. These two positions are reflections of each other across the hyperplane defined by the $K$ predecessors.

If we choose the "wrong" candidate position at step $i$, it will affect not only $x_i$, but also all subsequent positions $x_{i+1}, \dots, x_n$ because they are built sequentially relative to $x_i$. A **partial reflection operator** $g_i$ is a mathematical operator that reflects the entire "tail" of a realization (all vertices from $i$ to $n$) across the hyperplane defined by $x_{i-K}, \dots, x_{i-1}$.

Applying a partial reflection is highly powerful because it preserves all the [discretization edge](discretization-edge.md) distances (both within the reflected tail, and across the boundary from the predecessors to $x_i$). However, it *changes* the distances between the reflected tail and the unreflected head (vertices $1, \dots, i-K-1$). Thus, partial reflections are used to search for realizations that satisfy [pruning edges](pruning-edge.md), which bridge across these boundaries, without breaking any of the local discretization constraints.

## Formal definition

Let $x = (x_1, \dots, x_n)$ be a possible realization in $\mathbb{R}^K$ satisfying all discretization constraints. For each $i > K$, let $R^i_x(y)$ represent the reflection of a point $y \in \mathbb{R}^K$ through the hyperplane defined by the positions of the predecessors $x_{i-K}, \dots, x_{i-1}$:

$$
R^i_x(y) = (I - 2p_i p_i^T)(y - x_{i-1}) + x_{i-1}
$$

Where:
* $p_i \in \mathbb{R}^K$ is the unit normal vector to the hyperplane defined by $x_{i-K}, \dots, x_{i-1}$.
* $I$ is the identity matrix in $\mathbb{R}^{K \times K}$.

The **partial reflection operator** $g_i(x)$ is defined as:

$$
g_i(x) = (x_1, \dots, x_{i-1}, R^i_x(x_i), R^i_x(x_{i+1}), \dots, R^i_x(x_n))
$$

Key properties of $g_i$:
1. $g_i$ preserves all discretization distances: $g_i(x)$ remains a possible realization satisfying all discretization constraints.
2. Compositions of partial reflections satisfy a commutativity-like property across different hyperplanes (Proposition 1 in the reference).
3. Any possible realization $x' \in \hat{X}$ can be generated from $x \in \hat{X}$ by a unique composition:
   $$x' = g_{K+1}^{s_1} \circ g_{K+2}^{s_2} \circ \dots \circ g_n^{s_{n-K}}(x)$$
   Where $s \in \{0, 1\}^{n-K}$ is a binary vector representing the sequence of reflection choices.
