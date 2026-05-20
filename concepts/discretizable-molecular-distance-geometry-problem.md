---
term: Discretizable Molecular Distance Geometry Problem (DMDGP)
tags: [geometry, combinatorial-optimization]
sources: [goncalves2021new]
---

## Intuition

The Discretizable Molecular Distance Geometry Problem (DMDGP) is a subclass of the Distance Geometry Problem (DGP) where the search space can be discretized and represented as a binary tree. This transformation converts a difficult continuous global optimization problem into a discrete search problem that can be systematically explored.

Discretization is achieved by ordering the vertices of the graph such that each vertex (starting from the $K+1$-th) has known distances to its $K$ immediate predecessors. These distances correspond to [discretization edges](discretization-edge.md).
* For $K=1$, the next point is constrained by a distance to 1 predecessor, leading to 2 possible positions on a line.
* For $K=2$ (in a 2D plane), the next point's distance from 2 predecessors is determined by the intersection of 2 circles, which yields at most 2 candidate positions.
* For $K=3$ (in a 3D space), the next point's distance from 3 predecessors is determined by the intersection of 3 spheres, which yields at most 2 candidate positions.

If the predecessors are not collinear (or in general, do not lie in a lower-dimensional affine subspace), we are guaranteed to have exactly 2 candidate positions at each step. By building these positions sequentially, we construct a binary tree of candidate realizations of size $2^{|V|-K}$. Any additional distance constraints, called [pruning edges](pruning-edge.md), can then be used to prune branches of this tree.

## Formal definition

A DGP instance $(G, K)$ with $G = (V, E, d)$ is a $^K$DMDGP if there exists a vertex order $v_1, \dots, v_n \in V$ such that:

1. **Initial Clique**: The induced subgraph $G[\{v_1, \dots, v_K\}]$ is a clique (all pairwise distances are known, which allows fixing the coordinate system of the first $K$ vertices).
2. **Discretization Order**: For every $i > K$, $v_i$ is adjacent to $v_{i-1}, \dots, v_{i-K}$. The edges $\{v_{i-j}, v_i\}$ for $j = 1, \dots, K$ are called [discretization edges](discretization-edge.md) ($E_D$).
3. **Non-collinearity**: The Cayley-Menger determinant of the $K$ predecessors is non-zero:
   $$CM(v_{i-1}, \dots, v_{i-K}) \neq 0$$
   This condition ensures that the realizations of $v_{i-1}, \dots, v_{i-K}$ span an affine subspace of dimension $K-1$ (a hyperplane in $\mathbb{R}^K$), guaranteeing that the intersection of the $K$ spheres has exactly two distinct candidate positions $\{x_i^+, x_i^-\}$.
