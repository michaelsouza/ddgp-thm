---
term: Discretizable Distance Geometry Problem (DDGP)
tags: [geometry, combinatorial-optimization]
sources: [liberti2011number]
---

## Intuition

The Discretizable Distance Geometry Problem (DDGP) is a generalization of the [Discretizable Molecular Distance Geometry Problem (DMDGP)](discretizable-molecular-distance-geometry-problem.md). In both problems, the goal is to find a coordinate embedding of a graph in $\mathbb{R}^K$ by exploiting sequential vertex ordering, which allows the search space to be represented as a discrete binary tree and solved via the [Branch-and-Prune (BP) Algorithm](branch-and-prune.md).

The key difference lies in the choice of predecessors used to locate the next vertex $v$:
* In the **DMDGP**, vertex $v$ must have known distances to its $K$ **immediate** predecessors in the vertex order (i.e., $v_{i-1}, \dots, v_{i-K}$).
* In the **DDGP**, vertex $v$ only needs to have known distances to **any** $K$ predecessors that appear earlier in the vertex order, not necessarily the immediate ones.

While the DDGP encompasses a larger set of instances, this relaxation has a significant theoretical cost. In the DMDGP, the use of immediate predecessors ensures that partial reflections across hyperplanes map valid embeddings to valid embeddings, guaranteeing that the number of solutions is always a power of two with probability 1. In the DDGP, because the hyperplanes are defined by non-adjacent predecessors, reflections do not necessarily preserve validity, meaning the "power of two" solution count property does not generally hold.

## Formal definition

Let $G = (V, E)$ be an undirected graph, $d: E \to \mathbb{R}_{+}$ be an edge weight function, $K > 0$ be the dimension of the embedding space, and $<$ be a total order on $V$. Let $\gamma(v) = \{u \in V \mid u < v\}$ be the set of predecessors of $v$, and let $\rho(v)$ be the rank of $v$ in $<$.

A DGP instance $(G, d, K)$ is in the **DDGP** (or $\text{DDGP}_K$) if there exists a total order $<$ and a subset of initial vertices $V_0 = \{v \in V \mid \rho(v) \le K\}$ such that:

1. **Initial Clique**: The induced subgraph $G[V_0]$ is a clique $\mathbf{K}_K$, and a valid partial embedding $\bar{x}: V_0 \to \mathbb{R}^K$ exists.
2. **Discretization Condition**: For all $v \in V \setminus V_0$:
   $$|N(v) \cap \gamma(v)| \ge K$$
   That is, each vertex $v$ of rank $> K$ is adjacent to at least $K$ of its predecessors.
3. **Simplex Condition**: For each $v \in V \setminus V_0$, there exists a subset $U_v \subseteq N(v) \cap \gamma(v)$ of cardinality $K$ that forms a clique $\mathbf{K}_K$ satisfying [Strict Simplex Inequalities](strict-simplex-inequalities.md):
   $$\Delta_{K-1}(U_v) > 0$$

The decision problem of the DDGP is to determine if there exists a valid embedding $x: V \to \mathbb{R}^K$ that extends the partial embedding $\bar{x}$.
