---
term: Branch-and-Prune (BP) Algorithm
tags: [algorithms, combinatorial-optimization]
sources: [goncalves2021new]
---

## Intuition

The Branch-and-Prune (BP) algorithm is the standard method for solving Discretizable Molecular Distance Geometry Problems (DMDGP). Because the DMDGP search space is discrete, it can be represented as a binary tree of candidate positions. The BP algorithm systematically explores this tree using a depth-first search (DFS).

At each level of the tree (representing a vertex $v_i$), the algorithm:
1. **Branches**: Generates exactly two candidate positions $\{x_i^+, x_i^-\}$ using $K$-lateration (circle/sphere intersection).
2. **Prunes**: Checks if the candidate positions satisfy any known distances associated with [pruning edges](pruning-edge.md) (in the DMDGP, edges $\{v_h, v_i\}$ where $h < i - K$). If a candidate position violates a distance constraint within a given tolerance $\varepsilon$, the algorithm prunes that entire branch of the search tree immediately, avoiding the need to explore its descendants.

If a leaf node is successfully reached (level $n$), a valid full realization of the graph has been found. The algorithm can then backtrack to find alternative solutions, or terminate. BP is highly efficient because pruning edges occur frequently in real-world instances like protein conformations, dramatically cutting down the exponential size of the tree.

## Formal definition

Let $(G, K)$ be a DMDGP instance with $G = (V, E, d)$, and let $E = E_D \cup E_P$ be the partition of edges into [discretization edges](discretization-edge.md) $E_D$ and [pruning edges](pruning-edge.md) $E_P$. 

The BP algorithm can be recursively defined for step $i$ (with $i > K$, after fixing the initial clique $x_1, \dots, x_K$):

1. **Base Case**: If $i > n$, return the valid realization $x = (x_1, \dots, x_n)$.
2. **Lateration Step**: Find the two solutions $\{x_i^+, x_i^-\} \subset \mathbb{R}^K$ to the system:
   $$\| x_\ell - x_i \|^2 = d_{\ell,i}^2, \quad \forall \ell \in \{i-K, \dots, i-1\}$$
3. **Feasibility and Recursive Call**:
   * For the first candidate position $x_i = x_i^+$:
     * Check if it is feasible under all pruning constraints ending at $i$:
       $$\forall h < i \text{ such that } \{v_h, v_i\} \in E_P, \quad \left| \| x_h - x_i^+ \| - d_{hi} \right| < \varepsilon$$
     * If feasible, recursively call $BP(i+1, n, G, x)$.
   * For the second candidate position $x_i = x_i^-$:
     * Check if it is feasible:
       $$\forall h < i \text{ such that } \{v_h, v_i\} \in E_P, \quad \left| \| x_h - x_i^- \| - d_{hi} \right| < \varepsilon$$
     * If feasible, recursively call $BP(i+1, n, G, x)$.
