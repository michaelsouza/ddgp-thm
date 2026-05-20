---
term: Symmetry-based Build-up (SBBU)
tags: [algorithms, optimization]
sources: [goncalves2021new]
---

## Intuition

The Symmetry-based Build-up (SBBU) algorithm is a novel algorithm for solving the DMDGP proposed by Gonçalves et al. (2021). 

While the classic Branch-and-Prune (BP) algorithm uses standard depth-first search (DFS) with backtracking, which can become extremely slow if [pruning edges](pruning-edge.md) occur late or are highly sparse, SBBU takes a different approach by exploiting the mathematical symmetry of the discretizable space.

Instead of searching the tree, SBBU solves the DMDGP as a sequence of **nested subproblems**, each spanned by a single pruning edge.
1. The algorithm first orders the pruning edges in a specific way: sorting in increasing order of the target vertex index, followed by a decreasing order of the source vertex index.
2. It starts with an initial possible realization that satisfies only the [discretization edges](discretization-edge.md).
3. For each pruning edge $\{i, j\}$ in the sequence, the algorithm identifies the set of "necessary symmetry vertices" ($S^{ij}$) that control the relative positions between $i$ and $j$.
4. It then performs an exhaustive search over the local symmetry combination, applying compositions of partial reflections *only* to the candidate position of $x_j$ until the distance constraint $\|x_i - x_j\| = d_{ij}$ is satisfied.
5. Once the correct local reflection choice is found, it updates the realization. Because of the symmetry properties, this update is guaranteed to not violate any of the previous subproblem constraints.

By using this build-up approach, SBBU avoids the deep backtracking typical of DFS in sparse graphs, resulting in substantial speedups on protein-like instances.

## Formal definition

Given a feasible DMDGP $(G, K)$ with $G = (V, E, d)$ and $E = E_D \cup E_P$, where $E_D$ are discretization edges and $E_P$ are pruning edges, SBBU operates as follows:

1. **Pruning Edge Ordering**: Sort the pruning edges $e_1, \dots, e_m \in E_P$ such that $\{u, w\} < \{i, j\}$ if $w < j$ or ($w = j$ and $u > i$).
2. **Subproblem Definition**: Define $P^{ij} := \{e' \in E_P \mid e' < \{i, j\}\}$. The subproblem $(G^{ij}, K)$ includes discretization edges $E_D$ and the subset of pruning edges $P^{ij} \cup \{i, j\}$.
3. **Symmetry Set for Subproblem**: Define the set of necessary symmetry vertices $S^{ij}$ for subproblem $(G^{ij}, K)$ as:
   $$S^{ij} = \left\{v_\ell \in \{v_{i+K+1}, \dots, v_j\} \;\middle|\; \not\exists \{u, w\} \in P^{ij} \text{ such that } u+K < \ell \le w \right\}$$
4. **Iterative Reflection Correction**: Let $x$ be the current possible realization. Search for the local binary decision vector $\bar{s} \in \tilde{B}^{ij} = \{\bar{s} \in \{0, 1\}^{j-i-K} \mid \bar{s}_\ell = 0 \text{ if } v_{i+K+\ell} \notin S^{ij}\}$ such that the reflected position $x_j'$ satisfies the new constraint:
   $$x_j' = R(x_j, \bar{s}) = (R_x^{i+K+1})^{\bar{s}_1} \circ \dots \circ (R_x^j)^{\bar{s}_{j-i-K}}(x_j)$$
   $$\|x_i - x_j'\| = d_{ij}$$
5. **Update**: Update the realization and proceed to the next pruning edge in the sorted sequence.
