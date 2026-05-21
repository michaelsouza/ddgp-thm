---
term: DMDGP Symmetries
tags: [geometry, group-theory]
sources: [goncalves2021new]
---

## Intuition

A fascinating property of the Discretizable Molecular Distance Geometry Problem (DMDGP) is that its search space contains a high degree of symmetry. Because each candidate position $x_i$ for a vertex $v_i$ is determined by $K$-lateration, it has precisely two options ($x_i^+$ or $x_i^-$) that are reflections of each other across a local hyperplane. In these notes, the first $K+1$ vertices are fixed as an initial rigid simplex, so the trivial global reflection of this simplex is not counted as a solution-generating decision.

This local reflection structure forms a group of symmetry operations. Crucially, this means that if we find *just one* valid realization that satisfies all distance constraints (both [discretization edges](discretization-edge.md) and [pruning edges](pruning-edge.md)), we can generate all other valid solutions simply by applying specific combinations of partial reflections. The vertices whose reflections do not violate any pruning edge constraints are called **symmetry vertices** ($S$). 

Knowing this structure allows us to:
1. Predict the exact number of solutions a priori from the graph topology.
2. Speed up search algorithms by focusing only on finding a single initial solution, rather than searching for all of them, since all others are easily derived.

## Formal definition

Let $(G, K)$ be a feasible $^K$DMDGP with a vertex order $v_1, \dots, v_n$. With the first $K+1$ vertices fixed, the set of nontrivial **symmetry vertices** $S \subseteq \{v_{K+2}, \dots, v_n\}$ is defined as:

$$
S := \left\{v_\ell \in \{v_{K+2}, \dots, v_n\} \;\middle|\; \not\exists \{v_i, v_j\} \in E \text{ with } i + K < \ell \le j\right\}
$$

Where:
* The reflection at $v_{K+1}$ is the global reflection of the initial simplex and is removed by fixing the first $K+1$ vertices.
* If a vertex $v_\ell \in S$, it means there is no [pruning edge](pruning-edge.md) that spans across the local reflection hyperplane at $\ell$ (i.e., no edge goes from a vertex $v_i$ before $\ell-K$ to a vertex $v_j$ at or after $\ell$).

Under standard non-degeneracy conditions, with probability 1, the total number of incongruent valid realizations of a feasible DMDGP instance is:

$$
|X| = 2^{|S|}
$$

Any solution $x(s') \in X$ can be generated from an initial valid solution $x(s) \in X$ by setting $s'_\ell = s_\ell$ for all non-seed vertices where the corresponding $v_\ell \notin S$, and freely toggling the binary choice for indices in $S$.
