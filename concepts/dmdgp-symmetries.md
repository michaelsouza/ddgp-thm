---
term: DMDGP Symmetries
tags: [geometry, group-theory]
sources: [goncalves2021new]
---

## Intuition

A fascinating property of the Discretizable Molecular Distance Geometry Problem (DMDGP) is that its search space contains a high degree of symmetry. Because each candidate position $x_i$ for a vertex $v_i$ is determined by $K$-lateration, it has precisely two options ($x_i^+$ or $x_i^-$) that are reflections of each other across a local hyperplane.

This local reflection structure forms a group of symmetry operations. Crucially, this means that if we find *just one* valid realization that satisfies all distance constraints (both discretization and pruning), we can generate all other valid solutions simply by applying specific combinations of partial reflections. The vertices whose reflections do not violate any pruning edge constraints are called **symmetry vertices** ($S$). 

Knowing this structure allows us to:
1. Predict the exact number of solutions a priori from the graph topology.
2. Speed up search algorithms by focusing only on finding a single initial solution, rather than searching for all of them, since all others are easily derived.

## Formal definition

Let $(G, K)$ be a feasible $^K$DMDGP with a vertex order $v_1, \dots, v_n$. The set of **symmetry vertices** $S \subseteq V$ is defined as:

$$
S := \left\{v_\ell \in V \;\middle|\; \not\exists \{v_i, v_j\} \in E \text{ with } i + K < \ell \le j\right\}
$$

Where:
* $v_{K+1}$ is always in $S$ because the first $K$ vertices define a global symmetry hyperplane.
* If a vertex $v_\ell \in S$ (for $\ell > K + 1$), it means there is no "pruning edge" that spans across the local reflection hyperplane at $\ell$ (i.e., no edge goes from a vertex $v_i$ before $\ell-K$ to a vertex $v_j$ at or after $\ell$).

Under standard non-degeneracy conditions, with probability 1, the total number of incongruent valid realizations of a feasible DMDGP instance is:

$$
|X| = 2^{|S|}
$$

Any solution $x(s') \in X$ can be generated from an initial valid solution $x(s) \in X$ by setting $s'_\ell = s_\ell$ for all $\ell$ where $v_{K+\ell} \notin S$, and freely toggling the binary choice for indices in $S$.
