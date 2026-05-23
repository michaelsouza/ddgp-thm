---
term: Active Edge Set
tags: [geometry, combinatorial-optimization]
sources: []
---

## Intuition

The active edge set is the set of distance constraints currently being enforced
in a local DDGP subproblem.

In the rank-count theory, base edges and pruning edges are treated uniformly as
active constraints. This makes it possible to ask which branch-mask operations
preserve all currently relevant distances, regardless of whether those distances
come from the discretization construction or from additional pruning.

## Formal definition

Let $L_P$ be the [local vertex set](local-vertex-set.md) associated with a
pruning edge set $P$, and let $E_0$ be the base edge set induced by the DDGP
predecessor sets. The active edge set is usually

$$
F=E_0[L_P]\cup P.
$$

More generally, $F\subseteq E[L_P]$ is any chosen set of local distance
constraints whose preservation is being tested. The generic rank-count theorem
assumes

$$
E_0[L_P]\subseteq F,
$$

so that every local base edge needed by the Branch-and-Prune sign cover remains
enforced.
