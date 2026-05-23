---
term: Local Active Subproblem
tags: [geometry, combinatorial-optimization]
sources: []
---

## Intuition

A local active subproblem is the part of a DDGP instance needed to evaluate a
chosen set of pruning constraints, together with the active edges imposed on
that part.

The rank-count theory works locally because a pruning edge only depends on the
vertices required to construct its endpoints. For multiple pruning edges, the
local subproblem is the union of those dependency regions.

## Formal definition

Given a pruning edge set $P$, define the [local vertex set](local-vertex-set.md)

$$
L_P=
\bigcup_{\{u,v\}\in P}
\left(\operatorname{Fix}_U(u)\cup\operatorname{Fix}_U(v)\right),
$$

and the [local binary decision set](local-binary-decision-set.md)

$$
B_P=L_P\setminus V_0.
$$

A local active subproblem is the pair

$$
(L_P,F),
$$

where $F\subseteq E[L_P]$ is the [active edge set](active-edge-set.md), usually

$$
F=E_0[L_P]\cup P.
$$

The associated solution set is the [local solution code](local-solution-code.md)
$\Xi_F\subseteq\mathbb F_2^{B_P}$.
