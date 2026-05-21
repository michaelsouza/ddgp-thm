---
term: DDGP Labelled-Violation Rank Count
tags: [algorithms, geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The DDGP labelled-violation rank count is a candidate SBBU-like counting rule for local DDGP subproblems.

Instead of enumerating the binary Branch-and-Prune tree, it builds graph-derived reflection masks and computes one rank over $\mathbb{F}_2$.

The construction replaces the DMDGP idea of interval partial reflections by dependency-graph masks. A mask is allowed locally when it preserves the induced base edges. Pruning edges then impose labelled linear violations. The local solution count is predicted by the dimension of the mask combinations whose labelled violations cancel.

## Formal definition

Let a DDGP instance in $\mathbb{R}^K$ have fixed initial simplex

$$
V_0=\{v_1,\dots,v_{K+1}\},
$$

predecessor sets $U_i$, base edge set $E_0$, and pruning edge set $P$.

Define the local vertex set

$$
L_P
=\bigcup_{\{u,v\}\in P}
\left(\operatorname{Fix}_U(u)\cup\operatorname{Fix}_U(v)\right),
$$

and the local decision set

$$
B_P=L_P\setminus V_0.
$$

Build a family $\mathcal{G}$ of graph generators

$$
g=(m_g,C_g),
$$

where $m_g\in\mathbb{F}_2^{B_P}$ is a dependency-cone or base-clique closure mask, and $C_g$ is the mirror clique associated with that mask.

Keep only generators that are [mirror-compatible](mirror-compatible-graph-mask.md) with the induced local base edge set

$$
E_0[L_P].
$$

For each remaining generator, define its [labelled violation kernel](labelled-violation-kernel.md) signature against the pruning edge set $P$. This gives a violation vector

$$
\nu_g\in\mathbb{F}_2^{\mathcal{L}},
$$

where labels have the form $(e,C_g)$ with $e\in P$.

Let $M$ be the matrix whose columns are the mask vectors $m_g$, and let $V$ be the matrix whose columns are the labelled violation vectors $\nu_g$. Then the predicted local valid rank is

$$
r_P
=
\operatorname{rank}_{\mathbb{F}_2}
\begin{bmatrix}
M\\
V
\end{bmatrix}
-
\operatorname{rank}_{\mathbb{F}_2}(V).
$$

The predicted number of local solutions is

$$
|\Xi_P|=2^{r_P}.
$$

## Algorithm

1. Compute all discretization dependency sets $\operatorname{Fix}_U(i)$.
2. Build $L_P$ and $B_P$.
3. Build dependency-cone masks and base-clique closure masks on $B_P$.
4. Keep only masks that are mirror-compatible with $E_0[L_P]$.
5. Label every pruning-edge violation by both the pruning edge and the mirror clique.
6. Compute $r_P=\operatorname{rank}([M;V])-\operatorname{rank}(V)$ over $\mathbb{F}_2$.
7. Return $2^{r_P}$.

This is a graph-combinatorial prediction. It does not require enumerating local realizations, although enumeration remains useful for validating the conjecture on small instances.

## Status

This is currently a conjectural counting rule, not a proved theorem.

Experiments support it for generic DDGP and DMDGP instances in dimensions $K=2$ and $K=3$, with one or multiple pruning edges. The main known limitation is nongeneric geometry: accidental mirror coincidences can create extra valid realizations that are invisible to the purely combinatorial mirror-compatibility test.
