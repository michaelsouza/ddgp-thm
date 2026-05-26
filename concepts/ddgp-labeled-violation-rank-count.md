---
term: DDGP Labeled-Violation Rank Count
tags: [algorithms, geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The DDGP labeled-violation rank count is an SBBU-like counting theorem for
generic local DDGP subproblems.

Instead of enumerating the binary Branch-and-Prune tree, it builds graph-derived reflection masks and computes one rank over $\mathbb{F}_2$.

The construction replaces the DMDGP idea of interval partial reflections by
dependency-graph generators. Base edges and pruning edges are treated uniformly
as active edge constraints. The local solution count is the dimension of the
feasible branch shift space represented by generator combinations with zero net
labeled violation.

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

Let the active edge set be

$$
F=E_0[L_P]\cup P.
$$

Build a family $\mathcal{G}$ of graph generators

$$
g=(m_g,C_g),
$$

where $m_g\in\mathbb{F}_2^{B_P}$ is the mask of a cone generator or base
generator, and $C_g$ is the mirror clique associated with that mask.

For each generator, define its [labeled violation kernel](labeled-violation-kernel.md) signature against the active edge set $F$. This gives a violation vector

$$
\nu_g\in\mathbb{F}_2^{\mathcal{L}},
$$

where labels have the form $(e,C_g)$ with $e\in F$.

Let $M$ be the matrix whose columns are the mask vectors $m_g$, and let $V$ be the matrix whose columns are the labeled violation vectors $\nu_g$. Then the local solution rank is

$$
r_F
=
\operatorname{rank}_{\mathbb{F}_2}
\begin{bmatrix}
M\\
V
\end{bmatrix}
-
\operatorname{rank}_{\mathbb{F}_2}(V).
$$

The generic number of local solutions is

$$
|\Xi_F|=2^{r_F}.
$$

Presentation ambiguities are handled by the [labeled obstruction quotient](labeled-obstruction-quotient.md). The rank formula computes the dimension of the feasible branch shift space without explicitly constructing that quotient.

## Algorithm

1. Compute all discretization dependency sets $\operatorname{Fix}_U(i)$.
2. Build $L_P$ and $B_P$.
3. Build the active edge set $F=E_0[L_P]\cup P$.
4. Build cone generators and base generators on $B_P$.
5. Label every active-edge violation by both the edge and the mirror clique.
6. Compute $r_F=\operatorname{rank}([M;V])-\operatorname{rank}(V)$ over $\mathbb{F}_2$.
7. Return $2^{r_F}$.

This is a graph-combinatorial count. It does not require enumerating local
realizations, although enumeration remains useful for testing implementations on
small instances.

## Status

The Feasible Branch Shift Theorem states that, for $K\ge2$, generic active-edge
symmetries are exactly the feasible branch shifts represented by
zero-net-violation graph-generator presentations. The rank-count formula is the
linear-algebra corollary obtained by computing the dimension of that shift
space.

Small exact-enumeration examples are useful for checking the implementation and
for showing the matrices used by the count. The main known limitation is
nongeneric geometry: accidental mirror coincidences can create extra valid
realizations that are invisible to the purely combinatorial labeled test.

A proof formulation is maintained in the research note [Labeled-Violation Rank
Count: Generic Proof Formulation](../research/labeled-violation-rank-proof.md).
The proof uses the [labeled obstruction quotient](labeled-obstruction-quotient.md)
and a bracket-separation argument on a stretched moment-curve specialization.
