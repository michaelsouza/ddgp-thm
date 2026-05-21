---
term: Local Solution Code
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The local solution code is the set of binary branch strings that realize a local DDGP subproblem.

It is useful when studying the number of local solutions induced by a [pruning edge](pruning-edge.md). Instead of counting coordinate realizations directly, we encode each local realization by the binary choices made by the [Branch-and-Prune](branch-and-prune.md) construction after the initial rigid simplex has been fixed.

This makes it possible to test whether local DDGP solution sets have algebraic structure similar to the symmetry structure of the [DMDGP](dmdgp-symmetries.md).

## Definition

Let $G=(V,E,d)$ be a DDGP instance in $\mathbb{R}^K$ with predecessor sets $U_i$ and fixed initial simplex

$$
V_0=\{v_1,\dots,v_{K+1}\}.
$$

For a pruning edge $e=\{u,v\}$, define the local vertex set

$$
L_e=\operatorname{Fix}_U(u)\cup\operatorname{Fix}_U(v).
$$

The local decision set is

$$
B_e=L_e\setminus V_0.
$$

Each local realization can be encoded by a binary vector

$$
s\in \mathbb{F}_2^{B_e},
$$

where each coordinate records which of the two possible positions was selected when placing the corresponding non-seed vertex.

Given a set of active local edges $F\subseteq E$, the local solution code is

$$
\Xi(L_e,F)
=
\{s\in\mathbb{F}_2^{B_e}
\mid
s \text{ realizes all distance constraints in } F\}.
$$

For the local base subproblem, $F$ contains the induced base edges. For the pruned local subproblem, $F$ also contains the pruning edge $e$.

## Stabilizer

The XOR stabilizer of a local solution code is

$$
\operatorname{Stab}(\Xi)
=
\{g\in\mathbb{F}_2^{B_e}
\mid
\Xi\oplus g=\Xi\}.
$$

If $\Xi$ is an affine subspace of $\mathbb{F}_2^{B_e}$, then $\Xi$ is a single coset of its stabilizer and

$$
|\Xi| = |\operatorname{Stab}(\Xi)| = 2^{\dim \operatorname{Stab}(\Xi)}.
$$

This is a candidate local analogue of the DMDGP symmetry count.

Current experiments test whether this stabilizer can be generated from graph-derived masks, especially dependency-cone masks and [base-clique closure masks](base-clique-closure-mask.md).

A graph-derived mask is an individual symmetry candidate only when it is [mirror-compatible](mirror-compatible-graph-mask.md) with the relevant edge constraints.

For multiple active edge constraints, even individually non-compatible generators can be combined through a [labelled violation kernel](labelled-violation-kernel.md). When this kernel satisfies the [labelled presentation property](labelled-presentation-property.md), the local solution code is exactly one affine coset of the projected kernel. The resulting rank formula is recorded as the [DDGP labelled-violation rank count](ddgp-labelled-violation-rank-count.md).
