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

## Feasible branch shifts

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

This is the local analogue of the DMDGP symmetry count.

In the generic rank-count theory, this stabilizer is the feasible branch shift
space generated from graph-derived generators, especially [cone
generators](cone-generator.md) and [base generators](base-generator.md), through
the labeled violation kernel.

A graph-derived mask is an individual symmetry candidate only when it is [mirror-compatible](mirror-compatible-graph-mask.md) with the relevant edge constraints.

For multiple active edge constraints, even individually non-compatible
generators can be combined through a [labeled violation
kernel](labeled-violation-kernel.md). When this kernel satisfies the
[labeled presentation property](labeled-presentation-property.md), the local
solution code is exactly one affine coset of the feasible branch shift space.
The resulting rank formula is recorded as the [DDGP labeled-violation rank
count](ddgp-labeled-violation-rank-count.md).
