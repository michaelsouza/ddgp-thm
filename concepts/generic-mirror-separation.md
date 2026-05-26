---
term: Generic Mirror Separation
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

Generic mirror separation is the genericity condition that different
mirror-induced edge violations do not cancel accidentally.

For one graph mask, the geometric rule is simple: reflecting one endpoint of an edge through a mirror hyperplane preserves the edge length only when the other endpoint lies on that hyperplane. In DDGP language, this forced case is recorded by [mirror-compatible graph masks](mirror-compatible-graph-mask.md).

With several graph masks, cancellations can occur, but generic cancellations
are labeled: a violation of active edge $e$ through mirror clique $C$ can
cancel only with another violation carrying the same pair $(e,C)$.

In the current proof formulation, this is the geometric content behind the
[labeled presentation property](labeled-presentation-property.md).

## Formal definition

Let $\mathcal{G}$ be a set of graph generators

$$
g=(m_g,C_g),
$$

where $m_g$ is a binary support mask and $C_g$ is the mirror clique.

Let $F$ be an active edge set and let

$$
\nu(g)\in\mathbb{F}_2^{\mathcal{L}}
$$

be the labeled violation vector used in the [labeled violation kernel](labeled-violation-kernel.md), with labels of the form $(e,C_g)$.

A local active subproblem satisfies generic mirror separation if every branch
shift with nonzero labeled obstruction changes at least one active-edge length
outside a proper algebraic exceptional set. Equivalently, the labeled
graph-generator presentation is exact: relative to any reference solution
$s^\ast\in\Xi_F$,

$$
\Xi_F=s^\ast\oplus M(\ker\nu).
$$

A necessary geometric part of this property is that active-edge preservation is equivalent to the existence of a labeled-zero presentation of the same branch mask:

$$
\exists \alpha\in\mathbb F_2^{\mathcal G}
\quad
M\alpha=s\oplus t
\quad\text{and}\quad
\sum_g \alpha_g\nu(g)=0.
$$

Equivalently, the only generic edge-violation cancellations are the labeled
cancellations captured by the kernel of $\nu$.

This excludes accidental cases where a vertex outside a mirror clique nevertheless lies on the mirror hyperplane, or where reflections through different mirror hyperplanes preserve a pruning distance for special coordinate reasons.

## Proof role

Generic mirror separation is not an additional algorithmic step. It is the
genericity statement needed in the proof of the [DDGP labeled-violation rank
count](ddgp-labeled-violation-rank-count.md).

The algebraic part is finite: for every branch difference with nonzero
canonical obstruction, the accidental preservation of all active edges is
contained in a proper algebraic subset. The geometric part is edge-sign
separation: different labels $(e,C)$ must correspond to generically distinct
edge-sign characters after quotienting out pure presentation relations.

The current proof formulation establishes this for $K\ge2$ using the
[labeled obstruction quotient](labeled-obstruction-quotient.md) and a
bracket-separation argument. Small exact-enumeration examples are useful for
checking implementations and illustrating the rank count.
