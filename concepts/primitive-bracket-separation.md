---
term: Primitive Bracket Separation
tags: [algebra, geometry]
sources: []
---

## Intuition

Primitive bracket separation is the generic non-cancellation principle used in
the rank-count proof.

It says that different labeled obstructions have distinct leading algebraic
signatures after a suitable specialization. Therefore a nonzero labeled
obstruction class cannot vanish identically through accidental cancellation.

## Formal definition

For a labeled obstruction $(e=\{a,b\},C)$, define the primitive obstruction
monomial

$$
\mu_{e,C}=[C,a][C,b].
$$

For $K\ge2$, specialize generic coordinates to a stretched moment curve

$$
x_i(t)=(t^{w_i},t^{2w_i},\dots,t^{Kw_i}),
\qquad w_i=N^i,
$$

with $N$ sufficiently large. Under this specialization, distinct labels
$(e,C)$ have distinct leading exponents. Hence any nonempty sum of primitive
obstruction monomials has a unique leading term and is not the zero polynomial.

This separates nonzero classes in the [labeled obstruction
quotient](labeled-obstruction-quotient.md).
