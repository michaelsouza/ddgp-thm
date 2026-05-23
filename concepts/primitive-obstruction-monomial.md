---
term: Primitive Obstruction Monomial
tags: [algebra, geometry]
sources: []
---

## Intuition

A primitive obstruction monomial is the basic algebraic witness that an active
edge changes length under reflection through a mirror clique.

It combines the signed heights of both edge endpoints over the same mirror
hyperplane. If either endpoint lies on the mirror, the monomial vanishes and the
edge is mirror-compatible.

## Formal definition

Let $e=\{a,b\}$ be an edge and let $C$ be a $K$-vertex mirror clique. Using the
[oriented bracket](oriented-bracket.md), define

$$
\mu_{e,C}=[C,a][C,b].
$$

If $a\in C$ or $b\in C$, then $\mu_{e,C}=0$. Otherwise it is generically
nonzero.

In the rank-count proof, primitive obstruction monomials are used to show that
distinct labelled violations $(e,C)$ do not cancel outside a proper algebraic
exceptional set.
