---
term: Oriented Bracket
tags: [algebra, geometry, linear-algebra]
sources: []
---

## Intuition

An oriented bracket is a determinant that measures the signed height of a vertex
relative to a mirror clique.

The rank-count proof uses brackets to separate labeled obstructions
algebraically. If two labeled violations have different primitive bracket
monomials, they cannot cancel generically.

## Formal definition

Let $C=\{c_1,\dots,c_K\}$ be a $K$-clique in $\mathbb R^K$, and let

$$
\widehat{x}_i=(1,x_i)^\top
$$

be the homogeneous coordinate vector. For a vertex $z$, the oriented bracket is

$$
[C,z]
=
\det
\begin{bmatrix}
\widehat{x}_{c_1}&\widehat{x}_{c_2}&\cdots&\widehat{x}_{c_K}&\widehat{x}_z
\end{bmatrix}.
$$

It is proportional to the signed height of $z$ over
$H_C=\operatorname{aff}(C)$. In particular, $[C,z]=0$ exactly when $z$ lies in
the affine span of $C$.
