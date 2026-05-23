---
term: Rank-Nullity Projected Kernel
tags: [linear-algebra, geometry]
sources: []
---

## Intuition

The rank-nullity projected kernel is the linear-algebra identity behind the
DDGP labelled-violation count.

It computes the dimension of zero-obstruction branch masks without explicitly
constructing the labelled obstruction quotient. This is why the final counting
algorithm only needs ranks of binary matrices.

## Formal definition

Let

$$
T:\mathbb F_2^{\mathcal G}\to
\mathbb F_2^{B_P}\oplus\mathbb F_2^{\mathcal L_F},
\qquad
T(\alpha)=(M\alpha,V_F\alpha).
$$

Let $\pi$ be the projection from $\operatorname{im}T$ onto the violation
coordinates. Then

$$
\operatorname{im}\pi=\operatorname{im}V_F
$$

and

$$
\ker\pi
=
\{(M\alpha,0)\mid V_F\alpha=0\}
=
\mathcal K_F\oplus\{0\}.
$$

Rank-nullity applied to $\pi$ gives

$$
\dim\mathcal K_F
=
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F).
$$
