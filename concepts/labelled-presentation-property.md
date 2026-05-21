---
term: Labelled Presentation Property
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The labelled presentation property is the exact geometric statement needed to
turn the DDGP labelled-violation rank count into a theorem.

The rank formula builds graph generators, records which active edges they
violate, and labels each violation by the mirror clique that caused it. A branch
difference is predicted to preserve all active edges when it can be written as
a combination of graph generators whose labelled violations cancel.

The labelled presentation property says that this prediction is exact:
graph-visible labelled cancellations are precisely the generic geometric
cancellations of active-edge distances.

## Formal definition

Let $B$ be the local binary decision set and let $\mathcal G$ be a family of
graph generators

$$
g=(m_g,C_g),
$$

where $m_g\in\mathbb F_2^B$ is the branch mask and $C_g$ is the mirror clique.

Let $F$ be an active edge set. The labelled violation map is

$$
V_F:\mathbb F_2^{\mathcal G}\to \mathbb F_2^{\mathcal L_F},
$$

and the mask map is

$$
M:\mathbb F_2^{\mathcal G}\to \mathbb F_2^B.
$$

The lifted labelled span is

$$
Z_F
=
\operatorname{span}_{\mathbb F_2}
\{(m_g,\nu_F(g))\mid g\in\mathcal G\}
\subseteq
\mathbb F_2^B\oplus\mathbb F_2^{\mathcal L_F}.
$$

The zero-label projected kernel is

$$
\mathcal K_F
=
\{h\in\mathbb F_2^B\mid (h,0)\in Z_F\}
=
M(\ker V_F).
$$

Equivalently, let $Q_F$ be the [labelled obstruction
quotient](labelled-obstruction-quotient.md), and let

$$
\omega_F:\mathbb F_2^B\to Q_F
$$

be the canonical obstruction map. Then

$$
\mathcal K_F=\ker\omega_F.
$$

The graph-generator family has the labelled presentation property for $F$ when,
for every feasible local solution $s_0\in\Xi_F$,

$$
\Xi_F=s_0\oplus\mathcal K_F.
$$

Equivalently, a branch difference $h$ preserves the active local edge set
generically if and only if it has at least one graph-generator presentation

$$
h=M\alpha
$$

with zero total labelled violation

$$
V_F\alpha=0.
$$

## Role in the DDGP count

If the labelled presentation property holds, then the local solution code is an
affine subspace and

$$
|\Xi_F|=|\mathcal K_F|.
$$

The projected kernel rank is

$$
\dim \mathcal K_F
=
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F).
$$

Therefore the local count is

$$
|\Xi_F|
=
2^{
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F)
}.
$$

The remaining geometric burden is to prove generic edge-sign separation: outside
a proper algebraic exceptional set, active-edge distance changes with different
labels cannot cancel unless the labels cancel in $\mathbb F_2^{\mathcal L_F}$.
