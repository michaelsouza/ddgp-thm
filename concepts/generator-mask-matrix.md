---
term: Generator Mask Matrix
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The generator mask matrix collects the branch masks of all graph generators.

It records only the binary effect of generator combinations on the BP sign
cover. The geometric obstruction data is stored separately in the
[labelled violation matrix](labelled-violation-matrix.md).

## Formal definition

Let $\mathcal G=\{g_1,\dots,g_m\}$ be a family of DDGP graph generators with
branch masks $m_{g_j}\in\mathbb F_2^{B_P}$. The generator mask matrix is

$$
M=
\begin{bmatrix}
| & | & & |\\
m_{g_1} & m_{g_2} & \cdots & m_{g_m}\\
| & | & & |
\end{bmatrix}.
$$

It represents the linear mask map

$$
M:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{B_P},
$$

where

$$
M\alpha
=
\bigoplus_{g\in\operatorname{supp}(\alpha)}m_g.
$$

This matrix is the top block in the rank-count matrix
$\begin{bmatrix}M\\V_F\end{bmatrix}$.
