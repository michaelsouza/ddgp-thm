---
term: Cayley-Menger Determinant
tags: [geometry, linear-algebra]
sources: [goncalves2021new, liberti2011number]
---

## Intuition

The Cayley-Menger determinant is a highly useful algebraic tool in distance geometry. It allows one to calculate the volume of a simplex (such as a triangle, tetrahedron, etc.) in a given dimension using *only* the pairwise distances between its vertices, without knowing their coordinates.

In the context of the DMDGP, we want to ensure that a set of $K$ points in $\mathbb{R}^K$ does not lie in a lower-dimensional subspace (e.g., three points in $\mathbb{R}^3$ are not collinear). If they were collinear, the sphere intersection for the next point would result in a whole circle of infinite possible points instead of exactly two discrete points, which would break the discretization of the search space. By checking that the Cayley-Menger determinant of the $K$ points is non-zero, we guarantee that their volume in $K-1$ dimensions is non-zero, meaning they span a proper hyperplane.

## Formal definition

Let $v_1, \dots, v_m$ be a set of $m$ vertices, and let $d_{ij}$ denote the Euclidean distance between $v_i$ and $v_j$. The Cayley-Menger determinant $CM(v_1, \dots, v_m)$ is the determinant of the $(m+1) \times (m+1)$ matrix:

$$
CM(v_1, \dots, v_m) = \begin{vmatrix}
0 & 1 & 1 & 1 & \dots & 1 \\
1 & 0 & d_{12}^2 & d_{13}^2 & \dots & d_{1m}^2 \\
1 & d_{21}^2 & 0 & d_{23}^2 & \dots & d_{2m}^2 \\
1 & d_{31}^2 & d_{32}^2 & 0 & \dots & d_{3m}^2 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
1 & d_{m1}^2 & d_{m2}^2 & d_{m3}^2 & \dots & 0
\end{vmatrix}
$$

The volume $V$ of the $(m-1)$-dimensional simplex spanned by the realizations $x_1, \dots, x_m \in \mathbb{R}^K$ of these vertices satisfies:

$$
V^2 = \frac{(-1)^m}{2^{m-1} ((m-1)!)^2} CM(v_1, \dots, v_m)
$$

For a set of $K$ points in a $^K$DMDGP, we require $CM(v_{i-1}, \dots, v_{i-K}) \neq 0$, which ensures that their $(K-1)$-dimensional volume is strictly positive. This requirement is formally stated as the [Strict Simplex Inequalities](strict-simplex-inequalities.md).

