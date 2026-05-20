---
term: Strict Simplex Inequalities
tags: [geometry, mathematics]
sources: [liberti2011number]
---

## Intuition

Strict simplex inequalities are geometric conditions that ensure a set of points forms a non-degenerate simplex in a given dimension. For example:
* In 2D, 3 points form a non-degenerate triangle if they satisfy the strict triangle inequalities (which means they are not collinear and have a strictly positive area).
* In 3D, 4 points form a non-degenerate tetrahedron if they satisfy the strict tetrahedral inequalities (meaning they are not coplanar and have a strictly positive volume).

In distance geometry, especially for problems like the [DMDGP](discretizable-molecular-distance-geometry-problem.md) and [DDGP](discretizable-distance-geometry-problem.md), ensuring that the $K$ predecessors used for embedding the next vertex satisfy these inequalities is absolutely critical. According to the [Sphere Intersection Property (SIP)](sphere-intersection-property.md), the next vertex is placed at the intersection of $K$ spheres. If the $K$ centers of these spheres were degenerate (i.e. lying in a lower-dimensional affine subspace), their intersection would yield either an infinite number of points or a single tangent point with probability zero. By enforcing strict simplex inequalities, we guarantee that the centers span a proper hyperplane, which ensures that the intersection yields exactly two distinct candidate positions with probability 1.

## Formal definition

Let $U = \{x_1, \dots, x_{K+1}\}$ be a set of $K+1$ points in $\mathbb{R}^K$. The volume $\Delta_K(U)$ of the $K$-simplex on $U$ is given by the Cayley-Menger formula:

$$
\Delta_K(U) = \sqrt{\frac{(-1)^{K+1}}{2^K (K!)^2} |D'|}
$$

Where:
* $D$ is the $(K+1) \times (K+1)$ symmetric distance matrix whose $(i,j)$-th component is $\|x_i - x_j\|^2$.
* $D'$ is the $(K+2) \times (K+2)$ matrix obtained by bordering $D$ with a left column $(0, 1, \dots, 1)^T$ and a top row $(0, 1, \dots, 1)$ (see [Cayley-Menger Determinant](cayley-menger-determinant.md) for details).

The **Strict Simplex Inequalities** require that the volume of the simplex is strictly positive:

$$
\Delta_K(U) > 0
$$

Which is equivalent to requiring that the Cayley-Menger determinant $|D'| \neq 0$ and has the appropriate sign $(-1)^{K+1}$.

In a $^K\text{DMDGP}$ or $^K\text{DDGP}$, for each vertex $v$ with rank $> K$, the $K$ predecessors $U_v$ must satisfy:

$$
\Delta_{K-1}(U_v) > 0
$$

Which guarantees that the $K$ predecessor points span a $(K-1)$-dimensional hyperplane in $\mathbb{R}^K$.
