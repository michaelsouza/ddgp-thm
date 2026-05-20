---
term: K-lateration
tags: [geometry, algorithms]
sources: [goncalves2021new]
---

## Intuition

Lateration is the process of determining the position of an object by measuring its distance from several reference points with known coordinates.
* **1-lateration**: Intersection of 1 sphere in $\mathbb{R}^1$ (2 points).
* **2-lateration (Trilateration in 2D)**: Intersection of 2 circles in $\mathbb{R}^2$ (yields at most 2 points). This is commonly used in land surveying and local positioning systems.
* **3-lateration (Trilateration in 3D)**: Intersection of 3 spheres in $\mathbb{R}^3$ (yields at most 2 points). This is the key mechanism behind GPS (Global Positioning System) receivers.

In $K$-lateration, we generalize this to $K$ dimensions: we find the coordinates of a point $x_i \in \mathbb{R}^K$ by intersecting $K$ hyperspheres centered at $K$ known predecessor positions. When the centers are not in a lower-dimensional affine subspace, the intersection yields exactly two candidate positions, which can be computed analytically or algebraically.

## Formal definition

Given $K$ reference points $x_{i-K}, \dots, x_{i-1} \in \mathbb{R}^K$ whose coordinates are already determined, and exact target distances $d_{\ell, i} > 0$ for each reference point $\ell \in \{i-K, \dots, i-1\}$, $K$-lateration is the problem of finding the candidate realizations $\{x_i^+, x_i^-\} \subset \mathbb{R}^K$ for the $i$-th point that satisfy the system of quadratic equations:

$$
\|x_\ell - x_i\|^2 = d_{\ell, i}^2, \quad \forall \ell \in \{i-K, \dots, i-1\}
$$

Where $\|\cdot\|$ represents the Euclidean norm in $\mathbb{R}^K$. If the reference points do not lie in a lower-dimensional affine subspace (verified by $CM(v_{i-1}, \dots, v_{i-K}) \neq 0$), there are precisely two distinct solutions:

$$
x_i^+ \quad \text{and} \quad x_i^-
$$

Which are reflections of each other across the unique hyperplane defined by the reference points.
