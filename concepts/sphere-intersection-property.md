---
term: Sphere Intersection Property (SIP)
tags: [geometry, mathematics]
sources: [liberti2011number]
---

## Intuition

The Sphere Intersection Property (SIP) is the geometric foundation that enables the discretization of the Distance Geometry Problem (DGP) in $\mathbb{R}^K$. In the context of sequentially embedding a graph, the position of each new vertex $v$ is determined by its distances to a set of $K$ already embedded predecessors. Geometrically, finding this position corresponds to finding the intersection of $K$ spheres (or hyperspheres in higher dimensions) centered at those predecessors.

The SIP states that, under non-degenerate conditions (such as the centers of the spheres not lying in a lower-dimensional subspace), the intersection of $K$ hyperspheres in $\mathbb{R}^K$ generally consists of either $0$ or $2$ distinct points. Here, "generally" means that the set of configurations for which the intersection is a single point (tangential intersection) has Lebesgue measure zero in the space of all possible configurations.

This property is what allows search algorithms like [Branch-and-Prune (BP)](branch-and-prune.md) to formulate the search space as a discrete binary tree rather than a continuous space, converting a continuous optimization problem into a combinatorial search.

## Formal definition

Let $U_v = \{u_1, \dots, u_K\} \subset \mathbb{R}^K$ be the positions of $K$ already embedded vertices, and let $d_{u_i, v} > 0$ be the known distances from these vertices to a new vertex $v$. The position $z \in \mathbb{R}^K$ of $v$ must lie at the intersection of the $K$ hyperspheres:

$$
S_i = \{ z \in \mathbb{R}^K \mid \|z - x_{u_i}\| = d_{u_i, v} \}, \quad \forall i \in \{1, \dots, K\}
$$

According to the Sphere Intersection Property, if the points in $U_v$ are affinely independent (which is guaranteed if they satisfy the [Strict Simplex Inequalities](strict-simplex-inequalities.md)), then:

$$
\left| \bigcap_{i=1}^{K} S_i \right| \le 2
$$

Furthermore, with probability 1 (Lebesgue measure 1), if the intersection is non-empty, then:

$$
\left| \bigcap_{i=1}^{K} S_i \right| = 2
$$

This guarantees that the sequential placement of each vertex admits exactly two candidate coordinates $\{z', z''\} \subset \mathbb{R}^K$, which are reflections of each other across the hyperplane containing $U_v$.
