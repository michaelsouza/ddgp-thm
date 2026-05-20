---
term: Distance Geometry Problem (DGP)
tags: [geometry, mathematics]
sources: [goncalves2021new]
---

## Intuition

The Distance Geometry Problem (DGP) is a fundamental inverse problem: given a set of points and a subset of the pairwise Euclidean distances between them, the goal is to reconstruct the coordinates of the points in a target $K$-dimensional space.

Imagine you have a scaffolding or a protein structure, and you can only measure the distances between certain pairs of joints or atoms. The DGP asks if you can find a valid set of coordinates for all the joints/atoms that matches all of your measurements. If no such coordinates exist, the problem is infeasible. If multiple incongruent configurations satisfy all the distances, the solution is not unique.

This problem is strongly NP-hard in general, but arises naturally in many applications:
1. **Structural Biology**: Determining 3D molecular/protein conformations using distance constraints obtained from Nuclear Magnetic Resonance (NMR) experiments (usually $K = 3$).
2. **Sensor Network Localization**: Finding the positions of sensors in a network given partial distance measurements between neighboring sensors (usually $K = 2$ or $K = 3$).
3. **Data Science**: Dimensionality reduction and manifold learning (e.g., Multidimensional Scaling).

## Formal definition

Given a simple, undirected, weighted graph $G = (V, E, d)$, where the weight function $d : E \to (0, \infty)$ represents the known exact distances, and a positive integer $K$, the Distance Geometry Problem consists in finding a realization (or embedding) $x : V \to \mathbb{R}^K$ such that:

$$
\forall \{u, v\} \in E, \quad \|x_u - x_v\| = d_{uv}
$$

Where:
* $x_u = x(u) \in \mathbb{R}^K$ is the position assigned to vertex $u \in V$.
* $\|\cdot\|$ denotes the standard Euclidean norm in $\mathbb{R}^K$.
* $d_{uv} = d(\{u, v\})$ is the exact distance constraint between $u$ and $v$.
