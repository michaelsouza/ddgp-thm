---
term: Initial Rigid Simplex
tags: [geometry, combinatorial-optimization]
sources: []
---

## Intuition

The initial rigid simplex is the fixed coordinate seed used to remove global
Euclidean motions from a DDGP instance.

After this seed is fixed, the remaining vertices are placed by lateration and
carry binary Branch-and-Prune choices. In the counting convention used by the
rank-count report, branch decisions begin at vertex $v_{K+2}$.

## Formal definition

For a DDGP instance in $\mathbb R^K$ with order $v_1,\dots,v_n$, the initial
rigid simplex is

$$
V_0=\{v_1,\dots,v_{K+1}\}.
$$

The vertices in $V_0$ are assumed to be affinely independent and fixed in
coordinates. They do not contribute branch bits. The local binary decision set
therefore removes them:

$$
B_P=L_P\setminus V_0.
$$
