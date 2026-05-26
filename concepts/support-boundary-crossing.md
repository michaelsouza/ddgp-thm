---
term: Support Boundary Crossing
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

A support boundary crossing occurs when an active edge has exactly one endpoint
inside a branch mask's support.

Boundary crossings are where a reflection can break an edge length: one endpoint
moves and the other does not. The crossing is harmless only when the fixed
endpoint lies on the reflection hyperplane, which is encoded by membership in
the mirror clique.

## Formal definition

Let $g=(m_g,C_g)$ be a DDGP graph generator with support

$$
S_g=\operatorname{supp}(m_g).
$$

An active edge $e=\{a,b\}$ crosses the support boundary of $g$ when

$$
|e\cap S_g|=1.
$$

The crossing is a labeled violation when the fixed endpoint is not in the
mirror clique:

$$
e\setminus S_g\not\subseteq C_g.
$$

Otherwise the crossing is mirror-compatible, because the fixed endpoint lies on
$H_{C_g}=\operatorname{aff}(C_g)$.
