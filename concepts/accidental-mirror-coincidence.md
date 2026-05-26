---
term: Accidental Mirror Coincidence
tags: [geometry, algebra]
sources: []
---

## Intuition

An accidental mirror coincidence is a nongeneric geometric alignment that makes
a labeled violation harmless at one special coordinate realization.

The rank-count formula predicts generic behavior from graph data. If an endpoint
that should be off a reflection hyperplane happens to lie on it, a distance can
be preserved accidentally, producing more local solutions than the generic
rank-count prediction.

## Formal definition

Let $g=(m_g,C_g)$ be a graph generator with support $S_g$, and let
$e=\{a,b\}$ be an active edge with

$$
|e\cap S_g|=1.
$$

Generically, this crossing is a violation when the fixed endpoint is not in
$C_g$. An accidental mirror coincidence occurs at a specific realization when
that fixed endpoint is nevertheless geometrically contained in

$$
H_{C_g}=\operatorname{aff}(C_g).
$$

Equivalently, the relevant signed-height polynomial vanishes at that point even
though it is not identically zero. Such coincidences lie in a proper
[algebraic exceptional set](algebraic-exceptional-set.md).
