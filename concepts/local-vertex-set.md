---
term: Local Vertex Set
tags: [geometry, combinatorial-optimization]
sources: []
---

## Intuition

The local vertex set is the set of vertices whose positions are needed to test a
chosen set of pruning edges.

Because DDGP predecessor sets may be nonconsecutive, the vertices needed for an
endpoint are its dependency ancestors rather than every earlier vertex in the
order.

## Formal definition

For a pruning edge set $P$, the local vertex set is

$$
L_P
=
\bigcup_{\{u,v\}\in P}
\left(\operatorname{Fix}_U(u)\cup\operatorname{Fix}_U(v)\right),
$$

where $\operatorname{Fix}_U(i)$ is the
[discretization dependency set](discretization-dependency-set.md), or fixing
set, of vertex $i$.

For one pruning edge $e=\{u,v\}$ this specializes to

$$
L_e=\operatorname{Fix}_U(u)\cup\operatorname{Fix}_U(v).
$$
