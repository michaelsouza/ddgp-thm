---
term: Base-Admissible Block
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

A base-admissible block is a set of local branch decisions that can be reflected
as one block through a mirror clique without breaking the local base graph.

This concept is used in the rank-count proof to turn a zero-labelled generator
presentation into an actual Branch-and-Prune sign transformation. The base
edges are important: they ensure that the block is a valid operation on the BP
sign cover, not just an abstract reflection of an arbitrary vertex subset.

## Formal definition

Let $C$ be a [mirror clique](mirror-clique.md), let $S\subseteq B_P$ be a set of
local decision vertices, and let $E_0[L_P]$ be the induced local base edge set.
The set $S$ is $C$-admissible when:

1. every vertex in $S$ is downstream of $C$ in the DDGP order;
2. every local base edge crossing the boundary of $S$ has its fixed endpoint in
   $C$:

$$
|e\cap S|=1
\quad\Longrightarrow\quad
e\setminus S\subseteq C,
\qquad e\in E_0[L_P].
$$

If $S$ is $C$-admissible, reflecting exactly the vertices in $S$ through
$\operatorname{aff}(C)$ realizes the BP bit flip $s\mapsto s\oplus\chi_S$.
