---
term: Mirror-Compatible Graph Mask
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

A mirror-compatible graph mask is a binary mask whose support can be reflected through a predecessor clique without breaking a selected set of distance constraints.

This concept is a graph-theoretic approximation of a local reflection symmetry. The mask identifies which local decision vertices are toggled. The mirror clique identifies the vertices that define the reflection hyperplane. A distance edge crossing the boundary of the mask is preserved when the endpoint outside the mask lies on the mirror clique.

This is the DDGP analogue of the DMDGP fact that a partial reflection preserves all discretization edges but may violate pruning edges that cross the reflection boundary.

## Formal definition

Let $L_e$ be a local DDGP vertex set and let

$$
B_e=L_e\setminus V_0
$$

be its local binary decision set. A graph mask is a pair

$$
(m,C),
$$

where $m\in\mathbb{F}_2^{B_e}$ is a binary mask and $C$ is a predecessor clique acting as the mirror.

The support of $m$ is

$$
\operatorname{supp}(m)=\{z\in B_e \mid m_z=1\}.
$$

Given an edge set $F$, the mask $(m,C)$ is mirror-compatible with $F$ when every edge of $F$ crossing the support boundary has its fixed endpoint in the mirror clique:

$$
\forall \{a,b\}\in F,\quad
|\{a,b\}\cap \operatorname{supp}(m)|=1
\Longrightarrow
\{a,b\}\setminus \operatorname{supp}(m) \subseteq C.
$$

Equivalently:

- edges entirely inside the support are preserved because both endpoints are reflected;
- edges entirely outside the support are unchanged;
- crossing edges are preserved only when the non-reflected endpoint lies on the mirror.

## Use in local counting

For a pruning edge $e=\{u,v\}$, a graph mask can be tested against:

- the local base edge set, to decide whether it preserves the local base solution code;
- the single pruning edge $e$, to decide whether it survives the pruning constraint.

For one pruning edge, this individual mirror-compatibility test is often enough to predict the number of local valid solutions from the rank of the compatible graph-mask span.

This criterion is intended as a generic graph-theoretic test. Nongeneric coordinate coincidences can create additional symmetries, for example when an endpoint not belonging to the mirror clique nevertheless lies geometrically on the mirror hyperplane.

For multiple active edges, and for the current active-edge treatment of base and pruning edges together, individual mirror-compatibility is too restrictive. The corresponding linear object is the [labelled violation kernel](labelled-violation-kernel.md), where violations are labelled by both edge and mirror clique before cancellation.
