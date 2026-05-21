---
term: Endpoint-Separating Cone
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

An endpoint-separating cone is a dependency cone that contains exactly one endpoint of a pruning edge.

In the DMDGP, a pruning edge removes partial reflections whose reflection level lies between the two endpoints in the discretization order. In the DDGP, the discretization order is not enough because dependencies may be nonconsecutive. A dependency cone is the natural replacement for a linear interval.

## Definition

Let $D_U$ be the discretization dependency DAG of a DDGP instance. For a non-seed vertex $w$, define its forward dependency cone as

$$
\operatorname{Cone}_U(w)
=
\{z\in V \mid w \leadsto z \text{ in } D_U\}\cup\{w\}.
$$

For a pruning edge $e=\{u,v\}$, the cone of $w$ is endpoint-separating for $e$ when

$$
|\operatorname{Cone}_U(w)\cap\{u,v\}|=1.
$$

Such a cone is important because a reflection or binary toggle supported on that cone can move one endpoint of the pruning edge without moving the other.

## Hypothesis role

Endpoint-separating cones are candidate DDGP analogues of DMDGP reflection levels crossed by pruning edges.

They are not, by themselves, guaranteed to determine the number of local solutions. Experiments show that cone-based predictions can be strong approximations, but some DDGP instances have valid binary symmetries that are not simple dependency-cone toggles.
