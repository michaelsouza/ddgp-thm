---
term: Branch Mask
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

A branch mask records which local Branch-and-Prune bits are toggled by a
candidate symmetry or generator.

The mask is a binary object. It says which branch decisions change, but it does
not by itself record the mirror hyperplane that gives the operation its
geometric meaning. The rank-count theory therefore pairs each mask with a
[mirror clique](mirror-clique.md).

## Formal definition

Let $B_P$ be the [local binary decision set](local-binary-decision-set.md). A
branch mask is a vector

$$
m\in\mathbb F_2^{B_P}.
$$

Its support is

$$
\operatorname{supp}(m)=\{z\in B_P\mid m_z=1\}.
$$

Applying the mask to a BP sign string $s$ gives

$$
s\mapsto s\oplus m.
$$

In the rank-count construction, branch masks appear as columns of the
[generator mask matrix](generator-mask-matrix.md).
