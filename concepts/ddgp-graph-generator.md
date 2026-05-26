---
term: DDGP Graph Generator
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

A DDGP graph generator is a graph-derived candidate operation on the local BP
sign cover.

It has two pieces of information: a binary branch mask describing which local
decisions are toggled, and a mirror clique describing the reflection hyperplane
used to interpret the operation geometrically. Keeping both pieces is essential
because the same mask can have different geometric presentations.

## Formal definition

For a local DDGP subproblem with decision set $B_P$, a graph generator is a pair

$$
g=(m_g,C_g),
$$

where

$$
m_g\in\mathbb F_2^{B_P}
$$

is the [branch mask](branch-mask.md), and $C_g$ is the
[mirror clique](mirror-clique.md).

The rank-count theory uses two generator families:

- [cone generators](cone-generator.md);
- [base generators](base-generator.md).

The generator family is denoted by $\mathcal G$.
