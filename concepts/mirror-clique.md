---
term: Mirror Clique
tags: [geometry, group-theory]
sources: []
---

## Intuition

A mirror clique is the predecessor clique whose affine span acts as a reflection
hyperplane for a branch decision or block operation.

In the rank-count theory, the mirror clique is the geometric label attached to
a branch mask. The same binary mask can have different obstruction behavior
when interpreted through different mirror cliques.

## Formal definition

For a non-seed vertex $v_i$ in a DDGP instance, the predecessor set

$$
U_i=\{u_1,\dots,u_K\}
$$

is a $K$-vertex clique. Its affine span

$$
H_{U_i}=\operatorname{aff}(U_i)
$$

is the mirror hyperplane exchanging the two generic lateration positions of
$v_i$.

For a graph generator $g=(m_g,C_g)$, $C_g$ is the mirror clique and

$$
H_{C_g}=\operatorname{aff}(C_g)
$$

is the associated mirror hyperplane.
