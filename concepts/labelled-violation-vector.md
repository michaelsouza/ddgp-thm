---
term: Labelled Violation Vector
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

A labelled violation vector is the column of obstruction data associated with
one graph generator.

It records which active edges the generator violates, while keeping the mirror
clique label attached to each violation. These vectors become the columns of the
labelled violation matrix.

## Formal definition

Let $g=(m_g,C_g)$ be a DDGP graph generator and let $\mathcal L_F$ be the
[labelled violation set](labelled-violation-set.md). The labelled violation
vector of $g$ is

$$
\nu_F(g)\in\mathbb F_2^{\mathcal L_F}
$$

defined by

$$
\nu_F(g)_{(e,C)}=1
\quad\Longleftrightarrow\quad
C=C_g \text{ and } g\text{ violates }e.
$$

The [labelled violation matrix](labelled-violation-matrix.md) $V_F$ is obtained
by placing these vectors as columns.
