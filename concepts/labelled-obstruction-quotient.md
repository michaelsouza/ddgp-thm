---
term: Labelled Obstruction Quotient
tags: [linear-algebra, geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The labelled obstruction quotient is the presentation-independent version of
the labelled violation vector.

A branch difference can usually be written as a combination of graph generators
in more than one way. Different presentations may have different labelled
violation vectors, even though they produce the same binary branch mask. The
only labelled information that can be geometrically meaningful is therefore the
violation vector modulo the labels created by pure presentation relations.

This quotient is the correct object for the [labelled presentation
property](labelled-presentation-property.md).

## Formal definition

Let $\mathcal G$ be the graph-generator family, let

$$
M:\mathbb F_2^{\mathcal G}\to\mathbb F_2^B
$$

be the mask map, and let

$$
V_F:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{\mathcal L_F}
$$

be the labelled violation map for an active edge set $F$.

The pure presentation-label space is

$$
R_F=V_F(\ker M).
$$

It records labelled violations that can be introduced or removed by changing a
generator presentation without changing the total branch mask.

The labelled obstruction quotient is

$$
Q_F=\mathbb F_2^{\mathcal L_F}/R_F.
$$

For a branch difference $h\in\mathbb F_2^B$, choose any generator coefficient
vector $\alpha$ with

$$
M\alpha=h.
$$

The canonical labelled obstruction of $h$ is

$$
\omega_F(h)=V_F\alpha+R_F\in Q_F.
$$

This is well-defined because any two presentations of $h$ differ by an element
of $\ker M$.

The projected labelled kernel is exactly

$$
\ker\omega_F=M(\ker V_F).
$$

Thus, the DDGP rank-count theorem can be stated as:

$$
h \text{ preserves all active edges generically}
\quad\Longleftrightarrow\quad
\omega_F(h)=0.
$$

The rank formula computes the dimension of this kernel without explicitly
constructing $Q_F$.
