---
term: Labelled Violation Kernel
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The labelled violation kernel is a linear-algebraic object used to combine local graph-mask symmetries in the presence of multiple pruning edges.

For one pruning edge, it is often enough to keep graph masks that are individually [mirror-compatible](mirror-compatible-graph-mask.md) with that edge. With multiple pruning edges, this is too restrictive: two masks may each violate the same pruning edge, but their composition may cancel the violation.

The cancellation is not purely by pruning-edge parity. Experiments show that the violation should be labelled by both:

- the pruning edge being violated;
- the mirror clique of the graph mask causing the violation.

Thus, two violations cancel only when they have the same labelled type.

## Formal definition

Let $\mathcal{G}$ be a family of graph-derived masks

$$
g=(m_g,C_g),
$$

where $m_g\in\mathbb{F}_2^{B}$ is the binary mask and $C_g$ is its mirror clique. Assume each $g\in\mathcal{G}$ is mirror-compatible with the local base edge set.

Let $P$ be a set of pruning edges. Define the labelled violation set

$$
\mathcal{L}
=
\{(e,C_g) \mid g\in\mathcal{G},\ e\in P,\ g \text{ is not mirror-compatible with } e\}.
$$

The violation signature of a generator $g$ is the vector

$$
\nu(g)\in\mathbb{F}_2^{\mathcal{L}},
$$

where

$$
\nu(g)_{(e,C)}
=
1
\quad\Longleftrightarrow\quad
C=C_g
\text{ and }
g \text{ is not mirror-compatible with } e.
$$

Extend $\nu$ linearly to combinations of generators. The labelled violation kernel is

$$
\ker \nu
=
\{\alpha\in\mathbb{F}_2^{\mathcal{G}} \mid \nu(\alpha)=0\}.
$$

The corresponding mask space is the projection

$$
\mathcal{K}_P
=
\left\{
\bigoplus_{g\in\mathcal{G}}\alpha_g m_g
\ \middle|\
\alpha\in\ker\nu
\right\}
\subseteq \mathbb{F}_2^B.
$$

## Rank identity

Let $M$ be the matrix whose columns are the support masks $m_g$, and let $V$ be the matrix whose columns are the labelled violation vectors $\nu(g)$.

The rank of the projected kernel can be computed without enumerating the kernel:

$$
\operatorname{rank}(\mathcal{K}_P)
=
\operatorname{rank}
\begin{bmatrix}
M\\
V
\end{bmatrix}
-
\operatorname{rank}(V).
$$

This follows by viewing the column span of $\begin{bmatrix}M\\V\end{bmatrix}$ and taking its intersection with the subspace whose violation coordinates are zero.

## Counting hypothesis

For generic DDGP instances with a pruning-edge set $P$, experiments suggest that the local valid solution count is

$$
|\Xi_P|
=
2^{\operatorname{rank}(\mathcal{K}_P)}.
$$

This extends the one-pruning-edge mirror-compatible mask rule to multiple pruning edges. The corresponding algorithmic form is the [DDGP labelled-violation rank count](ddgp-labelled-violation-rank-count.md).
