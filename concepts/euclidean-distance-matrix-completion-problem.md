---
term: Euclidean Distance Matrix Completion Problem (EDMCP)
tags: [geometry, optimization]
sources: [liberti2011number]
---

## Intuition

The Euclidean Distance Matrix Completion Problem (EDMCP) is a classical problem in distance geometry that is closely related to the [Distance Geometry Problem (DGP)](distance-geometry-problem.md). In both problems, we are given a graph where some edges have weights representing Euclidean distances, and we want to find a corresponding set of coordinates in a Euclidean space.

The fundamental difference between DGP and EDMCP is how the dimension $K$ of the embedding space is treated:
* In the **DGP**, the dimension $K$ is part of the **input**. We are asked to embed the graph in a specific dimension (for example, $K=3$ for molecules or $K=2$ for sensor networks). This constraint makes the DGP strongly NP-hard even for small values of $K$.
* In the **EDMCP**, the dimension $K$ is part of the **output**. We are simply asked if there exists *any* dimension $K$ in which the graph can be embedded. Alternatively, we may want to find the minimum dimension $K$ that allows for a valid embedding.

If the dimension $K$ is not constrained, the problem can be solved in polynomial time using semidefinite programming (SDP) and analyzing the eigenvalues of the completed matrix. However, finding the minimum dimension $K$ or solving the problem when $K$ is fixed (which reduces to the DGP) remains computationally challenging.

## Formal definition

Let $G = (V, E)$ be a simple undirected graph, and let $d: E \to \mathbb{R}$ be a partial distance function. The **Euclidean Distance Matrix Completion Problem (EDMCP)** asks whether there exists a dimension $K \in \mathbb{N}$ and a set of points $x_1, \dots, x_n \in \mathbb{R}^K$ such that:

$$
\|x_i - x_j\| = d_{ij}, \quad \forall \{i,j\} \in E
$$

Equivalently, let $D \in \mathbb{R}^{n \times n}$ be a symmetric matrix with zero diagonal ($D_{ii} = 0$). We say $D$ is a Euclidean Distance Matrix (EDM) if there exist points $x_1, \dots, x_n$ in some Euclidean space such that $D_{ij} = \|x_i - x_j\|^2$. The EDMCP is the problem of completing a partially defined symmetric matrix $A$ (where only entries corresponding to edges in $E$ are given) into a full EDM:

$$
\text{Find } D \in \text{EDM} \quad \text{such that} \quad D_{ij} = d_{ij}^2, \quad \forall \{i,j\} \in E
$$

If a completion exists, the coordinates $x_i$ can be retrieved by finding a Gram matrix $B$ using the Schoenberg transform:

$$
B = -\frac{1}{2} J D J
$$

Where $J = I - \frac{1}{n}\mathbf{1}\mathbf{1}^T$ is the centering matrix. The matrix $D$ is a valid EDM if and only if $B$ is positive semidefinite ($\succeq 0$). The rank of $B$ then corresponds to the minimum dimension $K$ required for the embedding.
