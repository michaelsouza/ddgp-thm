Information Sciences 559 (2021) 1-7

ELSEVIER

Contents lists available at ScienceDirect

Information Sciences

journal homepage: www.elsevier.com/locate/ins

# A note on the Cayley-Menger determinant and the Molecular Distance Geometry Problem

![img-0.jpeg](img-0.jpeg)

Luiz Leduino de Salles Neto $^{a,\ast}$, Carlile Lavor $^{b}$, Weldon Lodwick $^{c}$

$^{a}$ Federal University of São Paulo, Brazil
$^{b}$ University of Campinas, Brazil
$^{c}$ University of Colorado Denver, United States

# ARTICLE INFO

Article history:

Received 17 January 2020

Received in revised form 20 September 2020

Accepted 27 December 2020

Available online 27 January 2021

Keywords:

Cayley-Menger determinant

Molecular Distance Geometry Problem

Nonlinear optimization.

# ABSTRACT

Using information from protein geometry and distance data provided by Nuclear Magnetic Resonance (NMR) experiments, the Molecular Distance Geometry Problem (MDGP) can be solved by a combinatorial approach, called Branch-and-Prune (BP). The primal version of BP algorithm seeks MDGP graph realizations, while the dual BP looks for completions of associated partial distance matrices. These two algorithms are very similar when distance values are precise. In the literature, there are some proposals for extending the primal BP to take care of NMR uncertainties. Using Cayley-Menger determinant, we present a global optimization approach that also allows the dual BP to deal with NMR uncertainties.

© 2021 Elsevier Inc. All rights reserved.

# 1. Introduction

The origin of Distance Geometry (DG) [25] is associated to Karl Menger, when he characterized geometric ideas using the concept of distance [26], and one of the most important applications of DG is in computational biology, where the 3D protein structure is calculated using distance data provided by Nuclear Magnetic Resonance (NMR) experiments [5].

The problem is called the Molecular Distance Geometry Problem (MDGP) [24], formally defined in terms of a graph  $G = (V, E, d)$ , where the sets  $V, E$  represent the atoms and the atomic pairs for which a distance is available, respectively (we say that  $G$  is simple if  $u \neq v$ , for  $\{u, v\} \in E$ , and undirected when  $\{u, v\} = \{v, u\}$ ):

Definition 1. Given a simple undirected graph  $G = (V, E, d)$  whose edges are weighted by  $d: E \to (0, \infty)$ , find a function  $x: V \to \mathbb{R}^3$  such that

$$
\forall \{u, v \} \in E, \| x _ {u} - x _ {v} \| _ {2} = d _ {u, v}, \tag {1}
$$

where  $x_{u} = x(u), x_{v} = x(v), d_{u,v} = d(\{u,v\})$

The Euclidean norm and the 3D space are the most common in molecular geometry, but for other applications, different norms and spaces  $\mathbb{R}^k$  ( $k \neq 3$ ) can be used [3,20,21,28].

The MDGP can be solved iteratively using a combinatorial method, called Branch-and-Prune (BP) [6,22], based on a vertex order [7,17,19]  $v_{1},\ldots ,v_{n}\in V$  such that (we denote  $x_{i}$  instead of  $x_{v_i}$  and  $d_{i,j}$  instead of  $d_{v_i v_j}$ ):

https://doi.org/10.1016/j.ins.2020.12.072

0020-0255/© 2021 Elsevier Inc. All rights reserved.

---

Luiz Leduino de Salles Neto, C. Lavor and W. Lodwick

Information Sciences 559 (2021) 1-7

1. For  $v_{1}, v_{2}, v_{3}$ , there exist  $x_{1}, x_{2}, x_{3} \in \mathbb{R}^{3}$  satisfying Eqs. (1);
2. For  $i &gt; 3$ , there exist  $v_{i-3}, v_{i-2}, v_{i-1}$  such that

$$
\left\{v _ {i - 3}, v _ {i} \right\}, \left\{v _ {i - 2}, v _ {i} \right\}, \left\{v _ {i - 1}, v _ {i} \right\} \in E; \tag {2}
$$

3. For  $i &gt; 3$

$$
d _ {i - 3, i - 2} + d _ {i - 2, i - 1} &gt; d _ {i - 3, i - 1}. \tag {3}
$$

MDGP instances with this order are called the Discretizable Molecular Distance Geometry Problem (DMDGP) [14,15].

Property 1 says that we can start with fixed positions for  $v_{1}$ ,  $v_{2}$ ,  $v_{3}$ . Property 2 allows us to define a system

$$
\left| \left| x _ {i} - x _ {i - 1} \right| \right| _ {2} = d _ {i - 1, i}, \tag {4}
$$

$$
\left| \left| x _ {i} - x _ {i - 2} \right| \right| _ {2} = d _ {i - 2, i},
$$

$$
\left| \left| x _ {i} - x _ {i - 3} \right| \right| _ {2} = d _ {i - 3, i},
$$

which provides up to two possible positions for  $v_{i}, i &gt; 3$  (from Property 3), given  $x_{i-1}, x_{i-2}, x_{i-3} \in \mathbb{R}^{3}$ . For each position for  $v_{i}$  we obtain another two for  $v_{i+1}$ , implying that the DMDGP search space can be represented by a binary tree [15].

For  $i &gt; 4$ , we have to check if there is an additional edge  $\{v_i, v_i\} \in E, j &lt; i - 3$ , implying another equation to the system related to  $v_i: ||x_i - x_j||_2 = d_{i,i}$  (such value  $d_{i,i}$  is called a pruning distance). If the points  $x_i, x_{i-1}, x_{i-2}, x_{i-3} \in \mathbb{R}^3$  are not in the same plane, we obtain a unique solution  $x_i^*$  for  $v_i$ , supposing  $||x_i^* - x_j||_2 = d_{i,i}$ . However, both positions for  $v_i$  may not satisfy the additional equation and must be pruned. In that case, we consider the other possible positions for  $v_{i-1}$  and repeat the procedure [15]. The search ends when a path from the root of the tree to a leaf node is found, in such a way that the positions for all the vertices in the DMDGP order satisfy Eqs. (1).

Distances  $d_{i-1,i}$  and  $d_{i-2,i}$  are known a priori, considered as precise values, since they are related to bond lengths and bond angles of a protein. However, if distances  $d_{i-3,i}$  are provided by NMR experiments, they contain uncertainties [4].

In [16], distances  $d_{i-3,i}$  were represented as interval distances  $\left[d_{i-3,i}, \overline{d}_{i-3,i}\right]$ , where  $d_{i-3,i} \leqslant d_{i-3,i} \leqslant \overline{d}_{i-3,i}$ , and an extension of the BP algorithm was proposed, called iBP, where the idea is to sample values from  $\left[d_{i-3,i}, \overline{d}_{i-3,i}\right]$  in order to solve the system (4) related to  $v_i$  [10,12]. Choosing many values, the search space increases exponentially, and for small samples, a solution may not be found [7,30].

For the branching phase of iBP, there are two recent and different approaches to avoid such sampling strategy: one is based on geometric algebra [1,2] and the other uses single-level interval analysis [9]. For the pruning phase, however, there is not a well established solution method. Some results are given in [11,18].

All strategies mentioned above for dealing with NMR uncertainties are given in terms of the atomic Cartesian coordinates. In [23], a dual BP algorithm is proposed to solve DMDGP instances based only on the calculation of (precise) distance values. Since the NMR data are provided as interval distances, the dual BP must be able to consider this kind of data. This is our motivation in this paper.

The next section explains all the details of our contribution and Section 3 provides some computational results. We conclude the paper with new research directions.

## 2. Cayley-Menger Determinant and the DMDGP

In the definition of the DMDGP, the strict triangular inequality (Property 3) is required to guarantee a finite search space, and since the bond angles in protein molecules are not linear, this is a realistic assumption [5].

The triangular inequality can be generalized using the so called Cayley-Menger determinant [8,27], defined for a set  $S_{n+1} = \{p_0, p_1, \ldots, p_n\}$  of points in  $\mathbb{R}^n$ , given by

$$
C M (S _ {n + 1}) = \left| \begin{array}{c c c c c} 0 &amp; 1 &amp; 1 &amp; 1 &amp; 1 \\ 1 &amp; 0 &amp; d _ {0, 1} ^ {2} &amp; \dots &amp; d _ {0, n} ^ {2} \\ 1 &amp; d _ {0, 1} ^ {2} &amp; 0 &amp; \dots &amp; d _ {1, n} ^ {2} \\ \vdots &amp; \vdots &amp; \vdots &amp; \ddots &amp; \vdots \\ 1 &amp; d _ {0, n} ^ {2} &amp; d _ {1, n} ^ {2} &amp; \dots &amp; 0 \end{array} \right|,
$$

where  $d_{ij} = ||p_i - p_j||_2, \forall i,j \in \{0,\dots,n\}$ .

It is possible to prove that five points in  $\mathbb{R}^3$  have zero Cayley-Menger determinant [29], and using this fact for the DMDGP order, considering  $S_{5} = \{x_{i - 4}, x_{i - 3}, x_{i - 2}, x_{i - 1}, x_{i}\}$ , we get

---

Luiz Leduino de Salles Neto, C. Lavor and W. Lodwick

Information Sciences 559 (2021) 1-7

$$
C M (S _ {5}) = 0 \Rightarrow \left| \begin{array}{c c c c c c} 0 &amp; 1 &amp; 1 &amp; 1 &amp; 1 &amp; 1 \\ 1 &amp; 0 &amp; d _ {i - 4, i - 3} ^ {2} &amp; d _ {i - 4, i - 2} ^ {2} &amp; d _ {i - 4, i - 1} ^ {2} &amp; \delta^ {2} \\ 1 &amp; d _ {i - 4, i - 3} ^ {3} &amp; 0 &amp; d _ {i - 3, i - 2} ^ {2} &amp; d _ {i - 3, i - 1} ^ {2} &amp; d _ {i - 3, i} ^ {2} \\ 1 &amp; d _ {i - 4, i - 2} ^ {2} &amp; d _ {i - 3, i - 2} ^ {2} &amp; 0 &amp; d _ {i - 2, i - 1} ^ {2} &amp; d _ {i - 2, i} ^ {2} \\ 1 &amp; d _ {i - 4, i - 1} ^ {2} &amp; d _ {i - 3, i - 1} ^ {2} &amp; d _ {i - 2, i - 1} ^ {2} &amp; 0 &amp; d _ {i - 1, i} ^ {2} \\ 1 &amp; \delta^ {2} &amp; d _ {i - 3, i} ^ {2} &amp; d _ {i - 2, i} ^ {2} &amp; d _ {i - 1, i} ^ {2} &amp; 0 \end{array} \right| = 0,
$$

resulting in a quadratic equation with two possible values for  $\delta^2$ , where  $\delta$  is the distance associated to the pair  $\{i - 4, i\}$ . This is the equation that must be solved, in the dual BP, replacing the primal Eq. (4).

When NMR uncertainties are taken into account, distance values  $d_{i-4,i-1}$  and  $d_{i-3,i}$  may be represented as interval distances. Instead of sampling values from intervals  $\left[\underline{d}_{i-4,i-1}, \overline{d}_{i-4,i-1}\right]$  and  $\left[\underline{d}_{i-3,i}, \overline{d}_{i-3,i}\right]$ , we propose a new way to calculate the possible values for  $\delta$ .

## 2.1. Calculating intervals for  $d_{i-4,i}$

Rewriting  $CM(S_5)$  in terms of  $\underline{d}_{i-4,i-1}$  and  $\underline{d}_{i-3,i}$ , with  $\lambda_1 \in \left[0, \overline{d}_{i-4,i-1} - \underline{d}_{i-4,i-1}\right]$  and  $\lambda_2 \in \left[0, \overline{d}_{i-3,i} - \underline{d}_{i-3,i}\right]$ , we have

$$
\left| \begin{array}{c c c c c c} 0 &amp; 1 &amp; 1 &amp; 1 &amp; 1 &amp; 1 \\ 1 &amp; 0 &amp; d _ {i - 4, i - 3} ^ {2} &amp; d _ {i - 4, i - 2} ^ {2} &amp; \left(\underline {{d}} _ {i - 4, i - 1} + \lambda_ {1}\right) ^ {2} &amp; \delta^ {2} \\ 1 &amp; d _ {i - 4, i - 3} ^ {2} &amp; 0 &amp; d _ {i - 3, i - 2} ^ {2} &amp; d _ {i - 3, i - 1} ^ {2} &amp; \left(\underline {{d}} _ {i - 3, i} + \lambda_ {2}\right) ^ {2} \\ 1 &amp; d _ {i - 4, i - 2} ^ {2} &amp; d _ {i - 3, i - 2} ^ {2} &amp; 0 &amp; d _ {i - 2, i - 1} ^ {2} &amp; d _ {i - 2, i} ^ {2} \\ 1 &amp; \left(\underline {{d}} _ {i - 4, i - 1} + \lambda_ {1}\right) ^ {2} &amp; d _ {i - 3, i - 1} ^ {2} &amp; d _ {i - 2, i - 1} ^ {2} &amp; 0 &amp; d _ {i - 1, i} ^ {2} \\ 1 &amp; \delta^ {2} &amp; \left(\underline {{d}} _ {i - 3, i} + \lambda_ {2}\right) ^ {2} &amp; d _ {i - 2, i} ^ {2} &amp; d _ {i - 1, i} ^ {2} &amp; 0 \end{array} \right| = 0.
$$

Solving this equation, we obtain the following functions:

$$
\begin{array}{l} \delta_ {1} ^ {2} \left(\lambda_ {1}, \lambda_ {2}\right) = \left[ c _ {1} + c _ {2} \lambda_ {1} + c _ {3} \lambda_ {1} ^ {2} + c _ {4} \lambda_ {2} + c _ {5} \lambda_ {1} \lambda_ {2} \right. \\ + c _ {6} \lambda_ {1} ^ {2} \lambda_ {2} + c _ {7} \lambda_ {2} ^ {2} + c _ {8} \lambda_ {1} ^ {2} \lambda_ {2} ^ {2} + c _ {9} \lambda_ {1} \lambda_ {2} ^ {2} \\ + \left(c _ {1 0} + c _ {1 1} \lambda_ {1} + c _ {1 2} \lambda_ {1} ^ {2} + c _ {1 3} \lambda_ {2} + c _ {1 4} \lambda_ {1} \lambda_ {2} \right. \\ + c _ {1 5} \lambda_ {1} ^ {2} \lambda_ {2} + c _ {1 6} \lambda_ {1} ^ {2} \lambda_ {2} ^ {2} \\ + c _ {1 7} \lambda_ {1} ^ {3} + c _ {1 8} \lambda_ {1} ^ {4} + c _ {1 9} \lambda_ {1} ^ {3} \lambda_ {2} + c _ {2 0} \lambda_ {1} ^ {4} \lambda_ {2} + c _ {2 1} \lambda_ {2} ^ {2} \\ + c _ {2 2} \lambda_ {1} \lambda_ {2} ^ {2} + + c _ {2 3} \lambda_ {1} ^ {3} \lambda_ {2} ^ {2} + c _ {2 4} \lambda_ {1} ^ {4} \lambda_ {2} ^ {2} \\ + c _ {2 5} \lambda_ {2} ^ {3} + c _ {2 6} \lambda_ {1} \lambda_ {2} ^ {3} + c _ {2 7} \lambda_ {1} ^ {2} \lambda_ {2} ^ {3} + c _ {2 8} \lambda_ {1} ^ {3} \lambda_ {2} ^ {3} + c _ {2 9} \lambda_ {1} ^ {4} \lambda_ {2} ^ {3} \\ \left. + c _ {3 0} \lambda_ {2} ^ {4} + c _ {3 1} \lambda_ {1} \lambda_ {2} ^ {4} + c _ {3 2} \lambda_ {1} ^ {3} \lambda_ {2} ^ {4} + c _ {3 3} \lambda_ {1} ^ {3} \lambda_ {2} ^ {4} + c _ {3 4} \lambda_ {1} ^ {4} \lambda_ {2} ^ {4}\right) ^ {1 / 2} ] \times c _ {3 5} \\ \end{array}
$$

and

$$
\begin{array}{l} \delta_ {2} ^ {2} \left(\lambda_ {1}, \lambda_ {2}\right) = \left[ c _ {1} + c _ {2} \lambda_ {1} + c _ {3} \lambda_ {1} ^ {2} + c _ {4} \lambda_ {2} + c _ {5} \lambda_ {1} \lambda_ {2} \right. \\ + c _ {6} \lambda_ {1} ^ {2} \lambda_ {2} + c _ {7} \lambda_ {2} ^ {2} + c _ {8} \lambda_ {1} ^ {2} \lambda_ {2} ^ {2} + c _ {9} \lambda_ {1} \lambda_ {2} ^ {2} \\ - \left(c _ {1 0} + c _ {1 1} \lambda_ {1} + c _ {1 2} \lambda_ {1} ^ {2} + c _ {1 3} \lambda_ {2} + c _ {1 4} \lambda_ {1} \lambda_ {2} \right. \\ + c _ {1 5} \lambda_ {1} ^ {2} \lambda_ {2} + c _ {1 6} \lambda_ {1} ^ {2} \lambda_ {2} ^ {2} \\ + c _ {1 7} \lambda_ {1} ^ {3} + c _ {1 8} \lambda_ {1} ^ {4} + c _ {1 9} \lambda_ {1} ^ {3} \lambda_ {2} + c _ {2 0} \lambda_ {1} ^ {4} \lambda_ {2} + c _ {2 1} \lambda_ {2} ^ {2} \\ + c _ {2 2} \lambda_ {1} \lambda_ {2} ^ {2} + + c _ {2 3} \lambda_ {1} ^ {3} \lambda_ {2} ^ {3} + c _ {2 4} \lambda_ {1} ^ {4} \lambda_ {2} ^ {3} \\ + c _ {2 5} \lambda_ {2} ^ {3} + c _ {2 6} \lambda_ {1} \lambda_ {2} ^ {3} + c _ {2 7} \lambda_ {1} ^ {2} \lambda_ {2} ^ {3} + c _ {2 8} \lambda_ {1} ^ {3} \lambda_ {2} ^ {3} + c _ {2 9} \lambda_ {1} ^ {4} \lambda_ {2} ^ {3} \\ \left. + c _ {3 0} \lambda_ {2} ^ {4} + c _ {3 1} \lambda_ {1} \lambda_ {2} ^ {4} + c _ {3 2} \lambda_ {1} ^ {2} \lambda_ {2} ^ {4} + c _ {3 3} \lambda_ {1} ^ {3} \lambda_ {2} ^ {4} + c _ {3 4} \lambda_ {1} ^ {4} \lambda_ {2} ^ {4}\right) ^ {1 / 2} ] \times c _ {3 5}, \\ \end{array}
$$

where the values for  $c_{i}$ ,  $i = 1, \dots, 35$ , are available at https://cutt.ly/BfPgzoe.

Since  $\delta_1^2 (\lambda_1,\lambda_2)$  and  $\delta_2^2 (\lambda_1,\lambda_2)$  are continuous functions, for  $\lambda_{1}\in \left[0,\overline{d}_{i - 4,i - 1} - \underline{d}_{i - 4,i - 1},\right]$  and  $\lambda_{2}\in \left[0,\overline{d}_{i - 3,i} - \underline{d}_{i - 3,i},\right],\delta_1^2 (\lambda_1,\lambda_2)$  and  $\delta_2^2 (\lambda_1,\lambda_2)$  result in other two intervals, implying that the possible values for  $d_{i - 4,i}$ $(\delta_t,t = 1,2)$  can then be calculated by solving the following optimization problems:

---

Luiz Leduino de Salles Neto, C. Lavor and W. Lodwick

Information Sciences 559 (2021) 1-7

$$
\min  _ {\lambda_ {1}, \lambda_ {2} \in \mathbb {R}} \quad \delta_ {i} ^ {2} \left(\lambda_ {1}, \lambda_ {2}\right) \tag {7}
$$

$$
\lambda_ {1} \in \left[ 0, \bar {d} _ {i - 4, i - 1} - \underline {{d}} _ {i - 4, i - 1} \right] \quad \lambda_ {2} \in \left[ 0, \bar {d} _ {i - 3, i} - \underline {{d}} _ {i - 3, i} \right]
$$

and

$$
\max  _ {\lambda_ {1}, \lambda_ {2} \in \mathbb {R}} \quad \delta_ {i} ^ {2} \left(\lambda_ {1}, \lambda_ {2}\right) \tag {8}
$$

$$
\lambda_ {1} \in \left[ 0, \bar {d} _ {i - 4, i - 1} - \underline {{d}} _ {i - 4, i - 1} \right] \quad \lambda_ {2} \in \left[ 0, \bar {d} _ {i - 3, i} - \underline {{d}} _ {i - 3, i} \right]
$$

## 3. Computational results

Since the proposed method can be used for any sequence of five connected atoms, we will generate random DMDGP instances with that size, according to the suggestions in [13].

By fixing the lengths of covalent bonds $(d_{i-1,i} = 1.0)$ and the values of covalent angles $(\theta_{i-2,i} = 2.0$ radians), instances with five atoms are defined by torsion angles $\omega_{i-4,i-1} \in [0,2\pi]$ and $\omega_{i-3,i} \in [0,2\pi]$. We randomly select a value $\omega$ from the set $\{60^{\circ}, 180^{\circ}, 300^{\circ}\}$, making a perturbation randomly selecting another value from $[0.9\omega, 1.1\omega]$.

In order to calculate the associated Cartesian coordinates of the points $(x_{i},y_{i},z_{i})\in \mathbb{R}^{3}$, $i = 1,\ldots ,5$, we use the following homogeneous matrices [15]:

$$
\left[ \begin{array}{c} x _ {1} \\ y _ {1} \\ z _ {1} \\ 1 \end{array} \right] = B _ {1} \left[ \begin{array}{c} 0 \\ 0 \\ 0 \\ 1 \end{array} \right], \quad \left[ \begin{array}{c} x _ {2} \\ y _ {2} \\ z _ {2} \\ 1 \end{array} \right] = B _ {1} B _ {2} \left[ \begin{array}{c} 0 \\ 0 \\ 0 \\ 1 \end{array} \right], \quad \left[ \begin{array}{c} x _ {3} \\ y _ {3} \\ z _ {3} \\ 1 \end{array} \right] = B _ {1} B _ {2} B _ {3} \left[ \begin{array}{c} 0 \\ 0 \\ 0 \\ 1 \end{array} \right]
$$

and

$$
\left[ \begin{array}{c} x _ {4} \\ y _ {4} \\ z _ {4} \\ 1 \end{array} \right] = B _ {1} B _ {2} B _ {3} B _ {4} \left[ \begin{array}{c} 0 \\ 0 \\ 0 \\ 1 \end{array} \right], \quad \left[ \begin{array}{c} x _ {5} \\ y _ {5} \\ z _ {5} \\ 1 \end{array} \right] = B _ {1} B _ {2} B _ {3} B _ {4} B _ {5} \left[ \begin{array}{c} 0 \\ 0 \\ 0 \\ 1 \end{array} \right],
$$

where

$$
B _ {1} = \left[ \begin{array}{c c c c} 1 &amp; 0 &amp; 0 &amp; 0 \\ 0 &amp; 1 &amp; 0 &amp; 0 \\ 0 &amp; 0 &amp; 1 &amp; 0 \\ 0 &amp; 0 &amp; 0 &amp; 1 \end{array} \right], \quad B _ {2} = \left[ \begin{array}{c c c c} - 1 &amp; 0 &amp; 0 &amp; - d _ {1, 2} \\ 0 &amp; 1 &amp; 0 &amp; 0 \\ 0 &amp; 0 &amp; - 1 &amp; 0 \\ 0 &amp; 0 &amp; 0 &amp; 1 \end{array} \right],
$$

$$
B _ {3} = \left[ \begin{array}{c c c c} - \cos \theta_ {1, 3} &amp; - \sin \theta_ {1, 3} &amp; 0 &amp; - d _ {2, 3} \cos \theta_ {1, 3} \\ \sin \theta_ {1, 3} &amp; - \cos \theta_ {1, 3} &amp; 0 &amp; d _ {2, 3} \sin \theta_ {1, 3} \\ 0 &amp; 0 &amp; 1 &amp; 0 \\ 0 &amp; 0 &amp; 0 &amp; 1 \end{array} \right],
$$

and

$$
B _ {i} = \left[ \begin{array}{c c c c} - \cos \theta_ {i - 2, i} &amp; - \sin \theta_ {i - 2, i} &amp; 0 &amp; - d _ {i - 1, i} \cos \theta_ {i - 2, i} \\ \sin \theta_ {i - 2, i} \cos \omega_ {i - 3, i} &amp; - \cos \theta_ {i - 2, i} \cos \omega_ {i - 3, i} &amp; - \sin \omega_ {i - 3, i} &amp; d _ {i - 1, i} \sin \theta_ {i - 2, i} \cos \omega_ {i - 3, i} \\ \sin \theta_ {i - 2, i} \sin \omega_ {i - 3, i} &amp; - \cos \theta_ {i - 2, i} \sin \omega_ {i - 3, i} &amp; \cos \omega_ {i - 3, i} &amp; d _ {i - 1, i} \sin \theta_ {i - 2, i} \sin \omega_ {i - 3, i} \\ 0 &amp; 0 &amp; 0 &amp; 1 \end{array} \right],
$$

for $i\in \{4,5\}$.

From the calculations above, we can obtain the distances between the atomic pairs $i - 4, i - 1$ and $i - 3, i$.

We will assume that these distances are equal to $\underline{d}_{i-4,i-1}$ and $\underline{d}_{i-3,i}$, respectively, and uncertainties will be given by $\lambda_1, \lambda_2 \in \{0.05, 0.10, 0.15, 0.20\}$.

To solve all the instances, we used the software AMPL, with the solver Baron 17.4.1, on a Lenovo notebook with 6 GB RAM and Intel celeron 1.6 GHz.

---

Luiz Leduino de Salles Neto, C. Lavor and W. Lodwick

Information Sciences 559 (2021) 1-7

All generated data and the AMPL model file are available at https://cutt.ly/BfPgzoe. Tables 1-4 provide the results.

# 4. Conclusions

The dual BP algorithm, proposed in [23], works in the distance space. This means that while the primal BP seeks graph realizations, the dual BP seeks completions of partial distance matrices, making realizations and distance matrix completions dual to each other.

When distances are precise values, both algorithms are similar. However, to deal with interval distances, the dual BP appears to be more suitable, since the input data of the problem are distance values. For the first time, this paper presents some new results on this research direction.

The intervals obtained for  $d_{i-4,i}$  are the mathematically possible values, calculated in terms of the distances

$$
d _ {i - 4, i - 3}, d _ {i - 4, i - 2}, d _ {i - 3, i - 2}, d _ {i - 3, i - 1}, d _ {i - 2, i - 1}, d _ {i - 2, i}, d _ {i - 1, i},
$$

and

$$
\left[ \underline {{d}} _ {i - 3, i}, \overline {{d}} _ {i - 3, i} \right], \left[ \underline {{d}} _ {i - 4, i - 1}, \overline {{d}} _ {i - 4, i - 1} \right],
$$

which can be very useful to validate NMR data.

Table 1 Computational time and solutions obtained for  $\lambda_1\leqslant 0.05$  and  $\lambda_{2}\leqslant 0.05$

|  ωi-4,i-1 | ωi-3,i | di-4,i-1 | di-3,i | di-4,i | Time (s)  |
| --- | --- | --- | --- | --- | --- |
|  300.08 | 327.63 | 2.0449 | 1.9010 | [1.7675, 1.8732] ∪ [2.1577, 2.3683] | 0.41  |
|  177.78 | 191.20 | 2.5813 | 2.5754 | [3.3584, 3.5107] | 0.25  |
|  188.12 | 321.31 | 2.5783 | 1.9287 | [2.8652, 2.9541] | 0.36  |
|  302.14 | 62.63 | 2.0325 | 2.0617 | [1.8973, 2.0202] ∪ [2.4671, 2.6322] | 0.27  |
|  196.03 | 323.16 | 2.5690 | 1.9202 | [2.8223, 2.9615] | 0.20  |

Table 2 Computational time and solutions obtained for  $\lambda_1\leqslant 0.10$  and  $\lambda_{2}\leqslant 0.10$

|  ωi-4,i-1 | ωi-3,i | di-4,i-1 | di-3,i | di-4,i | Time (s)  |
| --- | --- | --- | --- | --- | --- |
|  325.72 | 180.93 | 1.9090 | 2.5815 | [2.8798, 2.9561] | 0.28  |
|  167.72 | 295.05 | 2.5742 | 2.0761 | [2.9319, 3.1130] | 0.30  |
|  63.77 | 180.52 | 2.0663 | 2.5815 | [2.9867, 3.0640] | 0.48  |
|  171.15 | 163.74 | 2.577 | 2.5687 | [3.3363, 3.6520] | 0.32  |
|  183.58 | 302.11 | 2.5809 | 2.0326 | [2.9501, 3.0507] | 0.29  |

Table 3 Computational time and solutions obtained for  $\lambda_1\leqslant 0.15$  and  $\lambda_{2}\leqslant 0.15$

|  ωi-4,i-1 | ωi-3,i | di-4,i-1 | di-3,i | di-4,i | Time (s)  |
| --- | --- | --- | --- | --- | --- |
|  282.06 | 56.52 | 2.1599 | 2.0245 | [2.0308, 2.4121] ∪ [2.619, 3.0026] | 0.32  |
|  176.98 | 288.46 | 2.5811 | 2.1182 | [3.0109, 3.1453] | 0.26  |
|  308.00 | 62.73 | 1.9982 | 2.0623 | [1.8619, 2.2347] ∪ [2.4048, 2.8633] | 0.37  |
|  192.77 | 57.69 | 2.5736 | 2.0315 | [2.9009, 3.1184] | 0.26  |
|  175.27 | 62.56 | 2.5804 | 2.0610 | [2.9635, 3.1101] | 0.24  |

Table 4 Computational time and solutions obtained for  $\lambda_1\leqslant 0.20$  and  $\lambda_{2}\leqslant 0.20$

|  ωi-4,i-1 | ωi-3,i | di-4,i-1 | di-3,i | di-4,i | Time (s)  |
| --- | --- | --- | --- | --- | --- |
|  162.37 | 322.81 | 2.5664 | 1.9218 | [2.8139, 3.0860] | 0.24  |
|  58.16 | 327.80 | 2.0344 | 1.9004 | [1.7496, 2.8284] | 0.28  |
|  57.84 | 291.12 | 2.0325 | 2.1012 | [1.9507, 3.0548] | 0.27  |
|  307.47 | 197.10 | 2.0012 | 2.5673 | [2.8568, 3.1468] | 0.26  |
|  55.26 | 290.70 | 2.0171 | 2.1038 | [1.9394, 3.0419] | 0.24  |

---

Luiz Leduino de Salles Neto, C. Lavor and W. Lodwick

Information Sciences 559 (2021) 1-7

Differently from the "classical" way to deal with interval distances, where the dual BP should sample values from the intervals  $\left[d_{i-3,i}, \overline{d}_{i-3,i}\right]$  and  $\left[d_{i-4,i-1}, \overline{d}_{i-4,i-1}\right]$  in order to calculate  $d_{i-4,i}$ , the new approach proposed here avoids sampling process and allows us to obtain all the possible values for  $d_{i-4,i}$ .

If an interval distance  $d_{i-4,i}$  is provided by NMR data and the sampled value is outside this interval, another samples must be selected from  $\left[d_{i-3,i}, \overline{d}_{i-3,i}\right]$  and  $\left[d_{i-4,i-1}, \overline{d}_{i-4,i-1}\right]$ . Using the proposed approach, we just intersect all the related intervals. For example, considering the last line of Table 1 and supposing that NMR provides  $d_{i-4,i} \in [2.91, 3.55]$ , the new interval for  $d_{i-4,i}$  would be

$$
d _ {i - 4, i} \in [ 2. 8 2 2 3, 2. 9 6 1 5 ] \cap [ 2. 9 1, 3. 5 5 ] = [ 2. 9 1, 2. 9 6 ].
$$

Computational experiments indicate that torsion angle values are related to the fact that distance values for  $d_{i-4,i}$  can be given by an interval or an union of intervals. Our further study will investigate this in a more detailed way.

# CRediT authorship contribution statement

Luiz Leduino de Salles Neto: Investigation, Software, Data curation. Carlile Lavor: Methodology, Validation, Writing - review &amp; editing. Weldon Lodwick: Conceptualization, Writing - original draft.

# Declaration of Competing Interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

# Acknowledgements

We would like to thank the Brazilian research agencies CNPq and FAPESP for their financial support and the valuable comments made by the reviewers.

# Appendix A. Supplementary data

Supplementary data associated with this article can be found, in the online version, at https://doi.org/10.1016/j.ins.2020.12.072.

# References

[1] R. Alves, C. Lavor, Geometric algebra to model uncertainties in the discretizable molecular distance geometry problem, Advances in Applied Clifford Algebra 27 (2017) 439-452.
[2] R. Alves, C. Lavor, C. Souza, M. Souza, Clifford algebra and discretizable distance geometry, Mathematical Methods in the Applied Sciences 41 (2018) 3999-4346.
[3] A. Baez-Sanchez, C. Lavor, On the estimation of unknown distances for a class of Euclidean distance matrix completion problems with interval data, Linear Algebra and its Applications 592 (2020) 287-305.
[4] S. Billinge, P. Duxbury, D. Gonçalves, C. Lavor, A. Mucherino, Assigned and unassigned distance geometry: applications to biological molecules and nanostructures, 4OR, 14 (2016) 337-376.
[5] S. Billinge, P. Duxbury, D. Gonçalves, C. Lavor, A. Mucherino, Recent results on assigned and unassigned distance geometry with applications to proteinmolecules and nanostructures, Annals of Operations Research 271 (2018) 161-203.
[6] R. Carvalho, C. Lavor, F. Protti, Extending the geometric build-up algorithm for the molecular distance geometry problem, Information Processing Letters 108 (2008) 234-237.
[7] A. Cassioli, O. Gunluk, C. Lavor, L. Liberti, Discretization vertex orders in distance geometry, Discrete Applied Mathematics 197 (2015) 27-41.
[8] A. Cayley, A theorem in the geometry of position, Cambridge Mathematical Journal 2 (1941) 267-271.
[9] T. Costa, H. Bouwmeester, W. Lodwick, C. Lavor, Calculating the possible conformations arising from uncertainty in the molecular distance geometry problem using constraint interval analysis, Information Sciences 415-416 (2017) 41-52.
[10] C. Dambrosio, V. Ky, C. Lavor, L. Liberti, N. Maculan, New error measures and methods for realizing protein graphs from distance data, Discrete and Computational Geometry 57 (2017) 371-418.
[11] L. Dorst, Boolean combination of circular arcs using orthogonal spheres, Advances in Applied Clifford Algebra 29 (2019) 1-21.
[12] D. Gonçalves, A. Mucherino, C. Lavor, L. Liberti, Recent advances on the interval distance geometry problem, Journal of Global Optimization 69 (2017) 525-545.
[13] C. Lavor, On generating instances for the molecular distance geometry problem, in: L. Liberti, N. Maculan (Eds.), Global Optimization: From Theory to Implementation, Springer, Berlin, 2006, pp. 405-414.
[14] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, Recent advances on the discretizable molecular distance geometry problem, European Journal of Operational Research 219 (2012) 698-706.
[15] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, The discretizable molecular distance geometry problem, Computational Optimization and Applications 52 (2012) 115-146.
[16] C. Lavor, L. Liberti, A. Mucherino, The interval BP algorithm for the discretizable molecular distance geometry problem with interval data, Journal of Global Optimization 56 (2013) 855-871.
[17] C. Lavor, L. Liberti, B. Donald, B. Worley, B. Bardiaux, T. Malliavin, M. Nilges, Minimal NMR distance information for rigidity of protein graphs, Discrete Applied Mathematics 256 (2019) 91-104.
[18] C. Lavor, R. Alves, Oriented conformal geometric algebra and the molecular distance geometry problem, Advances in Applied Clifford Algebra 29 (2019) 1-19.

---

Luiz Leduino de Salles Neto, C. Lavor and W. Lodwick

Information Sciences 559 (2021) 1-7

[19] C. Lavor, M. Souza, L. Mariano, L. Liberti, On the polinomiality of finding  $^8$ DMDGP re-orders, Discrete Applied Mathematics 267 (2019) 190-194.
[20] C. Lavor, Comments on distance geometry and data science, TOP 28 (2020) 340-345.
[21] L. Liberti, Distance geometry and data science, TOP 28 (2020) 271-339.
[22] L. Liberti, C. Lavor, N. Maculan, A branch-and-prune algorithm for the molecular distance geometry problem, International Transactions in Operational Research 15 (2008) 1-17.
[23] L. Liberti, C. Lavor, On a relationship between graph realizability and distance matrix completion, in: A. Migdalas, A. Sifaleras, K. Christos, J. Papathanasiou, E. Stiakakis (Eds.), Optimization Theory, Decision Making, and Operations Research Applications, Springer, 2013, pp. 39-48.
[24] L. Liberti, C. Lavor, N. Maculan, A. Mucherino, Euclidean distance geometry and applications, SIAM Review 56 (2014) 3-69.
[25] L. Liberti, C. Lavor, Six mathematical gems from the history of distance geometry, International Transactions in Operational Research 23 (2016) 897-920.
[26] K. Menger, Untersuchungen über allgemeine Metrik, Mathematische Annalen 100 (1928) 75-163.
[27] K. Menger, New foundation of Euclidean geometry, American Journal of Mathematics 53 (1931) 721-745.
[28] A. Mucherino, C. Lavor, L. Liberti, N. Maculan (Eds.), Distance Geometry: Theory, Methods, and Applications, Springer, New York, 2013.
[29] M. Sippl, H. Scheraga, Cayley-Menger coordinates, Proceedings of the National Academy of Sciences USA 83 (1986) 2283-2287.
[30] B. Worley, F. Delhommel, F. Cordier, T. Malliavin, B. Bardiaux, N. Wolff, M. Nilges, C. Lavor, L. Liberti, Tuning interval branch-and-prune for protein structure determination, Journal of Global Optimization 72 (2018) 109-127.