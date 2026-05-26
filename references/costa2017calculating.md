Information Sciences 415-416 (2017) 41-52

ELSEVIER

Contents lists available at ScienceDirect

Information Sciences

journal homepage: www.elsevier.com/locate/ins

CrossMark

# Calculating the possible conformations arising from uncertainty in the molecular distance geometry problem using constraint interval analysis

T.M. Costa a,\*, H. Bouwmeester b, W.A. Lodwick c, C. Lavor d

$^{a}$ Post-Doc and CNPq fellow: Department of Mathematics and Statistical Sciences, University of Colorado Denver, Denver, Colorado, USA
$^{b}$ Department of Mathematical and Computer Sciences, Metropolitan State University of Denver, Denver, Colorado, USA
$^{c}$ Department of Mathematical Statistical Science, University of Colorado Denver, Denver, Colorado, USA
$^{d}$ Department of Applied Mathematics, University of Campinas, Campinas, São Paulo, Brazil

# ARTICLE INFO

Article history:

Received 8 November 2016

Revised 2 June 2017

Accepted 16 June 2017

Available online 17 June 2017

Keywords:

Molecular distance geometry problems

Interval uncertainties

Constraint interval analysis

# ABSTRACT

The calculation of the 3D structure of a protein molecule is important because it is associated to its biological function. Nuclear Magnetic Resonance (NMR) experiments can provide distance information between atoms that are close enough in a given protein and the problem is how to use these distances to determine the protein structure. Using the chemistry of proteins and supposing all the distances are precise values, it is possible to define an atomic order  $\nu_{1},\dots,\nu_{n}$ , such that the distances related to the pairs  $\{\nu_{i - 3},\nu_i\}$ ,  $\{\nu_{i - 2},\nu_i\}$ ,  $\{\nu_{i - 1},\nu_i\}$  are available, and solve the problem iteratively using a combinatorial method, called Branch-and-Prune (BP). However, due to uncertainty in NMR data, the distances associated with pairs  $\{\nu_{i - 3},\nu_i\}$  may not be precise, which implies that there are many difficulties in applying the BP algorithm to this scenario. The use of standard interval arithmetic can be directly applied to the algorithm, but it is known that it generates overestimations. This paper proposes a new methodology to compute possible conformations on the presence of uncertainties arising from NMR distance measurements using a constraint interval analysis approach. Some numerical examples are presented.

© 2017 Elsevier Inc. All rights reserved.

# 1. Introduction

The three-dimensional structure of molecules, for example protein molecules, is related to its physical-chemical properties [22]. A typical problem in molecular conformation is one in which some, but not all, distances between atoms of a molecule are given and the problem is to calculate its 3D structure. When a nuclear magnetic resonance (NMR) machine is used to measure the distances between atoms, some distances are more accurate than others. The more accurate measured distances are considered as real numbers, whereas the less accurate are considered here to be intervals that contain the accurate distance. In practice, the measurements of the distances between atoms that are close to each other are known

http://dx.doi.org/10.1016/j.ins.2017.06.015

0020-0255/© 2017 Elsevier Inc. All rights reserved.

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

more precisely (real-valued), whereas the measurements of the distances between atoms that are farther apart have error or cannot be measured.

Our problem arises from the general Distance Geometry Problem (DGP) [12,21], which is to determine the coordinates of a set of points in a given Euclidean space $\mathbb{R}^K$, where some of the distances are known. The DGP is often stated using a graph theoretic representation, which is the way we will state the problem. The vertices of the graph corresponding to our problem are the atoms and the pairs of atoms for which we know the distances are the weighted edges of the graph (the distances being the weights). The definition of the DGP is as follows.

**Definition 1.** Given an integer $K &gt; 0$ and a simple undirected graph $G = (V, E, d)$ whose edges are weighted by a nonnegative function $d: E \to (0, \infty)$, find a function $x: V \to \mathbb{R}^K$ such that

$$
\forall \{u, v \} \in E, \| x (u) - x (v) \| = d (u, v), \tag{1}
$$

where $\| \cdot \|$ is the Euclidean norm.

We are interested in the structural representation (its conformation) of a molecule, given distance information between some of its atoms. Therefore, the space of interest is $\mathbb{R}^3$ ($K = 3$), and our problem is to determine the coordinates in $\mathbb{R}^3$ of the molecule so that the distances between the computed coordinates (the function $x: V \to \mathbb{R}^3$ of Definition 1) match the given distances of our data set. It is known that the set of solutions of a DGP can be empty, finite, or uncountably infinite (see [12]). Thus, to deal with the calculations associated with molecular structures, it is useful to obtain conditions under which the problem has a finite number of solutions, when the solution set is not empty [19]. We consider the Discretizable Molecular Distance Geometry Problem (DMDGP) [4,6], which is a particular class of DGP with an extra assumption: the existence of an indexing of the vertices of its associated graph satisfying some properties [5]. Formally, the DMDGP is defined as follows.

**Definition 2.** Consider a graph $G = (V, E, d)$ of a DGP for $K = 3$ with a vertex order on $V$, denoted by $v_{1}, v_{2}, \ldots, v_{n}$, such that

1. For $v_{1}, v_{2}, v_{3}$, the coordinates $x(v_{1}), x(v_{2}), x(v_{3}) \in \mathbb{R}^{3}$ are known and satisfy Eq. (1);
2. For $i = 4, \ldots, n$, the set $\{v_{i-3}, v_{i-2}, v_{i-1}, v_i\}$ is a clique;
3. For $i = 3, \ldots, n$

$$
d \left(v _ {i - 2}, v _ {i}\right) &lt; d \left(v _ {i - 2}, v _ {i - 1}\right) + d \left(v _ {i - 1}, v _ {i}\right).
$$

The DMDGP is to find a function $x:V\to \mathbb{R}^3$ such that

$$
\forall \left\{v _ {i}, v _ {j} \right\} \in E, \| x \left(v _ {i}\right) - x \left(v _ {j}\right) \| = d \left(v _ {i}, v _ {j}\right). \tag{2}
$$

In what follows, $x_{i}$ will denote $x(v_{i})$ and $d_{i,j}$ will denote $d(v_{i}, v_{j})$.

We can interpret the input data of a DMDGP in the following way. The set of vertices, $V$, is the set of atoms of a molecule, and $d_{i,j}$ is the distance between the atoms $v_{i}$ and $v_{j}$. Thus, uncovering molecular structures is to construct all functions $x:V\to \mathbb{R}^3$ that satisfy Eq. (2). If we have all the distances between the vertex pairs, a unique structure can be obtained [12]. Given an incomplete set of distances, the possible conformations may not be unique, and all of them from partial data is sought.

Most DMDGP literature treats the problem from a deterministic viewpoint, that is, the distances are represented by real numbers. However, for a particular protein molecule whose distances are measured by a NMR experiment, it is known that the distance data given by these machines have measurement errors. Consequently, in order to more realistically model such a problem, we need to consider inaccuracies in the distance data. We will assume that the distance data from an NMR machine are represented by intervals of real numbers where, for atoms that are close together, the interval has zero width. That is, we only know the bounds on the machine errors, not their distribution within the bounds. For accurate measurements, the lower/upper bounds will be the same, the width of the interval is zero. Also, we will assume that:

1. An NMR machine measures distances of an actual existent molecule, which is obvious, but an important fact;
2. The NMR experiment does not change the molecular structure during its measurement;
3. The interval distances that arise are due to measurement errors;
4. Choosing one particular value for the distance within any one of the interval distances as "the" distance between two atoms within the molecule causes all atoms to adjust themselves as a rigid motion. This means that the position of all atoms are $100\%$ correlated, arising because it is an actual molecule that is being subjected to the NMR machine, which is our first point.

The introduction of interval inaccuracies to the DMDGP generates a new problem, called the Interval DMDGP (iDMDGP) [12].

**Definition 3.** Given a DMDGP graph $G = (V, E, d)$ with a vertex order on $V$, denoted by $v_{1}, v_{2}, \ldots, v_{n}$, we say that $G$ represents an instance of the iDMDGP if their weights are now interval distances and there exists a subset $E' \subset E$ with weights given by real numbers such that, for $i = 4, \ldots, n$

$$
\{\{v _ {i - 3}, v _ {i - 2} \}, \{v _ {i - 3}, v _ {i - 1} \}, \{v _ {i - 2}, v _ {i - 1} \}, \{v _ {i - 2}, v _ {i} \}, \{v _ {i - 1}, v _ {i} \} \} \subset E ^ {\prime}.
$$

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

In the iDMDGP, we want to find a function $x: V \to \mathbb{R}^3$ that satisfies

$$
\forall \{v_i, v_j\} \in E, \, \underline{d}(v_i, v_j) \leq \|x(v_i) - x(v_j)\| \leq \overline{d}(v_i, v_j),
$$

where $[\underline{d}(v_i, v_j), \overline{d}(v_i, v_j)]$ is the interval distance related to the pair $\{v_i, v_j\} \in E$ (for $\{v_i, v_j\} \in E'$, $\underline{d}(v_i, v_j) = \overline{d}(v_i, v_j)$).

The definition of the set $E'$ is justified by the fact that, in protein structure calculations, all the bond lengths and covalent bond angles are considered as real numbers [22]. Other pair of atoms can also be included in $E'$, depending on the precision of their distances, which could be real numbers. That is, we could have $\underline{d} = \overline{d}$ for precise data.

Two areas of mathematics which represent interval errors are interval analysis and generalized uncertainty theory [14]. The latter may require more information about the way the errors are distributed within the interval, but we are assuming we do not have distribution data. So our approach is an interval analysis one. However, since interval analysis can be considered as a particular case of generalized uncertainty theory (see [16]), what we develop herein can be considered as a first step toward solving problems with more generalized uncertainty types. Interval analysis is a natural way to deal with problems whose data inaccuracies arise from measurement errors such as those generated by an NMR machine.

The main goal of this study is to use a newer approach to interval analysis, called constraint interval analysis, to control the propagation of uncertainty when solving the problem. To this end, we develop an interval approach, which guarantees that the interval obtained from the calculations of the coordinates,

$$
[x]: V \to \mathbb{R}^3, \, [x] = \{x \mid \underline{x} \leq x \leq \overline{x}\},
$$

contains possible solutions $x: V \to \mathbb{R}^3$ of the associated iDMDGP.

We use $\mathbb{R}^3$ to denote $[x_1] \times [x_2] \times [x_3]$, $\underline{x}_i \leq x_i \leq \overline{x}_i$, $i = 1, 2, 3$, $x_i \in \mathbb{R}$, the space of three dimensional intervals. The approach developed herein uses a symbolic representation of an interval, constraint intervals [13,14], and it is this symbolic representation that is propagated in the calculations, and not the interval itself. When numerical values of intervals are needed, the possible real-valued coordinates are approximated. It is clear that if we wish to compute numerical intervals on a computer that guarantees that all solutions lie within the computed interval, enclosure methods (see [14]) are needed that take into account numerical method errors and roundoff errors. Two systems which are used for enclosures are:

1. C-XSC, http://www2.math.uni-wuppertal.de/~xsc/xsc/cxsc.html;
2. INTLAB, http://www.ti3.tu-harburg.de/rump/intlab/.

Enclosure methods are beyond the scope of this work.

This article is arranged as follows. Section 2 presents the justification for the use of constraint interval analysis in the iDMDGP and introduces a short discussion of constraint interval arithmetic as it relates to the iDMDGP. Section 3 presents a method for controlling the uncertainty in the iDMDGP based on a constraint interval analysis approach. Section 4 contains some comparative numerical examples that illustrate our approach and show the efficiency of our method by comparing it to the arithmetic used by Lavor et al. [7]. Lastly, Section 5 presents the conclusion and subsequent steps to our line of research. We emphasize that what is presented in this article is the theoretical basis for controlling interval error.

## 2. The role of intervals and constraint interval analysis in the iDMDGP

NMR errors are given as intervals, so that it is clear that when we calculate the effects of interval uncertainties beginning at a particular atom, we propagate these interval uncertainties through subsequent arithmetic computations. The control of the propagation of error is the focus of this article.

## 2.1. Interval computations

It is known that calculations based on standard (naive) interval arithmetic in the presence of dependencies add errors that are not present in the problem (see [17]). For example, suppose we have the following real-valued expression,

$$
x + 5 - x. \tag{3}
$$

Using the standard Moore interval arithmetic to evaluate the expression,

$$
([x] + 5) - [x], \tag{4}
$$

for $[x] = [1, 2]$, we have,

$$
\begin{aligned}
([1, 2] + 5) - [1, 2] &amp;= [6, 7] + [-2, -1] \tag{5} \\
&amp;= [4, 6].
\end{aligned}
$$

However, if we represent the interval as

$$
[x] = [1, 2] = 1 + \lambda_x, \lambda_x \in [0, 1], \tag{6}
$$

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

then,

$$
\begin{array}{l}
\left(\left(1 + \lambda_{x}\right) + 5\right) - \left(1 + \lambda_{x}\right) = 6 + \lambda_{x} - 1 - \lambda_{x} \tag{7} \\
= 5.
\end{array}
$$

That is, standard (naive) interval arithmetic includes errors that are not in the uncertainty of the problem which will always occur in the presence of dependencies. For our example, there is no uncertainty in the solution even though there is uncertainty in the variable $x$. For the non-naive standard interval arithmetic approach to obtain the correct evaluation of (3) to be 5 (see [18]), the interval [1,2] must be partitioned in such a way that the widths (diameters in $\mathbb{R}^n$) of the partitions converge to zero where one takes the union of the results from standard interval arithmetic applied to each partition.

For example, let $[1,2] = [1,\frac{3}{2}]\cup [\frac{3}{2},2]$. For $[1,\frac{3}{2}]$,

$$
\begin{array}{l}
\left(\left[1, \frac{3}{2}\right] + 5\right) - \left[1, \frac{3}{2}\right] = \left[6, \frac{13}{2}\right] + \left[-\frac{3}{2}, -1\right] \\
= \left[\frac{9}{2}, \frac{11}{2}\right].
\end{array}
$$

On $[\frac{3}{2}, 2]$, we have,

$$
\begin{array}{l}
\left(\left[ \frac{3}{2}, 2 \right] + 5\right) - \left[ \frac{3}{2}, 2 \right] = \left[ \frac{13}{2}, 7 \right] + \left[-2, -\frac{3}{2}\right] \\
= \left[\frac{9}{2}, \frac{11}{2}\right].
\end{array}
$$

The solution is obtained by the union of the two partitions, giving us

$$
\left[\frac{9}{2}, \frac{11}{2}\right] \cup \left[\frac{9}{2}, \frac{11}{2}\right] = \left[\frac{9}{2}, \frac{11}{2}\right],
$$

which is an improvement to the solution [4, 6] obtained without partition. Note that to obtain the correct solution 5, the interval [1, 2] needs an union of an infinite number of partitions whose diameters (widths) converge to zero after which the infinite union needs to be taken. Of course, in practice, a finite number of partitions is taken to obtain an approximate solution. Naive or non-naive standard interval arithmetic with dependencies such as (4) will need to be partitioned to reduce standard interval arithmetic's inherent overestimation (or other methods such as nesting calculations). In nesting calculations, dependencies would need to be known a-priori. On the other hand, the representation (6) requires, in general, a global optimization method, as we shall see, though for this simple example it does not. However, no a-priori knowledge is needed for representation (6), the dependencies are immediately taken into account by the representation itself as happened when (6) was used in the calculation of the simple example.

## 2.2. Constraint interval arithmetic

The basic concepts associated with the usual/standard interval arithmetic and analysis can be found in [17,18]. The interest here is in newer ways of representing and analyzing intervals. An interval, for this study, is represented by a continuous non-decreasing real function $f(\lambda_x):[0,1] \longrightarrow [\underline{x},\overline{x}]$, such that

$$
f(0) = \underline{x} \quad \text{and} \quad f(1) = \overline{x},
$$

which we call a generalized constraint interval (see [16]). We limit ourselves to a particular function $f(\lambda_x):[0,1] \longrightarrow [\underline{x},\overline{x}]$, given by

$$
f(\lambda_{x}) = \underline{x} + w_{x} \lambda_{x}, \quad w_{x} = \overline{x} - \underline{x}, \quad 0 \leq \lambda_{x} \leq 1, \tag{8}
$$

which is called a constraint interval representation. Lodwick, in [13] (see also [14]), introduced this representation and used (8) to define an arithmetic, called constraint interval arithmetic (CIA), given by

$$
\begin{array}{l}
[\underline{x}, \overline{x}] * [\underline{y}, \overline{y}] = \left[ \min_{0 \leq \lambda_{x}, \lambda_{y} \leq 1} \left\{ (\underline{x} + w_{x} \lambda_{x}) * (\underline{y} + w_{y} \lambda_{y}) \right\}, \right. \\
\left. \max_{0 \leq \lambda_{x}, \lambda_{y} \leq 1} \left\{ (\underline{x} + w_{x} \lambda_{x}) * (\underline{y} + w_{y} \lambda_{y}) \right\} \right], \tag{9}
\end{array}
$$

where $* \in \{+, -, \cdot, \div\}$. Distinct from the standard interval arithmetic and analysis found in [17,18], constraint interval analysis preserves the dependence or independence of the arithmetic operations in an explicit way through the parameters $\lambda$ defined on $[0, 1]$. In particular, $[x] - [x] = [0, 0]$ and $[x] \div [x] = 1$ (whenever $0 \notin [x]$). Constraint interval arithmetic is the implementation of the united extension (see [17]), whereas the standard interval analysis is not without taking an extra step of partitioning into subintervals, taking the limit as the diameters of the subintervals converge to zero, and taking the union of the limit over all these convergent partitions.

---

A subset of constraint interval arithmetic, called single level constraint interval arithmetic (SLCIA), is defined as$$\left\lbrack \underset{0 \leq \lambda \leq 1}{\underline{x}} \right\rbrack*\left\lbrack \underset{0 \leq \lambda \leq 1}{\underline{y}} \right\rbrack\left(\left(\underset{0 \leq \lambda \leq 1}{\underline{x}} + w_{x}\lambda \right)*\left(\underset{0 \leq \lambda \leq 1}{\underline{y}} + w_{y}\lambda \right) \right), \underset{0 \leq \lambda \leq 1}{\max}\left{ \left(\underset{0 \leq \lambda \leq 1}{\underline{x}} + w_{x}\lambda \right)*\left(\underset{0 \leq \lambda \leq 1}{\underline{y}} + w_{y}\lambda \right) \right\rbrack.$$

The theory and algebraic properties of SLCIA can be found in [2]. Single-level interval arithmetic uses a single $\lambda \in \left\lbrack {0,1} \right\rbrack$ to parametrize all intervals and so captures the notion of complete dependence.

### Intervals applied to the iDMDGP

We mentioned in the Introduction that a molecule whose pair-wise distances are measured by an NMR machine from a particular single “snapshot” does not undergo change to its structure while being measured and when one assigns to a particular pair of atoms (a particular distance within the bounded measured distance represented by the given NMR interval error), all other distances between atoms “adjust” themselves so that a chosen position of one atom affects all other atoms as in a rigid movement. This means that the position between atoms in a molecule are completely dependent on any one coordinate displacement. The mathematical representation of complete dependency across various interval distances is via a representation such as (6), called single-level constraint interval arithmetic (SLCIA) (see [2]), which is a subset of constraint interval arithmetic (CIA) (see [13], [14]). We will use SLCIA as the representation of intervals in our analysis of the iDMDGP. If one were to use standard interval analysis on the intervals, which assumes completely independent positions within each interval distance, then the errors grow in a similar way as they did with (4) at each encounter with interval distance calculations. Moreover, complete independence of intervals does not model the problem where the NMR machine has measured the distances between atoms of an actual molecule.

The use of a partitioning strategy for the iDMDGP as outlined by Moore (see [17], [18]) was investigated by Cassioli et. al. [1]. Each iteration of the BP algorithm has to be partitioned implying that the Moore partitioning strategy is intractable. Our method is different. We propagate the parameterization of the interval and not the interval itself. Then, after the end of the error propagation via SLCIA, all the coordinates of the atoms of the molecule will be a function of a single bounded parameter.

Independence of atom placement within the given interval distance data comes after the parameterized coordinates are obtained. Each calculated coordinate position is a function of a single parameter whose value can vary between 0 and 1. Given a SLCIA parameterized coordinate, the parameter of each coordinate position can independently vary within its constraint interval, [0, 1]. Thus, there are two steps to our method. The first step is to calculate the uncertainty in the coordinates of the atoms as a function of the same single parameter variable via single-level interval arithmetic. The second step is to compute all positions given the uncertainty in the position of each atom as all possible variations of the parameter within [0, 1] independently.

###### Remark 4.

When standard interval arithmetic and analysis is used in computations as in (5), errors not associated with the existent uncertainties are introduced and propagated. In the context of the iDMDGP, this error that is not a part of the problem uncertainty, but in the way the uncertainty is computed, is illustrated at the end of Example 6.

## Controlling uncertainty in the iDMDGP

There is a method for DMDGP, called branch-and-prune (BP), for finding all solutions up to rotations and translations [9], [10]. The main step of the algorithm is to solve a quadratic system to get the two possible positions for atom $i, i>3$, in terms of the positions of atoms $i-3, i-2, i-1$ and the distances $d_{i-1,i}, d_{i-2,i}, d_{i-3,i}$. This recursive procedure defines a binary tree containing all possible positions for each vertex $v_{i}$ on the respective layer. Each DMDGP solution can be represented as a path from the root to a leaf node of the tree.

Because of the properties of the DMDGP order, quadratic systems may be replaced by matrix multiplications, and the computational results presented in [3] demonstrate that the second approach guarantees more numerical stability in the BP algorithm.

We will consider a molecule as a chain of $n$ atoms, indexed by $1,\ldots,n$, described by its internal coordinates, given by the bond lengths $d_{i-1,i}$ (the Euclidean distance between $x_{i-1}$ and $x_{i}$), for $i=2,\ldots,n$, the bond angles $\theta_{i-2,i}$ (the angle defined by the atoms $i-2,i-1,i$), for $i=3,\ldots,n$, and the torsion angles $\omega_{i-3,i}$ (the angle between the normals through the planes defined by the atoms $i-3,i-2,i-1$ and $i-2,i-1,i$), for $i=4,\ldots,n$. The values $d_{i-1,i}, \theta_{i-2,i}, \cos(\omega_{i-3,i})$ can be calculated using the distances between the atoms $i-3,i-2,i-1,i$, for $i=4,\ldots,n$. Based on the positions for $x_{1},x_{2},x_{3},\ldots,x_{i-1}$, the two possible values for $x_{i}=(x_{i1},x_{i2},x_{i3})^{T}\in\mathbb{R}^{3}$ (related to the two values for $\sin(\omega_{i-3,i})=\pm\sqrt{1-\cos^{2}(\omega_{i-3,i})}$), can be obtained using the following matrix multiplications [8]:$$\left\lbrack \begin{matrix} x_{i1} x_{i2} x_{i3} 1 \end{matrix} \right\rbrack = B_{1}B_{2}\ldots B_{i}\left\lbrack \begin{matrix} 0 0 0 1 \end{matrix} \right\rbrack,\ \forall i = 1,\ldots,n,$$

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

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

for  $i = 4,\ldots ,n$

We will extend this matrix approach to DMDGP with interval distances using single-level interval analysis.

## 3.1. SLCIA matrix method

The interval matrix approach generates the coordinate values in  $\mathbb{R}^3 = (\mathbb{R}\times \mathbb{R}\times \mathbb{R})$  of the set atoms of a given molecule, where the uncertainties are represented as interval distances  $[d_{i - 3,i}]$ ,  $i = 4,\dots ,n$ .

Interval uncertainties of the matrix method appear in the torsion angles  $\omega_{i}$ , for each  $i = 4,\dots ,n$ . In this case, there are interval sines and cosines,  $\sin^I (\omega_i)$  and  $\cos^I (\omega_i)$ , that are input data to the problem. The formula used to calculate  $\cos^I (\omega_i)$  is an interval extension of the method of Lavor et al. [8] and is calculated by replacing  $d_{i - 3,i}$  by the interval  $[d_{i - 3,i}] = [\underline{d}_{i - 3,i},\overline{d}_{i - 3,i}]$ , which results in

$$
\cos^ {I} \left(\omega_ {i}\right) = \frac {2 \left(d _ {i - 2 , i - 1}\right) ^ {2} \left(\left(d _ {i - 3 , i - 2}\right) ^ {2} + \left(d _ {i - 2 , i}\right) ^ {2} - \left(\left[ d _ {i - 3 , i} \right]\right) ^ {2}\right) - \left(d _ {i - 3 , i - 2 , i - 1} ^ {(i)}\right) \left(d _ {i - 2 , i - 1 , i} ^ {(i)}\right)}{\sqrt {4 \left(d _ {i - 3 , i - 2}\right) ^ {2} \left(d _ {i - 2 , i - 1}\right) ^ {2} - \left(d _ {i - 3 , i - 2 , i - 1} ^ {(i)}\right) ^ {2}} \sqrt {4 \left(d _ {i - 2 , i - 1}\right) ^ {2} \left(d _ {i - 2 , i}\right) ^ {2} - \left(d _ {i - 2 , i - 1 , i} ^ {(i)}\right) ^ {2}}}, \tag {12}
$$

where

$$
d _ {i - 3, i - 2, i - 1} ^ {(i)} = d _ {i - 3, i - 2} ^ {2} + d _ {i - 2, i - 1} ^ {2} - d _ {i - 3, i - 1} ^ {2},
$$

$$
d _ {i - 2, i - 1, i} ^ {(i)} = d _ {i - 2, i - 1} ^ {2} + d _ {i - 2, i} ^ {2} - d _ {i - 1, i} ^ {2}.
$$

An example of a consistency check that needs to be done as the coordinates are computed comes as a result of computing the interval cosine of the torsion angle. Since (12) expresses the uncertainty contained in  $[d_{i-3,i}]$ , which is given by the error associated with an NMR machine, we need to check if this uncertainty is an overestimation that may have been caused by the machine misreading the distance. That is, we need to verify that all the values within the interval provided are in fact possible. If a given uncertainty in the distance is an overestimation, then it is possible to obtain  $\cos^I(\omega_i) = [a,b]$ , with  $a &lt; -1$  and/or  $b &gt; 1$ . In this case, the interval  $\cos^I(\omega_i) = [a,b]$  not only represents the possible values of the cosine of the torsion angle, but it also represents impossible values. Thus, we constrain  $\cos^I(\omega_i) \subseteq [-1,1]$ . If  $\cos^I(\omega_i)$  is overestimated, we can reduce the width of  $[d_{i-3,i}]$  so that  $[a,b] \subseteq [-1,1]$ .

The single-level constraint interval expression of (12) is

$$
\cos \left(\omega_ {i} (\lambda)\right) = \frac {2 \left(d _ {i - 2 , i - 1}\right) ^ {2} \left(\left(d _ {i - 3 , i - 2}\right) ^ {2} + \left(d _ {i - 2 , i}\right) ^ {2} - \left(\boldsymbol {d} _ {i - 3 , i} (\boldsymbol {\lambda})\right) ^ {2}\right) - \left(d _ {i - 3 , i - 2 , i - 1} ^ {(i)}\right) \left(d _ {i - 2 , i - 1 , i} ^ {(i)}\right)}{\sqrt {4 \left(d _ {i - 3 , i - 2}\right) ^ {2} \left(d _ {i - 2 , i - 1}\right) ^ {2} - \left(d _ {i - 3 , i - 2 , i - 1} ^ {(i)}\right) ^ {2}} \sqrt {4 \left(d _ {i - 2 , i - 1}\right) ^ {2} \left(d _ {i - 2 , i}\right) ^ {2} - \left(d _ {i - 2 , i - 1 , i} ^ {(i)}\right) ^ {2}}}, \tag {13}
$$

where  $d_{i - 3,i}(\lambda) = \underline{d}_{i - 3,i} + \lambda (\overline{d}_{i - 3,i} - \underline{d}_{i - 3,i})$ $0\leq \lambda \leq 1$  , implying

$$
\sin \left(\omega_ {i} (\lambda)\right) = \pm \sqrt {1 - \left(\cos \left(\omega_ {i} (\lambda)\right)\right) ^ {2}}. \tag {14}
$$

The following assumptions are made for the interval matrix approach for a given set of atoms  $V = \{v_{1}, v_{2}, \dots, v_{n}\}$  of a molecule:

1. The values of the distances  $d_{i-3,i-2}$ ,  $d_{i-3,i-1}$ ,  $d_{i-2,i-1}$ ,  $d_{i-2,i}$ ,  $d_{i-1,i}$  are real positive values and the distances  $[d_{i-3,i}]$ ,  $0 \leq \underline{d}_{i-3,i} &lt; \overline{d}_{i-3,i}$ , are intervals for all  $i = 4, \ldots, n$ ;
2. For each  $i = 3,\dots ,n$ , the angle defined by  $v_{i - 2},v_{i - 1},v_{i}$ , denoted by  $\theta_{i}$ , is in  $(0,\pi)$  and can be obtained accurately since these angles are computed by the Law of Cosines using precise distance data;

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

3. For each $i = 3, \ldots, n$, $\cos^i(\omega_i) \subseteq [-1, 1]$.

Using the matrices (11), the first three positions are given by

$$
x _ {1} = \left(x _ {1 1} x _ {1 2} x _ {1 3}\right) ^ {T} = \left[ \begin{array}{c} 0 \\ 0 \\ 0 \end{array} \right],
$$

$$
x _ {2} = \left(x _ {2 1} x _ {2 2} x _ {2 3}\right) ^ {T} = \left[ \begin{array}{c} - d _ {1, 2} \\ 0 \\ 0 \end{array} \right],
$$

$$
x _ {3} = \left(x _ {3 1} x _ {3 2} x _ {3 3}\right) ^ {T} = \left[ \begin{array}{c} - d _ {1, 2} + d _ {2, 3} \cos (\theta_ {3}) \\ d _ {2, 3} \sin (\theta_ {3}) \\ 0 \end{array} \right].
$$

For $i = 4,\dots ,n$ we have

$$
\left[ \begin{array}{c} x _ {i 1} (\lambda) \\ x _ {i 2} (\lambda) \\ x _ {i 3} (\lambda) \\ 1 \end{array} \right] = B _ {1} B _ {2} B _ {3} \dots B _ {i} (\lambda) \left[ \begin{array}{c} 0 \\ 0 \\ 0 \\ 1 \end{array} \right], \tag {15}
$$

where

$$
B _ {i} (\lambda) = \left[ \begin{array}{c c c c} - \cos (\theta_ {i}) &amp; - \sin (\theta_ {i}) &amp; 0 &amp; - d _ {i - 1, i} \cos (\theta_ {i}) \\ \sin (\theta_ {i}) \cos (\omega_ {i} (\lambda)) &amp; - \cos \theta_ {i} \cos (\omega_ {i} (\lambda)) &amp; - \sin (\omega_ {i} (\lambda)) &amp; d _ {i - 1, i} \sin (\theta_ {i}) \cos (\omega_ {i} (\lambda)) \\ \sin (\theta_ {i}) \sin (\omega_ {i} (\lambda)) &amp; - \cos (\theta_ {i}) \sin (\omega_ {i} (\lambda)) &amp; \cos (\omega_ {i} (\lambda)) &amp; d _ {i - 1, i} \sin (\theta_ {i}) \sin (\omega_ {i} (\lambda)) \\ 0 &amp; 0 &amp; 0 &amp; 1 \end{array} \right].
$$

That is, the first three entries of the last column of $B_{1}B_{2}B_{3}\ldots B_{i}(\lambda)$ are our coordinates.

Fixing a path in the BP tree, each coordinate is minimized and maximized with respect to the parameter $\lambda$ separately, coordinate by coordinate (the last coordinate, fixed at 1, is omitted):

$$
\left[ x _ {i} \right] = \left(\left[ x _ {i 1} \right], \left[ x _ {i 2} \right], \left[ x _ {i 3} \right]\right) ^ {T} = \tag {16}
$$

$$
\left(\left[ \min  _ {\lambda \in [ 0. 1 ]} x _ {i 1} (\lambda), \max  _ {\lambda \in [ 0. 1 ]} x _ {i 1} (\lambda) \right], \left[ \min  _ {\lambda \in [ 0. 1 ]} x _ {i 2} (\lambda), \max  _ {\lambda \in [ 0. 1 ]} x _ {i 2} (\lambda) \right], \left[ \min  _ {\lambda \in [ 0. 1 ]} x _ {i 3} (\lambda), \max  _ {\lambda \in [ 0. 1 ]} x _ {i 3} (\lambda) \right]\right).
$$

Note that we do not propagate the interval vector $[x_i]$ in the computations forward, only its symbolic single-level representation. At any point after the symbolic computation, we can compute the interval value of $x_{i}(\lambda)$ by taking the min/max over $\lambda \in [0,1]$, using (16).

## 4. Numerical examples

This section presents two numerical examples in which the proposed approach in this paper and the methodology iBP are compared. Example 5 computes the interval coordinates by means of the single-level constraint interval arithmetic and compares this with the iBP graphically. Example 6 computes the interval coordinates by means of both the single-level constraint interval arithmetic and standard interval arithmetic, and the comparison is numerically done through the width of each resulting interval coordinate.

Example 5. Consider an iDMDGP instance, with the following distances:

$$
d _ {1, 2} = 1; d _ {1, 3} = \sqrt {2}; d _ {1, 4} = [ 1. 6 8 2 1, 1. 7 8 2 1 ] = 1. 6 8 2 1 + (0. 1) \lambda ,
$$

$$
d _ {2, 3} = 1; d _ {2, 4} = \sqrt {2}; d _ {2, 5} = [ 1. 3 6 4 2, 1. 4 6 4 2 ] = 1. 3 6 4 2 + (0. 1) \lambda ,
$$

$$
d _ {3, 4} = 1; d _ {3, 5} = \sqrt {3}.
$$

From (12), it follows that

$$
\cos \left(\omega_ {4} (\lambda)\right) = \frac {2 (1) ^ {2} \left((1) ^ {2} + (\sqrt {2}) ^ {2} - (1 . 6 8 2 1 + (0 . 1) \lambda) ^ {2}\right) - (0) (2)}{\sqrt {4 (1) ^ {2} (1) ^ {2} - (0) ^ {2}} \sqrt {4 (1) ^ {2} (\sqrt {2}) ^ {2} - (2) ^ {2}}} \tag {17}
$$

$$
= - 0. 0 0 5 \lambda^ {2} - 0. 1 6 8 2 1 \lambda + 8. 5 2 7 0 \times 1 0 ^ {- 2} \tag {18}
$$

and

$$
\cos \left(\omega_ {5} (\lambda)\right) = \frac {2 (1) ^ {2} \left((1) ^ {2} + (\sqrt {3}) ^ {2} - (1 . 3 6 4 2 + (0 . 1) \lambda) ^ {2}\right) - (0) (2)}{\sqrt {4 (1) ^ {2} (1) ^ {2} - (0) ^ {2}} \sqrt {4 (1) ^ {2} (\sqrt {3}) ^ {2} - (2) ^ {2}}} \tag {19}
$$

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

$$
= -0.0025\sqrt{2}\lambda^{2} - 0.06821\sqrt{2}\lambda + 0.53474\sqrt{2}. \tag{20}
$$

Thus, two values for $\sin(\omega_4(\lambda))$ and $\sin(\omega_5(\lambda))$ are possible:

$$
\sin(\omega_4(\lambda)) = \pm\sqrt{1 - (-0.005\lambda^{2} - 0.16821\lambda + 8.5270 \times 10^{-2})^2} \tag{21}
$$

and

$$
\sin(\omega_5(\lambda)) = \pm\sqrt{1 - (-0.0025\sqrt{2}\lambda^{2} - 0.06821\sqrt{2}\lambda + 0.53474\sqrt{2})^2}. \tag{22}
$$

Doing the calculations, we obtain all the possible solutions to the problem, given by

$$
[x_4]^0 = \left[ \min_{0 \leq \lambda \leq 1} x_4^0(\lambda); \max_{0 \leq \lambda \leq 1} x_4^0(\lambda) \right] = \begin{bmatrix} [-1.0879, -0.91473] \\ [1, 1] \\ [-1, -0.99613] \end{bmatrix},
$$

$$
[x_4]^1 = \left[ \min_{0 \leq \lambda \leq 1} x_4^1(\lambda); \max_{0 \leq \lambda \leq 1} x_4^1(\lambda) \right] = \begin{bmatrix} [-1.0879, -0.91473] \\ [1, 1] \\ [0.99613, 1] \end{bmatrix},
$$

$$
[x_5]^10 = \left[ \min_{0 \leq \lambda \leq 1} x_5^{10}(\lambda); \max_{0 \leq \lambda \leq 1} x_5^{10}(\lambda) \right] = \begin{bmatrix} [-0.024967, 0.0072587] \\ [-0.06948, 0.07194] \\ [0.91746, 1.0900] \end{bmatrix},
$$

$$
[x_5]^11 = \left[ \min_{0 \leq \lambda \leq 1} x_5^{11}(\lambda); \max_{0 \leq \lambda \leq 1} x_5^{11}(\lambda) \right] = \begin{bmatrix} [-2.1509, -1.8367] \\ [-0.06948, 0.07194] \\ [0.90229, 1.0753] \end{bmatrix},
$$

$$
[x_5]^00 = \left[ \min_{0 \leq \lambda \leq 1} x_5^{00}(\lambda); \max_{0 \leq \lambda \leq 1} x_5^{00}(\lambda) \right] = \begin{bmatrix} [-2.1509, -1.8367] \\ [-0.06948, 0.07194] \\ [-1.0753, -0.90229] \end{bmatrix},
$$

$$
[x_5]^01 = \left[ \min_{0 \leq \lambda \leq 1} x_5^{01}(\lambda); \max_{0 \leq \lambda \leq 1} x_5^{01}(\lambda) \right] = \begin{bmatrix} [-0.024967, 0.0072587] \\ [-0.06948, 0.07194] \\ [-1.0900, -0.91746] \end{bmatrix}.
$$

Example 6. Consider now another $i$ DMDGP instance, with the following distance data:

$$
d_{1,2} = 1.5; \quad d_{1,3} = 2.441311; \quad [d_{1,4}] = 3.529049 + 0.2\lambda;
$$

$$
d_{2,3} = 1.486607; \quad d_{2,4} = 2.473863; \quad [d_{2,5}] = 3.879950 + 0.2\lambda;
$$

$$
d_{3,4} = 1.5; \quad d_{3,5} = 2.722132; \quad [d_{3,6}] = 2.631300 + 0.2\lambda;
$$

$$
d_{4,5} = 1.8; \quad d_{4,6} = 2.519921; \quad [d_{4,7}] = 3.551027 + 0.2\lambda;
$$

$$
d_{5,6} = 1.658312; \quad d_{5,7} = 2.3; \quad d_{6,7} = 1.489966.
$$

For this example, we have that $\cos^i(\omega_5) = [-1.0032, -0.66034]$. Imposing $\cos^i(\omega_5) = [-1, -0.66034]$ and using this result in (12), we obtain a reduced interval $[d_{2,5}] = 3.87995 + (0.1982)\lambda$. Applying the single level interval matrix multiplication approach, we obtain the following interval coordinates:

$$
[x_4]^0 = \begin{bmatrix} [-3.345269, -2.861396] \\ [1.512404, 1.685216] \\ [-1.194424, -0.653922] \end{bmatrix}; \quad [x_4]^1 = \begin{bmatrix} [-3.345269, -2.861396] \\ [1.512404, 1.685216] \\ [0.653922, 1.194424] \end{bmatrix};
$$

$$
[x_5]^00 = \begin{bmatrix} [-4.508006, -3.468350] \\ [2.456851, 3.3974009] \\ [-1.190682, 0.363596] \end{bmatrix}; \quad [x_5]^11 = \begin{bmatrix} [-4.508006, -3.468350] \\ [2.456851, 3.391340] \\ [-0.363596, 1.190682] \end{bmatrix};
$$

$$
[x_5]^01 = \begin{bmatrix} [-3.451132, -3.018725] \\ [2.853783, 3.439784] \\ [-1.856095, -1.166811] \end{bmatrix}; \quad [x_5]^10 = \begin{bmatrix} [-3.451132, -3.018725] \\ [2.853783, 3.439784] \\ [1.166811, 1.856095] \end{bmatrix};
$$

$$
[x_6]^000 = \begin{bmatrix} [-3.477389, -2.030772] \\ [3.772464, 3.993819] \\ [-0.545199, 0.452728] \end{bmatrix}; \quad [x_6]^{111} = \begin{bmatrix} [-3.477389, -2.030772] \\ [3.772464, 3.993819] \\ [-0.452728, 0.545199] \end{bmatrix};
$$

$$
[x_6]^{001} = \begin{bmatrix} [-4.414031, -3.643398] \\ [1.462924, 3.398071] \\ [0.480271, 1.689132] \end{bmatrix}; \quad [x_6]^{110} = \begin{bmatrix} [-4.414031, -3.643398] \\ [1.462924, 3.398071] \\ [-1.689132, -0.480271] \end{bmatrix};
$$

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

$$
\begin{array}{l}
\left[ x _ {6} \right] ^ {0 1 1} = \left[ \begin{array}{c} \{- 3. 6 4 3 3 9 8, - 2. 8 0 6 0 3 6 \} \\ [ 3. 3 9 8 0 7 1, 4. 0 1 5 4 9 8 ] \\ \{- 0. 7 4 5 0 7 0, 0. 4 8 0 2 7 1 \} \end{array} \right]; \left[ x _ {6} \right] ^ {1 0 0} = \left[ \begin{array}{c} \{- 3. 6 4 3 3 9 8, - 2. 8 0 6 0 3 6 \} \\ [ 3. 3 9 8 0 7 1, 4. 0 1 5 4 9 8 ] \\ \{- 0. 4 8 0 2 7 1, 0. 7 4 5 0 7 0 \} \end{array} \right]; \\
\left[ x _ {6} \right] ^ {0 1 0} = \left[ \begin{array}{l} \{- 2. 0 7 2 1 0 9, - 1. 3 6 4 7 9 5 \} \\ [ 2. 2 4 8 3 7 2, 3. 9 7 4 0 0 9 ] \\ \{- 2. 7 0 0 2 4 0, - 0. 5 4 5 1 9 9 \} \end{array} \right]; \left[ x _ {6} \right] ^ {1 0 1} = \left[ \begin{array}{l} \{- 2. 0 7 2 1 0 9, - 1. 3 6 4 7 9 5 \} \\ [ 2. 2 4 8 3 7 2, 3. 9 7 4 0 0 9 ] \\ [ 0. 5 4 5 1 9 9, 2. 7 0 0 2 4 0 \} \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 0 0 0} = \left[ \begin{array}{l} \{- 4. 5 7 4 0 5 1, - 1. 9 8 1 6 8 9 \} \\ [ 4. 6 8 0 7 3 1, 5. 2 4 0 1 7 9 ] \\ \{- 1. 5 3 5 0 6 0, 0. 9 4 6 8 2 6 \} \end{array} \right]; \left[ x _ {7} \right] ^ {1 1 1 1} = \left[ \begin{array}{l} \{- 4. 5 7 4 0 5 1, - 1. 9 8 1 6 8 9 \} \\ [ 4. 6 8 0 7 3 1, 5. 2 4 0 1 7 9 ] \\ \{- 0. 9 4 6 8 2 6, 1. 5 3 5 0 6 0 \} \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 0 0 1} = \left[ \begin{array}{l} \{- 4. 2 3 9 1 5 6, - 2. 7 3 1 2 4 2 \} \\ [ 4. 3 4 7 5 5 9, 4. 8 1 9 9 0 9 ] \\ [ 0. 4 5 4 2 6 9, 1. 5 8 0 4 4 1 ] \end{array} \right]; \left[ x _ {7} \right] ^ {1 1 1 0} = \left[ \begin{array}{l} \{- 4. 2 3 9 1 5 6, - 2. 7 3 1 2 4 2 \} \\ [ 4. 3 4 7 5 5 9, 4. 8 1 9 9 0 9 ] \\ \{- 1. 5 8 0 4 4 1, - 0. 4 5 4 2 6 9 \} \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 1 1 1} = \left[ \begin{array}{l} \{- 5. 0 1 0 8 6 7, - 3. 2 8 5 0 8 8 \} \\ [ 4. 0 3 0 3 2 3, 5. 2 5 0 7 7 6 ] \\ \{- 1. 6 2 3 3 9 9, 0. 3 9 1 1 9 5 \} \end{array} \right]; \left[ x _ {7} \right] ^ {1 0 0 0} = \left[ \begin{array}{l} \{- 5. 0 1 0 8 6 7, - 3. 2 8 5 0 8 8 \} \\ [ 4. 0 3 0 3 2 3, 5. 2 5 0 7 7 6 ] \\ \{- 0. 3 9 1 1 9 5, 1. 6 2 3 3 9 9 \} \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 0 1 1} = \left[ \begin{array}{l} \{- 5. 8 6 9 9 9 2, - 5. 0 1 0 8 6 7 \} \\ [ 1. 9 9 0 8 5 4, 4. 0 3 0 3 2 3 ] \\ [ 0. 3 9 1 1 9 5, 2. 3 8 6 8 4 1 ] \end{array} \right]; \left[ x _ {7} \right] ^ {1 1 0 0} = \left[ \begin{array}{l} \{- 5. 8 6 9 9 9 2, - 5. 0 1 0 8 6 7 \} \\ [ 1. 9 9 0 8 5 4, 4. 0 3 0 3 2 3 ] \\ \{- 2. 3 8 6 8 4 1, - 0. 3 9 1 1 9 5 \} \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 0 1 0} = \left[ \begin{array}{l} \{- 4. 9 1 1 7 1 3, - 2. 9 9 4 0 2 6 \} \\ [ 2. 4 3 3 6 3 4, 4. 7 2 6 0 5 8 ] \\ [ 0. 6 2 1 3 7 4, 2. 6 2 7 7 6 7 ] \end{array} \right]; \left[ x _ {7} \right] ^ {1 1 0 1} = \left[ \begin{array}{l} \{- 4. 9 1 1 7 1 3, - 2. 9 9 4 0 2 6 \} \\ [ 2. 4 3 3 6 3 4, 4. 7 2 6 0 5 8 ] \\ \{- 2. 6 2 7 7 6 7, - 0. 6 2 1 3 7 4 ] \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 1 0 0} = \left[ \begin{array}{l} \{- 2. 3 4 7 3 6 4, - 1. 3 5 3 6 1 5 \} \\ [ 3. 1 1 5 9 3 8, 5. 1 1 2 1 7 1 ] \\ \{- 3. 9 0 0 5 4 7, - 1. 5 3 5 0 6 0 ] \end{array} \right]; \left[ x _ {7} \right] ^ {1 0 1 1} = \left[ \begin{array}{l} \{- 2. 3 4 7 3 6 4, - 1. 3 5 3 6 1 5 \} \\ [ 3. 1 1 5 9 3 8, 5. 1 1 2 1 7 1 ] \\ \{- 3. 9 0 0 5 4 7, - 1. 5 3 5 0 6 0 ] \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 1 0 1} = \left[ \begin{array}{l} \{- 2. 7 3 1 2 4 2, - 1. 1 9 1 7 6 9 \} \\ [ 3. 5 2 1 9 3 8, 4. 9 6 1 2 6 1 ] \\ \{- 3. 4 5 4 0 7 6, 0. 4 5 4 2 6 9 ] \end{array} \right]; \left[ x _ {7} \right] ^ {1 0 1 0} = \left[ \begin{array}{l} \{- 2. 7 3 1 2 4 2, - 1. 1 9 1 7 6 9 \} \\ [ 3. 5 2 1 9 3 8, 4. 9 6 1 2 6 1 ] \\ \{- 0. 4 5 4 2 6 9, 3. 4 5 4 0 7 6 ] \end{array} \right]; \\
\left[ x _ {7} \right] ^ {0 1 1 0} = \left[ \begin{array}{l} \{- 2. 9 9 4 0 2 6, - 1. 9 9 0 6 0 5 \} \\ [ 4. 7 2 6 0 5 8, 5. 0 6 9 9 4 2 ] \\ \{- 1. 7 9 7 9 0 2, 0. 6 2 1 3 7 4 ] \end{array} \right]; \left[ x _ {7} \right] ^ {1 0 0 1} = \left[ \begin{array}{l} \{- 2. 9 9 4 0 2 6, - 1. 9 9 0 6 0 5 \} \\ [ 4. 7 2 6 0 5 8, 5. 0 6 9 9 4 2 ] \\ \{- 0. 6 2 1 3 7 4, 1. 7 9 7 9 0 2 ] \end{array} \right].
\end{array}
$$

Supposing that the molecule given to the NMR machine has the positions

$$
x _ {1} = (0, 0, 0) ^ {T}; x _ {2} = (- 1. 5, 0, 0) ^ {T}; x _ {3} = (- 2, 1. 4, 0) ^ {T}; x _ {4} = (- 3. 1, 1. 6, - 1) ^ {T};
$$

$$
x _ {5} = (- 4. 3, 2. 8, - 0. 4) ^ {T}; x _ {6} = (- 4. 4, 2. 1, 1. 1) ^ {T}; x _ {7} = (- 4. 6, 3. 4, 1. 8) ^ {T},
$$

the corresponding iDMDGP solution is the one given by  $[x_4]^0, [x_5]^{00}, [x_6]^{001}, [x_7]^{0010}$ .

Although this example deals with a small number of atoms, it shows the advantage and difference between SLCIA and standard interval arithmetic. More precisely, to represent the input data of this iDMDGP instance using standard interval arithmetic, we have

$$
d _ {1. 2} = 1. 5; d _ {1. 3} = 2. 4 4 1 3 1 1; [ d _ {1. 4} ] = [ 3. 5 2 9 0 4 9, 3. 7 2 9 0 4 9];
$$

$$
d _ {2. 3} = 1. 4 8 6 6 0 7; d _ {2. 4} = 2. 4 7 3 8 6 3; [ d _ {2. 5} ] = [ 3. 8 7 9 9 5 0, 4. 0 7 9 9 5];
$$

$$
d _ {3. 4} = 1. 5; d _ {3. 5} = 2. 7 2 2 1 3 2; [ d _ {3. 6} ] = [ 2. 6 3 1 3 0, 2. 8 3 1 3 0];
$$

$$
d _ {4. 5} = 1. 8 d _ {4. 6} = 2. 5 1 9 9 2 1; [ d _ {4. 7} ] = [ 3. 5 5 1 0 2 7, 3. 7 5 1 0 2 7];
$$

$$
d _ {5. 6} = 1. 6 5 8 3 1 2; d _ {5. 7} = 2. 3; d _ {6. 7} = 1. 4 8 9 9 6 6.
$$

By using (12), we have  $\cos^l (\omega_4) = [-0.884312, -0.515256]$  and  $\sin^l (\omega_4) = \pm [0.466896, 0.857036]$ . Considering  $\sin^l (\omega_4) = [0.466896, 0.857036]$ , we have

$$
\left[ \tilde {x} _ {4} \right] ^ {0} = \left[ \begin{array}{l} \left[ \tilde {x} _ {4 1} \right] ^ {0} \\ \left[ \tilde {x} _ {4 2} \right] ^ {0} \\ \left[ \tilde {x} _ {4 3} \right] ^ {0} \\ 1 \end{array} \right] = B _ {1} B _ {2} B _ {3} B _ {4} ^ {b} \left[ \begin{array}{l} 0 \\ 0 \\ 0 \\ 1 \end{array} \right],
$$

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

where

$$
\begin{array}{l}
B _ {4} ^ {k _ {s}} = \left[ \begin{array}{c c c c}
- \cos (\theta_ {4}) &amp; - \sin (\theta_ {4}) &amp; 0 &amp; - d _ {4} \cos (\theta_ {4}) \\
\sin (\theta_ {4}) \cos^ {t} (\omega_ {4}) &amp; - \cos (\theta_ {4}) \cos^ {t} (\omega_ {4}) &amp; - \sin^ {t} (\omega_ {4}) &amp; d _ {4} \sin (\theta_ {4}) \cos^ {t} (\omega_ {4}) \\
\sin (\theta_ {4}) \sin^ {t} (\omega_ {4}) &amp; - \cos (\theta_ {4}) \sin^ {t} (\omega_ {4}) &amp; \cos^ {t} (\omega_ {4}) &amp; d _ {4} \sin (\theta_ {4}) \sin^ {t} (\omega_ {4}) \\
0 &amp; 0 &amp; 0 &amp; 1
\end{array} \right] \\
= \left[ \begin{array}{c c c c}
0.372212 &amp; -0.928148 &amp; 0 &amp; 0.558318 \\
[-0.884312, -0.478234] &amp; [-0.329152, -191784] &amp; [-0.857036, -0.466896] &amp; [-1.231159, -0.717351] \\
0.433349, 0.795456] &amp; [0.173784, 0.318999] &amp; [-0.884312, -0.515256] &amp; [0.650023, 1.193184] \\
0 &amp; 0 &amp; 0 &amp; 1
\end{array} \right],
\end{array}
$$

which implies that

$$
[ \tilde {x} _ {4} ] ^ {0} = \left[ \begin{array}{l}
[ \tilde {x} _ {4 1} ] ^ {0} \\
[ \tilde {x} _ {4 2} ] ^ {0} \\
[ \tilde {x} _ {4 3} ] ^ {0}
\end{array} \right] = \left[ \begin{array}{l}
[-3.347215, -2.86341] \\
[1.511709, 1.684521] \\
[-1.193184, -0.650023]
\end{array} \right].
$$

Continuing in this manner, the associated iDMDGP solution that contains the original molecular structure is given by

$$
[ \tilde {x} _ {5} ] ^ {0 0} = \left[ \begin{array}{c}
[-5.163071, -3.234144] \\
[2.237908, 3.462696] \\
[-1.509219, 0.726549]
\end{array} \right],
$$

$$
[ \tilde {x} _ {6} ] ^ {0 0 1} = \left[ \begin{array}{c}
[-6.087774, -2.338009] \\
[0.826793, 3.828476] \\
[-1.212206, 3.456812]
\end{array} \right],
$$

$$
[ \tilde {x} _ {7} ] ^ {0 0 3 0} = \left[ \begin{array}{c}
[-8.007555, -1.031493] \\
[1.160911, 6.470073] \\
[-1.447948, 6.346757]
\end{array} \right].
$$

The width of  $[x_4]^0$  is a bit smaller than the width of  $[\tilde{x}_4]^0$ . However, for the next atoms, the difference between the interval coordinates given by SLCIA and Moore's arithmetic are quite different. If we look at  $x_5$ , we have

$$
[ x _ {5} ] ^ {0 0} = \left[ \begin{array}{c}
[-4.508003, -3.468348] \\
[2.456854, 3.391343] \\
[-1.190683, 0.363595]
\end{array} \right], [ \tilde {x} _ {5} ] ^ {0 0} = \left[ \begin{array}{c}
[-5.163071, -3.234144] \\
[2.237908, 3.462696] \\
[-1.509219, 0.726549]
\end{array} \right],
$$

where

$$
\operatorname {width} [ x _ {5} ] ^ {0 0} = \left[ \begin{array}{c}
1.0397 \\
0.93449 \\
1.5543
\end{array} \right], \quad \operatorname {width} [ \tilde {x} _ {5} ] ^ {0 0} = \left[ \begin{array}{c}
1.9289 \\
1.2248 \\
2.2358
\end{array} \right].
$$

Actually, the increasing difference in standard interval arithmetic is proportional to the number of atoms, whereas the widths of SLCIA is related to the actual error in the data.

The visualization of Example 5, Fig. 1, illustrates what the SLCIA and the standard interval arithmetic uncertainty propagation looks like. The SLCIA is a subset of CIA which in turn is a subset of the standard interval arithmetic (see [2]). Here subset means that when calculating with intervals, the results are subsets, that is, the interval result obtained by SLCIA is a subset of CIA. The interval result of CIA is a subset of the interval result using the standard interval arithmetic. Moreover, when there are no repetitions of variables in a series of calculations (that is, each variable in a series of calculations is independent of any other), CIA and Moore arithmetic will produce the same resultant interval (see [14,15]). Moreover, since both standard interval arithmetic and SLCIA are global min/max, they provide outer boxes that contain the true solutions as can be seen in Fig. 1.

The outer red boxes are obtained using the standard Moore interval arithmetic on the coordinates, that is, the outer red box was obtained via a min/max considering each coordinate as independent of any other. The inner blue boxes are obtained using SLCIA via Eq. (16). We also take random samples represented as points within the two boxes, which are the corresponding actual conformation coordinate values given uncertainty in the distance data. These are the red and blue "swaths" of Fig. 1. We consider the distance data from the Example 5, but with a 0.1 error for the distances  $d_{i-3,i}$ . In particular, the single level symbolic solutions are the blue "swaths" within the min/max SLCIA inner blue box. The red "swaths" are the coordinates calculated as being independent and are found within the outer red box. Both "swaths" were obtained by random sampling and computing the coordinates deterministically.

## 5. Conclusion

The BP algorithm is an efficient method to solve the DMDGP, specially when symmetric properties due to the DMDGP order are considered [11,20]. However, when the DMDGP is associated to protein conformation, where distance data are

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

![img-0.jpeg](img-0.jpeg)
Fig. 1. Five atom molecule example.

provided by an NMR machine, an extension of BP, called iBP [7], must be applied to take care of the uncertainties of the NMR data. In this case, the iBP cannot guarantee anymore that a solution will be found because the way it considers the uncertainties of the problem.

We propose a method to control the uncertainty in the iDMDGP (DMDGP with interval distances), using single-level interval analysis applied to the matrix multiplication used in the BP iteration. Differently from the iBP, our approach maintains the error propagation under control, avoiding sampling values from the interval distances, as done by the iBP. This is the main advantage of the proposed methodology in this article.

We believe that the first step to solve iDMDGP instances has been done, but a lot of work still remains to be done, mainly related to the development of efficient ways to solve the global optimization problems associated to the single-level interval analysis and to the incorporation to our method other interval distances provided by the NMR experiments that can reduce the search space of an iDMDGP instance. This is the next step of our research.

# Acknowledgments

The authors are grateful to the reviewers and editor for the constructive criticisms regarding this article. The authors also wish to thank CNPq and FAPESP for partially supporting this research and of course, for the support from our respective universities.

# References

[1] A. Cassioli, B. Bordeaux, G. Bouvier, A. Mucherino, R. Alves, L. Liberti, M. Nilges, C. Lavor, T. Malliavin, An algorithm to enumerate all possible protein conformations verifying a set of distance constraints, BMC Bioinf. 16 (2015) 16-23.
[2] Y. Chalco-Cano, W.A. Lodwick, B. Bede, Single level constraint interval arithmetic, Fuzzy Set Syst. 256 (2015) 146-168.
[3] C. Lavor, L. Liberti, N. Maculan, Computational experience with themolecular distance geometry problem, in: J. Printer (Ed.), Global Optimization: Scientific and Engineering Case Studies, Springer, Berlin, 2006, 213-225.
[4] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, Recent advances on the discretizable molecular distance geometry problem, Eur. J. Oper. Res. 219 (2012) 698-706.
[5] C. Lavor, J. Lee, A. Lee-St John, L. Liberti, A. Mucherino, M. Sviridenko, Discretization orders for distance geometry problems, Optim. Lett. 6 (2012) 783-796.
[6] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, The discretizable molecular distance geometry problem, Comput. Optim. Appl. 52 (2012) 115-146.
[7] C. Lavor, L. Liberti, A. Mucherino, The interval BP algorithm for the discretizable molecular distance geometry problem with interval data, J. Global Optim. 56 (2013) 855-871.
[8] C. Lavor, R. Alves, W. Figueiredo, A. Petraglia, N. Maculan, Clifford algebra and the discretizable molecular distance geometry problem, Adv. Appl. Clifford Algebra 25 (2015) 925-942.
[9] L. Liberti, C. Lavor, N. Maculan, A branch-and-prune algorithm for the molecular distance geometry problem, Int. Trans. Oper. Res. 15 (2008) 1-17.
[10] L. Liberti, C. Lavor, A. Mucherino, N. Maculan, Molecular distance geometry methods: from continuous to discrete, Int. Trans. Oper. Res. 18 (2010) 33-51.
[11] L. Liberti, B. Masson, J. Lee, C. Lavor, A. Mucherino, On the number of realizations of certain Henneberg graphs arising in protein conformation, Discrete Appl. Math. 165 (2014) 213-232.
[12] L. Liberti, C. Lavor, N. Maculan, A. Mucherino, Euclidean distance geometry and applications, SIAM Rev. 56 (2014) 3-69.
[13] W.A. Lodwick, Constrained interval arithmetic, CCM Report, 138, 1999.
[14] W.A. Lodwick, Interval and fuzzy analysis: An unified approach., in: P.W. Hawkes (Ed.), Advances in Imagining and Electronic Physics, 148, Elsevier, 2007, pp. 75-192.
[15] W.A. Lodwick, O. Jenkins, Constrained intervals and interval spaces, Soft comput. 17 (2013) 1393-1402.

---

T.M. Costa et al./Information Sciences 415-416 (2017) 41-52

[16] W.A. Lodwick, D. Dubois, Interval linear systems as a necessary step in fuzzy linear systems, Fuzzy Set Syst. 274 (2015) 227-251.
[17] R.E. Moore, Interval Analysis, Prentice-Hall, Englewood Cliffs, New Jersey, 1966.
[18] R.E. Moore, R.B. Kearfott, M.J. Cloud, Introduction to Interval Analysis, Society for Industrial and Applied Mathematics, Philadelphia, PA, USA, 2009.
[19] A. Mucherino, C. Lavor, L. Liberti, The discretizable distance geometry problem, Optim. Lett. 6 (2012) 1671-1686.
[20] A. Mucherino, C. Lavor, L. Liberti, Exploiting symmetry properties of the discretizable molecular distance geometry problem, J. Bioinform. Comput. Biol. 10 (2012) 1242009(1-15).
[21] A. Mucherino, C. Lavor, L. Liberti, N. Maculan, Distance Geometry: Theory, Methods, and Applications, Springer, New York, 2013.
[22] T. Schlick, Molecular Modelling and Simulation: An Interdisciplinary Guide, Springer, New York, 2002.