Linear Algebra and its Applications 592 (2020) 287-305

![img-0.jpeg](img-0.jpeg)
ELSEVIER

# Contents lists available at ScienceDirect

# Linear Algebra and its Applications

www.elsevier.com/locate/laa

![img-1.jpeg](img-1.jpeg)

# An dres David Baez Sánchez $^{a,\ast}$ , Carlile Lavor $^{b}$

a Mathematics Department, Federal University of Technology, Av. Sete de Setembro, 3165, 80230-901, Curitiba, Paraná, Brazil
b Applied Mathematics Department, University of Campinas, R. Sergio Buarque de Holanda, 651, 13083-859, Campinas, São Paulo, Brazil

# ARTICLE INFO

Article history:

Received 24 September 2019

Accepted 27 January 2020

Available online 28 January 2020

Submitted by R. Brualdi

MSC:

15A83

Keywords:

Interval Euclidean distance matrix

completion problems

Discretizable molecular distance

geometry problems

Positive semidefinite matrix

# ABSTRACT

We consider some Euclidean distance matrix completion problems whose structure is inspired by molecular conformation problems. Some matrix distances are given precisely or in terms of intervals and other values are unknown. We present theoretical results that can be useful in methods to estimate the possible values of the unknown matrix distances.

© 2020 Elsevier Inc. All rights reserved.

# 1. Introduction

Using nuclear magnetic resonance (NMR) spectroscopy to estimate distances between nuclei of atoms on proteins, some of them remain unknown due to experimental limita

https://doi.org/10.1016/j.laa.2020.01.036

0024-3795/© 2020 Elsevier Inc. All rights reserved.

---

tions. Additionally, because of the non-static behavior of molecules, some measurements may not be given precisely *[1]*. The problem, one of the main applications of Distance Geometry *[2, 3, 4, 5]*, is how to use this distance information to calculate the positions of the protein atoms *[6, 7, 8, 9]*.

A matrix $D\in\mathbb{R}^{n\times n}$, whose elements are given by

$d_{ij}=||x_{i}-x_{j}||_{2}^{2},$

for points $x_{1},...,x_{n}$ in some Euclidean space, is called an Euclidean Distance Matrix (EDM). This is one of the fundamental mathematical concepts in Distance Geometry *[10]*. If $k$ is the dimension of the affine span of $x_{1},...,x_{n}$, we say that $D$ has embedding dimension $k$ or that $D$ is an EDM with a $k$-dimensional realization *[4, 11, 12]*.

In the context of NMR and protein calculations, where we do not have all the distances, we want to complete the associated matrix considering possible 3-dimensional realizations.

EDM completion problems have been extensively studied considering precise *[13, 14, 15, 16, 17]*, and noise data *[18, 19]*.

In protein structure calculations, we may assume that the distances related to bond lengths and bond angles are represented by real numbers *[4]*. Using special atomic orderings for proteins and NMR distance information *[20, 21, 22, 23, 24, 25]*, the following values can be known a priori ($n$ is the number of atoms):

- $d_{1,2},\ldots,d_{n-1,n}$;
- $d_{1,3},\ldots,d_{n-2,n}$;
- $d_{1,4},\ldots,d_{n-3,n}$.

The first two lists of distances can be represented by real numbers, since they are associated to bond lengths and bond angles. The third list, however, may be given by intervals of real numbers, due to imprecision of NMR data *[23]*.

For each instance of the problem, additional distance values $d_{ij}$, $|j-i|>3$, can be provided by NMR if the corresponding atoms are close enough *[26]*. When the distances are precisely given, the problem can be solved by a combinatorial approach *[27, 28]*. The problem is NP-hard in general and the great challenge is to propose efficient methods to solve instances considering the uncertainties of NMR data *[29, 30, 31, 32, 33, 34, 35, 36]*.

In this paper, we present new theoretical results that can be useful to estimate the unknown values of distance matrices, mainly in problems related to protein conformation.

The atomic orderings mentioned above will motivate us to consider a special class of completion problems, described in Section 2, and to define the $\mathbf{d}_{1n}-$EDM completion problem, given in Section 3 (for precise distances) and in Section 4 (for interval distances). In Section 5, we consider $\mathbf{d}_{1n}-$EDM completion problems with fixed embedding dimensions. Section 6 presents our main results and Section 7 provides some examples. Final remarks are presented in Section 8.

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

## 2. A class of EDM completion problems

To motivate the class of EDM completion problems that we will consider, let us start with an example.

**Example 1.** Consider the following EDM completion problem with $n = 6$:

$$
D = \left( \begin{array}{cccc}
0 &amp; d_{12} &amp; d_{13} &amp; [d_{14}] &amp; \mathbf{d}_{15} &amp; \mathbf{d}_{16} \\
d_{12} &amp; 0 &amp; d_{23} &amp; d_{24} &amp; [d_{25}] &amp; \mathbf{d}_{26} \\
d_{13} &amp; d_{23} &amp; 0 &amp; d_{34} &amp; d_{35} &amp; [d_{36}] \\
[d_{14}] &amp; d_{24} &amp; d_{34} &amp; 0 &amp; d_{45} &amp; d_{46} \\
\mathbf{d}_{15} &amp; [d_{25}] &amp; d_{35} &amp; d_{45} &amp; 0 &amp; d_{56} \\
\mathbf{d}_{16} &amp; \mathbf{d}_{26} &amp; [d_{36}] &amp; d_{46} &amp; d_{56} &amp; 0
\end{array} \right).
\tag{1}
$$

Some values in $D$ are precisely known, some $d_{ij}$ are given by intervals (represented in brackets) and others are unknown (represented in bold types). We want to determine the possible values for the unknown distances in the matrix $D = (d_{ij})$, such that the completed matrix corresponds to an EDM with a 3-dimensional realization.

Note that in completion problem (1), if $|j - i| \leq 3$, the distance $d_{ij}$ is given in terms of an interval or a real number. This is the kind of EDM completion problem we are interested in, where there is positive constant $\rho$, $2 &lt; \rho \leq n$, such that if $|j - i| \leq \rho$, the distance $d_{ij}$ is given in terms of an interval or a real number.

If the only unknown value is $\mathbf{d}_{1n}$, the problem is called the $\mathbf{d}_{1n}$-EDM completion problem. To illustrate the relevance of $\mathbf{d}_{1n}$-EDM completion problems, consider again Example 1. Note that if one knows how to solve $\mathbf{d}_{1n}$-EDM completion problems with interval data for $n = 5$, then it would be possible to estimate a set of values for $\mathbf{d}_{15}$ using the known distances in the square sub-matrix $D_{15} = (d_{ij})$, $i,j = 1,2,3,4,5$. Similarly, we could estimate a set of values for $\mathbf{d}_{26}$, using the corresponding square sub-matrix $D_{26} = (d_{ij})$, $i,j = 2,3,4,5,6$. In the same way, if one is able to solve $\mathbf{d}_{1n}$-EDM completion problems with interval data for $n = 6$, then we can obtain an estimate $\mathbf{d}_{16}$, using all known distances, and the values previously obtained for $\mathbf{d}_{15}$ and $\mathbf{d}_{26}$.

In the two following sections we discuss $\mathbf{d}_{1n}$-EDM completion problems with precise and imprecise data.

## 3. $\mathbf{d}_{1n}$-EDM completion problems with precise data

Consider the following $\mathbf{d}_{1n}$-EDM completion problem given by the matrix

$$
D = \left( \begin{array}{cccc}
0 &amp; d_{12} &amp; \dots &amp; d_{1(n-1)} &amp; \mathbf{d}_{1n} \\
d_{12} &amp; 0 &amp; \dots &amp; d_{2(n-1)} &amp; d_{2n} \\
\vdots &amp; \vdots &amp; \vdots &amp; \vdots &amp; \vdots \\
\mathbf{d}_{1n} &amp; d_{2n} &amp; \dots &amp; d_{n(n-1)} &amp; 0
\end{array} \right).
\tag{2}
$$

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

If $\mathbf{d}_{1n}$ is a value such that $D$ is an EDM, it is necessary that all known distances actually correspond to inter-points distances in some Euclidean space. From now on, we assume that $D_{1(n-1)}$ is an EDM, where $D_{1(n-1)}$ is the leading $n-1$ principal sub-matrix of $D$ obtained by eliminating the last row and last column of $D$.

Let $[\mathbf{d}_{1n}]$ be the set of possible values for $\mathbf{d}_{1n}$ such that $D$ is an EDM. The set $[\mathbf{d}_{1n}]$ can be completely determined from the values of $D$ as established in the following theorem.

**Theorem 1.** Consider a $\mathbf{d}_{1n}$-EDM completion problem given by the matrix $D$ in (2) and let $X = (x_{ij})$ be the matrix defined as

$$
x_{ij} = d_{1i} + d_{1j} - d_{ij},\ i,j = 1,2, \dots, n - 1. \tag{3}
$$

Let $u, v$ be column vectors defined as

$$
u_1 = 0 = v_1 \text{ and } u_i = 1,\ v_i = d_{1i} - d_{in},\ \text{for } i = 2, \dots, n - 1. \tag{4}
$$

Also, consider $a = u^T X^+ u$, $b = 2u^T X^+ v$, $c = v^T X^+ v$ and $r_1, r_2$ the solutions of the quadratic equation

$$
0 = a \omega^2 + (b - 2) \omega + c = f(\omega),
$$

where $X^{+}$ denotes the pseudo-inverse of $X$. Then,

(i) if $r_1, r_2$ are real numbers and $XX^+ u - u = XX^+ v - v = 0$, then $[\mathbf{d}_{1n}] = [r_1, r_2]$;
(ii) if $XX^{+}u - u \neq 0$ and $\omega$ is a solution of $-\omega (XX^{+}u - u) = XX^{+}v - v$, with $\omega \geq 0$ and $f(\omega) \leq 0$, then $[\mathbf{d}_{1n}] = \{\omega\}$;
(iii) $[\mathbf{d}_{1n}]$ is an empty set, otherwise.

With some mild additional hypotheses, the set $[\mathbf{d}_{1n}]$ can also be obtained by a different approach. In [37], the author considered the equivalent one missing entry EDM completion problem represented by

$$
\tilde{D} = \left( \begin{array}{ccc} 0 &amp; d^T &amp; \alpha \\ d &amp; D &amp; c \\ \alpha &amp; c^T &amp; 0 \end{array} \right), \tag{5}
$$

where $\alpha$ is unknown, $D$ is an $m \times m$ EDM with embedding dimension $k$, and $D_1 = \left( \begin{array}{cc} 0 &amp; d^T \\ d &amp; D \end{array} \right)$ and $D_2 = \left( \begin{array}{cc} D &amp; c \\ c^T &amp; 0 \end{array} \right)$ are both $(m + 1) \times (m + 1)$ EDM. Theorem 7.1 in [37] established that $\tilde{D}$ is an EDM if and only if $\alpha_l \leq \alpha \leq \alpha_u$, where $\alpha_l = f(d,c) - \sqrt{f(d,d)f(c,c)}$ and $\alpha_u = f(d,c) + \sqrt{f(d,d)f(c,c)}$, with $f(x,y) = \dfrac{e^T(x + y)}{n} + \dfrac{e^T De}{n^2} - \dfrac{e^{T} De}{n^2}$.

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

$\frac{1}{2}\left(x - \frac{De}{n}\right)^T B^+ \left(y - \frac{De}{n}\right)$ and $B = -JDJ/2$, for $J = (I_m - \frac{ee^T}{n})$ and $e$ the vector of all 1's in $\mathbb{R}^m$.

Note that, in terms of notation used in [37], the vectors $u$ and $v$ in (4) can be written as $u = \begin{bmatrix} 0 \\ e \end{bmatrix}$ and $v = \begin{bmatrix} 0 \\ d - c \end{bmatrix}$, and $[\mathbf{d}_{1n}] = [\alpha_l, \alpha_u]$.

For the sake of completeness, in the next section we develop our own approach on the completion problem (2). This also will be useful to prove Theorem 2, which characterizes the values in $[\mathbf{d}_{1n}]$ corresponding to embedding of minimal dimension.

Theorem 2. Consider a $\mathbf{d}_{1n}$-EDM completion problem with precise data and assume that the solution set $[\mathbf{d}_{1n}]$ is not empty. Then $d \in [\mathbf{d}_{1n}]$ corresponds to a realization of minimal embedding dimension if and only if $d$ is an extreme point of $[\mathbf{d}_{1n}]$.

## 3.1. Proofs of Theorems 1 and 2

Let us start by recalling a fundamental result on Euclidean distance matrices and positive semidefinite (PSD) matrices.

Proposition 3 (Gower [38], [39]). Let $D$ be a $n \times n$ symmetric, non-negative matrix with zero diagonal and $\mathbf{e}$ the $n$-dimensional vector with all entries equal to 1. Then $D$ is an EDM if and only if $F = -\frac{1}{2} L^T DL$ is a PSD matrix, where $L = (I - \mathbf{se}^T)$ for any $\mathbf{s}$ such that $\mathbf{s}^T e = 1$ and $\mathbf{s}^T D \neq \mathbf{0}$. In this case, the matrix $F$ is called a Gram matrix of $D$ and the rank of $F$ is equal to the embedding dimension of $D$.

Since the constant $\frac{1}{2}$ is not relevant to determine the positive definiteness or the rank of matrix $F$, we consider just the matrix $-L^T DL$. Next results describe the structure of the matrix $-L^T DL$, when $\mathbf{s}$ is equal to $e_1$, the first canonical vector of $\mathbb{R}^n$.

Lemma 4. Let $B = (b_{ij})$ be a $n \times n$ symmetric matrix with zero diagonal and $A = -L^T BL$, with $L = (I - \mathbf{se}^T)$ and $\mathbf{s} = e_1$. Then, $A = (a_{ij})$ satisfies that $a_{ij} = b_{1i} + b_{1j} - b_{ij}$, for $i,j = 1,\dots,n$ or, equivalently,

$$
A = - L ^ {T} B L = \left( \begin{array}{cc} X &amp; y \\ y ^ {T} &amp; 2 b _ {1 n} \end{array} \right),
$$

where $X = (x_{ij})$, $x_{ij} = b_{1i} + b_{1j} - b_{ij}$, and $y = (y_i)$, $y_i = b_{1i} + b_{1n} - b_{in}$, for $i,j = 1,\ldots ,n - 1$.

Proof. Let $A_{i}, A_{.j}$ be the row $i$ and column $j$ in a given matrix $A$, respectively. Considering $L = (I - \mathbf{se}^T) = (L_{ij})$, we have that

$$
L _ {. j} = e _ {j} - e _ {1} \tag {6}
$$

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

and, therefore, row $i$ in $L^T$ is given by

$$
L_{i.}^{T} = e_{i}^{T} - e_{1}^{T}. \tag{7}
$$

When considering the term $a_{ij}$ from the matrix $A = -L^T BL$, we have that

$$
a_{ij} = -L_{i.}^{T}(BL)_{,j} = -L_{i.}^{T}BL_{,j}, \tag{8}
$$

hence, expressions (6) and (7) lead us to

$$
\begin{array}{l}
a_{ij} = -L_{i.}^{T}BL_{,j} \\
\quad = -(e_{i}^{T} - e_{1}^{T})B(e_{j} - e_{1}) \\
\quad = -(e_{i}^{T} - e_{1}^{T})(B_{,j} - B_{,1}) \\
\quad = -(e_{i}^{T}B_{,j} - e_{i}^{T}B_{,1} - e_{1}^{T}B_{,j} + e_{1}^{T}B_{,1}) \\
\end{array}
$$

$$
a_{ij} = -(b_{ij} - b_{i1} - b_{1j} + b_{11}).
$$

From the fact that $B$ is symmetric with zero diagonal, last equality is equivalent to

$$
a_{ij} = b_{1i} + b_{1j} - b_{ij}. \quad \square
$$

Combining the results from Proposition 3 and Lemma 4, we obtain the next lemma.

**Lemma 5.** $D$ is an EDM if and only if the matrix

$$
- L^{T}DL = \left( \begin{array}{cc} X &amp; y \\ y^{T} &amp; 2\mathbf{d}_{1n} \end{array} \right) \tag{9}
$$

is PSD, where

$$
X = (x_{ij}), \ x_{ij} = d_{1i} + d_{1j} - d_{ij}, \tag{10}
$$

for $i,j = 1,2,\ldots,n-1$, and $y \in \mathbb{R}^{n-1}$ is defined as

$$
y_{i} = d_{1i} + d_{1n} - d_{in}, \tag{11}
$$

for $i = 1,\ldots,n-1$.

In order to determine conditions on $\mathbf{d}_{1n}$ such that the matrix in (9) is PSD, it will be useful to consider the following result.

**Theorem 6** (Albert [40]). Let

$$
S = \left( \begin{array}{cc} S_{11} &amp; S_{12} \\ S_{12}^{T} &amp; S_{22} \end{array} \right),
$$

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

where $S_{11}$ and $S_{22}$ are symmetric matrices. $S$ is PSD if and only if

(i) $S_{11}$ is PSD,
(ii) $S_{22} - S_{12}^T S_{11}^+ S_{12}$ is PSD, and
(iii) $S_{11}S_{11}^{+}S_{12} = S_{12}$

where $S_{11}^{+}$ is the pseudo-inverse of $S_{11}$.

Applying this result to the matrix (9), we obtain immediately the following proposition.

**Proposition 7.** Let $X$ and $y$ be as defined in (10) and (11). $D$ is an EDM if and only if

(i) $X$ is PSD,
(ii) $2\mathbf{d}_{1n} - y^T X^+ y \geq 0$,
(iii) $XX^{+}y = y$

Conditions in Proposition 7 can be expressed explicitly in terms of $\mathbf{d}_{1n}$, as presented in the following result.

**Proposition 8.** Let $X$ be defined as in (10). $D$ is an EDM if and only if $\omega = \mathbf{d}_{1n}$ satisfies

(i) $a\omega^2 + (b - 2)\omega + c \leq 0$,
(ii) $-\omega (XX^{+}u - u) = XX^{+}v - v$

with $a = u^T X^+ u$, $b = 2u^T X^+ v$, and $c = v^T X^+ v$, where $u, v \in \mathbb{R}^{n-1}$ are defined as

$$
u_1 = 0 \text{ and } u_i = 1, i = 2, \dots, n - 1, \tag{12}
$$

$$
v_1 = 0 \text{ and } v_i = d_{1i} - d_{in}, i = 2, \dots, n - 1. \tag{13}
$$

**Proof.** By the assumptions considered in the beginning of Section 3, we have that the sub-matrix $D_{1(n-1)} = (d_{ij})$, $i,j = 1,\dots,n-1$, is an EDM. Therefore, Proposition 3 ensures that the matrix $X$ is PSD, because $X$ is a positive multiple of the Gram matrix associated to $D_{1(n-1)}$.

Note that condition (ii) in Proposition 7 is a real inequality on the variable $\mathbf{d}_{1n}$. In fact, considering $\omega = \mathbf{d}_{1n}$, $y$ as in (11) and $u, v$ defined as in (12) and 13, we have that $y = \omega u + v$, which implies that

$$
\begin{array}{l}
y^T X^+ y = (\omega u + v)^T X^+ (\omega u + v) \tag{14} \\
= \omega^2 (u^T X^+ u) + \omega (2 u^T X^+ v) + v^T X^+ v \tag{15} \\
= a \omega^2 + b \omega + c, \tag{16}
\end{array}
$$

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

with $a = u^T X^+ u$, $b = 2u^T X^+ v$, and $c = v^T X^+ v$. Hence, inequality (ii) in Proposition 7 can be rewritten in terms of $\omega = \mathbf{d}_{1n}$ as

$$
a \omega^ {2} + b \omega + c - 2 \omega \leq 0. \tag {17}
$$

Finally, considering the condition $XX^{+}y = y$ in Proposition 7 and that $y = \omega u + v$, we obtain

$$
\begin{array}{l} X X ^ {+} y = y \\ X X ^ {+} (\omega u + v) = \omega u + v \\ - \omega (X X ^ {+} u - u) = (X X ^ {+} v - v). \quad \square \\ \end{array}
$$

Proposition 8 implies that the set $[\mathbf{d}_{1n}]$ of possible completion values for $\mathbf{d}_{1n}$, such that $D$ is an EDM, can be obtained by determining all values $\omega$ satisfying both conditions in Proposition 8. This is the main line of reasoning in the proof of Theorem 1 presented below.

Proof of Theorem 1. Consider first that $r_1, r_2$ are real numbers and that $XX^+ u - u = 0 = XX^+ v - v$. In this case, condition (ii) in Proposition 8 is trivially satisfied for all $\omega$. Thus, we are only concerned with condition (i) in Proposition 8.

Since $X^{+}$ is PSD if and only if $X$ is PSD, then $y^{T}X^{+}y\geq 0$, for all $y$ and also $a,c\geq 0$. From expression (14), we obtain therefore that, $a\omega^2 +b\omega +c\geq 0$, for all $\omega$, implying that $b^{2} - 4ac\leq 0$. On the other hand, we must have $(b - 2)^{2} - 4ac = b^{2} - 4ac - 4b + 4\geq 0$ because $a\omega^2 +(b - 2)\omega +c$ has real roots (precisely $r_1$ and $r_2$), concluding that $-4b + 4\geq 0$ which is equivalent to $b\leq 1$. Under this condition, it is straightforward to show that $r_1$ and $r_2$ must be non-negative numbers.

The condition $XX^{+}u - u = 0$ implies also that $a &gt; 0$. In fact, if $a = u^T X^+ u = 0$, we would obtain $X^{+}u = 0$ (because $X^{+}$ is PSD), and therefore, $u = XX^{+}u = 0$, which is an absurd by the definition of $u$. Hence, $f(\omega) \leq 0$ if and only if $\omega = \mathbf{d}_{1n}$ belongs to the interval $[r_1, r_2]$.

If $r_1, r_2$ are not real numbers, $f(\omega) &gt; 0$ for all $\omega$, due to the continuity of $f$ and the fact that $f(0) = c \geq 0$. Therefore, there is no $\omega$ satisfying $f(\omega) \leq 0$ (condition (i) in Proposition 8), which implies that $[\mathbf{d}_{1n}]$ is empty.

Considering now that $XX^{+}u - u \neq 0$, there is at least a non-zero component, let us say component $j$, which implies that the unique possible value for $\omega$, such that $-\omega (XX^{+}u - u) = XX^{+}v - v$, is given by $\omega = \frac{(XX^{+}v - v)_j}{(XX^{+}u - u)_j}$. If $\omega$ is non-negative and also satisfies $f(\omega) \leq 0$, then $[\mathbf{d}_{1n}] = \{\omega\}$. Otherwise, condition (ii) in Proposition 8 is not satisfied and $[\mathbf{d}_{1n}]$ is empty. Finally, if $XX^{+}u - u = 0$ and $XX^{+}v - v \neq 0$, condition (ii) in Proposition 8 can not be satisfied and the set $[\mathbf{d}_{1n}]$ is empty.

Before considering Theorem 2, recall that, according to Proposition 3, the embedding dimension of $D$ is equal to the rank of $-L^T DL$. In order to determine the relationship

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

between this rank and $\omega$, we will use a result concerning the rank of a block matrix. In what follows, $r(A)$ denotes the rank of the matrix $A$ and $\mathcal{R}(A)$ denotes the space generated by its columns.

**Lemma 9** (Marsaglia and Styan [41], Meyer [42]). Let $A, B, C, D$ be given matrices with dimensions such that the block matrix

$$
\left( \begin{array}{cc} A &amp; B \\ C &amp; D \end{array} \right)
$$

is well defined. If $\mathcal{R}(B) \subset \mathcal{R}(A)$ and $\mathcal{R}(C^T) \subset \mathcal{R}(A^T)$, then

$$
r \left( \begin{array}{cc} A &amp; B \\ C &amp; D \end{array} \right) = r(A) + r(D - CA^+ B).
$$

Lemma 9 is used to obtain the next proposition.

**Proposition 10.** Let $\mathbf{d}_{1n}$ be such that $D$ is an EDM (i.e. the conditions in Proposition 7 are satisfied with $X$ and $y$ defined by expressions (10) and (11), respectively). Then,

$$
r(-L^T D L) = r(X) + r(2 \mathbf{d}_{1n} - y^T X^+ y).
$$

**Proof.** If the conditions in Proposition 7 are satisfied, then in particular $y$ belongs to $\mathcal{R}(X)$, so we can apply Lemma 9 to the block matrix (9) to obtain the desired result. $\square$

**Remark 1.** Since the expression $\mathbf{d}_{1n} - y^T X^+ y$ is a scalar, Proposition 10 implies that the embedding dimension of the completed matrix $D$ is either equals to $r(X) + 1$ or to $r(X)$. Since $X$ is a positive multiple of the Gram matrix of the leading sub-matrix $D_{1(n-1)}$, from Proposition 3 it follows that $r(X)$ is equal to the embedding dimension of $D_{1(n-1)}$. Hence, $D_{1(n-1)}$ and $D$ have the same embedding dimensions if and only if $2\mathbf{d}_{1n} - y^T X^+ y = 0$. Note also that the expression $2\mathbf{d}_{1n} - y^T X^+ y$ has already appeared in Proposition 7, where we have proved its relationship with the quadratic polynomial $a\omega^2 + (b - 2)\omega + c$ with $\omega = \mathbf{d}_{1n}$. These comments lead us to the following corollaries.

**Corollary 11.** Let $X, u, v, a, b, c$ as defined in Theorem 1 and consider $d \in [\mathbf{d}_{1n}]$. The embedding dimension of the completed matrix $D$ is equal to $r(X) + h(d)$, where

$$
h(\omega) = \left\{ \begin{array}{ll} 0, &amp; \text{if } a\omega^2 + (b - 2)\omega + c = 0 \\ 1, &amp; \text{otherwise.} \end{array} \right.
$$

**Corollary 12.** Assume that $D_{1(n-1)}$ is an EDM with embedding dimension $k$. Let $X, u, v, a, b$ and $c$ as defined in Theorem 1. Then, $\mathbf{d}_{1n}$ is such that the embedding dimension of the completed matrix $D$ is equal to $k$, if and only if $\omega = \mathbf{d}_{1n}$ is a solution of $a\omega^2 + (b - 2)\omega + c = 0$.

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

Now, we finally prove Theorem 2.

Proof of Theorem 2. By Theorem 1 we have that if  $[\mathbf{d}_{1n}]$  is not empty, then  $[\mathbf{d}_{1n}]$  is a singleton or  $[\mathbf{d}_{1n}] = [r_1,r_2]$ , where  $r_1,r_2$  are the solutions of the quadratic equation  $a\omega^2 +(b - 2)\omega +c = 0$ .

If  $[\mathbf{d}_{1n}]$  is a singleton, let us say  $[\mathbf{d}_{1n}] = \{d\}$ , this means that there is only one possible completion for matrix  $D$  and therefore only one possible embedding dimension. This unique embedding dimension is of course minimal and  $d$  is clearly an extreme point of  $[\mathbf{d}_{1n}] = \{d\}$ .

By Corollary 12, there are at most two possible values for the embedding dimension of the completed matrix  $D$ ; either  $r(X)$  or  $r(X) + 1$ . Also,  $d = r_1$  and  $d = r_2$  are the only values that could complete the matrix with an embedding dimension equal to  $r(X)$  (which would be minimal). Hence, if  $[\mathbf{d}_{1n}] = [r_1,r_2], d \in [\mathbf{d}_{1n}]$  corresponds to realizations of minimal embedding dimension if and only if  $d = r_1$  or  $d = r_2$ , i.e., if and only if  $d$  is one of the extreme points of  $[\mathbf{d}_{1n}]$ .

Remark 2. In general, for EDM completion problems with more than one unknown value, it can be proved, from spectrahedra theory or directly from properties of positive semidefinite matrices (see [43] or [44]), that completion values corresponding to minimal embedding dimension must be extreme points of the set of possible completion values. The reciprocal may not be true in general, but Theorem 2 establishes both implications when there is only one unknown distance.

In the next section, we consider  $\mathbf{d}_{1n}$ -EDM completion problems with imprecise data.

# 4.  $\mathbf{d}_{1n}$ -EDM completion problems with imprecise data

If the value  $\mathbf{d}_{1n}$  is unknown, but there is available information (precise or in terms of intervals) about all other entries in the associated matrix, we represent this situation as follows:

$$
D (x, \mathbf {d} _ {\mathbf {1 n}}) = \left( \begin{array}{c c} &amp; \mathbf {d} _ {\mathbf {1 n}} \\ &amp; A (x) \\ \mathbf {d} _ {\mathbf {1 n}} &amp; \end{array} \right), x \in I, \tag {18}
$$

for some set  $I$ . The gray area  $A(x)$  represents the known data (precise or imprecise) depending on some variable  $x$ . The notation  $D(x, \mathbf{d}_{1n})$  indicates the two different components in the distance matrix: the known part depending on  $x$  and the unknown part  $\mathbf{d}_{1n}$ . Note that this framework does not consider any specific dimension for  $x$  or  $I$ , neither some specific positions for the imprecise terms in  $A$ , which allows to consider different kinds of imprecision located at different positions in the matrix.

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

Example 2. Consider the following completion problem given by

$$
D = \left( \begin{array}{cccc}
0 &amp; d_{12} &amp; [d_{13}] &amp; [d_{14}] &amp; \mathbf{d}_{15} \\
d_{12} &amp; 0 &amp; [d_{23}] &amp; d_{24} &amp; [d_{25}] \\
[d_{13}] &amp; [d_{23}] &amp; 0 &amp; d_{34} &amp; d_{35} \\
[d_{14}] &amp; d_{24} &amp; d_{34} &amp; 0 &amp; d_{45} \\
\mathbf{d}_{15} &amp; [d_{25}] &amp; d_{35} &amp; d_{45} &amp; 0
\end{array} \right), \tag{19}
$$

which can also be represented as

$$
D(x, \mathbf{d}_{15}) = \left( \begin{array}{cccc}
0 &amp; d_{12} &amp; x_1 &amp; x_3 &amp; \mathbf{d}_{15} \\
d_{12} &amp; 0 &amp; x_2 &amp; d_{24} &amp; x_4 \\
x_1 &amp; x_2 &amp; 0 &amp; d_{34} &amp; d_{35} \\
x_3 &amp; d_{24} &amp; d_{34} &amp; 0 &amp; d_{45} \\
\mathbf{d}_{15} &amp; x_4 &amp; d_{35} &amp; d_{45} &amp; 0
\end{array} \right),
$$

with $x = (x_1, x_2, x_3, x_4) \in I = [d_{13}] \times [d_{23}] \times [d_{14}] \times [d_{25}]$.

In order to distinguish two different completion problems obtained by fixing a value $x$ or by considering all the possible values of $x$, we introduce the following definitions.

For each $x \in I$, let $[\mathbf{d}_{\mathbf{1n}}(x)]$ be defined as

$$
[\mathbf{d}_{\mathbf{1n}}(x)] = \{\mathbf{d} \in \mathbb{R} : D(x, \mathbf{d}) \text{ is an EDM}\}, \tag{20}
$$

and let $[\mathbf{d}_{\mathbf{1n}}]$ be defined as

$$
[\mathbf{d}_{\mathbf{1n}}] = \{\mathbf{d} \in \mathbb{R} : D(x, \mathbf{d}) \text{ is an EDM for some } x \in I\}. \tag{21}
$$

When the data are precisely defined, the set $I$ is just a point and both sets are the same.

# 5. The set $[\mathbf{d}_{\mathbf{1n}}]_k$

In some applications, we may be interested in completion values corresponding to a specific embedding dimension (dimension three in the case of protein conformation problems), which induce us to consider the following definition:

$$
[\mathbf{d}_{\mathbf{1n}}]_k = \{\mathbf{d} \in \mathbb{R} : D(x, \mathbf{d}) \text{ is a } k\text{-dimensional EDM for some } x \in I\}.
$$

Depending on the value of $k$, this set can take very different forms as illustrated in the next example.

---

A.D. Baez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

![img-2.jpeg](img-2.jpeg)
Fig. 1. Realization for matrix $D(x_{1},x_{2})$ given in Example 3.

Example 3. Consider the following completion problem, given by

$$
D \left(x _ {1}, x _ {2}\right) = \left( \begin{array}{c c c c} 0 &amp; 4 &amp; x _ {1} &amp; \mathbf {d} _ {1 4} \\ 4 &amp; 0 &amp; 1 &amp; 4 \\ x _ {1} &amp; 1 &amp; 0 &amp; x _ {2} \\ \mathbf {d} _ {1 4} &amp; 4 &amp; x _ {2} &amp; 0 \end{array} \right), \tag {22}
$$

with $(x_{1},x_{2})\in [d_{13}]\times [d_{34}]$. Since the second row (column) is completely known, all distances to point $p_2$ are well established. We can associate this distance matrix with a 3-dimensional realization where points $p_1$ and $p_4$ lie on a sphere with center at $p_2$ and radius 2, and $p_3$ lies on a sphere with center at $p_2$ and radius 1, as illustrated in Fig. 1.

Without additional information, $x_{1}$ and $x_{2}$ are constraint to take values in the interval [1,9] and $\mathbf{d}_{14}$ takes values in the interval [0, 16]. Because the value 0 is incompatible with a 3-dimensional realization (three points will correspond to a two-dimensional realization at most), we have that $[\mathbf{d}_{14}]_3 = (0,16)$. On the other hand, the whole interval [0, 16] is compatible with two-dimensional realizations, implying that $[\mathbf{d}_{14}]_2 = [0,16]$.

If $x_{1} = x_{2} = 4$, the points $p_{1}, p_{2}, p_{3}$ and $p_{2}, p_{3}, p_{4}$ form isosceles triangles as illustrated in Fig. 2. In this case, $\mathbf{d}_{14}$ takes values in the interval [0, 15]. The value $\mathbf{d}_{14} = 15$ corresponds to a situation where the points define a parallelogram with four sides equal to 2 and minor diagonal equal to 1. The situation $\mathbf{d}_{14} = 0$, when $p_{1}$ and $p_{4}$ coincide, can be interpreted as a limit situation of the parallelogram folding through the segment joining $p_{2}$ and $p_{3}$. Thus, we have that $[\mathbf{d}_{14}]_{2} = \{0, 15\}$ and $[\mathbf{d}_{14}]_{3} = (0, 15)$.

Finally, consider that $x_{1} = 4$ and $x_{2} \in [3.61, 4]$. In this case, $\mathbf{d}_{14}$ takes again values in the interval [0, 15], but there are several distances compatible with 2-dimensional realizations: some close to 15 (when the polygon is extended) and some close to 0 (when the polygon is folded). The set $[d_{14}]_2$ will take the form of a disjoint union of intervals and no longer a unique connected interval. A better estimation of this set will be obtained in Example 4, as a consequence of the result presented in the following section.

---

A.D. Baez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

![img-3.jpeg](img-3.jpeg)
Fig. 2. Realization for matrix  $D(4,4)$  in Example 3.

# 6. Estimation of the set  $[\mathbf{d}_{1n}]_k$  for minimal embedding dimension  $k$

Consider the EDM completion problem related to (2). If we assume that the known distance values already correspond to a 3-dimensional realization, it is impossible to complete the matrix and obtain an embedding dimension less than 3. Hence, if there is a completion with embedding dimension 3, this is in fact the minimal embedding dimension. This is the reason why we are interested in completion values corresponding to the minimal embedding dimension.

In this section, we consider the completion problem related to (18) and develop an estimation of the set  $[\mathbf{d}_{1n}]_k$ , when  $k$  is the minimal embedding dimension.

Before establishing our main result (Theorem 13), let us introduce some functions and related optimization problems. For a fixed  $x \in I$ ,  $l(x)$  and  $u(x)$  are defined, respectively, as the lower and upper bounds for the interval  $[\mathbf{d}_{1n}(x)]$ , i.e.

$$
[ \mathbf {d} _ {1 n} (x) ] = [ l (x), u (x) ].
$$

We consider also the following optimization problems:

$$
\begin{array}{l l} \min  (\max) &amp; l (x) \\ \text {s . t .} x \in I &amp; \text {a n d} \end{array} \quad \begin{array}{l l} \min  (\max) &amp; u (x) \\ \text {s . t .} x \in I &amp; \text {a n d} \end{array} . \tag {23}
$$

Let  $l_{min}, l_{max}$  and  $u_{min}, u_{max}$  be the corresponding optimal values of the optimization problems defined in (23). The following theorem shows how one can use these values to estimate the set  $[\mathbf{d}_{1n}]_k$ , for the minimal embedding dimension  $k$ .

Theorem 13. Let  $I$  be a set in some Euclidean space and assume that, for each  $x \in I$ , the leading  $n - 1$  principal square sub-matrix of  $A(x)$  is an EDM with embedding dimension  $k$  and that there exists at least one possible value for  $\mathbf{d}_{1n}$  such that  $D(x, \mathbf{d}_{1n})$  is an EDM with embedding dimension  $k$ . Assume also that  $l_{min}, l_{max}$  and  $u_{min}, u_{max}$  are attained in  $I$ . Then,

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

$$
[ \mathbf {d} _ {\mathbf {1 n}} ] _ {k} \subset [ l _ {m i n}, l _ {m a x} ] \cup [ u _ {m i n}, u _ {m a x} ]
$$

and

$$
\left\{l _ {m i n}, l _ {m a x}, u _ {m i n}, u _ {m a x} \right\} \subset \left[ \mathbf {d} _ {\mathbf {1 n}} \right] _ {k}.
$$

Proof. From the hypothesis, for each value $x$, the leading $n - 1$ principal square submatrix of $A(x)$ is an EDM with embedding dimension $k$, which implies that any possible completion value $\mathbf{d}_{1n}$ must be associated to a realization with embedding dimension at least equal to $k$. Since we are assuming that there exists at least one possible value for $\mathbf{d}_{1n}$ such that $D(x, \mathbf{d}_{\mathbf{1n}})$ is an EDM with embedding dimension $k$, then, for each $x$, the minimal embedding dimension of a possible completion $D(x, \mathbf{d}_{\mathbf{1n}})$ is $k$.

Theorem 2 establishes that completion values corresponding to minimal embedding dimension must be extreme points of the set of possible completion values $[\mathbf{d}_{1n}(x)]$ and that all extreme values are associated with matrix completions having minimal embedding dimension. Hence, for all $x$, $l(x)$ and $u(x)$ are values associated to matrix completions with minimal embedding dimension $k$.

If $d \in [\mathbf{d}_{\mathbf{1n}}]_k$, then, $D(x,d)$ is a $k$-dimensional EDM for some $x \in I$. By Theorem 2, this implies that $d$ is an extreme point of $[\mathbf{d}_{\mathbf{1n}}(x)]$ and, therefore, $d$ is either equal to $l(x)$ or equal to $u(x)$. If $d = l(x)$, $l_{min} \leq d \leq l_{max}$ by definition of $l_{min}$ and $l_{max}$. Similarly, if $d = u(x)$, $u_{min} \leq d \leq u_{max}$. This proves the first statement of the proposition.

From the hypothesis and the definition of $l_{min}$, we have that $l_{min} = l(x^{*})$, for some $x^{*} \in I$, implying that $l_{min}$ is the left extreme of $[\mathbf{d}_{\mathbf{1n}}(x^{*})]$. By Theorem 2, $D(x^{*}, l_{min})$ has minimal embedding dimension, which in this case is equal $k$, implying that $l_{min} \in [\mathbf{d}_{\mathbf{1n}}]_k$. Since the reasoning is analogous for $l_{max}, u_{min}$ and $u_{max}$, we conclude that

$$
\left\{l _ {m i n}, l _ {m a x}, u _ {m i n}, u _ {m a x} \right\} \subset [ \mathbf {d} _ {\mathbf {1 n}} ] _ {k}. \quad \square
$$

Based on the previous result, we propose an estimation for the set $[\mathbf{d}_{\mathbf{1n}}]_k$, given by

$$
\left[ \mathbf {d} _ {\mathbf {1 n}} \right] _ {k} \approx \left[ l _ {\min }, l _ {\max } \right] \cup \left[ u _ {\min }, u _ {\max } \right]. \tag {24}
$$

## 7. Examples

In this section we consider some examples to illustrate numerically the kind of results obtained when using the expression (24). The values of objective functions $l(x)$ and $u(x)$ were computed implementing the results described in Theorem 1 using the open software GNU-Octave, and the optimal values of $l(x)$ and $u(x)$ were estimated using the numerical global optimization tools also available in GNU-Octave.

---

A.D. Baez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

Example 4. Consider the completion problem given by

$$
D (x) = \left( \begin{array}{c c c c} 0 &amp; 4 &amp; 4 &amp; \mathbf {d} _ {1 4} \\ 4 &amp; 0 &amp; 1 &amp; 4 \\ 4 &amp; 1 &amp; 0 &amp; x \\ \mathbf {d} _ {1 4} &amp; 4 &amp; x &amp; 0 \end{array} \right),
$$

with  $x \in [3.61, 4]$ . For each  $x$ , the matrix  $\left( \begin{array}{ccc}0 &amp; 4 &amp; 4\\ 4 &amp; 0 &amp; 1\\ 4 &amp; 1 &amp; 0 \end{array} \right)$  corresponds to a realization with embedding dimension 2 (three non-colinear points). The optimal values of  $l(x)$  and  $u(x)$ , subject to  $x \in [3.61, 4]$ , were computed as  $l_{min} = 0$ ,  $l_{max} = 0.041763$ ,  $u_{min} = 14.56796$ , and  $u_{max} = 15$ . Thus, according to Theorem 13, we can estimate  $[\mathbf{d}_{14}]_2$  as  $[\mathbf{d}_{14}]_2 \approx [0, 0.041763] \cup [14.56796, 15]$ . This result improves the initial naive guess obtained in Example 3.

Example 5. Consider the completion problem corresponding to

$$
D (x _ {1}, x _ {2}) = \left( \begin{array}{c c c c c} 0 &amp; 2. 3 3 8 4 5 &amp; 6. 4 4 9 5 7 &amp; x _ {1} &amp; \mathbf {d} _ {1 5} \\ 2. 3 3 8 4 5 &amp; 0 &amp; 2. 3 2 9 8 9 &amp; 5. 4 4 8 9 6 &amp; x _ {2} \\ 6. 4 4 9 5 7 &amp; 2. 3 2 9 8 9 &amp; 0 &amp; 2. 2 6 2 6 2 &amp; 5. 4 6 4 8 4 \\ x _ {1} &amp; 5. 4 4 8 9 6 &amp; 2. 2 6 2 6 2 &amp; 0 &amp; 2. 3 7 4 0 6 \\ \mathbf {d} _ {1 5} &amp; x _ {2} &amp; 5. 4 6 4 8 4 &amp; 2. 3 7 4 0 6 &amp; 0 \end{array} \right),
$$

with  $x_{1} \in [13.71887, 13.91887]$  and  $x_{2} \in [5.19658, 5.39679]$ . Assuming that, for all values of  $x_{1}$  and  $x_{2}$ , the distances correspond to some 3-dimensional realization, the optimal values of  $l(x_{1}, x_{2})$  and  $u(x_{1}, x_{2})$  were estimated as  $l_{min} = 12.37632$ ,  $l_{max} = 13.17182$ ,  $u_{min} = 14.39595$ , and  $u_{max} = 14.81403$ . So, the estimation of  $[\mathbf{d}_{15}]_3$  is given by  $[\mathbf{d}_{15}]_3 \approx [12.37632, 13.17182] \cup [14.39595, 14.81403]$ .

Example 6. Consider the following completion problem given by

$$
\left( \begin{array}{c c c c c c} 0 &amp; 2. 3 3 8 4 5 &amp; 6. 4 4 9 5 7 &amp; x _ {1} &amp; \mathbf {d} _ {1 5} &amp; \mathbf {d} _ {1 6} \\ 2. 3 3 8 4 5 &amp; 0 &amp; 2. 3 2 9 8 9 &amp; 5. 4 4 8 9 6 &amp; x _ {2} &amp; \mathbf {d} _ {2 6} \\ 6. 4 4 9 5 7 &amp; 2. 3 2 9 8 9 &amp; 0 &amp; 2. 2 6 2 6 2 &amp; 5. 4 6 4 8 4 &amp; x _ {3} \\ x _ {1} &amp; 5. 4 4 8 9 6 &amp; 2. 2 6 2 6 2 &amp; 0 &amp; 2. 3 7 4 0 6 &amp; 9. 8 0 4 4 1 \\ \mathbf {d} _ {1 5} &amp; x _ {2} &amp; 5. 4 6 4 8 4 &amp; 2. 3 7 4 0 6 &amp; 0 &amp; 5. 9 9 2 2 1 \\ \mathbf {d} _ {1 6} &amp; \mathbf {d} _ {2 6} &amp; x _ {3} &amp; 9. 8 0 4 4 1 &amp; 5. 9 9 2 2 1 &amp; 0 \end{array} \right),
$$

with  $x_{1}\in [13.71887,13.91887]$ ,  $x_{2}\in [5.19658,5.39679]$ , and  $x_{3}\in [17.10319,17.30310]$ .

Since the completion problem associated with the sub-matrix  $D_{15}$  is the same problem considered in Example 5, we already have an estimation for  $[\mathbf{d}_{15}]_3$ , given by  $[\mathbf{d}_{15}]_3 \approx$

---

A.D. Baez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

[12.37632, 13.17182]  $\cup$  [14.39595, 14.81403]. Following a similar process for  $\mathbf{d}_{26}$ , we obtain  $[\mathbf{d}_{26}]_3 \approx [13.37584, 13.81981] \cup [20.74529, 21.46747]$ .

Now, we can use these estimations as new constraints in order to obtain an estimation for  $[d_{16}]_3$ , by considering the completion problem given by

|  0 | 2.33845 | 6.44957 | x1 | x4 | d16  |
| --- | --- | --- | --- | --- | --- |
|  2.33845 | 0 | 2.32989 | 5.44896 | x2 | x5  |
|  6.44957 | 2.32989 | 0 | 2.26262 | 5.46484 | x3  |
|  x1 | 5.44896 | 2.26262 | 0 | 2.37406 | 9.80441  |
|  x4 | x2 | 5.46484 | 2.37406 | 0 | 5.99221  |
|  d16 | x5 | x3 | 9.80441 | 5.99221 | 0  |

taking  $x_{1}, x_{2}, x_{3}$  as before,  $x_{4} \in [12.37632, 13.17182] \cup [14.39595, 14.81403]$  and  $x_{5} \in [13.37584, 13.81981] \cup [20.74529, 21.46747]$ . Considering all the possible situations, we obtained the following estimations for  $[\mathbf{d}_{16}]_3$ :

$[\mathbf{d}_{16}]_3 \approx [24.69692, 26.40726] \cup [24.69692, 26.40726]$ , if  $x_4 \in [12.37632, 13.17182]$  and  $x_5 \in [13.37584, 13.81981]$ ;

$[\mathbf{d}_{16}]_3 \approx [30.15657, 33.75378] \cup [30.15657, 33.71428]$ , if  $x_4 \in [12.37632, 13.17182]$  and  $x_5 \in [20.74529, 21.46747]$ ;

$[\mathbf{d}_{16}]_3 \approx [34.33843, 37.70451] \cup [35.68868, 40.15630]$ , if  $x_4 \in [14.39595, 14.81403]$  and  $x_5 \in [20.74529, 21.46747]$ ;

$[\mathbf{d}_{16}]_3 \approx [22.61954, 25.46312] \cup [23.83490, 26.02530]$ , if  $x_4 \in [14.39595, 14.81403]$  and  $x_5 \in [13.37584, 13.81981]$ .

Considering the overlapping between all sets, we finally obtain a general estimation given by

$[\mathbf{d}_{16}]_3 \approx [22.61954, 26.40726] \cup [30.15657, 33.75378] \cup [34.33843, 40.15630]$ .

# 8. Final comments

The main contribution of this paper is a theoretical development of a method for estimation values for the set  $[\mathbf{d}_{1n}]_k$ , where  $k$  is minimum. In particular, this can be very useful in combinatorial methods applied to NMR protein calculations that must deal with uncertainties of NMR data [45].

In order to illustrate this possibility, let us suppose that in Example 6 an additional distance information is provided, given by  $\mathbf{d}_{16} \in [28.09, 32.49]$ . This new information can be used to improve the estimation for  $\mathbf{d}_{16}$ , because from the second partial estimation previously obtained, we could conclude now that  $[\mathbf{d}_{16}]_3 \approx [30.15657, 32.49]$ . In addition

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

to that, this new information will allow us to improve the estimations for $\mathbf{d}_{15}$ and $\mathbf{d}_{26}$ resulting in $[\mathbf{d}_{15}]_3 \approx [12.37632, 13.17182]$ and $[\mathbf{d}_{26}]_3 \approx [20.74529, 21.46747]$.

There are some results for taking care of uncertainties in the branching phase of combinatorial methods like BP algorithm [46–49], but it remains a big challenge to do the same for the pruning phase of BP. The interval reductions discussed above can open new research directions, mainly related to the dual BP algorithm proposed in [50], considering just precise distance values.

## Declaration of competing interest

The authors declare that they have no conflict of interest.

## Acknowledgements

The authors would like to thank the Brazilian research agencies CNPq and FAPESP, for their financial support, and the anonymous referee for valuable suggestions.

## References

[1] K. Wüthrich, Protein structure determination in solution by nuclear magnetic resonance spectroscopy, Science 243 (1989) 45–50, https://doi.org/10.1126/science.2911719.
[2] L. Liberti, C. Lavor, Six mathematical gems from the history of distance geometry, Int. Trans. Oper. Res. 23 (2016) 897–920, https://doi.org/10.1111/itor.12170.
[3] L. Liberti, C. Lavor, Euclidian Distance Geometry: An Introduction, Springer, New York, 2017.
[4] L. Liberti, C. Lavor, N. Maculan, A. Mucherino, Euclidean distance geometry and applications, SIAM Rev. 56 (2014) 3–69, https://doi.org/10.1137/120875909.
[5] A. Mucherino, C. Lavor, L. Liberti, N. Maculan, Distance Geometry: Theory, Methods, and Applications, Springer, New York, 2013.
[6] S. Billinge, P. Duxbury, D. Gonçalves, C. Lavor, A. Mucherino, Assigned and unassigned distance geometry: applications to biological molecules and nanostructures, 4OR 14 (2016) 337–376, https://doi.org/10.1007/s10288-016-0314-2.
[7] S. Billinge, P. Duxbury, D. Gonçalves, C. Lavor, A. Mucherino, Recent results on assigned and unassigned distance geometry with applications to protein molecules and nanostructures, Ann. Oper. Res. 271 (2018) 161–203, https://doi.org/10.1007/s10479-018-2989-6.
[8] G. Crippen, T. Havel, Distance Geometry and Molecular Conformation, Wiley, New York, 1988.
[9] C. Lavor, L. Liberti, W. Lodwick, T.M. da Costa, An Introduction to Distance Geometry Applied to Molecular Geometry, Springer Briefs, Springer, New York, 2017.
[10] L. Blumenthal, Theory and Applications of Distance Geometry, Oxford University Press, Oxford, 1953.
[11] J. Alencar, C. Lavor, L. Liberti, Realizing Euclidean distance matrices by sphere intersection, Discrete Appl. Math. 256 (2019) 5–10, https://doi.org/10.1016/j.dam.2018.06.003.
[12] N. Moreira, L. Duarte, C. Lavor, C. Torezzan, A novel low-rank matrix completion approach to estimate missing entries in Euclidean distance matrix, Comput. Appl. Math. 37 (2018) 4989–4999, https://doi.org/10.1007/s40314-018-0613-7.
[13] A. Alfakih, A. Khandani, H. Wolkowicz, Solving Euclidean distance matrix completion problems via semidefinite programming, Comput. Optim. Appl. 12 (1999) 13–30, https://doi.org/10.1023/A:1008655427845.
[14] M. Bakonyi, C. Johnson, The Euclidean distance matrix completion problem, SIAM J. Matrix Anal. Appl. 16 (2) (1995) 646–654, https://doi.org/10.1137/S0895479893249757.
[15] I. Dokmanić, R. Parhizkar, J. Ranieri, M. Vetterli, Euclidean distance matrices: essential theory, algorithms, and applications, IEEE Signal Process. Mag. 32 (2015) 12–30, https://doi.org/10.1109/msp.2015.2398954.

---

A.D. Baez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

[16] H. Fang, D.P. O'Leary, Euclidean distance matrix completion problems, Optim. Methods Softw. 27 (4–5) (2012) 695–717, https://doi.org/10.1080/10556788.2011.643888.

[17] M. Laurent, A tour d'horizon on positive semidefinite and Euclidean distance matrix completion problems, in: Fields Institute Communications, vol. 18, American Mathematical Society, United States, 1998, pp. 51–76.

[18] E.J. Candès, Y. Plan, Matrix completion with noise, Proc. IEEE 98 (6) (2010) 925–936, https://doi.org/10.1109/JPROC.2009.2035722.

[19] E.J. Candès, B. Recht, Exact matrix completion via convex optimization, Found. Comput. Math. 9 (6) (2009) 717, https://doi.org/10.1007/s10208-009-9045-5.

[20] A. Cassioli, O. Gunluk, C. Lavor, L. Liberti, Discretization vertex orders in distance geometry, Discrete Appl. Math. 197 (2015) 27–41, https://doi.org/10.1016/j.dam.2014.08.035.

[21] V. Costa, A. Mucherino, C. Lavor, A. Cassioli, L. Carvalho, N. Maculan, Discretization orders for protein side chains, J. Global Optim. 60 (2014) 333–349, https://doi.org/10.1007/s10898-013-0135-1.

[22] C. Lavor, J. Lee, A.L.-S. John, L. Liberti, A. Mucherino, M. Sviridenko, Discretization orders for distance geometry problems, Optim. Lett. 6 (2012) 783–796, https://doi.org/10.1007/s11590-011-0302-6.

[23] C. Lavor, L. Liberti, B. Donald, B. Worley, B. Bardiaux, T. Malliavin, M. Nilges, Minimal nmr distance information for rigidity of protein graphs, Discrete Appl. Math. 256 (2019) 91–104, https://doi.org/10.1016/j.dam.2018.03.071.

[24] C. Lavor, A. Mucherino, L. Liberti, N. Maculan, On the computation of protein backbones by using artificial backbones of hydrogens, J. Global Optim. 50 (2011) 329–344, https://doi.org/10.1007/s10898-010-9584-y.

[25] C. Lavor, M. Souza, L. Mariano, L. Liberti, On the polinomiality of finding $K$ dmdgp re-orders, Discrete Appl. Math. 267 (2019) 190–194, https://doi.org/10.1016/j.dam.2019.07.021.

[26] B. Donald, Algorithms in Structural Molecular Biology, MIT Press, Boston, 2011.

[27] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, The discretizable molecular distance geometry problem, Comput. Optim. Appl. 52 (2012) 115–146, https://doi.org/10.1007/s10589-011-9402-6.

[28] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, Recent advances on the discretizable molecular distance geometry problem, European J. Oper. Res. 219 (2012) 698–706, https://doi.org/10.1016/j.ejor.2011.11.007.

[29] A. Cassioli, B. Bordiaux, G. Bouvier, A. Mucherino, R. Alves, L. Liberti, M. Nilges, C. Lavor, T. Malliavin, An algorithm to enumerate all possible protein conformations verifying a set of distance constraints, BMC Bioinform. 16 (2015) 16–23, https://doi.org/10.1186/s12859-015-0451-1.

[30] T. Costa, H. Bouwmeester, W. Lodwick, C. Lavor, Calculating the possible conformations arising from uncertainty in the molecular distance geometry problem using constraint interval analysis, Inform. Sci. 415–416 (2017) 41–52, https://doi.org/10.1016/j.ins.2017.06.015.

[31] C. Dambrosio, V. Ky, C. Lavor, L. Liberti, N. Maculan, New error measures and methods for realizing protein graphs from distance data, Discrete Comput. Geom. 57 (2017) 371–418, https://doi.org/10.1007/s00454-016-9846-7.

[32] D. Gonçalves, A. Mucherino, C. Lavor, L. Liberti, Recent advances on the interval distance geometry problem, J. Global Optim. 69 (2017) 525–545, https://doi.org/10.1007/s10898-016-0493-6.

[33] C. Lavor, L. Liberti, A. Mucherino, The interval branch-and-prune algorithm for the discretizable molecular distance geometry problem with inexact distances, J. Global Optim. 56 (2013) 855–871, https://doi.org/10.1007/s10898-011-9799-6.

[34] L. Liberti, C. Lavor, Open research areas in distance geometry, in: Open Problems in Optimization and Data Analysis, in: Springer Optimization and Its Applications, Springer, Cham, 2018, pp. 183–223.

[35] M. Souza, C. Lavor, A. Muritiba, N. Maculan, Solving the molecular distance geometry problem with inaccurate distance data, BMC Bioinform. 14 (2013) S71–S76, https://doi.org/10.1186/1471-2105-14-S9-S7.

[36] B. Worley, F. Delhommel, F. Cordier, T. Malliavin, B. Bardiaux, N. Wolff, M. Nilges, C. Lavor, L. Liberti, Tuning interval branch-and-prune for protein structure determination, J. Global Optim. 72 (2018) 109–127, https://doi.org/10.1007/s10898-018-0635-0.

[37] A. Alfakih, Euclidean Distance Matrices and Their Applications in Rigidity Theory, Springer, 2018.

[38] J. Gower, Euclidean distance geometry, Math. Sci. 7 (1) (1982) 1–14.

[39] J. Gower, Properties of Euclidean and non-Euclidean distance matrices, Linear Algebra Appl. 67 (1985) 81–97, https://doi.org/10.1016/0024-3795(85)90187-9.

---

A.D. Báez Sánchez, C. Lavor / Linear Algebra and its Applications 592 (2020) 287-305

[40] A. Albert, Conditions for positive and nonnegative definiteness in terms of pseudoinverses, SIAM J. Appl. Math. 17 (2) (1969) 434–440, https://doi.org/10.1137/0117041.

[41] G. Marsaglia, G.P.H. Styan, Equalities and inequalities for ranks of matrices, Linear Multilinear Algebra 2 (1974) 269–292, https://doi.org/10.1080/03081087408817070.

[42] C. Meyer, Generalized inverses and ranks of block matrices, SIAM J. Appl. Math. 25 (1973) 597–602, https://doi.org/10.1137/0125057.

[43] L. Chua, D. Plaumann, R. Sinn, C. Vinzant, Ch. Gram spectrahedra, in: Contemporary Mathematics, vol. 697, 2017, pp. 81–105.

[44] M. Ramana, A. Goldman, Some geometric results in semidefinite programming, J. Global Optim. 7 (1995) 33–50, https://doi.org/10.1007/BF01100204.

[45] T.E. Malliavin, A. Mucherino, C. Lavor, L. Liberti, Systematic exploration of protein conformational space using a distance geometry approach, J. Chem. Inf. Model. 59 (10) (2019) 4486–4503, https://doi.org/10.1021/acs.jcim.9b00215.

[46] R. Alves, C. Lavor, Geometric algebra to model uncertainties in the discretizable molecular distance geometry problem, Adv. Appl. Clifford Algebr. 27 (2017) 439–452, https://doi.org/10.1007/s00006-016-0653-2.

[47] R. Alves, C. Lavor, C. Souza, M. Souza, Clifford algebra and discretizable distance geometry, Math. Methods Appl. Sci. 41 (2018) 3999–4346, https://doi.org/10.1002/mma.4422.

[48] C. Lavor, S. Xambó-Descamps, I. Zaplana, A Geometric Algebra Invitation to Space-Time Physics, Robotics and Molecular Geometry, SpringerBriefs in Mathematics, Springer International Publishing, 2018.

[49] C. Lavor, R. Alves, Oriented conformal geometric algebra and the molecular distance geometry problem, Adv. Appl. Clifford Algebr. 29 (2019) 1–19, https://doi.org/10.1007/s00006-018-0925-0.

[50] L. Liberti, C. Lavor, On a relationship between graph realizability and distance matrix completion, in: A. Migdalas, A. Sifaleras, K. Christos, J. Papathanasiou, E. Stiakakis (Eds.), Optimization Theory, Decision Making, and Operations Research Applications, Springer New York, New York, NY, 2013, pp. 39–48.