Discrete Applied Mathematics 267 (2019) 190-194

ELSEVIER

Contents lists available at ScienceDirect

Discrete Applied Mathematics

journal homepage: www.elsevier.com/locate/dam

SURGEON GENERAL MATHEMATICS

# On the polynomiality of finding $^K$ DMDGP re-orders

Carlile Lavor $^{a,1}$, Michael Souza $^{b,2}$, Luiz Mariano Carvalho $^{c}$, Leo Liberti $^{d,*3}$

$^{a}$ University of Campinas (IMECC-UNICAMP), 13081-970, Campinas SP, Brazil
$^{b}$ Federal University of Ceará (DEMA-UFC), 60020-181, Fortaleza CE, Brazil
$^{c}$ Rio de Janeiro State University (IME-UERJ), 20559-900, Rio de Janeiro RJ, Brazil
$^{d}$ LIX CNRS, École Polytechnique, 91128 Palaiseau, France

# ARTICLE INFO

Article history:

Received 13 July 2018

Received in revised form 10 June 2019

Accepted 19 July 2019

Available online 2 August 2019

Keywords:

Distance Geometry

Discretizable Molecular Distance Geometry

Problem

Vertex orders

# ABSTRACT

In Cassioli et al. (2015), the complexity of finding $^K$ DMDGP re-orders was stated to be NP-complete by inclusion, which fails to provide a complete picture. In this paper we show that this problem is indeed NP-complete for $K = 1$, but it is in $\mathbf{P}$ for each fixed $K \geq 2$.

© 2019 Elsevier B.V. All rights reserved.

# 1. Introduction

The fundamental problem in Distance Geometry (DG) is the DG Problem (DGP): given an integer $K &gt; 0$ and a simple, undirected, non-negatively edge-weighted graph $G = (V, E, d)$, with $d: E \to \mathbb{R}_+$, find positions in $\mathbb{R}^K$ for each vertex such that each edge, drawn as a segment, has length equal to the weight [22,23,26]. The set of positions of all the vertices in $V$ is called a realization of $G$. Many variants replace equality with inequalities to address data measurement error and noise [1,2,6,10,14,15,21,33,34]. The DGP has applications to many fields of science and engineering, including clock synchronization protocols, sensor network localization, robotics, nanostructures, and protein structure determination [3,4,9,11,18,31].

Most of the solution methods for the DGP on arbitrary graphs consist of search algorithms in continuous space [27], but if an appropriate vertex order is given, the DGP solution space becomes discrete [13,16,29].

Our motivation is based on a particular subclass of DGP graphs which arise in the modelling of proteins, where $K = 3$ [8,32]. More specifically, we look at graphs which are models of protein backbones. The edge set of these graphs, related to the pairs of atoms whose distances are known theoretically and experimentally, contains:

1.  a Hamiltonian path which is used to number the vertices from $v_{1}$ to $v_{n}$, where $n = |V|$;
2.  a clique defined on the first three vertices $v_{1}, v_{2}, v_{3}$ (in general, up to $v_{K}$);

https://doi.org/10.1016/j.dam.2019.07.021

0166-218X/© 2019 Elsevier B.V. All rights reserved.

---

C. Lavor, M. Souza, L.M. Carvalho et al. / Discrete Applied Mathematics 267 (2019) 190-194

(3) for each vertex $v_{i}$ in the Hamiltonian path, $i &gt; 3$, edges connecting $v_{i}$ to its three immediate predecessors, namely vertices $v_{i-1}, v_{i-2}, v_{i-3}$ (in general, for each $i &gt; K$, edges connecting $v_{i}$ to its $K$ immediate predecessors);
(4) possibly other edges not described by (1)-(3).

Assuming that the strict triangular inequality is satisfied for $v_{i-3}, v_{i-2}, v_{i-1}, i &gt; 3$, these graphs are rigid in $\mathbb{R}^3$ (a similar generalized statement holds for $K$ and strict simplex inequalities). The associated DGP subclass is called the Discretizable Molecular Distance Geometry Problem (DMDGP) [19,20]. In [28], it is proved that the number of its realizations is almost always a power of two, implying interesting symmetric properties in the solution space [12,24,30], where a Branch-and-Prune algorithm can be applied to find all these realizations [5,25].

We now look at an auxiliary problem: given an arbitrary graph, does its edge set satisfy the conditions (1)-(4) given above? In other words, can we find a vertex order such that (1)-(3) hold? The problem, informally known as "DMDGP order", was named Contiguous Trilateration Ordering Problem (CTOP) in [7], where it was shown (by reduction from Hamiltonian Path) that CTOP is NP-complete.

In order to address this computational complexity issue arising in protein structure determination applications of the DGP, "hand-crafted" vertex orders were proposed based on the atomic sequences in protein backbones [17,21]. The novelty introduced by such orders is that they allow the repetition of vertices in the order. More precisely, the vertex of rank $i$ in the order (for $i &gt; 3$) also has the possibility of being the same vertex having rank less than $i - 2$. Such orders are called repetition orders (or re-orders).

We generalize the fact that proteins exist in at most three spatial dimensions by replacing the number three with a general integer $K &gt; 0$. In this setting, the problem of finding a re-order relative to $K$ in an arbitrary graph is called the Re-Order Problem (ReOP). In [7, §3.3] it was stated that "because a re-order which never repeats any vertices is a $^K$ DMDGP order, the ReOP is also NP-complete for any fixed $K$ by restriction to the CTOP". This sentence is wrong: the correct sentence is "the ReOP is NP-complete by inclusion of the case when $K = 1$". In finding a re-order, since we are free to repeat vertex $v_i$ at rank $i + K$, it follows that re-orders do not need to be Hamiltonian paths. The ReOP is NP-complete for $K = 1$ because repeating vertex $v_i$ at rank $i + 1$ is equivalent to "failing to progress in the construction of the order". Therefore, for $K = 1$, the reduction from Hamiltonian Path to ReOP is trivial. This implies that ReOP is NP-complete for $K = 1$, and that it is NP-complete by inclusion of this case. For higher values of $K$, however, it turns out that re-orders can be found in polynomial time, including the case of interest for proteins $K = 3$.

## 2. Finding a $^K$ DMDGP re-order in polynomial time for $K &gt; 1$

We start by giving the precise definition of a $^K$ DMDGP re-order for a graph $G$ and the definition of an auxiliary $K$-clique incidence graph $G_K$ derived from $G$. We consider $K &gt; 1$.

**Definition 1.** Given a simple undirected graph $G = (V, E)$, a $^K$ DMDGP re-order for $G$ is a surjective function $r: \{1, \ldots, m\} \to V$, with length $m \in \mathbb{N}$ bounded by a polynomial in $|V|$, such that (we simplify the notation using $r_i$ instead of $r(i)$):

(1) $\{r_1, \ldots, r_K\}$ defines a $K$-clique in $G$;
(2) For $i \in \{K + 1, \ldots, m\}$, $\{r_{i - K + 1}, r_i\}$, $\{r_{i - K + 2}, r_i\}$, $\ldots$, $\{r_{i - 1}, r_i\} \in E$;
(3) For $i \in \{K + 1, \ldots, m\}$, $\{r_{i - K}, r_i\} \in E$ or $r_{i - K} = r_i$.

**Definition 2.** Given a simple undirected graph $G = (V, E)$, the graph $G_K = (V_K, E_K)$ is defined as follows:

(1) $u \in V_K$ iff $u$ is a set defined by the vertices of a $K$-clique in $G$;
(2) $\{u, v\} \in E_K$ iff $u \cup v$ is a set given by the vertices of a $(K + 1)$-clique in $G$.

In Fig. 1, we give an example of $G$ and $G_K$, for $K = 2$. The next theorem identifies $^K$ DMDGP re-orders for $G$ with particular walks in $G_K$.

**Theorem 1.** Given a simple undirected graph $G = (V, E)$, there is a $^K$ DMDGP re-order $r$ in $G$ if and only if there exists a walk $\gamma = (\gamma_1, \ldots, \gamma_{m - K + 1})$ in $G_K = (V_K, E_K)$ such that $\gamma_i \in V_K$ for all $i \leq m - K + 1$ and $\bigcup_{i \leq m - K + 1} \gamma_i = V$.

**Proof.** Given a $^K$ DMDGP re-order $r$ in $G$, with length $m \in \mathbb{N}$, we define $\gamma_i$, for $i \in \{1, \dots, m - K + 1\}$, in the following way:

$$
\gamma_1 = \{r_1, r_2, \dots, r_K\},
$$

$$
\gamma_2 = \{r_2, r_3, \dots, r_{K+1}\},
$$

$$
\vdots
$$

$$
\gamma_{m - K + 1} = \{r_{m - K + 1}, r_{m - K + 2}, \dots, r_m\}.
$$

---

C. Lavor, M. Souza, L.M. Carvalho et al. / Discrete Applied Mathematics 267 (2019) 190-194

![img-0.jpeg](img-0.jpeg)

![img-1.jpeg](img-1.jpeg)
Fig. 1. An example of  $G$  and  $G_{2}$ .

Since  $\gamma_{i} = \{r_{i}, r_{i+1}, \ldots, r_{i+K-1}\}$  defines a K-clique in  $G$  (by definition of  $^K$  DMDGP re-order), it follows that  $\gamma_{i} \in V_{K}$  for each  $i \in \{1, \ldots, m-K+1\}$ . Moreover, since  $\gamma_{i} \neq \gamma_{i+1}$ , by definition of  $G_{K}$  we have  $\{\gamma_{i}, \gamma_{i+1}\} \in E_{K}$ , which implies that  $\gamma$  is a walk in  $G_{K}$  such that  $\bigcup_{i} \gamma_{i} = V$  (because  $r$  is a surjection).

Conversely, let us suppose that there exists a walk  $\gamma$  of length  $p$  in  $G_{K}$  with  $\bigcup_{i\leq p}\gamma_{i} = V$ . Consider an arbitrary ordering on  $\gamma_{1}$  given by  $\gamma_{1} = (\gamma_{1}^{1},\gamma_{1}^{2},\ldots ,\gamma_{1}^{K})$ . For  $i\in \{2,\dots ,p\}$  and  $j\in \{1,\dots ,K\}$  define:

$$
\gamma_ {i} ^ {j} = \left\{ \begin{array}{l l} \gamma_ {i - 1} ^ {j} &amp; \text {i f} \gamma_ {i - 1} ^ {j} \in \gamma_ {i} \\ \gamma_ {i} ^ {*} &amp; \text {o t h e r w i s e}, \end{array} \right.
$$

where  $\gamma_{i}^{*}$  is the unique element of  $\gamma_{i}\setminus \gamma_{i - 1}$ . We claim that the function  $r:\{1,\ldots ,pK\} \to V$  defined by

$$
r _ {i} = \gamma_ {\{i / K \}} ^ {i - (\{i / K \} - 1) K}
$$

is a  $^K$  DMDGP re-order. We will show it satisfies the requirements of the definition:

(1)  $\{r_1, \ldots, r_K\}$  is a K-clique in  $G$ : by definition  $\{r_1, \ldots, r_K\} = \{\gamma_1^1, \gamma_1^2, \ldots, \gamma_1^K\} = \gamma_1$ , where  $\gamma_1$  is a K-clique in  $G$  by definition of  $G_K$ .
(2) For  $\ell \in \{K + 1,\dots ,pK\}$ ,  $\{r_{\ell -K + 1},r_{\ell}\}$ ,  $\{r_{\ell -K + 2},r_{\ell}\}$ , ...,  $\{r_{\ell -1},r_{\ell}\} \in E$ : by definition of  $r$ ,  $\{r_{\ell -K + 1},r_{\ell -K + 2},\ldots ,r_{\ell}\} \subset \gamma_{[\ell /K] - 1} \cup \gamma_{[\ell /K]}$ , for  $\ell \in \{K + 1,\dots ,pK\}$ , which implies that  $\{r_{\ell -j},r_{\ell}\} \in E$  for  $j \in \{1,\dots ,K - 1\}$  by definition of  $G_K$ .
(3) For  $\ell \in \{K + 1,\dots ,pK\}$ ,  $\{r_{\ell -K},r_{\ell}\} \in E$  or  $r_{\ell -K} = r_{\ell}$ : again by definition of  $r$  and  $G_{K}$ , for  $\ell = K + 1,\ldots ,pK$ , we have that  $r_{\ell -K} = r_{\ell}$  if  $r_{\ell}\in \gamma_{[\ell /K] - 1}\cap \gamma_{[\ell /K]}$ , and  $\{r_{\ell -K},r_{\ell}\} \in E$  if  $r_{\ell -K}\in \gamma_{[\ell /K] - 1}$  and  $r_{\ell}\in \gamma_{[\ell /K]}$ .

This concludes the proof.

The next result states the relationship between  $^K$  DMDGP re-orders in  $G$  and the connectivity of  $G_K$ .

Corollary 1. There is a  $^K$  DMDGP re-order in a simple undirected graph  $G = (V, E)$  if and only if there exists a connected component of  $G_K = (V_K, E_K)$  whose vertex union is  $V$ .

Proof. This follows immediately from Theorem 1.  $\square$

The last two results ensure that  $^K$  DMDGP re-orders can be found in polynomial time.

Lemma 1. Given a simple undirected graph  $G = (V, E)$ , the graph  $G_K = (V_K, E_K)$  can be constructed in  $O(K^3 | V|^{K+1})$  steps.

Proof. Constructing  $G_{K}$  requires constructing its edge set  $E_{K}$ . Each vertex of each  $(K + 1)$ -clique in  $G$  yields an incidence of two  $K$ -cliques in  $G$ . By definition, this means that there are  $K + 1$  edges in  $E_{K}$  for each  $(K + 1)$ -clique in  $G$ . In the worst case, the number of  $(K + 1)$ -cliques in  $G$  is  $\binom{|V|}{K + 1} = O(|V|^{K + 1})$ , which yields an  $O(K|V|^{K + 1})$  worst-case time complexity bound. Furthermore, considering the time for verifying whether a set of vertices is a clique is  $O(K^2)$ , we conclude that the worst case time complexity for constructing  $G_{K}$  is  $O(K^3|V|^{K + 1})$ .

Theorem 2. A  $^K$  DMDGP re-order for a simple undirected graph  $G = (V, E)$  can be found in  $O(|V|^{2K})$  steps.

Proof. Using DFS or BFS, we find connected components in  $O(n + m)$  steps for a graph with  $n$  vertices and  $m$  edges; the worst case is on dense graphs where  $m$  is  $O(n^2)$ .  $G_K$  has at worst  $\binom{|V|}{K}$  vertices and  $O\left(\binom{|V|}{K}^2\right)$  edges. As mentioned above,  $\binom{|V|}{K}$  is dominated by  $O(|V|^K)$ , so this step of the algorithm yields a  $O(|V|^{2K})$  time bound. Finding a spanning

---

C. Lavor, M. Souza, L.M. Carvalho et al. / Discrete Applied Mathematics 267 (2019) 190-194

tree requires $O(m + n\log n)$ steps using Prim's algorithm with Fibonacci heaps; in the worst case the $O(m) = O(n^{2})$ term dominates, so that also yields a $O(|V|^{2K})$ time bound. A walk is obtained from a tree by replacing each edge with two antiparallel arcs, and then exploring the corresponding Eulerian cycle, which can be done in $O(n) = O(|V|^{K})$ steps. Constructing $G_{K}$ from $G$ requires $O(K^3 |V|^{K + 1})$ steps by Lemma 1. The whole algorithm has a worst-case time complexity $O(K^3 |V|^{K + 1} + 2|V|^{2K} + |V|^{K}) = O(|V|^{2K})$, since $K &lt; |V|$. □

In conclusion, we have the following theorem which summarizes our discussion.

Theorem 3. Given a simple undirected graph $G = (V, E)$, the problem of finding a $^K$ DMDGP re-order is NP-complete for $K = 1$, and in $\mathbf{P}$ for each fixed $K \geq 2$.

## Acknowledgements

We thank the Brazilian research agencies CNPq and FAPESP for support. We thank A. Mucherino, D. Gonçalves, and R. Marques for useful discussions, as well as two anonymous referees for valuable suggestions. This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement n. 764759 ETN "MINOA"

## References

[1] R. Alves, C. Lavor, Geometric algebra to model uncertainties in the discretizable molecular distance geometry problem, Adv. Appl. Clifford Algebr. 27 (2017) 439-452.
[2] R. Alves, C. Lavor, C. Souza, M. Souza, Clifford algebra and discretizable distance geometry, Math. Methods Appl. Sci. 41 (2018) 3999-4346.
[3] S. Billinge, P. Duxbury, D. Gonçalves, C. Lavor, A. Mucherino, Assigned and unassigned distance geometry: applications to biological molecules and nanostructures, 40R 14 (2016) 337-376.
[4] S. Billinge, P. Duxbury, D. Gonçalves, C. Lavor, A. Mucherino, Recent results on assigned and unassigned distance geometry with applications to protein molecules and nanostructures, Ann. Oper. Res. 271 (2018) 161-203.
[5] R. Carvalho, C. Lavor, F. Protti, Extending the geometric build-up algorithm for the molecular distance geometry problem, Inform. Process. Lett. 108 (2008) 234-237.
[6] A. Cassioli, B. Bordiaux, G. Bouvier, A. Mucherino, R. Alves, L. Liberti, M. Nilges, C. Lavor, T. Malliavin, An algorithm to enumerate all possible protein conformations verifying a set of distance constraints, BMC Bioinformatics 16 (2015) 16-23.
[7] A. Cassioli, O. Gunluk, C. Lavor, L. Liberti, Discretization vertex orders in distance geometry, Discrete Appl. Math. 197 (2015) 27-41.
[8] V. Costa, A. Mucherino, C. Lavor, A. Cassioli, L. Carvalho, N. Maculan, Discretization orders for protein side chains, J. Global Optim. 60 (2014) 333-349.
[9] G. Crippen, T. Havel, Distance Geometry and Molecular Conformation, Wiley, New York, 1988.
[10] C. Dambrosio, V. Ky, C. Lavor, L. Liberti, N. Maculan, New error measures and methods for realizing protein graphs from distance data, Discrete Comput. Geom. 57 (2017) 371-418.
[11] B. Donald, Algorithms in Structural Molecular Biology, MIT Press, Boston, 2011.
[12] F. Fidalgo, D. Gonçalves, C. Lavor, L. Liberti, A. Mucherino, A symmetry-based splitting strategy for discretizable distance geometry problems, J. Global Optim. 71 (2018) 717-733.
[13] D. Gonçalves, A. Mucherino, Discretization orders and efficient computation of Cartesian coordinates for distance geometry, Optim. Lett. 8 (2014) 2111-2125.
[14] D. Gonçalves, A. Mucherino, C. Lavor, L. Liberti, Recent advances on the interval distance geometry problem, J. Global Optim. 69 (2017) 525-545.
[15] C. Lavor, R. Alves, Oriented conformal geometric algebra and the molecular distance geometry problem, Adv. Appl. Clifford Algebr. 29 (2019) 1-19.
[16] C. Lavor, J. Lee, A. Lee-St. John, L. Liberti, A. Mucherino, M. Sviridenko, Discretization orders for distance geometry problems, Optim. Lett. 6 (2012) 783-796.
[17] C. Lavor, L. Liberti, B. Donald, B. Worley, B. Bardiaux, T. Malliavin, M. Nilges, Minimal NMR distance information for rigidity of protein graphs, Discrete Appl. Math. 256 (2019) 91-104.
[18] C. Lavor, L. Liberti, W. Lodwick, T. Mendonça da Costa, An Introduction to Distance Geometry Applied to Molecular Geometry, in: SpringerBriefs, Springer, New York, 2017.
[19] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, The discretizable molecular distance geometry problem, Comput. Optim. Appl. 52 (2012) 115-146.
[20] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, Recent advances on the discretizable molecular distance geometry problem, European J. Oper. Res. 219 (2012) 698-706.
[21] C. Lavor, L. Liberti, A. Mucherino, The interval branch-and-prune algorithm for the discretizable molecular distance geometry problem with inexact distances, J. Global Optim. 56 (2013) 855-871.
[22] L. Liberti, C. Lavor, Six mathematical gems from the history of distance geometry, Int. Trans. Oper. Res. 23 (2016) 897-920.
[23] L. Liberti, C. Lavor, Euclidean Distance Geometry: An Introduction, Springer, New York, 2017.
[24] L. Liberti, C. Lavor, J. Alencar, G. Resende, Counting the Number of Solutions of  $^K$ DMDGP Instances, in: Lecture Notes in Computer Science, vol. 8085, 2013, pp. 224-230.
[25] L. Liberti, C. Lavor, N. Maculan, A branch-and-prune algorithm for the molecular distance geometry problem, Int. Trans. Oper. Res. 15 (2008) 1-17.
[26] L. Liberti, C. Lavor, N. Maculan, A. Mucherino, Euclidean distance geometry and applications, SIAM Rev. 56 (2014) 3-69.
[27] L. Liberti, C. Lavor, A. Mucherino, N. Maculan, Molecular distance geometry methods: from continuous to discrete, Int. Trans. Oper. Res. 18 (2010) 33-51.
[28] L. Liberti, B. Masson, J. Lee, C. Lavor, A. Mucherino, On the number of realizations of certain Henneberg graphs arising in protein conformation, Discrete Appl. Math. 165 (2014) 213-232.
[29] A. Mucherino, C. Lavor, L. Liberti, The discretizable distance geometry problem, Optim. Lett. 6 (2012) 1671-1686.
[30] A. Mucherino, C. Lavor, L. Liberti, Exploiting symmetry properties of the discretizable molecular distance geometry problem, J. Bioinform. Comput. Biol. 10 (2012) 1242009, (1-15).
[31] A. Mucherino, C. Lavor, L. Liberti, N. Maculan (Eds.), Distance Geometry: Theory, Methods and Applications, Springer, New York, 2013.

---

C. Lavor, M. Souza, L.M. Carvalho et al. / Discrete Applied Mathematics 267 (2019) 190-194

[32] S. Sallaume, S. Martins, L. Ochi, W. Gramacho, C. Lavor, L. Liberti, A discrete search algorithm for finding the structure of protein backbones and side chains, Int. J. Bioinform. Res. Appl. 9 (2013) 261-270.
[33] M. Souza, C. Lavor, A. Muritiba, N. Maculan, Solving the molecular distance geometry problem with inaccurate distance data, BMC Bioinformatics 14 (2013) S71-S76.
[34] M. Souza, A. Xavier, C. Lavor, N. Maculan, Hyperbolic smoothing and penalty techniques applied to molecular structure determination, Oper. Res. Lett. 39 (2011) 461-465.