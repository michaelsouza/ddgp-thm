Computational and Applied Mathematics (2021) 40:98

https://doi.org/10.1007/s40314-021-01479-6

Check for updates

# On the optimality of finding DMDGP symmetries

Carlile Lavor $^{1}$  · Andrés Oliveira $^{1}$  · Wagner Rocha $^{1}$  · Michael Souza $^{2}$

Received: 21 July 2020 / Revised: 5 March 2021 / Accepted: 8 March 2021 / Published online: 29 March 2021

© SBMAC - Sociedade Brasileira de Matemática Aplicada e Computacional 2021

# Abstract

The Discretizable Molecular Distance Geometry Problem (DMDGP) is a subclass of the Distance Geometry Problem, which aims to embed a weighted simple undirected graph in a Euclidean space, such that the distances between the points correspond to the values given by the weighted edges in the graph. The search space of the DMDGP is combinatorial, based on a total vertex order that implies symmetry properties related to partial reflections around planes defined by the Cartesian coordinates of three immediate and consecutive vertices that precede the so-called symmetry vertices. Since these symmetries allow us to know a priori the cardinality of the solution set and to calculate all the DMDGP solutions, given one of them, it would be relevant to identify these symmetries efficiently. Exploiting mathematical properties of the vertices associated with these symmetries, we present an optimal algorithm that finds such vertices.

Keywords DMDGP  $\cdot$  Branch-and-prune  $\cdot$  DMDGP symmetries  $\cdot$  Protein geometry

Mathematics Subject Classification 51K05  $\cdot$  68R05  $\cdot$  05C90

# 1 The Discretizable Distance Geometry Problem (DMDGP)

The main problem of the Euclidean Distance Geometry (Liberti et al. 2014; Liberti and Lavor 2016; Mucherino et al. 2013), called the Distance Geometry Problem (DGP), is given as follows (recent literature itens about the DGP can be found in Billinge et al. (2016), Billinge et al. (2018), Lavor et al. (2017), Liberti and Lavor (2017), Liberti and Lavor (2018), among others):

Definition 1 Given a simple undirected graph  $G = (V, E, d)$  whose edges are weighted by a function  $d: E \to [0, \infty)$  and a positive integer  $K$ , find a function  $x: V \to \mathbb{R}^K$  such that

$$
\forall \{u, v \} \in E, \| x _ {u} - x _ {v} \| = d _ {u v}, \tag {1}
$$

Springer

JBMC

---

Page 2 of 10

C. Lavor et al.

where $x_{u} = x(u)$, $x_{v} = x(v)$, $d_{uv} = d(\{u,v\})$, and $\| x_u - x_v\|$ is the Euclidean distance between $x_{u}$ and $x_{v}$.

For a general function $d: E \to [0, \infty)$, it may happen that $G = (V, E, d)$ cannot be embedded in $\mathbb{R}^K$, for some positive integer $K$, satisfying the system (1) (for example, this happens if there is a clique in $G$ such that the related triangular inequality is not satisfied). However, in DGP applications, the function $d$ is related to some existing frameworks and Eq. (1) have a solution. Hence, throughout the paper, we will assume that the DGP solution set is nonempty.

The DGP can be formulated as a global optimization problem,

$$
\min_{x_1,\ldots,x_n\in \mathbb{R}^K}\sum_{\{u,v\} \in E}\left(\| x_u - x_v\|^2 -d_{uv}^2\right)^2,
$$

where $|V| = n$, but results presented in Lavor et al. (2006) indicate that continuous optimization algorithms do not scale well for medium or large instances. A survey on different methods to the DGP is given in Liberti et al. (2010).

The DMDGP (Lavor et al. 2012a,b) is a combinatorial version of the DGP (Cassioli et al. 2015; Malliavin et al. 2019), first proposed in the context of quantum computing (Lavor et al. 2005; Marquezino et al. 2019) and related to the problem of calculating 3D protein structures ($K = 3$) using distance information provided by Nuclear Magnetic Resonance (NMR) experiments (Crippen and Havel 1988; Lavor et al. 2019; Wüthrich 1989).

The formal definition of the DMDGP (Lavor et al. 2012a) is the following:

**Definition 2** Given a DGP graph $G = (V, E, d)$ and a vertex order $v_1, \ldots, v_n$ such that

- $v_1, v_2, v_3$ can be fixed in $\mathbb{R}^3$ satisfying (1);
- $\forall i &gt; 3$, the set $\{v_{i-3}, v_{i-2}, v_{i-1}, v_i\}$ is a clique with

$$
d_{i-3,i-2} + d_{i-2,i-1} &gt; d_{i-3,i-1},
$$

find a function $x: V \to \mathbb{R}^3$ such that

$$
\forall \{u,v\} \in E, \| x_u - x_v \| = d_{uv}.
$$

The key point in the DMDGP definition is the vertex order (Lavor et al. 2012). For a general DGP graph, DMDGP orders may not exist and the problem to find them is NP-complete (Cassioli et al. 2015). However, for DGP graphs related to protein molecules, a DMDGP order can be obtained a priori (Lavor et al. 2019), based on information provided by protein geometry. Allowing some repetitions in the vertex order, the problem can be solved in polynomial time (Lavor et al. 2019).

The DMDGP order implies symmetry properties related to partial reflections around planes defined by the Cartesian coordinates of three consecutive vertices (see Fig. 1). In Fig. 1, there is a symmetry at vertex $v_9$, given by the plane defined by the positions of vertices $v_6, v_7, v_8$, generating two solutions.

Since the DMDGP symmetries allow us to know a priori the cardinality of the solution set and to calculate all the DMDGP solutions, given one of them (Liberti et al. 2013; Mucherino et al. 2012) [(which can be obtained by any algorithm proposed in the literature to solve the DMDGP (Liberti et al. 2014)], it would be relevant to identify the vertices related to these symmetries in an efficient way. This paper presents an optimal algorithm to do that.

The next section characterizes mathematically the DMDGP symmetries and Section 3 describes the main contribution of the paper. Section 4 provides the final comments.

Springer

---

On the optimality of finding DMDGP symmetries

![img-0.jpeg](img-0.jpeg)
Fig. 1 Two DMDGP solutions with symmetry at vertex  $v_{9}$

# 2 BP algorithm and DMDGP symmetries

In the DMDGP definition, the positions for the first three vertices guarantee that the solution set will contain just incongruent realizations [i.e., the DMDGP solution set does not contain solutions obtained by rotating and translating other solutions (Liberti et al. 2014)] and the strictness of the triangular inequality prevents an uncountable quantity of solutions (Lavor et al. 2012a).

The DMDGP vertex order induces that the search space can be modeled as a graph in a binary-tree shape which will be explored in a depth-first search. In the first three layers of such tree, there is only one node each, as it represents the oneness for the positions of each of the first three vertices. From the fourth vertex on, we have the representations of all the possible positions for the vertices  $v_{i}$ ,  $i = 4, \ldots, n$ , with 2 possibilities for  $v_{4}, 4$  for  $v_{5}, 8$  for  $v_{6}$ , and so on (two possible positions for  $v_{i}$  as the intersection of three spheres with centers in one of the possible choices  $x_{i-1}, x_{i-2}, x_{i-3} \in \mathbb{R}^{3}$  for the associated vertices).

By the DMDGP definition and the non-empty solution set assumption, the sphere intersection associated to any three consecutive vertices always results in two possibilities with probability one. The intersection can give just one point when the distance  $d_{i-3,i}$ ,  $i = 4, \ldots, n$ , is related to torsion angles 0 or  $\pi$  radians (angles defined by the planes given by  $v_{i-3}$ ,  $v_{i-2}$ ,  $v_{i-1}$  and  $v_{i-2}$ ,  $v_{i-1}$ ,  $v_i$ ), which occurs with probability zero.

In Carvalho et al. (2008), Liberti et al. (2008), a Branch-and-Prune (BP) method is presented for finding all incongruent solutions. The BP algorithm can be exponential in the worst case, which is consistent with the fact that the DMDGP is an NP-hard problem (Lavor et al. 2012a).

Additional distance information related to the pairs  $\{v_{j}, v_{i}\}$ ,  $j &lt; i - 3$ ,  $i = 5, \ldots, n$ , can be used to reduce the search space by pruning infeasible positions in the tree and the search ends when a path from the root of the tree to a leaf node is found by the BP algorithm, such that the positions relative to vertices in the path satisfy the DGP equations (1).

As we already mentioned in the Introduction, DMDGP symmetries are related to partial reflections around planes defined by the Cartesian coordinates of three consecutive vertices. The reflection is "partial" because just part of the structure is reflected, starting from the symmetry vertex, as illustrated in Fig. 1.

The identification of the symmetry vertices is given by the symmetry set (Mucherino et al. 2012), defined by

$$
S = \left\{v \in V: \not \exists \{u, w \} \in E \text {s u c h t h a t} u + 3 &lt;   v \leq w \text {a n d} v &gt; 3 \right\}, \tag {2}
$$

Springer

---

Page 4 of 10

C. Lavor et al.

![img-1.jpeg](img-1.jpeg)
Fig. 2 Symmetric solutions in the BP tree with  $S = \{v_4, v_5, v_8, v_9\}$

for a given DMDGP instance  $G = (V, E, d)$ , where we denote by  $u + 3$  the third vertex after  $u$ .

The intuition behind the definition of the  $S$  is that, if  $v \in S$ , the partial reflection around the plane defined by the positions of  $v - 1$ ,  $v - 2$ ,  $v - 3$  gives another solution, since there is no edge  $\{u, w\}$ ,  $u + 3 &lt; v \leq w$ , that implies infeasibility of the reflected positions. For example, consider the two possible positions in  $\mathbb{R}^3$  for  $v_4$ :  $x_4^0$  and  $x_4^1$ . Since  $v_4 \in S$ , for all DMDGP instances, for any solution found in the left subtree of the search space, having the node  $x_4^0$  as its root, there exists another one that is symmetric to the plane defined by the points  $x_1$ ,  $x_2$ ,  $x_3 \in \mathbb{R}^3$  (see Fig. 2).

Considering a small DMDGP instance given by

```latex
$V = \{v_{1},v_{2},v_{3},v_{4},v_{5},v_{6},v_{7},v_{8},v_{9},v_{10}\} ,$ $E = \{\{v_1,v_2\} ,\{v_1,v_3\} ,\{v_1,v_4\} ,$ $\{v_{2},v_{3}\} ,\{v_{2},v_{4}\} ,\{v_{2},v_{5}\} ,\{v_{2},v_{6}\} ,\{v_{2},v_{7}\} ,$ $\{v_{3},v_{4}\} ,\{v_{3},v_{5}\} ,\{v_{3},v_{6}\} ,$ $\{v_4,v_5\} ,\{v_4,v_6\} ,\{v_4,v_7\} ,$ $\{v_5,v_6\} ,\{v_5,v_7\} ,\{v_5,v_8\} ,$ $\{v_6,v_7\} ,\{v_6,v_8\} ,\{v_6,v_9\} ,\{v_6,v_{10}\}$ $\{v_7,v_8\} ,\{v_7,v_9\} ,\{v_7,v_{10}\}$ $\{v_8,v_9\} ,\{v_8,v_{10}\} ,$ $\{v_9,v_{10}\} \}$
```

we obtain

$$
S = \left\{v _ {4}, v _ {5}, v _ {8}, v _ {9} \right\}.
$$

In Liberti et al. (2013, 2014), it is proved that the cardinality of the DMDGP solution set is  $2^{|S|}$  (this occurs with probability 1, for the same reason mentioned in the third paragraph of this Section), which implies that there are  $2^4$  DMDGP solutions for the example above.

Representing a DMDGP solution by a sequence of zeros and ones, according to the positions obtained by the related sphere intersection, and supposing that the first solution found by BP is given by (remember that the positions for  $v_{1}$ ,  $v_{2}$ ,  $v_{3}$  are fixed)

$$
s _ {1} = (0, 0, 1, 0, 0, 0, 1),
$$

Springer

---

On the optimality of finding DMDGP symmetries

Page 5 of 10

the other solutions are

$s_2 = (0, 0, 1, 0, 0, 1, 0), \quad s_3 = (0, 0, 1, 0, 1, 0, 1), \quad s_4 = (0, 0, 1, 0, 1, 1, 0),$

$s_5 = (0,1,0,1,0,0,1),\quad s_6 = (0,1,0,1,0,1,0),\quad s_7 = (0,1,0,1,1,0,1),$

$s_8 = (0,1,0,1,1,1,0),\quad s_9 = (1,0,1,0,0,0,1),\quad s_{10} = (1,0,1,0,0,1,0),$

$s_{11} = (1,0,1,0,1,0,1),\quad s_{12} = (1,0,1,0,1,1,0),\quad s_{13} = (1,1,0,1,0,0,1),$

$s_{14} = (1,1,0,1,0,1,0),\quad s_{15} = (1,1,0,1,1,0,1),\quad s_{16} = (1,1,0,1,1,1,0).$

Since there is a symmetry at vertex  $v_{9}$ , the solution  $s_{2}$  is obtained reflecting the positions for vertices  $v_{9}$  and  $v_{10}$  (around the plane defined by  $v_{6}, v_{7}, v_{8}$ ) and maintaining the positions for the previous vertices (see Fig. 1). The other solutions are generated using the same idea (see Fig. 2).

In addition to the properties mentioned above, DMDGP symmetries also suggest other ways to explore the DMDGP search space (Fidalgo et al. 2018; Nucci et al. 2013).

# 3 Finding DMDGP symmetries

The edges of a DMDGP graph  $G = (V, E, d)$  can be divided in two disjoint sets,

$$
E = E _ {d} \cup E _ {p},
$$

where

$$
E _ {d} = \left\{\left\{v _ {i}, v _ {j} \right\} \in E: | j - i | \leq 3 \right\}
$$

is the branching set and

$$
E _ {p} = E \backslash E _ {d}
$$

is the pruning set (Lavor et al. 2012a).

The branching edges "model" the search space as a binary tree and the pruning edges "indicate" a feasible path. To find a solution, we go down the tree from the root until we reach a leaf node, performing all feasibility tests for the edges in  $E_{p}$  along the way.

From the definition of the set  $S$  (2), there is a method to calculate  $S$ , in  $O(|V||E_p|)$  steps that verifies if there are edges associated with  $S$  (Mucherino et al. 2012).

A more efficient method can be defined using the "complement" of  $S$ ,

$$
\bar {S} = \left\{v \in V: \exists \{u, w \} \in E \text { such that } u + 3 &lt;   v \leq w \right\}, \tag {3}
$$

implying that

$$
S \cup \bar {S} = V \backslash \{v _ {1}, v _ {2}, v _ {3} \} \tag {4}
$$

and

$$
S \cap \overline {{S}} = \emptyset .
$$

From the definition of  $\overline{S}$  (3), we see that only edges in  $E_{p}$  are used to calculate  $\overline{S}$ . So, let us define, for  $i = 1, \ldots, n - 4$ ,

$$
V _ {i} = \left\{j: \left\{v _ {i}, v _ {j} \right\} \in E _ {p} \right\}
$$

and

$$
m _ {i} = \left\{ \begin{array}{l l} \max  V _ {i} &amp; \text { if } V _ {i} \neq \emptyset , \\ 0 &amp; \text { otherwise} \end{array} \right. \tag {5}
$$

Springer

JBMC

---

Page 6 of 10

C. Lavor et al.

where $\max V_{i} = \max \{j:j\in V_{i}\}$.

Now, we define another set $\mathcal{E}_p$ associated to values $m_i\neq 0$:

1. Set $\mathcal{E}_p\coloneqq \emptyset$
2. For $i = 1, \dots, n - 4$, do
- If $m_i \neq 0$, then $\mathcal{E}_p := \mathcal{E}_p \cup \left\{\{v_i, v_{m_i}\}\right\}$.

**Proposition 1** The set $\overline{S}$ can be rewritten as

$$
\bar {S} = \bigcup_ {\left\{v _ {i}, v _ {m _ {i}} \right\} \in \mathcal {E} _ {p}} \left[ v _ {i + 4}, v _ {m _ {i}} \right], \tag {6}
$$

where

$$
[ a, b ] = \{v \in V: a \leq v \leq b \}
$$

and “$\leq$” is induced by the vertex order.

**Proof** First, we prove that

$$
\overline {{S}} \subset \bigcup_ {\{v _ {i}, v _ {m _ {i}} \} \in \mathcal {E} _ {p}} [ v _ {i + 4}, v _ {m _ {i}} ].
$$

If $v_{k} \in \overline{S}$, there exists $\{v_{i}, v_{j}\} \in E_{p}$ such that $i + 3 &lt; k \leq j$, which implies that $V_{i} = \{j : \{v_{i}, v_{j}\} \in E_{p}\}$ is not empty. Since $m_{i} = \max V_{i}$, we have $i + 4 \leq k \leq m_{i}$ and $v_{k} \in [v_{i + 4}, v_{m_{i}}]$, implying that $v_{k} \in \bigcup_{\{v_{i}, v_{m_{i}}\} \in \mathcal{E}_{p}} [v_{i + 4}, v_{m_{i}}]$.

Now, let us prove that

$$
\bigcup_ {\{v _ {i}, v _ {m _ {i}} \} \in \mathcal {E} _ {p}} [ v _ {i + 4}, v _ {m _ {i}} ] \subset \overline {{\overline {{S}}}}.
$$

If $v_{k} \in \bigcup_{\{v_{i}, v_{m_{i}}\} \in \mathcal{E}_{p}} [v_{i + 4}, v_{m_{i}}]$, there exists $i$ such that $\{v_{i}, v_{m_{i}}\} \in \mathcal{E}_{p}$ and $v_{k} \in [v_{i + 4}, v_{m_{i}}]$.

Since $\mathcal{E}_p\subset E_p$, there exists $\{v_i,v_j\} \in E_p$, with $i + 3 &lt; k\leq j$, which means that $v_{k}\in \overline{S}$.

To illustrate these new concepts, let us consider the same example given in Section 2. That is,

$$
\begin{array}{l}
E _ {d} = \left\{\left\{v _ {1}, v _ {2} \right\}, \left\{v _ {1}, v _ {3} \right\}, \left\{v _ {1}, v _ {4} \right\}, \right. \\
\left\{v _ {2}, v _ {3} \right\}, \left\{v _ {2}, v _ {4} \right\}, \left\{v _ {2}, v _ {5} \right\}, \\
\left\{v _ {3}, v _ {4} \right\}, \left\{v _ {3}, v _ {5} \right\}, \left\{v _ {3}, v _ {6} \right\}, \\
\left\{v _ {4}, v _ {5} \right\}, \left\{v _ {4}, v _ {6} \right\}, \left\{v _ {4}, v _ {7} \right\}, \\
\left\{v _ {5}, v _ {6} \right\}, \left\{v _ {5}, v _ {7} \right\}, \left\{v _ {5}, v _ {8} \right\}, \\
\left\{v _ {6}, v _ {7} \right\}, \left\{v _ {6}, v _ {8} \right\}, \left\{v _ {6}, v _ {9} \right\}, \\
\left\{v _ {7}, v _ {8} \right\}, \left\{v _ {7}, v _ {9} \right\}, \left\{v _ {7}, v _ {1 0} \right\}, \\
\left\{v _ {8}, v _ {9} \right\}, \left\{v _ {8}, v _ {1 0} \right\}, \\
\left\{v _ {9}, v _ {1 0} \right\} \},
\end{array}
$$

and

$$
E _ {p} = \left\{\left\{v _ {2}, v _ {6} \right\}, \left\{v _ {2}, v _ {7} \right\}, \left\{v _ {6}, v _ {1 0} \right\} \right\}.
$$

Springer

JBMC

---

On the optimality of finding DMDGP symmetries

Page 7 of 10

Calculating  $m_{i}$ , for  $i = 1, \ldots, 6$ , we obtain

$$
m _ {1} = 0, m _ {2} = 7, m _ {3} = 0, m _ {4} = 0, m _ {5} = 0, m _ {6} = 1 0,
$$

which implies that

$$
\mathcal {E} _ {p} = \left\{\left[ v _ {2}, v _ {7} \right], \left\{v _ {6}, v _ {1 0} \right\} \right\}
$$

and

$$
\begin{array}{l} \overline {{S}} = \bigcup_ {\left[ v _ {i}, v _ {m _ {i}} \right] \in \mathcal {E} _ {p}} \left[ v _ {i + 4}, v _ {m _ {i}} \right], \\ = \left[ v _ {6}, v _ {7} \right] \cup \left[ v _ {1 0}, v _ {1 0} \right] \\ = \left\{v _ {6}, v _ {7}, v _ {1 0} \right\}. \\ \end{array}
$$

The algorithm for finding the set  $S$  is based on the following lemmas.

Lemma 1 For a given DMDGP instance  $G = (V, E, d)$ , the set  $\mathcal{E}_p$  is calculated in  $O(|V| + |E_p|)$  steps.

Proof Note that we can initialize  $m_{i} = 0$ , for  $i = 1, \dots, n - 4$ , and in one sweep of the edges in  $E_{p}$ , we update the values for  $m_{i}$ . After that, we end up with the final values for  $m_{i}$  and the set  $\mathcal{E}_p$  in  $O(|V| + |E_p|)$  steps.

Lemma 2 For a given DMDGP instance  $G = (V, E, d)$  and  $\mathcal{E}_p$ , the set  $\overline{S}$  is calculated in  $O(|V|)$  steps.

Proof This follows from (6) and from the fact that if  $\{v_{i}, v_{m_{i}}\}, \{v_{j}, v_{m_{j}}\} \in \mathcal{E}_{p}$ , with  $v_{i} &lt; v_{j}$  and  $v_{j+3} \leq v_{m_{i}}$ , then

$$
\left[ v _ {i + 4}, v _ {m _ {i}} \right] \cup \left[ v _ {j + 4}, v _ {m _ {j}} \right] = \left[ v _ {i + 4}, \max  \left\{v _ {m _ {i}}, v _ {m _ {j}} \right\} \right].
$$

Lemma 3 For a given DMDGP instance  $G = (V, E, d)$ ,  $\mathcal{E}_p$  and  $\overline{S}$ , the set  $S$  is calculated in  $O(|V|)$  steps.

Proof As the set  $\overline{S}$  is ordered, this follows directly from (4).

A pseudo-code of the method is given in Algorithm 1. Since it is based on the calculation of the sets  $\mathcal{E}_p$ ,  $\overline{S}$  and  $S$ , in this order, the next result follows directly from the above lemmas.

Theorem 1 Given a DMDGP instance  $G = (V, E, d)$ , where  $E = E_d \cup E_p$ , the associated symmetry set  $S$  can be found in  $O(|V| + |E_p|)$  steps.

Corollary 1 The Algorithm 1 for calculating the symmetry set of a DMDGP instance is optimal.

Proof From the definition of the set  $E_{d}$ , we obtain that

$$
| E | = | E _ {d} | + | E _ {p} | = \Omega (| V | + | E _ {p} |)
$$

and, from Theorem 1, we get

$$
| E | + | S | = \Omega (| V | + | E _ {p} |) + O (| V | + | E _ {p} |).
$$

In order to represent a DMDGP graph  $G = (V, E, d)$  and find its symmetry set  $S$ , we have at least to read the set  $E$  and give as output the set  $S$ . Thus, to calculate  $S$ , we need  $\Omega(|E| + |S|)$  steps, which can be computed in  $\Omega(|V| + |E_p|)$  steps.

Springer

JBMC

---

Page 8 of 10

C. Lavor et al.

Algorithm 1: Optimal way to identify the set  $S$ .
1 Function OptimalWayToIdentifyS(V, Ep, K) //  $K \in \{2, 3, \ldots\}$
2  $\mathcal{E}_p \gets \text{Identify}\mathcal{E}_p(E_p, |V|, K)$ ;
3  $\overline{S} \gets \text{Identify}\overline{S}(\mathcal{E}_p, K)$ ;
4  $S \gets (V \setminus \{v_1, v_2, v_3, \ldots, v_K\}) \setminus \overline{S}$ ;
5 return S;
6 end
7 Function Identify $\mathcal{E}_p(E_p, n, K)$  //  $\{v_i, v_j\} \in E_p$  such that  $i &lt; j$
8  $M \gets \{m_1, m_2, \ldots, m_{n-(K+1)}\} \mid m_i = 0 \forall i$ ;
9 foreach  $\{v_i, v_j\} \in E_p$  do
10  $|m_i \gets \max \{m_i, j\}$ ;
11 end
12  $\mathcal{E}_p \gets \emptyset$ ;
13 foreach  $m_i \in M$  do
14 if  $(m_i \neq 0)$  then
15  $\mathcal{E}_p \gets \mathcal{E}_p \cup \{\{v_i, v_{m_i}\}\}$ ;
16 end
17 end
18 return  $\mathcal{E}_p$ ;
19 end
20 Function Identify $\overline{S}(\mathcal{E}_p, K)$  //  $\{v_i, v_{m_i}\} \in \mathcal{E}_p$  with  $i &lt; m_i$  and
21  $\{v_i, v_{m_i}\}, \{v_j, v_{m_j}\} \in \mathcal{E}_p$  then  $v_i &lt; v_j \Leftrightarrow i &lt; j$
22 if  $\mathcal{E}_p \neq \emptyset$  then
23  $i_1 \gets \min \{i : \{v_i, v_{m_i}\} \in \mathcal{E}_p\}$ ;
24  $\ell \gets i_1; r \gets m_i$ ;
25 foreach  $\{v_i, v_{m_i}\} \in \mathcal{E}_p$  do
26 if  $(r \geq i + K)$  then
27  $|r \gets \max\{m_i, r\}$ ;
28 else
29  $\overline{S} \gets \overline{S} \cup [v_{\ell+(K+1)}, v_r]$ ;
30  $\ell \gets i; r \gets m_i$ ;
31 end
32 end
33  $\overline{S} \gets \overline{S} \cup [v_{\ell+(K+1)}, v_r]$ ;
34 end
35 return  $\overline{S}$ ;
36 end

# 4 Final comments

We proposed an optimal algorithm for finding symmetry sets associated with DMDGP instances. We were motivated for distance geometry problems related to molecular geometry  $(K = 3)$ , but note that the dimension  $K$  of the realization space is also given as input, which means that the algorithm solves  $^K\mathrm{DMDGP}$  instances (Liberti et al. 2013).

Acknowledgements We would like to thank the Brazilian research agencies, CNPq and FAPESP, and the careful reading and important comments made by the reviewers.

Springer

JBMC

---

On the optimality of finding DMDGP symmetries

Page 9 of 10

# Declarations

Conflict of interest The authors declare that they have no conflict of interest.

# References

Billinge S, Duxbury P, Gonçalves D, Lavor C, Mucherino A (2016) Assigned and unassigned distance geometry: applications to biological molecules and nanostructures. 4OR 14:337-376
Billinge S, Duxbury P, Gonçalves D, Lavor C, Mucherino A (2018) Recent results on assigned and unassigned distance geometry with applications to protein molecules and nanostructures. Ann Oper Res 271:161-203
Carvalho R, Lavor C, Protti F (2008) Extending the geometric build-up algorithm for the molecular distance geometry problem. Inf Process Lett 108:234-237
Cassioli A, Gunluk O, Lavor C, Liberti L (2015) Discretization vertex orders in distance geometry. Discret Appl Math 197:27-41
Cassioli A, Bordiaux B, Bouvier G, Mucherino A, Alves R, Liberti L, Nilges M, Lavor C, Malliavin T (2015) An algorithm to enumerate all possible protein conformations verifying a set of distance constraints. BMC Bioinform 16:16-23
Crippen G, Havel T (1988) Distance geometry and molecular conformation. Wiley, Oxford
Fidalgo F, Gonçalves D, Lavor C, Liberti L, Mucherino A (2018) A symmetry-based splitting strategy for discretizable distance geometry problems. J Global Optim 71:717-733
Lavor C, Liberti L, Maculan N (2005) Grover's algorithm applied to the molecular distance geometry problem. In: Proceedings of the 7th Brazilian congress of neural networks
Lavor C, Liberti L, Maculan N (2006) Computational experience with the molecular distance geometry problem. In: Pintér J (ed) Global optimization: scientific and engineering case studies. Springer, Berlin, pp 213-225
Lavor C, Lee J, Lee-St JA, Liberti L, Mucherino A, Sviridenko M (2012) Discretization orders for distance geometry problems. Optim Lett 6:783-796
Lavor C, Liberti L, Maculan N, Mucherino A (2012) The discretizable molecular distance geometry problem. Comput Optim Appl 52:115-146
Lavor C, Liberti L, Maculan N, Mucherino A (2012) Recent advances on the discretizable molecular distance geometry problem. Eur J Oper Res 219:698-706
Lavor C, Liberti L, Lodwick W, da Costa TM (2017) An introduction to distance geometry applied to molecular geometry. Springer, Berlin
Lavor C, Liberti L, Donald B, Worley B, Bardiaux B, Malliavin T, Nilges M (2019) Minimal NMR distance information for rigidity of protein graphs. Discret Appl Math 256:91-104
Lavor C, Souza M, Mariano L, Liberti L (2019) On the polinomiality of finding  $^K$  DMDGP re-orders. Discret Appl Math 267:190-194
Liberti L, Lavor C (2016) Six mathematical gems from the history of distance geometry. Int Trans Oper Res 23:897-920
Liberti L, Lavor C (2017) Euclidean distance geometry: an introduction. Springer, Berlin
Liberti L, Lavor C (2018) Open research areas in distance geometry. In: Migalas A, Pardalos P (eds) Open problems in optimization and data analysis. Springer, Berlin, pp 183-223
Liberti L, Lavor C, Maculan N (2008) A branch-and-prune algorithm for the molecular distance geometry problem. Int Trans Oper Res 15:1-17
Liberti L, Lavor C, Mucherino A, Maculan N (2010) Molecular distance geometry methods: from continuous to discrete. Int Trans Oper Res 18:33-51
Liberti L, Lavor C, Alencar J, Resende G (2013) Counting the number of solutions of  $^K$  DMDGP instances. Lect Notes Comput Sci 8085:224-230
Liberti L, Lavor C, Maculan N, Mucherino A (2014) Euclidean distance geometry and applications. SIAM Rev 56:3-69
Liberti L, Masson B, Lee J, Lavor C, Mucherino A (2014) On the number of realizations of certain Henneberg graphs arising in protein conformation. Discret Appl Math 165:213-232
Malliavin T, Mucherino A, Lavor C, Liberti L (2019) Systematic exploration of protein conformational space using a distance geometry approach. J Chem Inf Model 59:4486-4503
Marquezino F, Portugal R, Lavor C (2019) A primer on quantum computing. Springer, Berlin
Mucherino A, Lavor C, Liberti L (2012) Exploiting symmetry properties of the discretizable molecular distance geometry problem. J Bioinform Comput Biol 10:1242009

Springer

JB/AC

---

Page 10 of 10

C. Lavor et al.

Mucherino A, Lavor C, Liberti L, Maculan N (eds) (2013) Distance geometry: theory, methods, and applications. Springer, Berlin
Nucci P, Nogueira L, Lavor C (2013) Solving the discretizable molecular distance geometry problem by multiple realization trees. In: Mucherino A, Lavor C, Liberti L, Maculan N (eds) Distance geometry: theory, methods, and applications. Springer, New York, pp 161-176
Wüthrich K (1989) Protein structure determination in solution by nuclear magnetic resonance spectroscopy. Science 243:45-50
Publisher's Note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.

Springer