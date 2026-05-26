J Glob Optim (2014) 60:333-349

DOI 10.1007/s10898-013-0135-1

# Discretization orders for protein side chains

Virginia Costa · Antonio Mucherino · Carlile Lavor · Andrea Cassioli · Luiz M. Carvalho · Nelson Maculan

Received: 15 January 2013 / Accepted: 18 December 2013 / Published online: 8 January 2014

© Springer Science+Business Media New York 2014

Abstract Proteins are important molecules that are widely studied in biology. Since their three-dimensional conformations can give clues about their function, an optimal methodology for the identification of such conformations has been researched for many years. Experiments of Nuclear Magnetic Resonance (NMR) are able to estimate distances between some pairs of atoms forming the protein, and the problem of identifying the possible conformations satisfying the available distance constraints is known in the scientific literature as the Molecular Distance Geometry Problem (MDGP). When some particular assumptions are satisfied, MDGP instances can be discretized, and solved by employing an ad-hoc algorithm, named the interval Branch &amp; Prune. When dealing with molecules such as proteins, whose chemical structure is known, a priori information can be exploited for generating atomic orderings that allow for the discretization. In previous publications, we presented a handcrafted order for the protein backbones. In this work, we propose 20 new orders for the 20 side chains that

V. Costa · N. Maculan

COPPE, Federal University of Rio de Janeiro, Rio de Janeiro, RJ, Brazil

e-mail: virscosta@cos.ufrj.br

N. Maculan

e-mail: maculan@cos.ufrj.br

A. Mucherino (☑)

IRISA, University of Rennes 1, Rennes, France

e-mail: Antonio.Mucherino@irisa.fr

C. Lavor

IMECC-UNICAMP, University of Campinas, Campinas, SP, Brazil

e-mail: clavor@ime.unicamp.br

A. Cassioli

LIX, École Polytechnique, Palaiseau, France

e-mail: cassioli@lix.polytechnique.fr

L. M. Carvalho

IME, State University of Rio de Janeiro, Rio de Janeiro, RJ, Brazil

e-mail: luizmc@ime.uerj.br

Springer

---

can be present in proteins. Computational experiments on artificial and real instances from NMR show the usefulness of the proposed orders.

## Keywords

Discretization orders Distance geometry Branch-and-prune Combinatorial Optimization Proteins

## Introduction

The Molecular Distance Geometry Problem (MDGP) consists in finding the suitable conformations for a certain molecule which satisfy a set of constraints based on some distances between pairs of its atoms [4,14]. When the distance information is given through a list of lower and upper bounds on the distances, i.e. by a list of real intervals, the problem is also referred to as interval MDGP (iMDGP) [12]. Over the years, its solution has been attempted by formulating global optimization problems in continuous spaces (see for example [13,21,22]), where a penalty function is generally employed in order to measure the violation of the distance constraints for given molecular conformations.

In this work, we consider a particular subset of problems belonging to the iMDGP class: the interval Discretizable Molecular Distance Geometry Problem (iDMDGP) [11]. The word “Discretizable” refers to the fact that the search conformational space can be reduced to a discrete set. More precisely, the iDMDGP is a combinatorial problem, whose search domain is a tree. This problem is NP-hard [10,19]. When the discretization is possible (the instance belongs to the iDMDGP class), we can employ an exact algorithm for its solution: the interval BP (iBP) algorithm [11]. This algorithm is potentially able to enumerate the entire solution set.

Approximations of the distances between pairs of atoms of a molecule can be found by applying experiments of Nuclear Magnetic Resonance (NMR) [7]. This technique can provide, in general, upper bounds for some pairwise distances between atoms, while the corresponding lower bounds can be identified from the so-called van der Waals radii of the considered atoms: if they are not chemically bonded, they should not be closer than a certain threshold. In most cases, a list of intervals related to pairs of hydrogen atoms can be defined.

Nowadays, NMR is the second most used technique for the identification of protein conformations, which are very important molecules performing several fundamental functions in living beings. Due to this great interest, three-dimensional conformations of proteins found by the scientific community are generally stored in a web database named the Protein Data Bank (PDB) [1]. To date, however, in most cases, distance geometry problems related to NMR experiments are formulated as continuous global optimization problems, whose solution is generally attempted by a method whose optimization core is represented by the meta-heuristic Simulated Annealing (SA) [16].

The main discretization assumption is related to the existence of an atomic order (that can be different from the natural ordering of the atoms in a molecule) such that, for each atom v, there are at least three preceding reference atoms u_{1}, u_{2} and u_{3} with known distance from v. This assumption gives us the possibility to compute some atomic positions for v by intersecting three Euclidean objects in the three-dimensional space. When the distance between two atoms, say u_{1} and v, is precise, the considered Euclidean object is a sphere centered in u_{1} and having radius d(u_{1}, v). If this distance is instead imprecise, i.e. it is represented by an interval [d(u_{1}, v), d̂(u_{1}, v)], then the Euclidean object is a spherical shell centered in u_{1} having minimum radius d̂(u_{1}, v) and maximum radius d̂(u_{1}, v). When, for a certain v, two reference distances are precise and only one is represented by an interval,

---

J Glob Optim (2014) 60:333-349

then a subset of atomic positions for $v$ can be computed by intersecting two spheres with one spherical shell. This intersection generally provides at most two disjoint curves in the three-dimensional space, from which sample positions can be taken. When the three distances are all exact, then the spherical shell collapses into a sphere, and the curves degenerate and define single points [10].

Since $i$ DMDGP instances could be non-discretizable in the order given to their atoms, other different orders may need to be searched for performing the discretization and applying the $i$ BP algorithm. Information about bond lengths and bond angles can be exploited for adding distance constraints in $i$ MDGP instances. In [11], some of the authors of this paper presented an artificial ordering for the atoms of protein backbones, where only information about bond lengths and angles was exploited. We point out that, as a consequence, this order is independent from NMR data.

In this paper, we present 20 new orderings for the side chains of the 20 amino acids that can be involved in the protein synthesis. All these orders satisfy the assumptions for the $i$ DMDGP, and they are based on the information about bond lengths and angles. As they are independent from NMR data, these orders can be reused for discretizing every instance concerning protein conformations: single orders for each amino acid can be attached to each other as the protein chain does. Some preliminary studies on side chain discretization orders were presented in [2,3,18].

The rest of the paper is organized as follows. In Sect. 2, we briefly discuss the discretization assumptions and present the $i$ BP algorithm. In Sect. 3, the new discretization orders for the 20 side chains are proposed. Section 4 gives some implementation detail, while Sect. 5 presents some experiments. Finally, Sect. 6 ends the paper with conclusions and work in progress.

## 2 The discretization and the $i$ BP algorithm

Let $G = (V, E, d)$ be a weighted undirected graph. In this paper, with a little abuse of notations, we will refer to $i \in V$ as a vertex of a graph $G$, as well as an atom of the corresponding molecule. We suppose that $d$ is the function $d : (i, j) \in E \longrightarrow d_{ij} \in \mathbb{R}_+$, which associates, to the edge $(i, j)$, the known distance $d_{ij}$ between two atoms $i$ and $j$. We remark that the distance information may not be precise, so that the distance $d_{ij}$ could actually be represented by a suitable interval $[l_{ij}, u_{ij}]$.

**Definition 2.1** The graph $G$ represents an instance of the $i$ DMDGP if there exists a “discretization order” for its vertices, i.e. if there exists an order for which the following assumptions are satisfied:

- $V_0 = (1,2,3) \subset V$ is a clique, and $\forall (i,j) \in E[V_0], l_{ij} = u_{ij}$;
- $\forall i \in V : i &gt; 3$, we have:

- $\{(i - 3, i), (i - 2, i), (i - 1, i)\} \subset E$,
- $l_{i-2,i} = u_{i-2,i}$ and $l_{i-1,i} = u_{i-1,i}$,
- $l_{i-2,i} &lt; l_{i-2,i-1} + l_{i-1,i}$.

Since the first clique only contains exact distances, we can define the coordinate system for the molecular conformation by fixing the coordinates of its first 3 atoms. This way, we avoid considering solutions that can be obtained from other solutions by translations and/or rotations. Moreover, the assumptions ensure that, for each atom $i$ having rank greater than 3, there are at least three reference distances, that are related to the three edges $(i - 3,i)$

Springer

---

J Glob Optim (2014) 60:333-349

$(i - 2, i)$, $(i - 1, i)$. As a consequence, the possible positions for the atom $i$ can be identified by intersecting, in the three-dimensional space, the three Euclidean objects that are associated to these three edges. Since the distances related to $(i - 2, i)$ and $(i - 1, i)$ are always exact, the Euclidean object related to each of these edges is a sphere. This may not be true however for the distance related to the edge $(i - 3, i)$ and, in such a case, the corresponding Euclidean object is a spherical shell. Because of the strict triangular inequality, the intersection of two spheres with one spherical shell gives at most 2 disjoint curves, with probability 1. In order to discretize, a finite number $D$ of points on these two curves can be chosen (see [11] for more details about the discretization of the curves).

For a given graph $G$, the edge set $E$ can be divided in two main parts. The subset $E_{d} = \{(j,i) \in E : i - 3 \leq j &lt; i\}$ contains all the edges that are related to the discretization distances, i.e. to the distances that are required by the $i$ DMDGP assumptions and that are involved in the intersections. These are actually the distances that allow us to discretize the problem. The complement $E_{p} = E \setminus E_{d}$ contains instead some additional distances that can be used for verifying the feasibility of computed atomic positions (obtained with the intersections). We refer to such distances as pruning distances. In our orders, pruning distances are supposed to be obtained by NMR experiments; they generally refer to pairs of hydrogen atoms and their upper bounds are generally smaller than 4–5 Å.

By considering the edges in $E_{d}$, a tree, which represents the search domain for the $i$ DMDGP, can be defined. This tree provides, layer by layer, the possible positions for given atoms of the molecule. The total number of branches this tree has can be huge for medium-sized molecules, but the pruning distances in $E_{p}$ can help in identifying infeasible branches that can be removed from the tree. In this way, an exhaustive search on the remaining branches may be possible. In such a case, all solutions to the $i$ DMDGP could be enumerated. This is an important consequence of the discretization.

The $i$ BP algorithm recursively calls itself for exploring new tree branches starting from the current node of the tree. By considering the available discretization distances, candidate atomic positions for the current atom are computed by performing the above mentioned intersections. The two distances $d_{i-1,i}$ and $d_{i-2,i}$ are always exact by hypothesis, and therefore only the edge $(i-3,i)$ may be associated to either an exact, or to an interval distance. We also remark that the distance $d_{i-3,i}$ may be zero: in this case, the two vertices $i-3$ and $i$ make actually reference to the same atom. In our orders, in fact, we allow for vertex repetitions, because they make it possible to have the necessary distance information near the atom that needs such an information [11]. The distances $d_{i-1,i}$ and $d_{i-2,i}$ can never be zero, because this would go against the strict triangular inequalities. There are three situations that we can occur while performing the intersections:

1. $l_{i-3,i} = u_{i-3,i} = 0 \Longrightarrow$ the current atom is a repeated atom: since it already appeared beforehand in the order on $V$, it can only be placed in the same position as its previous copy; there is no branching on the tree;
2. $l_{i-3,i} = u_{i-3,i} \neq 0 \Longrightarrow$ the three reference distances are all exact: there are generally only two candidate positions that can be obtained by intersecting the three spheres;
3. $l_{i-3,i} \neq u_{i-3,i} \Longrightarrow$ the candidate positions lie on two disjoint curves of the three-dimensional space: $D$ equidistant positions on each curve are chosen in order to discretize.

A sketch of the $i$ BP algorithm is given in Algorithm 1. It essentially requires 4 input arguments: the index $j$ of the current atom (in the discretization order $r$) to be placed, the set of distances $d$ (which can be either exact or represented by intervals), and the number $D$ of sample distances used for discretizing interval distances. More details about this algorithm can be found in [11].

Springer

---

J Glob Optim (2014) 60:333-349

|  Algorithm 1 The iBP algorithm.  |   |
| --- | --- |
|  1: | iBP(j,r,d,D)  |
|  2: | if (rj is a repeated atom) then  |
|  3: | iBP(j+1,r,d,D);  |
|  4: | else  |
|  5: | if (d(rj-3,rj) is exact) then  |
|  6: | b=2;  |
|  7: | else  |
|  8: | b=2D;  |
|  9: | end if  |
|  10: | for k∈{1,...,b} do  |
|  11: | compute the k-th atomic position xRjfor the rj-th atom;  |
|  12: | check the feasibility of position xRjusing pruning devices;  |
|  13: | if (xRjis feasible) then  |
|  14: | if (j=|r|) then  |
|  15: | a solution x is found, print it;  |
|  16: | else  |
|  17: | iBP(j+1,r,d,D);  |
|  18: | end if  |
|  19: | end if  |
|  20: | end for  |
|  21: | end if  |

We remark that, rather than explicitly computing two curves in the three-dimensional space, every time the reference distance  $(i - 3,i)$  is an interval, we select  $D$  sample distances  $\hat{d}$  from the interval  $[l_{i - 3,i},u_{i - 3,i}]$  and we compute  $2\times D$  atomic positions by performing several sphere intersections. For each chosen distance  $\hat{d}$ , we define the sphere centered in  $x_{i - 3}$  and with radius  $\hat{d}$ , and we intersect it with other two spheres  $S_{i - 1}$  and  $S_{i - 2}$ , related to the edges  $(i - 2,i)$  and  $(i - 1,i)$ . For each  $\hat{d}$ , two atomic positions are computed.

As already mentioned, fundamental in  $i\mathrm{BP}$  is its pruning phase: as soon as a new candidate position is computed, its feasibility can be verified by employing some pruning devices. The simplest (and probably most efficient) pruning device is the one that exploits the distances in  $E_{p}$ : if one of such distances is not satisfied by one of the candidate positions (for a given tolerance  $\varepsilon_{d}$ ), then the current branch can be pruned. We remark that, in case this information is not already included in the available data, the lower bounds for all known interval distances can be refined (increased) by exploiting the van der Waals (vdW) radii associated to each atom type. For each pair of nonbonded atoms  $(i,j)$ , their relative distance  $d_{ij}$  must exceed the sum of their respective vdW radii  $r_i^v, r_j^v$ , i.e:

$$
d _ {i j} \geq r _ {i} ^ {v} + r _ {j} ^ {v}. \tag {1}
$$

In the practice, this bound is generally relaxed by taking a fraction (we will consider the  $85\%$ ) of each radius. We refer to [8,17] for the vdW radius values.

Other pruning devices can be conceived and included in the iBP algorithm [15]. For example, the secondary structure of a protein gives an important structural information of the molecule. Secondary structures are local organizations of the amino acids forming the protein into a regular structure. In proteins,  $\alpha$ -helices and  $\beta$ -sheets can be found. The former is characterized by a sequence of amino acids that are folded in a regular helix whose pace and radius depend on a strong hydrogen bond between a hydrogen atom and an oxygen atom of two non-consecutive amino acids. Moreover, some torsion angles defined on the protein backbone must attain values belonging to predetermined ranges. For more details

Springer

---

J Glob Optim (2014) 60:333-349

about protein secondary structures, the reader is referred to [6]; for more details about the implementation of secondary structure-based pruning device, the reader is referred to [15].

Conservative values for the backbone torsion angles can be obtained from standard Ramachandran plots [20] and NMR measures [5]. Improper angles (as for example those forcing the existence of the peptide plane) can also be exploited for verifying the feasibility of computed atomic positions. Finally, the protein chirality [4] can also be exploited for pruning purposes.

Even if there are various pruning devices that might be included in  $i\mathrm{BP}$ , in this paper our main aim is to improve the efficiency of the pruning device based on the verification of the distances in  $E_{p}$ , by including the information about protein side chains. When this is the case, NMR distances between pairs of hydrogens belonging to side chains are included in  $E_{p}$ . As it is well known, the subset of atoms which is common for each amino acid takes part to the construction of the backbone of the protein, while the side chains are subsets of atoms, varying amino acid per amino acid, that are chemically bonded to the backbone. Some of these side chains are hydrophobic, i.e. they do not react with the protein solvent (water), and they tend therefore to be in the interior of the molecule. This concentration of side chains generates the so-called hydrophobic core of the protein. Since it is compact, distances between several hydrogen atoms in these side chains are likely to be detected by NMR, and therefore they are likely to be included in the subset  $E_{p}$ .

# 3 Discretization orders for side chains

A handcrafted order for the discretization of protein backbones was previously proposed in [11]. This order uses information about the chemical structure of the protein for allowing for the discretization, while NMR data are supposed to be exploited for pruning purposes.  $r_{\mathrm{PB}}$  is the order for a backbone formed by  $p$  amino acids (a drawing depicting the order for  $p = 3$  is in Fig. 1):

![img-0.jpeg](img-0.jpeg)

![img-1.jpeg](img-1.jpeg)
Fig. 1 The handcrafted order for a small protein backbone containing 3 amino acids

Springer

---

J Glob Optim (2014) 60:333-349

In $r_{\mathrm{PB}}$, the superscripts identify the amino acid label where each atom belongs to. This order specifies a particular atomic sequence for the first, the second and the last amino acid, as well as for the generic amino acid, the one having $i$, with $2 &lt; i &lt; p$, as a superscript. We remark that $H_N^0$ represents the second hydrogen atom bonded to the nitrogen of the first amino acid, while $O^{p+1}$ represents the second oxygen bonded to carbon $C$ of the last amino acid.

The theoretical concept related to this backbone order is the repetition order (re-order) [11]. A re-order builds a larger graph whose structure is derived from the original graph $G$ associated to an instance with the aim of satisfying the discretization assumptions.

Let $G = (V, E, d)$ be a graph associated to an instance of the $i$ DMDGP. The set of edges $E$, already partitioned in $E_d$ and $E_p$ in Sect. 2, can also be partitioned into those edges $(i, j) \in E'$ for which $d_{ij}$ is a real nonnegative number, and those edges $(i, j) \in E''$ for which $d_{ij}$ is represented by a suitable interval. Let $V' = V \cup \{0\}$.

Definition 3.1 A repetition order (re-order) is a sequence $r: \mathbb{N} \to V'$ with length $|r| \in \mathbb{N}$ (for which $r_i = 0$ for all $i &gt; |r|$) such that:

- $G[\{r_1, r_2, r_3\}]$ is a clique
- for all $i \in \{4, \ldots, |r|\}$ the sets $\{r_{i-2}, r_i\}$, $\{r_{i-1}, r_i\}$ are edges in $E'$;
- for all $i \in \{4, \ldots, |r|\}$ the set $\{r_{i-3}, r_i\}$ is either a singleton (i.e. $r_{i-3} = r_i$) or an edge in $E' \cup E''$.

Re-orders are generally discretization orders, because they can satisfy the assumptions of the $i$ DMDGP. The only discretization assumption that may not be satisfied by re-orders is actually the one related to the strict triangular inequality (see Definition 2.1). However, we remark that the subset of re-orders for which the strict triangular inequalities are not satisfied has Lebesgue measure equal to zero [9]. As a consequence, we can state that, with probability 1, all re-orders satisfy the discretization assumptions in Definition 2.1.

Table 1 contains order sets $r_{SC}$ for the 20 amino acids involved in the protein synthesis (SC refers to the well-known 3-letter labels for the 20 amino acids), including their side chains. We provide the orders for a given amino acid (backbone and side chain) by supposing that it is the first one in the amino acid sequence (there are always two hydrogens bonded to the nitrogen $N$). In case the amino acid is placed somewhere else in the sequence, the given ordering can be simply adapted (the order in Fig. 1 can be used as a reference). Since all the orders refer to the first amino acid in the protein, we use the superscripts for discriminating among atoms of the same type that appear more than once in the same amino acid. For example, in $r_{ARG}$ (Table 1), there are four nitrogens atoms $N^1$, $N^2$, $N^3$ and $N^4$. For each carbon, besides the indices, we also provide the commonly used Greek letter. Each hydrogen is labeled with the same Greek letter as for the corresponding carbon atom. Moreover, all hydrogens have a second index that gives information about the atom that is bonded to this hydrogen. In the sequence $r_{ARG}$, for example, there is the atom $H_N^{3,2}$: the two superscripts indicate that this is the second hydrogen bonded to the third nitrogen $N^3$. In Figs. 2, 3, and 4, we present the proposed discretization orders in a graphical representation. We have the following theorem.

Theorem 3.2 If edges related to pairs of atoms separated by at most 2 covalent bonds correspond to exact distances, then all $r_{SC}$ are re-orders.

Proof In the hypothesis of the theorem, all distances between pairs of bonded atoms (1 covalent bond), as well as all distances between atoms separated by only 2 covalent bonds,

Springer

---

J Glob Optim (2014) 60:333-349

Table 1 The new orders for the 20 amino acids involved in the protein synthesis

|  SC | Re-order  |
| --- | --- |
|  GLY | N1,H1,N1,H1,2,C1,N1,H1,1,C1,H1,2,C1,C1  |
|  ALA | N1,H1,N1,H1,2,C1,N1,H1,1,C1,H1,2,C1,H1,3,C1,C1,C1  |
|  PRO | N1,H1,N1,C1,H1,1,H1,2,C1,C1,C1,H1,2,C1,C1,H1,2,C1,C1,H1,2,C1,C1,H1,2,C1,C1,C1  |
|  LEU | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,H1,1,C1,C1,H1,1,C1,C1,H1,2,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  VAL | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,H1,1,C1,C1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  ILE | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,H1,1,C1,C1,C1,S1,H1,1,C1,C1,S1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  MET | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,H1,1,C1,C1,S1,H1,1,C1,C1,S1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  PHE | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  TRP | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  GLN | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  TYR | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  ASN | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  SER | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  CYS | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  THR | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  ARG | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  HIS | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  LYS | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  GLU | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |
|  ASP | N1,H1,N1,H1,2,C1,N1,H1,1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1  |

Springer

---

J Glob Optim (2014) 60:333-349

![img-2.jpeg](img-2.jpeg)
Fig. 2 Handcrafted orders for the amino acids LYS, GLU, ASP

are exact. This introduces in the instances at least 2 exact distances per atom. Moreover, torsion angles can be defined on quadruplets of bonded consecutive atoms: the minimal and maximal extensions for such torsion angles allow us to define an interval for the distance between the first and the last atom of the quadruplet (by the hypothesis, all other distances in the quadruplet are exact). This introduces interval distances in our instances. By inspection of all provided orders and by considering these exact and interval distances, the thesis of the theorem can be verified.

# 4 Generation of the proposed orders

The orders for the side chains presented in Sect. 3 only exploit a priori information derived from the chemical structure of the amino acids. For this reason,  $i$  DMDGP instances can be generated in two main steps. First, from the amino acid sequence of the protein, the corresponding orders in Sect. 3 can be combined in order to represent the whole molecular structure. The instance generated in this way corresponds to an  $i$  DMDGP instance without pruning distances based on NMR. The second step consists therefore in adding NMR distances to the instance, for allowing the selection of the branches of the tree where there are feasible solutions. Note that some of the intervals given by NMR distances may be already included in the instance (some distances between hydrogen atoms may be induced from the chemical structure of protein): in such a case, the two available intervals can be intersected.

In order to generate our instances, we consider a matrix  $A$  describing the existing covalent bonds between pairs of vertices in the given order. Such a matrix is very important for identifying the necessary distances between pairs of vertices in the order  $r$ . The  $(i,j)$  element of the squared matrix  $A$  is 1 when there is a covalent bond between the vertices  $i$  and  $j$ , and 0 otherwise. This matrix not only provides a map of all covalent bonds, but it is also necessary for finding the distances between vertices that are separated by a few other chemical bonds. Figure 5 shows two examples of this matrix.

Our procedure for generating  $i$  DMDGP instances (Algorithm 2) has two main nested loops enumerating all possible pairs  $(i, j)$  of vertices in the given order, with  $j &gt; i$ . For each

Springer

---

J Glob Optim (2014) 60:333-349

![img-3.jpeg](img-3.jpeg)
Fig. 3 Handcrafted orders for the amino acids ALA, PRO, GLY, LEU, VAL, ILE, MET, PHE, TRP

Springer

---

J Glob Optim (2014) 60:333-349

![img-4.jpeg](img-4.jpeg)
Fig. 4 Handcrafted orders for the amino acids GLN, TYR, ASN, SER, CYS, THR, ARG, HIS

$i &gt; 3$ , we know that the order  $r$  has been constructed so that the distances related to the pairs  $(i - 3, i)$ ,  $(i - 2, i)$  and  $(i - 1, i)$  are all available. Such distances can be either bond lengths (the element  $(i, j)$  of  $A$  is 1), zero ( $i$  and  $j$  represent the same atom), or obtained

Springer

---

J Glob Optim (2014) 60:333-349

![img-5.jpeg](img-5.jpeg)
(a)
Fig. 5 Two examples of the covalent bond matrix  $A$ . Each element  $(i, j)$  represents a covalent bond between the atom  $i$  and the atom  $j$ . a The matrix  $A$  for a Lysine (LYS), where  $nz$  is the number of covalent bonds. b The matrix  $A$  for the sequence LEU-GLU-LYS-VAL

![img-6.jpeg](img-6.jpeg)
(b)

by considering bond lengths and bond angles among different vertices. Only these distances are considered: any other distances related to pairs  $(i,j)$ , with  $j &gt; i + 3$ , are added to the instance only if they correspond to 0, i.e. when  $i$  and  $j$  represent the same atom.

Algorithm 2 Algorithm for generating iDMDGP instances
1: Read the sequence  $r_{SC}$ .
2: Read the matrix  $A$  that indicates all covalent bonds among the elements of  $r_{SC}$ .
3: Define distances in the initial clique  $(r_{SC}^{1}, r_{SC}^{2}, r_{SC}^{3})$ .
4: for all  $(i \in r_{SC}$  and  $i &gt; 3)$  do
5: for all  $(k \in \{3, 2, 1\})$  do
6: nBonds  $\leftarrow$  number of covalent bonds between the  $r_{SC}^{i-k}$  and  $r_{SC}^{i}$ .
7: if (nBonds = 0) then
8:  $r_{SC}^{i-k}$  and  $r_{SC}^{i}$  are the same atom.
9:  $d(r_{SC}^{i-k}, r_{SC}^{i}) \gets 0$
10: end if
11: if (nBonds = 1) then
12:  $d(r_{SC}^{i-k}, r_{SC}^{i})$  is the length of the covalent bond between  $r_{SC}^{i-k}$  and  $r_{SC}^{i}$ .
13: end if
14: if (nBonds = 2) then
15: Find the atom  $X$  bonded to both  $r_{SC}^{i-k}$  and  $r_{SC}^{i}$  (by using  $A$ ).
16: Compute  $d(r_{SC}^{i-k}, r_{SC}^{i})$  by inverting the cosine law and using the bond angle  $ang(r_{SC}^{i-k}, \hat{X}, r_{SC}^{i})$  and the bond lengths  $d(r_{SC}^{i-k}, X)$  and  $d(X, r_{SC}^{i})$ .
17: end if
18: if (nBonds = 3) then
19: The distance  $d(r_{SC}^{i-k}, r_{SC}^{i})$  is within an interval  $[d^L, d^U]$ .
20: Find the atoms  $X$  and  $Y$  forming, with  $r_{SC}^{i-k}$  and  $r_{SC}^{i}$ , a torsion angle  $\omega$  (by using  $A$ ).
21: Calculate  $d^L$  and  $d^U$  by computing the minimal and the maximal extension for  $\omega$ .
22: end if
23: end for
24: end for

Springer

---

In our implementation, we exploit information in the matrix A for finding out how many chemical bonds separate two vertices. There are two interesting situations: when i and j are separated by 2 bonds, and when they are separated by 3 bonds. In the first case, if k is the vertex connected to both i and j, the distance between them can be computed by considering the bond lengths (i, k) and (k, j), as well as the bond angle formed by the three atoms. The distance related to this pair (i, j) is evidently exact. In the case of 3 bonds, we can compute the lower and the upper bound for the distance between i and j, given by the minimal and maximal extension of the torsion angle that they can define. In this case, the distance related to (i, j) is an interval.

This procedure, described in details in Algorithm 2, is able to generate an instance of the iDMDGP for a given discretization order r. This order can be related to the protein backbone, for example, as well as to side chains, or it could combine both. This procedure provides an instance without NMR distances. In order to simulate them, we randomly choose a branch of the tree related to the instance without pruning distances, we compute all the distances between pairs of hydrogens belonging to this branch, and we consider in the instance all those distances that are shorter than 5Å. In order to simulate NMR data, which are generally imprecise, we randomly generate an interval, which contains the computed distance between the two hydrogens, and add this interval in the instance.

## Computational experiments

Computational experiments where we consider the 20 side chains that can be found in proteins are presented in this section. We consider two sets of instances. In the first one, the instances are artificially generated by applying the procedure detailed in Sect. 4. The sequences of amino acids are randomly generated in a way that each side chain is considered at least once. We remark that preliminary computational experiments presented in [2] already showed that the inclusion of the side chains in iDMDGP instances can help iBP in focusing the search on the feasible parts of the search domain. Our codes were partially written in C and C++.

Table 2 shows some experiments on the first set of considered instances. All instances are formed by 4 randomly chosen amino acids. The parameters of the iBP algorithm (Algorithm 1) for these experiments are D = 8 and ϵ_{d} = 0.1 (D is the discretization factor for the interval distances, and ϵ_{d} is the pruning distance tolerance, which aims to compensate some numerical errors that might propagate). In the table, we provide the total number of vertices forming our instances, the number of solutions found by iBP, as well as the CPU time (in seconds). The experiments show that all instances can be discretized, so that the iBP algorithm can be applied for their solution. In these experiments, iBP was able to produce a large number of solutions in a short time interval (only the pruning device based on the distances in E_{p} is here considered). This phenomenon can be alleviated by using the additional pruning devices discussed in Sect. 2.

Before considering our second set of instances, we firstly present an experiment performed on a single amino acid (the lysine) for which we enumerate the whole search tree, i.e. all the possible conformations for a given tolerance ϵ_{d} and discretization factor D. We set D = 4 and ϵ_{d} = 0.075. In this experiment, all pruning devices mentioned in Sect. 2 are used, and the bound fraction related to Eq. (1), for the vdW radii, is set to 0.90. We found a total number of 2981 solutions (see Fig. 6).

Our second set consists of 3 instances containing real NMR data. In this case, the instances are generated by applying Algorithm 2, and by adding real NMR distances to the obtained instance (to be included in the set E_{p}). We consider the three following molecules:

---

J Glob Optim (2014) 60:333-349

Table 2 Number of solutions for random sequences of four amino acids

|  Sequence | Number of vertices | Number of solutions | CPU time  |
| --- | --- | --- | --- |
|  ALA-THR-CYS-ILE | 94 | 68394 | 375.83  |
|  ASP-TYR-CYS-ALA | 94 | 315624 | 302.12  |
|  VAL-ALA-ILE-GLY | 86 | 99350 | 300.84  |
|  TYR-ASP-PHE-ASN | 117 | 71556 | 302.52  |
|  CYS-ASP-VAL-LYS | 106 | 8288 | 387.15  |
|  GLN-SER-LEU-ASN | 105 | 91238 | 301.23  |
|  LEU-GLU-LYS-VAL | 94 | 71682 | 335.20  |
|  ILE-ALA-SER-LYS | 106 | 10392 | 38.49  |
|  TRP-SER-CYS-GLY | 88 | 30078 | 301.12  |
|  PHE-ALA-GLY-LEU | 93 | 111168 | 301.21  |
|  PHE-VAL-ILE-VAL | 121 | 188240 | 326.82  |
|  MET-ALA-LEU-GLY | 91 | 17182 | 300.85  |
|  GLU-GLU-ASN-GLY | 85 | 42834 | 303.02  |
|  GLN-SER-LEU-LYS | 121 | 27412 | 303.13  |
|  CYS-SER-LYS-MET | 109 | 112320 | 301.85  |

![img-7.jpeg](img-7.jpeg)
Fig. 6 Six super-imposed conformations for lysine, chosen as representative of the 2981 found solutions. Different colors distinguish the shown conformations. One of them shows the atom types

2JMY A small protein composed by 15 amino acids and 282 atoms. It forms a single  $\alpha$ -helix (Fig. 7). Backbone atoms, as well as all side chain atoms, are included.
2KXA A protein formed by two helices linked by a flexible coil. It is formed by 24 amino acids and 332 atoms (Fig. 7). All backbone and side chain atoms are included.
2KSL A protein formed by 51 amino acids and 763 atoms. It is folded in four helices. This protein has a peculiarity: it forms the so-called disulfide bridge, i.e. a strong bond

Springer

---

J Glob Optim (2014) 60:333-349

![img-8.jpeg](img-8.jpeg)
Fig. 7 On the left, the conformation of protein 2JMY; on the right, the conformation of 2KXA

![img-9.jpeg](img-9.jpeg)
Fig. 8 Solution for the protein 2KSL. Note the presence of the disulfide bridge (in thick black) between two  $S_{y}$  atoms belonging to cysteines

between sulfides belonging to pairs of cysteines (Fig. 8). We consider the first 36 residues (264 atoms), so that only one disulfide bridge is present. We consider here all backbone atoms, the Carbons  $C_{\beta}$  for each side chain, and the cysteines.

Together with the considered values for  $D$  and  $\epsilon_d$ , Table 3 reports the number of nodes generated during the execution of the iBP algorithm (including duplicates), as well as the number of accepted nodes. We also report the CPU time necessary for finding the first solution for these three instances. All pruning devices in Sect. 2 are used in the experiments. They show that, despite a considerable number of repetitions in the side chains orders, it is possible to solve NMR instances by using this discrete approach. The increased complexity of the pruning phase poses several implementation and performance issues: the presented results might be improved by using future versions of iBP.

# 6 Conclusions

We provided 20 special orderings for the 20 different side chains that can be present in proteins, which allow us to discretize all iMDGP instances related to proteins. Only information about the chemical structure of the amino acids was exploited during the conception of such

Springer

---

J Glob Optim (2014) 60:333-349

Table 3 Experiments on instances containing real NMR data

|   | D | εd | # Generated nodes | # Accepted nodes | Time (s)  |
| --- | --- | --- | --- | --- | --- |
|  2JMY | 4 | 0.1 | 3791493170 | 1616994331 | 6665.55  |
|  2KXA | 5 | 0.075 | 3479815027 | 1122340143 | 7530.35  |
|  2KSL | 5 | 0.1 | 56544801 | 26462289 | 139.79  |

orders, while NMR distances are used only for pruning purposes. Preliminary computational experiments illustrated the applicability and the usefulness of the proposed orders.

This work represents an important step ahead for the long-term project of developing an exact algorithm to be able to provide all the possible molecular structures for a protein from NMR data. In [15], we had already shown that this was possible when considering protein backbones. This paper, through computational experiments, showed the same result for the side chains that can be present in protein conformations.

Acknowledgments We are thankful to Thérèse Malliavin and Douglas Gonçalves for the fruitful comments on this paper, as well as to the anonymous referees. We also wish to thank CAPES, that funded a 4-month visit to Rennes for Virginia Costa: part of this work was performed during such a visit. We are also thankful to the French Embassy in São Paulo and UNICAMP (which funded a 2-month chaire in UNICAMP for Antonio Mucherino), to the Brazilian research agencies FAPESP and CNPq, and to the French National Research Agency (ANR).

# References

1. Berman, H., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T., Weissig, H., Shindyalov, I., Bourne, P.: The protein data bank. Nucleic Acid Res. 28, 235-242 (2000)
2. Costa, V., Mucherino, A., Lavor, C., Carvalho, L.M., Maculan, N.: On Suitable Orders for Discretizing Molecular Distance Geometry Problems related to Protein Side Chains. In: IEEE Conference Proceedings, Federated Conference on Computer Science and Information Systems (FedCSIS12), pp. 397-402, Wroclaw, Poland, 2012. Workshop on Computational, Optimization (WCO12)
3. Costa, V., Mucherino, A., Carvalho, L.M., Maculan, N.: On the Discretization of iDMDGP instances regarding Protein Side Chains with rings. In: Proceedings of Distance Geometry and Applications (DGA13), Manaus, Amazonas, Brazil, pp. 99-102 (2013)
4. Crippen, G., Havel, T.: Distance Geometry and Molecular Conformation. Wiley, New York (1988)
5. Donald, B.: Algorithms in Structural Molecular Biology. MIT Press, Boston (2011)
6. Grishaev, A., Bax, A.: An empirical backbone-backbone hydrogen-bonding potential in proteins and its applications to nmr structure refinement and validation. J. Am. Chem. Soc. 126, 7281-7292 (2004)
7. Havel, T.: Distance geometry. In: Grant, D.M., Harris, R.K. (eds.) Encyclopedia of Nuclear Magnetic Resonance, pp. 1701-1710. Wiley, New York (1995)
8. Honig, B., Nicholls, A.: Classical electrostatics in biology and chemistry. Science 268, 1144-1149 (1995)
9. Lavor, C., Lee, J., Lee-St. John, A., Liberti, L., Mucherino, A., Sviridenko, M.: Discretization orders for distance geometry problems. Optim. Lett. 6(4), 783-796 (2012)
10. Lavor, C., Liberti, L., Maculan, N., Mucherino, A.: The discretizable molecular distance geometry problem. Comput. Optim. Appl. 52, 115-146 (2012)
11. Lavor, C., Liberti, L., Mucherino, A.: The interval Branch-and-Prune algorithm for the discretizable molecular distance geometry problem with inexact distances. J. Global Optim. 56, 855-871 (2013)
12. Liberti, L., Lavor, C., Maculan, N., Mucherino, A.: Euclidean distance geometry and applications. SIAM Rev. 56(1) (2014)
13. Liberti, L., Lavor, C., Mucherino, A., Maculan, N.: Molecular distance geometry methods: from continuous to discrete. Int. Trans. Oper. Res. 18, 33-51 (2010)
14. Mucherino, A., Lavor, C., Liberti, L., Maculan, N. (eds.): Distance Geometry: Theory, Methods and Applications, 410 pp. Springer, Berlin (2013)
15. Mucherino, A., Lavor, C., Malliavin, T., Liberti, L., Nilges, M., Maculan, N.: Influence of pruning devices on the solution of molecular distance geometry problems. In: Pardalos, P.M., Rebennack, S. (eds.) Lecture

Springer

---

J Glob Optim (2014) 60:333–349

Notes in Computer Science 6630, Proceedings of the 10th International Symposium on Experimental Algorithms (SEA11), Crete, Greece, pp. 206–217 (2011)

16. Nilges, M., Clore, G., Gronenborn, A.: Determination of three-dimensional structures of proteins from interproton distance data by hybrid distance geometry-dynamical simulated annealing calculations. Fed. Eur. Biochem. Soc. 229, 317–324 (1988)

17. Rocchia, W., Alexov, E., Honig, B.: Extending the applicability of the nonlinear Poisson–Boltzmann equation: multiple dielectric constants and multivalent ions. J. Phys. Chem. B 105, 6507–6514 (2001)

18. Sallaume, S., Martins, S., Ochi, L., Gramacho, W., Lavor, C., Liberti, L.: A discrete search algorithm for finding the structure of protein backbones and side chains. Int. J. Bioinf. Res. Appl. 9, 261–270 (2013)

19. Saxe, J.: Embeddability of weighted graphs in $k$-space is strongly NP-hard. In: Proceedings of 17th Allerton Conference in Communications, Control and Computing, pp. 480–489 (1979)

20. Schlick, T.: Molecular Modelling and Simulation: An Interdisciplinary Guide. Springer, New York (2002)

21. Souza, M., Lavor, C., Muritiba, A., Maculan, N.: Solving the molecular distance geometry problem with inaccurate distance data. BMC Bioinf. 14(Suppl 9), S7 (2013)

22. Souza, M., Xavier, A., Lavor, C., Maculan, N.: Hyperbolic smoothing and penalty techniques applied to molecular structure determination. Oper. Res. Lett. 39, 461–465 (2011)

Springer