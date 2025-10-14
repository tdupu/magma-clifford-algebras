# About magma-clifford-algebras

This is a repository for working with Clifford algebras, orders $O$ in Clifford algebras, Clifford-Bianchi groups $SL_2(O)$, and modular symbols for hyperbolic $n$-space.

Many of the functions are documented in this pre-print: 
[The Basic Theory Of Clifford-Bianchi Groups](https://arxiv.org/abs/2407.19122)

```
@article{DHIL2024,
  title={The Basic Theory of Clifford-Bianchi Groups for Hyperbolic n-Space},
  author={Dupuy, Taylor and Hilado, Anton and Ingalls, Colin and Logan, Adam},
  journal={arXiv preprint arXiv:2407.19122},
  year={2024}
}
```

The package was created by Taylor Dupuy and Adam Logan starting in 2023. The storage functions are based off earlier work of David Roe and Edgar Costa for the database of finite groups in the LMFDB.


**Installation of Magma**

Magma is publicly available throught the [Simons Agreement](http://magma.maths.usyd.edu.au/magma/simons_details).

Sage notebook users may prefer to use a [Jupyter](https://jupyter.org/install) notebook with a [magma kernel](https://github.com/edgarcosta/magma_kernel).



**Running the Code**

There are comments in the code for how to run things but these need to be updated. It would be a helpful contribution to update the examples. 


# Contents

There is a section on the [Intrinsics](#intrinsics) of the package.
1. [Attaching The Package](#attaching-the-package)
2. [Basic Computations With Clifford Algebras](#basic-computations-with-clifford-algebras)
3. [Order in Clifford Algebra](#orders-in-clifford-algebras)
4. [Clifford Units](#units)
5. [Using Clifford-Euclideanity](#division-algorithms)
6. [Closest Vectors](#closest-vectors)
7. [Cusps and Hyperbolic Space](#cusps-and-hyperbolic-space)
8. [Polyhedra and Polytopes](#working-with-polyhedron-and-polytopes)
9. [Crystallographic Groups: Fundamental Domains For Stabilizer of Infinity](#fundamental-domains-for-stabilizers-of-infinity)
10. [Fundamental Domains](#fundamental-domains)
11. [Polyhedra and Polytopes part 2](#polyhedra-and-polytopes-part-2)

There is a section on [Magma Fixes](#magma-fixes). General magma enhancements dealing with existing magma types. 

1. [Coset Tables and Presentations](#coset-tables-and-presentations)
2. [Canonical Cosets](#canonical-cosets)
3. [Transfer Maps](#transfer-maps)
4. [Abelianization of Group Homomorphisms](#abelianization-of-group-homomorphisms)
   

There is a section on [Writing to Disc and Loading from Disc](#writing-to-disc-and-loading-from-disc). This is mainly for a future application to an LMFDB database.
1. [Unix Remote Computation Cheat Sheet](#remote-computations)
2. [Saving Hecke Matrices](#saving-hecke-matrices-to-disc)
3. [Loading Hecke Matrices](#loading-hecke-matrices-from-disc)

# Intrinsics


### Attaching The Package

```
Z:=Integers();
Q:=Rationals();
AttachSpec("hecke.spec");
```

If computations to or from the disc are going to be made then we need to set the folder where the `DATA` will appear. 

```
SetLMFDBRootFolder(".");
```

There are several verbose flags which go up to 6 which can be helpful. You can set them by 
```
SetVerbose("Clifford",2);
```

### Basic Computations With Clifford Algebras

Here is a minimal introduction. 



You load the code by "attaching". This uses the magma intrinsics we have defined. These are different from magma functions in that they are global and the Magma types (which are essentially python classes) are defined. 

```
>Z:=Integers();
>Q:=Rationals();
>AttachSpec("cliff.spec");
```

If computations to or from the disc are going to be made then we need to set the folder where the `DATA` will appear. 

```
SetLMFDBRootFolder(".");
```

There are several verbose flags which go up to 6 which can be helpful. You can set them by 
```
SetVerbose("Clifford",2);
```


We have some functionality for making clifford algebras with our CliffAlg and CliffAlgElt type. We have orders with elements CliffOrd, CliffOrdElt and we have GrpPSL2Cliff types (these types don't have elements yet and matrices are implemented as 2x2 things with Magma's built-in AlgClffElt elements as entries so that Magma's coersion works better. Currently coercion via the ISA program isn't available to people developing packages but they promise to implement it in the future. )

```
d:=4;
cl := clifford_algebra(d,Q); //there are other ways to specify >clifford algebras in the code
O:=get_maximal_orders(cl)[1];
PSL2O:=clifford_bianchi_group(O) //this will take a while
```

We can extract magmas FPGrp using fp_group. We also have some functionality for converting between matrix elements and finite group elements. This is mostly used in our for modular symbol algorithms which we are currently developing. 
```
>fp_group(PSL2O);
```
```
LARGE OUTPUT 
```

One can extract the lattice of clifford vectors as a Magma Lat object. There is also some functionality for getting the basis as CliffOrdElt types. 
```
vec(O);
```
```
Lattice of rank 4 and degree 8
Determinant: 1/4
Factored Determinant: 1/2^2
Minimum: 1
Basis:
(1 1 1 0 1 0 0 0)
(0 2 0 0 0 0 0 0)
(0 0 2 0 0 0 0 0)
(0 0 0 0 2 0 0 0)
Basis Denominator: 2
```

We have some ways of working with explicit elements. 
```
I:=get_I(cl); //cl<I> :=clifford_algebra(d,Q) is not implemented. This gives you the generators as an algebra.
a:=2+2*I[2]/2;
b:=closest_vector(O,a);
Eltseq(a)[1]; // Eltseq returns a magma sequence
star(a);
prime(a);
bar(a);
nm(a); //returns a clifford vector 
Type(b);
Type(a);
```
```
2
(2 0 1 0 0 0 0 0)
( 2  0 -1  0  0  0  0  0)
( 2  0 -1  0  0  0  0  0)
(5 0 0 0 0 0 0 0)
CliffOrdElt
CliffAlgElt
```

We also make matrices. We can take pseudodeterminants and inverses. As stated previously we are using matrices with `AlgClffElt` entries so if you want to get the entry of a matrix and treat it like a `CliffAlgElt` you will need to coerce. 
```
mat:=Matrix(2,2, [cl| 1,2*I[1],0,1]);
sig:=Matrix(2,2,[cl|0,-1,1,0]);
pdet(mat);
mat^(-1);
```
```
(1 0 0 0 0 0 0 0)
[[1 0 0 0 0 0 0 0] [ 0 -2  0  0  0  0  0  0]]
[[0 0 0 0 0 0 0 0] [1 0 0 0 0 0 0 0]]
```

We also have some functionality for working with cusps. One can use a `*` notation or an operator notation for actions on cusps.
```
>cusp_set:=cusps(cl);
>my_cusp:=cusp_set!a; // the exclaimation point is magma's version of coercion
>infty:=infinity(cusp_set);
>zero:=zero(cusp_set);
>mat*my_cusp;
>mat(infty);
>sig*infty eq zero; //there is some functionality for equality and comparison le,lt etc
(2 2 1 0 0 0 0 0)
infinity
true
```

There are a bunch of other standard functions. There are the conjugations `prime`,`star`,`bar`; testing for membership of the Clifford vectors `is_vec`; the norm `nm` which is $x \bar{x}$.

There is also some functionality for extending Clifford algebras. 
This is useful for actions on hyperbolic space. 
The intrinsic `extend` will add an extra generator $j$ with $j^2=-1$. 
```
cl<[I]> := clifford_algebra(Q,[-1,-3]);
cl_big:=extend(cl); cl_big;
is_extension(cl_big,cl);
a:=cl.2;
b:=cl_big.1;
cl_big!Eltseq(cl_big,a) eq cl_big!a;
a+b;
a*b;
a-b;
```
```
Clifford algebra over Rational Field with structure constants [ -1, -3, -1 ]
true
true
(1 1 0 0 0 0 0 0)
(0 1 0 0 0 0 0 0)
(-1  1  0  0  0  0  0  0)
```


### Orders in Clifford Algebras

There is some functionality for computing maximal orders in a given rational clifford algebra. 
```
cl<[I]>:=clifford_algebra(Q,[-1,-1,-1]);
ord := maximal_orders(cl)[1];
```

There is functionality for working with the underlying 2^3 dimensional lattice of the order and the 4-dimensional lattice of Clifford vectors. The function `clifford_vector` will take a lattice element and convert it back into an order element. 
```
time lattice(ord);
time vec(ord);
```
```
Lattice of rank 8 and degree 8
Determinant: 1/4096
Factored Determinant: 1/2^12
Basis:
( 1  1  1 -1  1  1 -1  1)
( 1 -1  1  1  1 -1 -1 -1)
( 1 -1  1  1 -1  1  1  1)
( 1 -1  1 -1  1  1  1  1)
( 1 -1  1 -1 -1 -1 -1 -1)
( 1 -1 -1 -1 -1  1 -1 -1)
( 0  2  0  0  0  0  2  0)
( 0  0  2  0  0  2  0  0)
Basis Denominator: 4
Time: 0.000
Lattice of rank 4 and degree 8
Determinant: 1/4
Factored Determinant: 1/2^2
Minimum: 1
Basis:
( 1  1 -1  0 -1  0  0  0)
( 1 -1  1  0  1  0  0  0)
( 0  2  0  0  0  0  0  0)
( 0  0  2  0  0  0  0  0)
Basis Denominator: 2
Time: 0.000
```

The coverting radius of the Clifford vectors determines if an order is norm-Euclidean or not. This determines if division algorithms are available for the order.  In this example the covering radius is less than 1 so the order is norm-Euclidean.
```
time CoveringRadius(vec(ord));
```
```
1/2
Time: 0.000
```

For working with Hecke operators or $\Gamma_0$ subgroups it is often important to enumerate elements of the Clifford group which have a prime norm.
```
time primes_up_to_sample(ord,20);
```
```
{@
( 1/2 -5/2 -3/2    0 -3/2    0    0    0),
( 1/2 -7/2 -5/2    0 -1/2    0    0    0),
( 1/2 -5/2 -1/2    0 -1/2    0    0    0),
( 1/2 -7/2 -3/2    0 -3/2    0    0    0),
( 1/2 -3/2 -3/2    0 -1/2    0    0    0),
( 1/2 -3/2 -1/2    0 -1/2    0    0    0),
( 1 -1  0  0  0  0  0  0),
( 1/2 -7/2 -1/2    0 -1/2    0    0    0)
@}
Time: 0.010
```

We also have options for working with integral lattices. Here there is a second output and `(1/den)*lat` will be the original lattice output.
```
lat,den:=vec(O:small:=true,integral:=true);
```

### Units

Important in working with Clifford algebras are the units and Clifford units. The Clifford units are actually more important. Over a field of characteristic zero the Clifford units are the subgroup of the group of units generated by non-zero Clifford vectors. Elements of this group which are members of an order whose inverse is also a member of the order are Clifford group elements of the order. 

The intrinsic `units` is a simplified version of the `cliff_units` function.

```
cl:=clifford_algebra(Q,[-1,-1,-1]);
time O:=get_maximal_orders(cl)[1];
time c_units:=units(O);
```
```
Time: 0.700
Time: 1.060
```
The function is slower for d = 4 but the only one that has a chance for d = 6 or larger, etc.
--The `bound` argument specifies the maximum number of generators to use in the search for short vectors.
--The `find_all` argument specifies whether to find all units of the order, or just some.

The function works by finding all the short vectors in a lattice that is associated with the order, subject to the condition that the sum of any two short vectors is also a short vector. It then checks which of these short vectors are Clifford group elements, and returns those elements as the units of the order. The function returns a list of units of the order, a list of all units of the order (if find_all is true), and the matrix group generated by the known units.

Note that finding the order of the matrix group directly may be slow and it is recommended to reduce it modulo a small prime and use the CompositionTree functionality instead.

For each unit $u$ we can consider the representation $\pi_u(x) = uxu^*$ where $x$ is a Clifford vector. 
We can write this in is standard basis. Magma uses a right-multiplication-by-matrices convention so we have a version of this map that multiplies on the left and one that multiplies on the right. 
```
u1:=c_units[10];
u2:=c_units[7];
u1*u2 eq u2*u1;
pi_repn(u1)*pi_repn(u2) eq pi_repn(u1*u2);
pi_repn_right(u2)*pi_repn_right(u1) eq pi_repn_right(u1*u2);
```
```
false
true
true
```




### Division Algorithms

Orders in Clifford algebras which are norm-Euclidean have Euclidean algorithms which are implemented. These algorithms are central to `closest_vector` our modular symbol algorithms 

```
d := 3;
C:=clifford_algebra(d);
mygens:=generators(C);
ord := Order(Z,mygens); //Lipschitz Order
mos := maximal_orders(mygens); //There is a unique maximal order containing the Lipschitz order
O:=mos[1];
c:=(1+I[1]+2*I[2]+20*I[3]);
b:=c;
a:=10;
clifford_gcd(a,b);
```

ord := get_maximal_orders(cl)[1];
time lattice(ord);
time vec(ord);
time CoveringRadius(vec(ord));
time primes_up_to_sample(ord,20);
time irreds:=primes_of_size(ord,3);

### Cusps and Hyperbolic Space

```
cl<[I]> := clifford_algebra(Q,[-1,-3]);
cusp_set:=cusps(cl);
H:=hyperbolic_space_with_cusps(cl);
pt:=elt(H,[1,1,1,1]);
pt eq H!Eltseq(pt);
pt eq H!vector(pt);
pt eq H!clifford_vector(pt);
```
```
true
true
true
true
```

We can also action on points of hyperbolic space by elements of SL2(cl).
```
inv:=Matrix(2,2,[cl|0,-1,1,0]);
Cusps:=cusps(cl);
pt:=H!Cusps!cl!1;
pti:=H!Cusps!I[1];
pti eq inv*pti;
pt eq inv*inv*pt;
```
```
true
true
```

The action as currently written works through `maps_of_affine_space(H,g)` which converts a given matrix in SL2(cl) into a rational morphism of schemes from `affine_space(H)` to itself. This is eventually use to act on hyperbolic polyhedra. The method here is to act on the underlying schemes of halfspaces and then taking the equations of those to get the translated halfspaces. This should eventually be rewritten to work through the big clifford algebra.

### Closest Vectors  

For an already instantiated Clifford-Euclidean order there is support for find closest vectors. This is used in continued fraction algorithms which help build the modular symbols. We also have ways of converting between these vectors and 2x2 matrices with `AlgClffElt` entries.

```
a:=I[1]+I[2]/4;
alpha:=closest_vector(O,a);
alpha_quo:=vec_to_psl2quo(PSL2O,alpha);
alpha_matrix:=fp_element_to_matrix(PSL2O,alpha_quo);
translation_matrix(alpha) eq alpha_matrix;
```
There was an issue at one point with negative exponents so we added this test case.

```
alpha := I[1]-2*I[2];
vec_to_psl2quo(PSL2O,alpha);
y[3]^-1 * y[2]^-1 * y[3]^2 * y[2]^-2 * y[3]^-1 // output varies
```
        
Here is another case where we had to add something to test zero.

```
x:=I[1]/5-I[2]/4;
z:=closest_vector(O,x);
z eq 0;
vec_to_psl2quo(PSL2O,z) eq Identity(fp_group(PSL2O));
```

### Basic Operations On Polyhedra

When computing fundamental domains we sometimes need to act on polyhedron in the space of Clifford vectors by element of the clifford group or by translation. We have implemented some functions for working with this. You can create polytopes and polyhedron from clifford vectors and you can act on them by exponentiating. The right action of some Clifford group element $u$ on a clifford vector $x$ is given $x\mapsto u^* x u$ and we have used the notation `x^u$ to do this. 
```
cl<[I]>:=clifford_algebra(Q,[-1,-1]);
v1:=1+I[1]+2*I[2];
v2:=I[2];
v3:=I[1];

v3^(v1*v2) eq (v3^v1)^v2;
```
```
true
```

We can also convert between Magma's `ModTupFldElt` for vectors and clifford vectors and act on those.
```
vector(v3)^(v1*v2) eq vector(v3^(v1*v2));
```
```
true
```

We also have support for polyhedron (Magma's `TorPol`) from Clifford vectors and actions on the polyhedron.
```
ph:=polyhedron([v1,v2,v3]);
pt:=polytope([v1,v2,v3]);
pt^v3;
```
```
2-dimensional polytope with 3 generators:
    (-1, -1,  2),
    ( 0,  0,  1),
    ( 0, -1,  0)
```
Note that vertices of Polyhedron are a `TorPolElt` type which are different from `ModTupFldElt` and action on those types is also supported.

### Fundamental Domains For Stabilizers Of Infinity

WARNING: `side_pairing_infty` was formerly named `side_pairing`.

In order to compute the stability of $\infty$ under $PSL_2(O)$ we first compute a fundamental domain for an extension of $\lbrace \pi_u : u \in O^{\times}\rbrace \subset PSL_2(O)$ by $inv$ where $inv(x)=-\bar{x}$. We work with `units_matrix_group(O:action:=right)` to get actions on polyhedra. Also there is an element `inv_matrix(cl)` if `cl:=clifford_algebra(O)` where `O` is the order in the Clifford algebra we are working with. The map $inv$ is important because it is the map $x\mapsto -x^{-1}$ when $\vert x \vert=1$.

```
cl<[I]> := clifford_algebra(Q,[-1,-3]);
O:=load_clifford_orders(cl)[1];
is_norm_euclidean(O);
cl<[I]>:=clifford_algebra(O); //we have to redo this when we load from file
inv:=func<x|inversion_inv(cl,x)>;
M:=inv_matrix(cl);
M in units_matrix_group(O:action:=`right`,include_inv:=true);
```
```
true
```

For computing fundamental domains for the stabilizer of infinity there are two options. 
We have `include_inv` which computes the fundamental domain $F$ for the extended group. Then the full fundamental domain is given by $F \cup inv(F)$. 
To compute the fundamental domain we take the Voronoi cell and then partition it into cones for each of the facets. 
When the `include_inv` option is included extra subdivisions of the facets are included so that we have no cones which have a non-trivial intersection its image under $inv$.
```
FD_inv:=fundamental_domain_infty(O:include_inv:=true,inv_stable:=false);
FD_full,FD_half:=fundamental_domain_infty(O:include_inv:=false,inv_stable:=true);
FD_inv eq FD_half;
```
```
true
```
We can also compute the side pairing. The intrinsic `side_pairing_infty` return an associative array which for each facet $S_1$ of the fundamental domain returns $(g,S_2)$. The $g$ is the "step" across the facet $S_1$. That is $g(F) \cap F= S_1$. The facet $S_2$ is what you get from pulling back the intersection $g(F) \cap F$ to the original fundamental domain $F$. We always have $S_1=g*S_2$ since $S_2 = (g^-1)*S_1$.
```
time side_pair:=side_pairing_infty(O,FD_full);
```
```
Time: 23.610
```
We can then use the side pairing to compute a presentation for the group. 
```
sides:=Sorted([k: k in Keys(side_pair)]);
F1:=sides[1];
g1,F1_prime:=Explode(side_pair[F1]);
g2,F1_prime_prime:=Explode(side_pair[F1_prime]);
F1_prime_prime eq F1;
```
```
true
```
Here is a call for the stabilizer of infinity:
```
stab,to_stab,gens:=stabilizer_infty(O,FD_full,side_pair);
```

We can see that the generators come from the side pairing:
```
#Generators(stab) eq #{side_pair[k][1] : k in Keys(side_pair)};
```
```
true
```
We can also check that the relations are satisfied
```
rels:=[Eltseq(rel[1]*rel[2]^-1) : rel in Relations(stab)];
{evaluate_relation(gens,rel) : rel in rels};
```
```
{
[[1 0 0 0] [0 0 0 0]]
[[0 0 0 0] [1 0 0 0]],

[[-1  0  0  0] [0 0 0 0]]
[[0 0 0 0] [-1  0  0  0]]
}
```

### Fundamental Domains
There is support for the fundamental domains for SL2(O). There is also support for the side_pairing of this fundamental domain and action of SL2 matrices with Clifford algebra entries on these polyhedra. 

```
Z:=Integers();
Q:=Rationals();
RR:=RealField();
AttachSpec("hecke.spec");
SetLMFDBRootFolder(".");
cl<[I]> := clifford_algebra(Q,[-1,-3]);
O:=load_clifford_orders(cl)[1];
cl<[I]>:=clifford_algebra(O);
FD:=fundamental_domain(O);
```



### Polyhedra and Polytopes part 2
Polyhedra are defined to be a finite intersections of halfspaces. Some polyhedra are polytopes which are the (hyperbolic) convex hull of finitely many points. As an example we look at the output of `fundamental_domain` which outputs a finite number of polyhedra which make up the polyhedral complex that is the fundamental domain. 
```
Z:=Integers();
Q:=Rationals();
RR:=RealField();
AttachSpec("hecke.spec");
SetLMFDBRootFolder(".");
cl<[I]> := clifford_algebra(Q,[-1,-3]);
O:=load_clifford_orders(cl)[1];
cl<[I]>:=clifford_algebra(O);
FD_SL2:=fundamental_domain(O);
pol:=FD_SL2[1];
pol2:=FD_SL2[2];
```

We can take all of the halfspaces that define a given polyhedron.
```
halfspaces(pol)
```
```
{
Hyperbolic halfspace defined by -$.1 + 2*$.2 - 3*$.3>=0,
Hyperbolic halfspace defined by -$.1 - 6*$.2 + $.3>=0,
Hyperbolic halfspace defined by $.1 - $.3>=0,
Hyperbolic halfspace defined by -$.1>=0,
Hyperbolic halfspace defined by $.1 + 3*$.3 + 1>=0,
Hyperbolic halfspace defined by $.1^2 + $.2^2 + 3*$.3^2 + $.4^2 - 1>=0
}
```

The cluster is the zero dimensional scheme determined by intersection of finitely many sides of the polyhedron (provided they exist). 
```
cluster(pol);
```
```
Scheme over Rational Field defined by
$.1 - 550/3841*$.2 - 47712/3841*$.3^2 - 17193/3841*$.3 + 648/3841*$.4^4 + 290/3841*$.4^2 - 938/3841,
$.1 - 12*$.3^3 - 19*$.3^2 - 5*$.3,
$.1*$.4^2 - 9/16*$.1 + 5/4*$.2 + 33/16*$.3 - 9/4*$.4^2 + 9/4,
$.1 - 72/121*$.2*$.4^2 + 10/11*$.2 - 1344/121*$.3^2 - 309/121*$.3 - 148/121*$.4^2 + 148/121,
$.1 + 20/7*$.2 - 48/7*$.3*$.4^2 + 81/7*$.3 - 52/7*$.4^2 + 52/7,
$.1^2 + $.1 - 9*$.3^2 - 3*$.3,
$.1*$.2 + 39/8*$.1 + 5/2*$.2 - 54*$.3^2 - 111/8*$.3 - 9/2*$.4^2 + 9/2,
$.1 - $.2^2 - 12*$.3^2 - 3*$.3 - $.4^2 + 1,
$.1*$.3 + 3*$.3^2 + $.3,
$.1 - 8/13*$.2*$.3 + 4/13*$.2 - 144/13*$.3^2 - 37/13*$.3 - 12/13*$.4^2 + 12/13
```

There is support for computing the vertices of a polyhedron.
```
vertices(pol);
```
```
[
(0 0 0 1),
(-1/4 -1/2 -1/4   r2),
(-1/4 -1/2 -1/4   r1),
(-1/4    0 -1/4   r4),
(   0 -1/2 -1/3  r10),
(-1/4    0 -1/4   r3),
(    0 -1/18  -1/3    r8),
(   0 -1/2 -1/3   r9),
(    0 -1/18  -1/3    r7),
infinity
]
```
Given the vertices of a polyhedron we can then throw away the halfspaces whose defining hyperplane doesn't contain any of the vertices.
```
prune(pol);
```
```
Hyperbolic convex polyhedra defined by -$.1 + 2*$.2 - 3*$.3>=0
-$.1 - 6*$.2 + $.3>=0
$.1 - $.3>=0
-$.1>=0
$.1 + 3*$.3 + 1>=0
$.1^2 + $.2^2 + 3*$.3^2 + $.4^2 - 1>=0
```

We can take the polytope of a polyhedron. 
```
polt:=polytope(pol);
```
```
Convex hyperbolic polytope of dimension 4. Convex hull of [
(0 0 0 1),
(-1/4 -1/2 -1/4   r2),
(-1/4 -1/2 -1/4   r1),
(-1/4    0 -1/4   r4),
(   0 -1/2 -1/3  r10),
(-1/4    0 -1/4   r3),
(    0 -1/18  -1/3    r8),
(   0 -1/2 -1/3   r9),
(    0 -1/18  -1/3    r7),
infinity
].
```

We also have some support for facets of polytopes. This is currently bugged. 


# Modular Symbols 

Modular symbols are a presentation of relative first homology of the locally symmetric space compactified using the Satake compactification. Chains are linear combinations of chains $\lbrace a,b \rbrace$ where $a$ and $b$ are cusps. These represent oriented geodesics starting at $a$ and ending at $b$.

For an element $M$ in $PSL_2(K)$, where $K$ is the rational clifford algebra, we can write $(M)=\lbrace M(\infty),M(0) \rbrace$ and it is a theorem that the coset representatives for $\Gamma$ inside $PSL_2(O)$ generate the group of modular symbols. The theorem is coded up here and we need to find a way of writing a general modular symbol as a linear combination of these "principal symbols" coming from the coset representatives.

There are a couple ways we store modular symbol elements and there are functions in here that convert between them. At the end of the day modular symbol elements will be encoded into `ModSymCliffElt` objects which will allow the user to ignore these distinctions. Internally, we use all of the various representations of these abstract objects and need to convert between them in order to get generators and relations.

At the end of the day `ModSymCliffElt` will consist of an associative array of `SetCspCliffElt` which is our type for cusps on which `AlgMatElt` consisting of two by two matrices is allowed to act.

TUPLE REPRESENTATION:
There are tuples or "abstract" modular symbols which consist of tuples of the form
    `<<c1,<a1,b1>>, <c2,<a2,b2>>,...,<cm,<am,bm>>>`
where `ci` is in the base ring and `ai`,`bi` are clifford vectors of the `ClffAlgElt` type or infinity.

ASSOCIATIVE ARRAY REPRESENTATION:
For each `<ai,bi>` we will have an associative array ms such that `ms[<ai,bi>] = ci`. This prevents us from making double assignments. Since associative arrays need to have a common type for keys this requires us to create the `SetCspCliffElt` Type.

FREE MODULE REPRESENTATION:
We also store modular symbol elements initially in a free module in bijection with a right transversal (choice of coset representatives) for `Gammaquo` inside of `PSL2Oquo`, our `FPGrp` presentations of $PSL2_(O)$ and a given finite index subgroup.

We eventually move on past the free module representation and to do that we need to compute the relations. The meat of the algorithm comes from a continued fraction decomposition.


# Magma fixes

Some stuff for `TransferMap`,`AbelianQuotient` (of a morphism), `AbelianGeneratorLifts`,`canonical_cosets`,`presentation`,`coset_table` have been implemented. Some are to replace bugs in Magma.  

### Coset Tables and Presentations

The intrinsic `coset_table`,`presentation`,`coset_table_matrix` are an implementation of the Reidemeister-Schreier algorithm. Magma's version of this seems buggy. 
```
test := func<g|forall{r: r in Eltseq(beta^-1*m*beta)|r in ord} where m is &*[psl`sl2quo_to_mat[G.i]: i in [0] cat Eltseq(G!g)]>;
// weirdly much faster when I don't also test m (not just a factor of 2)
ct,tv := coset_table(G,test:limit := 10^6);
assert IsSatisfied(Relations(G),ct);
h := hom<G->SymmetricGroup(#tv)|ct>;
gamma0 := sub<G|h>;
assert forall{x: x in Generators(gamma0)|test(x)};
//time Rewrite(G,~gamma0); // sadly slow, and no longer needed
// This presents the finite-index subgroup, such as Gamma_0(p).
pres,pres_to := presentation(G,ct,tv,test);

// again, faster without testing m
test_m := func<m|forall{r: r in Eltseq(beta^-1*m*beta)|r in ord}>;
ctm,tvm := coset_table_matrices(G,test_m,psl`sl2quo_to_mat,cl);
tvm_comb := [<&*[G|G.(i gt Ngens(G) select Ngens(G)-i else i): i in t[2]],t[2]>: t in tvm];
// The Reidemeister-Schreier algorithm wants things to be in an FP-group,
// not a group of matrices.  tvm_comb is a combinatorial version of the
// transversal calculated in terms of matrices.
pres,pres_to := presentation(ctm,tvm_comb);

assert Image(pres_to) eq gamma0;
// make sure that the relations are as expected
assert forall{x: x in Relations(pres)|&*[psl`sl2quo_to_mat[G.i]: i in [0] cat Eltseq(pres_to(LHS(x)*RHS(x)^-1))] in [ScalarMatrix(2,cl!i): i in [1,-1]]};
assert AQInvariants(gamma0) eq AQInvariants(pres);
```

### Canonical Cosets
```
S6:=SymmetricGroup(6);
a:=S6.1;
b:=S6.2;
H:=sub<S6|a^2,b>;
R,get_rep:=RightTransversal(S6,H);
#canonical_cosets(S6,G,get_rep,lex_comparison);
Index(S6,H);
```

### Abelianization of Group Homomorphisms
Given an `GrpFP` we can compute its abelianization using `AbelianQuotient`. Given a `Map` of `GrpFP` we can compute the induced map from one to the other. This uses `AbelianGeneratorLifts` which takes preimages of the generators of `AbelianQuotient` in the original group. 

### Transfer Maps
Given two finitely presented group $G$ and $H$ with $H$ finite index in $G$ there is a natural map $G^{ab} \to H^{ab}$. 

# Writing to Disc and Loading from Disc

There are many computations which are expensive and tedious to recompute. 
There include computing an order $O$ in large dimension, computing a presentation for a Clifford-Bianchi group $SL_2(O)$ of an order, computing a presentation of a $\Gamma_0(\alpha)$ subgroup of $SL_2(O)$, or computing the Hecke operators which act on the space of modular symbols. There are some implementations for reading and writing from text files.

### Remote Computations

This is a cheat sheet for some commands for doing remote computations. (Use `nice` in front of remote commands to not suck up resources from people working interactively.)

How to find the PID of a process:

```
ps -ef | grep "command name"
```

How to kill a remote process with process identification number PID (which is listed in the first column of the output above).
```
kill -9 PID
```



### Saving Hecke Matrices to Disc

In this note we assume that you have successfully stored orders, Clifford-Bianchi groups, and Gamma0 groups of interest. We will them compute torsion-free Hecke matrices associated to modular symbols of some group of interest. 

```
d:=3;
cl := clifford_algebra(d,Q);
O:=load_clifford_orders(cl)[1];
PSL2O:=load_clifford_bianchi_group(O);
cl:=clifford_algebra(PSL2O);
O:=clifford_order(PSL2O);
beta := Matrix(2,2,[cl|2*(1+cl.2+cl.3),0,0,1-cl.2-cl.3]);
T1:=Cputime();
Gamma0,mapp:=load_Gamma0_fp(PSL2O,beta);
Cputime(T1);
```

```
entering load_lmfdb_clifford_order 
entering load_lmfdb_clifford_order 
entering load_lmfdb_sl2 
entering load_lmfdb_clifford_order
./DATA/gamma0/3.1_1.m64.sxOhZ
./DATA/gamma0/3.1_1.m64.sxOhZ
13.070
```

We now create the modular symbols. 
```
psl_fp := fp_group(PSL2O);
SL2O_quo:=sl2_fp(PSL2O);
quot_map:=PSL2O`sl2_to_psl2_fp;
PGamma0:=quot_map(mapp(Gamma0));
assert Index(psl_fp,PGamma0) eq 480; // 3*4*40
```
There is no output from this cell. 

```
MS:=modular_symbols(PGamma0,PSL2O);

```

We select a group of 2x2 matrices in our clifford algebra which we will use to construct our Hecke operators from. 
```
primes:=primes_up_to(O,3);
alpha_list := [Matrix(2,2,[cl|pi,0,0,cl`prime(Eltseq(pi))]) : pi in primes];
```

We then precompute our Hecke matrices. They are then stored as an associative array attribute in the `ModSymCliff` object `MS`. Anytime they are called after this, unless an option for recomputing is given it will pass the Hecke matrix that is already in storage.
```
Talpha_tf_list:=[];
t1:=Cputime();
for alpha in alpha_list[1..2] do
  time Append(~Talpha_tf_list,hecke_matrix_tf(MS,alpha));
end for;
t2:=Cputime(t1);
```
```
//Large output suppressed
```
We then store our Hecke operators. 

```
time0:=Sprintf("%o",Cputime(t2));
write_lmfdb_hecke(MSp,beta:timing:=time0,overwrite:=true);
```
```
./DATA/gamma0_hecke/3.1_1.m64.sxOhZ/3.1_1.m64.sxOhZ.2FohuJ
[[[[[1,-1,0,0],1],[[0,0,0,0],1]],[[[0,0,0,0],1],[[1,1,0,0],1]]],[[[[1,-3,-1,-1],2],[[0,0,0,0],1]],[[[0,0,0,0],1],[[1,3,1,-1],2]]]]|3.1_1.m64.sxOhZ.2FohuJ|[[[0,27,0,8,-8,0,-1,0],[1,18,0,0,-8,0,-1,0],[0,0,9,-8,8,0,-17,18],[0,-18,0,1,8,0,1,0],[0,-18,0,0,9,0,0,0],[0,0,0,8,-8,0,26,-27],[0,0,0,8,-8,0,-1,0],[0,0,8,0,0,-1,-18,18]],[[27,0,0,8,-8,0,-18,0],[18,1,0,0,-8,0,-18,0],[0,0,9,-8,8,-18,18,0],[-18,0,0,1,8,0,18,0],[-18,0,0,0,9,0,18,0],[0,0,0,8,-8,27,-18,0],[0,0,0,8,-8,0,9,0],[0,0,8,0,0,-18,26,1]]]|86.410
```

### Loading Hecke Matrices From Disc

With the modular symbols object `MS` and the clifford algebra object `beta` instantiated as above and having written Hecke matrices to disc one can restore them to the associative array attribute of `MS` by using a load command 
```
load_hecke_matrices_tf(MSp,beta);
```
