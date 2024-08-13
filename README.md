# magma-clifford-algebras
A Magma package for Clifford Algebras extending magmas AlgClff and AlgClffElt types. 



**Algorithms**

Dupuy, Taylor, Anton Hilado, Colin Ingalls, and Adam Logan. "The Basic Theory of Clifford-Bianchi Groups for Hyperbolic n-Space." *arXiv preprint arXiv:2407.19122* (2024).

Here is the bibtex:

```
@article{DHIL2024,
  title={The Basic Theory of Clifford-Bianchi Groups for Hyperbolic n-Space},
  author={Dupuy, Taylor and Hilado, Anton and Ingalls, Colin and Logan, Adam},
  journal={arXiv preprint arXiv:2407.19122},
  year={2024}
}
```



**Installation of Magma**

Magma is publicly available throught the [Simons Agreement](http://magma.maths.usyd.edu.au/magma/simons_details).

Sage notebook users may prefer to use a [Jupyter](https://jupyter.org/install) notebook with a [magma kernel](https://github.com/edgarcosta/magma_kernel). Magma.



**Running the Code**

There are comments in the code for how to run things but these need to be updated. It would be a helpful contribution to update the examples. 



Here is a minimal introduction. 



You load the code by "attaching". This uses the magma intrinsics we have defined. These are different from magma functions in that they are global and the Magma types (which are essentially python classes) are defined. 

```
>Z:=Integers();
>Q:=Rationals();
>AttachSpec("cliff.spec");
```



We have some functionality for making clifford algebras with our CliffAlg and CliffAlgElt type. We have orders with elements CliffOrd, CliffOrdElt and we have GrpPSL2Cliff types (these types don't have elements yet and matrices are implemented as 2x2 things with Magma's built-in AlgClffElt elements as entries so that Magma's coersion works better. Currently coercion via the ISA program isn't available to people developing packages but they promise to implement it in the future. )

```
>d:=4;
>cl := clifford_algebra(d,Q); //there are other ways to specify >clifford algebras in the code
>O:=get_maximal_orders(cl)[1];
>PSL2O:=clifford_bianchi_group(O : verbose:=true, word_length:=8) //this will take a while
```

We can extract magmas FPGrp using fp_group. We also have some functionality for converting between matrix elements and finite group elements. This is mostly used in our for modular symbol algorithms which we are currently developing. 
```
>fp_group(PSL2O);
LARGE OUTPUT 
```

One can extract the lattice of clifford vectors as a Magma Lat object. There is also some functionality for getting the basis as CliffOrdElt types. 
```
>vec(O);
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
>I:=get_I(cl); //cl<I> :=clifford_algebra(d,Q) is not implemented. This gives you the generators as an algebra.
>a:=2+2*I[2]/2;
>b:=closest_vector(O,a);
>Eltseq(a)[1]; // Eltseq returns a magma sequence
>star(a);
>prime(a);
>bar(a);
>nm(a); //returns a clifford vector 
>Type(b);
>Type(a);
2
(2 0 1 0 0 0 0 0)
( 2  0 -1  0  0  0  0  0)
( 2  0 -1  0  0  0  0  0)
(5 0 0 0 0 0 0 0)
CliffOrdElt
CliffAlgElt
```

We also make matrices. We can take pseudodeterminants and inverses. As stated previously we are using matrices with AlgClffElt entries so if you want to get the entry of a matrix and treat it like a CliffAlgElt you will need to coerce. 
```
>mat:=Matrix(2,2, [cl| 1,2*I[1],0,1]);
>sig:=Matrix(2,2,[cl|0,-1,1,0]);
>pdet(mat);
>mat^(-1);
(1 0 0 0 0 0 0 0)
[[1 0 0 0 0 0 0 0] [ 0 -2  0  0  0  0  0  0]]
[[0 0 0 0 0 0 0 0] [1 0 0 0 0 0 0 0]]
```

We also have some functionality for working with cusps. One can use a * notation or an operator notation for actions on cusps.
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

