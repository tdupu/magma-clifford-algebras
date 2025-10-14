WARNING: This doesn't match what is currently being stored and needs to be redone. 

--Taylor, October 2025

## Orders

The Magma type is LMFDBCliffOrd

| Column                          | Type    | Example                                                      | Note                                                         |
| ------------------------------- | ------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| label                           | text    | 3.1_1.m64.vtir6                                              | order_label.element_label                                    |
| quadratic_form                  | jsonb   | [[1,0],[1]]                                                  | upper triangular matrix                                      |
| basis                           | jsonb   | [ <[ 1, 1, 1, 1 ], 2>, <[ 0, 1, 0, 0 ], 1>, <[ 0, 0, 1, 0 ], 1>, <[ 0, 0, 0, 1 ], 1> ] |                                                              |
| is_maximal*                     | boolean | true                                                         | 0 false, 1 true                                              |
| clifford_unit_group_generators* | jsonb   |                                                              | TODO                                                         |
| clifford_unit_vectors*          | jsonb   |                                                              | TODO                                                         |
| clifford_vector_basis*          | jsonb   |                                                              |                                                              |
| clifford_vector_lattice*        | string  |                                                              | This should return the LMFDB label of the lattice.           |
| dimension                       | integer | 4                                                            | TODO                                                         |
| conjugates*                     | jsonb   |                                                              | store the elements that act by conjugation to give this order |
| is_norm_euclidean               | boolean | true                                                         |                                                              |
| discriminant                    | integer | -64                                                          |                                                              |
| is_euclidean*                   | boolean | true                                                         |                                                              |

*not implemented

## Clifford Bianchi Groups

SL2(O) for O and order in a Clifford algebra. The magma objects are LMFDBGrpSL2Cliff.

The example is SL2(O) for O the Hurwitz order.

| Column            | Type  | Example                                                      | Notes                           |
| ----------------- | ----- | ------------------------------------------------------------ | ------------------------------- |
| label             | text  | 3.1_1.m64.vtir6                                              | LMFDB label of the SL2(O) group |
| order             | text  | 3.1_1.m64.vtir6                                              | LMFDB label of the order        |
| matrix_generators | jsonb | [ [ [ <[ 1, 0, 0, 0 ], 1>, <[ 0, 1, 0, 0 ], 1> ], [ <[ 0, 0, 0, 0 ], 1>, <[ 1, 0, 0, 0 ], 1> ] ], [ [ <[ -1, 1, 1, -1 ], 2>, <[ 0, 0, 0, 0 ], 1> ], [ <[ 0, 0, 0, 0 ], 1>, <[ -1, -1, -1, -1 ], 2> ] ], [ [ <[ 0, 0, 0, 0 ], 1>, <[ -1, 0, 0, 0 ], 1> ], [ <[ 1, 0, 0, 0 ], 1>, <[ 0, 0, 0, 0 ], 1> ] ] ] |                                 |
| relations         | jsonb | [ [ 3, 3, 3, 3 ], [ 2, 2, -3, -3, 2 ], [ -1, 3, 3, 1, 3, 3 ], [ 2, 3, -2, -3, 2, 3, -2, 3 ], [ 3, -2, 3, -1, 2, 3, -2, -3, -1, 2 ], [ -1, 3, 2, -3, -2, -1, 3, 2, 3, -2 ], [ -1, 3, -2, 3, 2, -1, 3, -2, -3, 2 ], [ 1, 1, -2, 1, 2, -1, -1, -2, -1, 2 ], [ -1, -2, 1, 2, 1, 3, -2, -3, 1, 3, 2, -3 ], [ -2, 3, -2, -3, -2, 3, -2, 3, -2, 3, -2, -3 ], [ -1, -2, 3, -2, 3, -2, 1, 2, -3, 2, -3, 2 ] ] |                                 |

*not implemented



## Gamma0 subgroups

The following is the scheme for storing information about Gamma0(alpha). The magma objects are LMFDBGrpSL2CliffGamma0. 

| Column                | Type    | Example               | Notes                                                        |
| --------------------- | ------- | --------------------- | ------------------------------------------------------------ |
| label                 | text    | 3.1_1.m64.vtir6.2LXCx | HashOfOrder_HashOfGL2Element , the hash is taken in base 62 to use the whole alphabet. for [1,0;0,2] |
| defining_element      | jsonb   |                       | for Gamma0 elements this is the matrix that defines it       |
| matrix_generators     | jsonb   |                       | [jsonify(matrix) : matrix in matgens]                        |
| fp_generators | jsonb | | This refers the the generators for an existing SL2O group. |
| relations             | jsonb   |                       | jsonify_relations(FPGroup)                                   |
| order                 |         | 3.1_1.m64.vtir6       |                                                              |
| level_norm | integer |  | pdet(alpha) |
| matrix_coset_representatives | jsonb   |                       | A sequence of matrices giving the coset representatives      |
| fp_coset_representatives | jsonb   |                       | A sequence of words (integers in the range [-n..n], where n is the number of generators) giving the coset representatives                   |
| index                 | integer |                       | the index inside the ambient group                           |
| torsion | jsonb |                       | a sequence of integers from Magmas invariants. |
| rank  | integer |                       | rank of the abelianization                                   |

*not implemented



## Low Index Subgroups

This if for caching low index subgroups of Clifford-Bianchi groups.

| Column                       | Type | Example | Notes |
| ---------------------------- | ---- | ------- | ----- |
| label                        | text |         |       |
| group_label                  |      |         |       |
| index                        |      |         |       |
| generators                   |      |         |       |
| relations                    |      |         |       |
| rank                         |      |         |       |
| fp_coset_representatives     |      |         |       |
| order                        |      |         |       |
| is_gamma0                    |      |         |       |
| matrix_coset_representatives |      |         |       |
| torsion                      |      |         |       |
|                              |      |         |       |
|                              |      |         |       |
|                              |      |         |       |



## Modular Forms

This is for caching modular forms. This is issue #99

| Column       | Type | Example | Notes                                               |
| ------------ | ---- | ------- | --------------------------------------------------- |
| label        | text |         | group_label followed number, sorted by coefficients |
| group_type   |      |         | gamma0 or low_index                                 |
| group_label  |      |         |                                                     |
| index        |      |         |                                                     |
| matrices     |      |         |                                                     |
| coefficients |      |         |                                                     |
|              |      |         |                                                     |
|              |      |         |                                                     |
|              |      |         |                                                     |

