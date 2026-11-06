# Minimum certified dominating sets in subdivided coronas

Verification code and recorded computational results by **Marek Gałązka and Hanna Wdowicka**.

This repository contains verification code and its recorded outputs. The proofs
in the manuscript are independent of these finite computational checks.

## Run

```bash
python -m pip install -r requirements.txt
python verify.py
```

Run without Python's -O flag: assertions are part of the verification.
The script overwrites `results/` next to itself. Keep a copy of the recorded outputs if comparing runs. No network access is used by the
verification itself. The dependency-installation step may require a network.

## Recorded environment

```text
Python 3.12.14
NumPy 2.3.5
SciPy 1.17.0
random seed 20261106
```

The requirements file pins the recorded NumPy and SciPy versions. A compatible
Python version is needed. Running on another platform may change the reported
elapsed time or the particular optimum found internally by the solver.

## Methods

1. Generate every simple base graph up to vertex renaming, for 1 <= n <= 6.
   The numbers of representatives are checked against 1, 2, 4, 11, 34, 156.
2. Construct the complete subdivided corona, with separately indexed original
   vertices, vertices next to new leaves, leaves, and base-edge subdivisions.
3. Solve two independent full-graph integer models: ordinary domination and
   certified domination. They use only the defining adjacency constraints.
4. Enumerate all base vertex subsets to evaluate the proposed formulas and
   construct every predicted minimum certified dominating set.
5. For 60 small test cases, enumerate actual full-graph vertex subsets in
   increasing size and compare the complete optimum families, not just values.
6. For bipartite bases, independently implement integral augmenting paths and
   check the bounded edge set, minimum-cut choice, and certified-set formula.
7. Check the degree criterion and the characterization by independent subsets.
8. Check the path, cycle and complete-base enumeration formulas up to n = 13.

## Independent certified-domination model

For each full-graph vertex v, x_v is 1 exactly when v is selected. A second
binary variable y_v distinguishes the alternative of at least two outside
neighbors. Write d_v for its degree and q_v = sum_{u adjacent to v} x_u.
The model minimizes sum_v x_v subject to:

```text
x_v + q_v >= 1
q_v + 2*y_v <= d_v
d_v*x_v - q_v - d_v*y_v <= 0
y_v <= x_v
```

When x_v=1 and y_v=0, q_v=d_v, giving zero outside neighbors. When x_v=1 and
y_v=1, q_v<=d_v-2, giving at least two. When x_v=0, the extra conditions impose
no restriction beyond domination. Thus this model does not assume the paper's
reduction or any forced vertex. Returned sets are also checked directly.

## Recorded results

```text
208 cases with one new leaf at each base vertex, n <= 6
1020 assignments r_v in {1,2,3} on representatives with n <= 4
60 additional seeded cases, 7 <= n <= 12, r_v in {1,2,3,4}
1288 total test cases; 2576 integer-model solves; all passed
764 bipartite cases checked by the flow implementation
60 full-vertex-subset enumerations of all minimum solutions
34 additional formula checks for paths, cycles and complete bases
```

The test groups overlap: 1288 counts test cases, not pairwise distinct graphs.
The weighted cases include every assignment to the fixed labeling of each
unlabelled base representative. Some assignments are equivalent under vertex
renaming. The 34 additional enumeration checks evaluate the base-subset
formula, rather than solving an independent full-graph model at these sizes.

## Files

```text
verify.py                  complete verification program
requirements.txt           pinned numerical-library versions
results/summary.json       totals, environment, seed, elapsed time
results/cases.json         inputs, optimum values, minimizing base subsets,
                          and bounded-edge certificates for all test cases
results/cases.csv          compact tabular summary
results/enumeration.json   family counts and optimum values
```

## Conventions

Base vertices are 0,...,n-1; edges are unordered pairs. A subset of base
vertices is stored as an integer mask: bit v is 1 precisely when v is selected.
In `construct`, base vertices come first. For each base vertex in order, its new support vertices and leaves are appended in alternating pairs; the base-edge subdivision vertices follow in the input edge order.

AI assistance was used in developing this code. No computational check is
presented as a proof of a universal mathematical statement.

## Citation and scope

Please cite this software by its exact Git commit when using the recorded checks.
Machine-readable citation metadata are provided in `CITATION.cff`.
The associated manuscript is **Minimum certified dominating sets in subdivided coronas**,
by Marek Gałązka and Hanna Wdowicka (2026). No publication DOI has been assigned.

The repository contains code, computational results, and reproduction instructions.
The manuscript and submission correspondence are maintained separately.
