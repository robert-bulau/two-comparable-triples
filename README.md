# Forcing graphs and cyclic completion energy

This repository accompanies the paper **“Forcing Graphs and Cyclic Completion
Energy in Two-Comparable Triple Systems”** by Bulau Robert Nicolae.  The paper
studies pairwise two-comparable subsets of `[n]^3` through monotone partial
matrices, oriented forcing graphs, cyclic completion-energy packing, exact
slack identities, and recursive pressure supports.

The main open problem—whether `G(n) = O(n^(2-delta))` for an absolute positive
`delta`—is not solved.  Every theorem in the paper has a mathematical proof
independent of the code.

## Repository layout

```text
forcing-structures-two-comparable/
├── .gitignore
├── README.md
├── requirements.txt
├── paper/
│   ├── manuscript.tex
│   └── manuscript.pdf
├── results/
│   └── expected_output.txt
└── src/
    └── verify_finite.py
```

No license is included because the author has not selected one.  Add a license
before public redistribution if desired.

## Requirements

- Python 3.10 or later for the finite verification.
- A standard LaTeX installation with `pdflatex` and the packages imported by
  `paper/manuscript.tex` for rebuilding the paper.

The verification uses only the Python standard library.  The
`requirements.txt` file is therefore intentionally dependency-free.

## Reproduce the finite verification

From the repository root, run:

```console
python src/verify_finite.py
```

The program examines all `4^9 = 262144` arrays in
`{empty,1,2,3}^{[3] x [3]}`.  It retains exactly 712 monotone partial matrices,
checks all three coordinate representations of each one, and then verifies
small parameter cases of both explicit constructions.  A successful run must
match `results/expected_output.txt` exactly.

## What the script checks

| Paper location | Computational check |
|---|---|
| Lemma 3.1 and Theorem 3.3 | Common orientation; forcing-graph simplicity; degree, chain-neighbourhood, colour-matching, rainbow-cycle, edge, and wedge assertions |
| Theorem 4.1 | Endpoint exclusion for every row subset in every cyclic representation |
| Lemma 5.1 and Theorem 5.2 | Explicit saturated-inversion images and their disjoint cyclic packing |
| Theorem 6.1 | Pair partition and exact slack identity using rational arithmetic |
| Lemma 7.1 and Theorem 7.3 | Triangular pressure holes, both recursive encodings, and the three projection bounds |
| Proposition 3.5 | Biclique construction for `(2,2)`, `(2,3)`, `(3,3)`, and `(3,4)` |
| Proposition 7.5 | Sharp pressure construction for `(1,1)`, `(2,2)`, `(3,3)`, and `(4,6)` |

The exhaustive order-three run proves only this finite statement.  It is a
consistency audit of the definitions and proofs, not a computer-assisted proof
of any theorem for general `n`.

## Rebuild the manuscript

From the repository root, run:

```console
cd paper
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

The second pass resolves cross-references.  The checked PDF in this repository
was built from the included source.

## Numerical statements

The decimal `0.434258...` in Corollary 6.2 is the rounded value of the exact
algebraic constant `(sqrt(13)-1)/6`; it is not obtained by numerical
optimization.  No floating-point optimization, random search, external solver,
or numerical root-finding is used or cited in the retained manuscript.
