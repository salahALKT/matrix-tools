# Matrices-tools

[![tests](https://github.com/salahALKT/Matrices-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/salahALKT/Matrices-tools/actions/workflows/ci.yml)

A C++ matrix library written from scratch, with no dependencies: it uses only the standard library. The
`Matrix<T>` template lives in `Matrix.h` / `Matrix.tpp`; compile `overloads.cpp` alongside your code.
It covers the core of a linear algebra course, from Gaussian elimination to the singular value decomposition,
and works with complex numbers throughout.

```cpp
#include "Matrix.h"
using C = std::complex<double>;

Matrix<C> M({{1. + C(0, 1), 2, 4},
             {2. + C(0, 2), 4, 8},
             {1,            0, 0}});

Matrix<C> P(1, 1, myZero()), Q(1, 1, myZero());
Matrix<C> D = M.SVD(P, Q);      // M = P * D * Q, with P and Q unitary
Matrix<C> K = M.kernel();       // columns form a basis of ker(M)
C det = M.determinant();
```

## What it can do

| Area | Methods |
|---|---|
| Basics | `+`, `-`, `*`, scaling, `transpose`, `conjugate`, `~` (conjugate transpose), `identity`, `is_hermitian` |
| Elimination | `gauss_elimination` (row echelon or fully reduced form), `inverse`, `determinant` |
| Subspaces | `kernel`, `image`, `orthogonal_base`, `complete_base` (extend to a basis of the whole space) |
| Orthogonality | `gram_schmidt` (with or without normalization) |
| Decompositions | `QRD` (QR decomposition), `QR_algo` (eigenvalues by the QR algorithm), `diagonal_base` (eigenvectors), `SVD` |

## Building and testing

Requires a C++17 compiler (`g++` or `clang++`) and `make`.

```bash
make test     # build and run the tests
make demo     # build the demo in main.cpp
./demo
```

The [tests](tests/test_matrix.cpp) check each operation by its defining property rather than against stored
answers: `A * A⁻¹ = I`, `A * ker(A) = 0`, `Q` unitary and `QR = A`, `P` and `Q` unitary and `PDQ = M`. They
cover square, singular, tall and wide matrices with complex entries, and run on every push through GitHub
Actions.

## Using other number types

`Matrix<T>` is a template. The scalar type `T` needs `+`, `*`, `==` and unary `-`, plus a few helper functions,
declared in [`overloads.h`](overloads.h):

| Function | Meaning | Needed for |
|---|---|---|
| `myZero()`, `myUnit()` | 0 and 1 | everything |
| `~x` | complex conjugate | conjugate transpose |
| `myTo_string(x)` | text form | printing |
| `myInverse(x)`, `isZero(x)` | 1/x, test for zero | elimination and everything built on it |
| `mySqrt(x)` | square root of the modulus | normalization, SVD |

These are provided for `std::complex<double>`, which also covers real matrices. Another type (exact rationals,
for example) needs its own versions of these functions.

## Limitations

- **Eigenvalues:** they come from 50 iterations of the unshifted QR algorithm. That works well when eigenvalues
  are well separated, but converges slowly when two have nearly the same modulus.
- **Tolerance:** values below about 1e-5 in modulus are treated as zero (`isZero`). Fine for well-conditioned
  matrices, not for ill-conditioned ones.
- **Speed:** written for clarity rather than speed. Matrices are `std::vector<std::vector<T>>` and many
  operations copy.

---

*Cleanup and presentation (this README, the tests, CI and several bug fixes) were done with the help of Claude
(Anthropic) in 2026. The library itself is my own work.*
