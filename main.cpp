// A short tour of the Matrix class on a small complex matrix.
// Build and run with: make demo && ./demo

#include "Matrix.h"
#include <iostream>

using C = std::complex<double>;

int main() {
    const C i(0, 1);
    Matrix<C> M({{1. + i, 2, 4},
                 {2. + 2. * i, 4, 8},   // twice the first row, so M has rank 2
                 {1, 0, 0}});

    std::cout << "M =\n" << M << "\n";

    std::cout << "det(M) = " << myTo_string(M.determinant()) << "\n\n";

    std::cout << "A basis of ker(M):\n" << M.kernel() << "\n";

    // Singular value decomposition: M = P D Q with P, Q unitary and D diagonal.
    Matrix<C> P(1, 1, myZero()), Q(1, 1, myZero());
    Matrix<C> D = M.SVD(P, Q);
    std::cout << "Singular values (D):\n" << D << "\n";
    std::cout << "P * D * Q (should equal M):\n" << P * D * Q << "\n";

    // An invertible matrix and its inverse.
    Matrix<C> A({{2, 1}, {1, 3. + i}});
    std::cout << "A =\n" << A << "\n";
    std::cout << "A^-1 =\n" << A.inverse() << "\n";
    std::cout << "A * A^-1 =\n" << A * A.inverse();
}
