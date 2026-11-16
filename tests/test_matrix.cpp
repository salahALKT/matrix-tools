// Tests for the Matrix class. Each test checks a defining property of the result
// (A * A^-1 = I, A * ker(A) = 0, Q unitary and QR = A, ...) rather than exact values.
//
// Build and run with: make test

#include "../Matrix.h"
#include "../overloads.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using C = std::complex<double>;
using M = Matrix<C>;

static int failures = 0;
static int checks = 0;

static void check(bool ok, std::string const& what) {
    ++checks;
    if (!ok) {
        ++failures;
        std::cout << "  FAIL: " << what << std::endl;
    }
}

// Largest entry-wise difference between two matrices of the same size.
static double max_diff(M const& A, M const& B) {
    if (A.numRows() != B.numRows() || A.numCols() != B.numCols()) return INFINITY;
    double d = 0;
    for (size_t i = 0; i < A.numRows(); ++i)
        for (size_t j = 0; j < A.numCols(); ++j)
            d = std::max(d, std::abs(A.loc(i, j) - B.loc(i, j)));
    return d;
}

static bool close(M const& A, M const& B, double tol = 1e-8) { return max_diff(A, B) < tol; }

static bool is_unitary(M const& U) { return close(U * (~U), M::identity(U.numRows())) && close((~U) * U, M::identity(U.numCols())); }

static bool is_upper_triangular(M const& R) {
    for (size_t i = 0; i < R.numRows(); ++i)
        for (size_t j = 0; j < i && j < R.numCols(); ++j)
            if (std::abs(R.loc(i, j)) > 1e-8) return false;
    return true;
}

static bool is_diagonal(M const& D) {
    for (size_t i = 0; i < D.numRows(); ++i)
        for (size_t j = 0; j < D.numCols(); ++j)
            if (i != j && std::abs(D.loc(i, j)) > 1e-8) return false;
    return true;
}

// Test matrices
static const C i1(0, 1);
static M square3() { return M({{C(2, 1), 1, 0}, {1, C(3, 0), i1}, {0, -i1, 4}}); }   // invertible, complex
static M singular3() { return M({{C(1, 1), 2, 4}, {C(2, 2), 4, 8}, {1, 0, 0}}); }   // rank 2 (row 2 = 2 * row 1)
static M tall() { return M({{1, 2}, {3, C(4, 1)}, {5, 6}}); }                        // 3 x 2
static M wide() { return M({{1, 0, 2}, {0, C(1, -1), 1}}); }                        // 2 x 3

static void test_arithmetic() {
    M A({{1, 2}, {3, 4}}), B({{5, 6}, {7, 8}});
    check(close(A * B, M({{19, 22}, {43, 50}})), "product of 2x2 matrices");
    check(close(A + B, M({{6, 8}, {10, 12}})), "sum");
    check(close(A - A, M(2, 2, myZero())), "A - A = 0");
    check(close(A * M::identity(2), A), "A * I = A");
    M Z = square3();
    check(close(~(~Z), Z), "conjugate transpose twice gives back A");
    check(std::abs((~Z).loc(0, 1) - std::conj(Z.loc(1, 0))) < 1e-12, "conjugate transpose entries");
}

static void test_inverse_and_determinant() {
    M A = square3();
    check(close(A * A.inverse(), M::identity(3)), "A * inverse(A) = I");
    check(close(A.inverse() * A, M::identity(3)), "inverse(A) * A = I");
    check(std::abs(M({{2, 1}, {1, 3}}).determinant() - C(5)) < 1e-9, "determinant of a 2x2");
    check(std::abs(M({{1, 2, 3}, {0, 1, 4}, {5, 6, 0}}).determinant() - C(1)) < 1e-9, "determinant of a 3x3");
    check(std::abs(singular3().determinant()) < 1e-9, "determinant of a singular matrix is 0");
}

static void test_kernel_and_image() {
    M S = singular3();
    M K = S.kernel();
    check(K.numCols() == 1, "kernel of a rank-2 3x3 matrix has dimension 1");
    check(close(S * K, M(3, K.numCols(), myZero())), "A * kernel(A) = 0");
    M W = wide();
    M KW = W.kernel();
    check(close(W * KW, M(2, KW.numCols(), myZero())), "kernel of a wide matrix");
}

static void test_gram_schmidt() {
    M G = square3().gram_schmidt(true);
    check(is_unitary(G), "Gram-Schmidt gives orthonormal columns");
}

static void test_qr() {
    for (M A : {square3(), M({{1, 2}, {3, 4}})}) {
        M Q(1, 1, myZero()), R(1, 1, myZero());
        A.QRD(Q, R);
        check(is_unitary(Q), "QR: Q is unitary");
        check(is_upper_triangular(R), "QR: R is upper triangular");
        check(close(Q * R, A), "QR: Q * R = A");
    }
}

static void test_svd() {
    std::vector<std::pair<std::string, M>> cases = {
        {"invertible 3x3", square3()},
        {"singular 3x3", singular3()},
        {"tall 3x2", tall()},
        {"wide 2x3", wide()},
    };
    for (auto const& [name, A] : cases) {
        M P(1, 1, myZero()), Q(1, 1, myZero());
        M D = A.SVD(P, Q);
        check(is_diagonal(D), "SVD (" + name + "): D is diagonal");
        check(is_unitary(P), "SVD (" + name + "): P is unitary");
        check(is_unitary(Q), "SVD (" + name + "): Q is unitary");
        check(close(P * D * Q, A), "SVD (" + name + "): P * D * Q = A");
    }
}

int main() {
    std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"arithmetic", test_arithmetic},
        {"inverse and determinant", test_inverse_and_determinant},
        {"kernel and image", test_kernel_and_image},
        {"Gram-Schmidt", test_gram_schmidt},
        {"QR decomposition", test_qr},
        {"SVD", test_svd},
    };
    for (auto const& [name, run] : tests) {
        int before = failures;
        try {
            run();
        } catch (std::exception const& e) {
            ++failures;
            std::cout << "  FAIL: exception: " << e.what() << "\n";
        }
        std::cout << (failures == before ? "ok    " : "FAILED ") << name << std::endl;
    }
    std::cout << "\n" << checks << " checks, " << failures << " failures\n";
    return failures == 0 ? 0 : 1;
}
