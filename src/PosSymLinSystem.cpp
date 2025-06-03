#include "PosSymLinSystem.h"

#include <cassert>
#include <cmath>

PosSymLinSystem::PosSymLinSystem(Matrix& A, Vector& b, bool copyData)
    : LinearSystem(A, b, copyData) {
    assert(IsSymmetric(*mpA));
}

bool PosSymLinSystem::IsSymmetric(const Matrix& A) const {
    if (A.GetNumberOfRows() != A.GetNumberOfColumns()) {
        return false;
    }
    for (int i = 0; i < A.GetNumberOfRows(); ++i) {
        for (int j = i + 1; j < A.GetNumberOfColumns(); ++j) {
            if (std::abs(A(i + 1, j + 1) - A(j + 1, i + 1)) > 1e-9) {
                return false;
            }
        }
    }
    return true;
}

Vector PosSymLinSystem::Solve() {
    assert(mpA != nullptr && mpb != nullptr);
    assert(mSize > 0);

    Matrix& A = *mpA;
    Vector& b = *mpb;

    Vector x(mSize);
    Vector r = b - (A * x);
    Vector p = r;
    double rsold = r.Norm(2) * r.Norm(2);

    int max_iterations = mSize * 2;
    double tolerance = 1e-9;

    for (int i = 0; i < max_iterations; ++i) {
        Vector Ap = A * p;
        double alpha_num = rsold;
        double alpha_den = 0.0;
        for(int k=0; k<p.GetSize(); ++k) alpha_den += p[k] * Ap[k];

        if (std::abs(alpha_den) < 1e-12) {
            std::cerr << "Warning (CG): Denominator for alpha is near zero. Matrix may not be positive definite or system is ill-conditioned." << std::endl;
            break;
        }
        double alpha = alpha_num / alpha_den;

        x = x + (p * alpha);
        r = r - (Ap * alpha);

        double rsnew = r.Norm(2) * r.Norm(2);

        if (std::sqrt(rsnew) < tolerance) {
            break;
        }

        double beta = rsnew / rsold;
        p = r + (p * beta);
        rsold = rsnew;
    }
    if (std::sqrt(rsold) >= tolerance) {
        std::cout << "Conjugate Gradient did not converge within " << max_iterations << " iterations." << std::endl;
        std::cout << "Final residual norm: " << std::sqrt(rsold) << std::endl;
    }

    return x;
}
