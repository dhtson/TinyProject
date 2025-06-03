#include "LinearSystem.h"

#include <cassert>
#include <cmath>

LinearSystem::LinearSystem() : mSize(0), mpA(nullptr), mpb(nullptr), mOwnsPointers(false) {}
LinearSystem::LinearSystem(const LinearSystem& other) : mSize(0), mpA(nullptr), mpb(nullptr), mOwnsPointers(false) {
    mSize = other.mSize;
    mOwnsPointers = other.mOwnsPointers;
    
    if (mOwnsPointers) {
        if (other.mpA) mpA = new Matrix(*(other.mpA)); else mpA = nullptr;
        if (other.mpb) mpb = new Vector(*(other.mpb)); else mpb = nullptr;
    } else {
        mpA = other.mpA;
        mpb = other.mpb;
    }
}

LinearSystem::LinearSystem(Matrix& A, Vector& b, bool copyData) {
    assert(A.GetNumberOfRows() == A.GetNumberOfColumns());
    assert(A.GetNumberOfRows() == b.GetSize());
    mSize = A.GetNumberOfRows();

    mOwnsPointers = copyData;
    if (copyData) {
        mpA = new Matrix(A);
        mpb = new Vector(b);
    } else {
        mpA = &A;
        mpb = &b;
    }
}

LinearSystem::~LinearSystem() {
    if (mOwnsPointers) {
        delete mpA;
        delete mpb;
    }
}

Matrix LinearSystem::GetMatrix() const {
    assert(mpA != nullptr);
    return *mpA;
}
Vector LinearSystem::GetRHSVector() const {
    assert(mpb != nullptr);
    return *mpb;
}

Vector LinearSystem::Solve() {
    assert(mpA != nullptr && mpb != nullptr);
    assert(mSize > 0);

    Matrix A_work = *mpA;
    Vector b_work = *mpb;

    for (int k = 0; k < mSize -1; ++k) {
        int max_row = k;
        double max_val = std::abs(A_work(k + 1, k + 1));

        for (int i = k + 1; i < mSize; ++i) {
            if (std::abs(A_work(i + 1, k + 1)) > max_val) {
                max_val = std::abs(A_work(i + 1, k + 1));
                max_row = i;
            }
        }

        if (max_row != k) {
            for(int col_idx = 0; col_idx < mSize; ++col_idx) {
                double temp = A_work(k+1, col_idx+1);
                A_work(k+1, col_idx+1) = A_work(max_row+1, col_idx+1);
                A_work(max_row+1, col_idx+1) = temp;
            }
            
            double temp_b = b_work(k+1);
            b_work(k+1) = b_work(max_row+1);
            b_work(max_row+1) = temp_b;
        }
        
        if (std::abs(A_work(k + 1, k + 1)) < 1e-9) {
            std::cerr << "Warning: Matrix is singular or nearly singular during Gaussian elimination." << std::endl;
            return Vector(mSize);
        }

        for (int i = k + 1; i < mSize; ++i) {
            double factor = A_work(i + 1, k + 1) / A_work(k + 1, k + 1);

            for (int j = k; j < mSize; ++j) {
                A_work(i + 1, j + 1) -= factor * A_work(k + 1, j + 1);
            }

            b_work(i + 1) -= factor * b_work(k + 1);
        }
    }

    Vector x(mSize);

    for (int i = mSize - 1; i >= 0; --i) {
        double sum_ax = 0.0;

        for (int j = i + 1; j < mSize; ++j) {
            sum_ax += A_work(i + 1, j + 1) * x(j + 1);
        }

        if (std::abs(A_work(i + 1, i + 1)) < 1e-9) {
            std::cerr << "Warning: Matrix is singular (division by zero in back substitution)." << std::endl;
            return Vector(mSize);
        }

        x(i + 1) = (b_work(i + 1) - sum_ax) / A_work(i + 1, i + 1);
    }

    return x;
}
