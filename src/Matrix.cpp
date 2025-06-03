#include "Matrix.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>

void Matrix::AllocateMemory() {
    if (mNumRows > 0 && mNumCols > 0) {
        mData = new double*[mNumRows];

        for (int i = 0; i < mNumRows; ++i) {
            mData[i] = new double[mNumCols];

            for (int j = 0; j < mNumCols; ++j) {
                mData[i][j] = 0.0;
            }
        }
    } else {
        mData = nullptr;
    }
}

void Matrix::DeallocateMemory() {
    if (mData) {
        for (int i = 0; i < mNumRows; ++i) {
            delete[] mData[i];
        }
        delete[] mData;
        mData = nullptr;
    }
}

Matrix::Matrix() : mNumRows(0), mNumCols(0), mData(nullptr) {}
Matrix::Matrix(int numRows, int numCols) : mNumRows(numRows), mNumCols(numCols) {
    assert(numRows >= 0 && numCols >= 0);
    AllocateMemory();
}

Matrix::Matrix(const Matrix& otherMatrix) {
    mNumRows = otherMatrix.mNumRows;
    mNumCols = otherMatrix.mNumCols;
    AllocateMemory();
    if (mData) {
        for (int i = 0; i < mNumRows; i++) {
            for (int j = 0; j < mNumCols; j++) {
                mData[i][j] = otherMatrix.mData[i][j];
            }
        }
    }
}

Matrix::~Matrix() {
    DeallocateMemory();
}

int Matrix::GetNumberOfRows() const { return mNumRows; }
int Matrix::GetNumberOfColumns() const { return mNumCols; }

Matrix& Matrix::operator=(const Matrix& otherMatrix) {
    if (this == &otherMatrix) {
        return *this;
    }
    DeallocateMemory();
    mNumRows = otherMatrix.mNumRows;
    mNumCols = otherMatrix.mNumCols;
    AllocateMemory();
    if (mData) {
        for (int i = 0; i < mNumRows; i++) {
            for (int j = 0; j < mNumCols; j++) {
                mData[i][j] = otherMatrix.mData[i][j];
            }
        }
    }
    return *this;
}

double& Matrix::operator()(int i, int j) {
    assert(i >= 1 && i <= mNumRows); 
    assert(j >= 1 && j <= mNumCols); 
    return mData[i - 1][j - 1];       
}

const double& Matrix::operator()(int i, int j) const {
    assert(i >= 1 && i <= mNumRows); 
    assert(j >= 1 && j <= mNumCols); 
    return mData[i - 1][j - 1];       
}

Matrix Matrix::operator+() const {
    return *this;
}

Matrix Matrix::operator-() const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = -mData[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator+(const Matrix& m1) const {
    assert(mNumRows == m1.mNumRows && mNumCols == m1.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] + m1.mData[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& m1) const {
    assert(mNumRows == m1.mNumRows && mNumCols == m1.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] - m1.mData[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& m1) const {
    assert(mNumCols == m1.mNumRows);
    Matrix result(mNumRows, m1.mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < m1.mNumCols; j++) {
            double sum = 0.0;
            for (int k = 0; k < mNumCols; k++) {
                sum += mData[i][k] * m1.mData[k][j];
            }
            result.mData[i][j] = sum;
        }
    }
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] * scalar;
        }
    }
    return result;
}

Vector Matrix::operator*(const Vector& v) const {
    assert(mNumCols == v.GetSize());
    Vector result(mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        double sum = 0.0;
        for (int j = 0; j < mNumCols; j++) {
            sum += mData[i][j] * v[j];
        }
        result[i] = sum;
    }
    return result;
}

Matrix operator*(double scalar, const Matrix& m) {
    return m * scalar;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    os << "[";
    for (int i = 0; i < m.mNumRows; ++i) {
        if (i > 0) os << " ";
        os << "[";
        for (int j = 0; j < m.mNumCols; ++j) {
            os << std::fixed << std::setprecision(4) << std::setw(8) << m.mData[i][j];
            if (j < m.mNumCols - 1) {
                os << ", ";
            }
        }
        os << "]";
        if (i < m.mNumRows - 1) {
            os << "\n";
        }
    }
    os << "]";
    return os;
}

Matrix Matrix::Transpose() const {
    Matrix result(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            result.mData[j][i] = mData[i][j];
        }
    }
    return result;
}

Matrix GetSubMatrix(const Matrix& mat, int p, int q) {
    int n = mat.GetNumberOfRows();
    Matrix temp(n - 1, n - 1);
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        if (row == p) continue;
        j = 0;
        for (int col = 0; col < n; col++) {
            if (col == q) continue;
            temp(i + 1, j + 1) = mat(row + 1, col + 1);
            j++;
        }
        i++;
    }
    return temp;
}

double Matrix::Determinant() const {
    assert(mNumRows == mNumCols);
    assert(mNumRows > 0);

    int n = mNumRows;
    double det = 0;

    if (n == 1) {
        return mData[0][0];
    }
    if (n == 2) {
        return mData[0][0] * mData[1][1] - mData[0][1] * mData[1][0];
    }

    Matrix subMatrix(n - 1, n - 1);
    int sign = 1;

    for (int f = 0; f < n; f++) {
        int subi = 0;
        for (int i = 1; i < n; i++) {
            int subj = 0;
            for (int j = 0; j < n; j++) {
                if (j == f) continue;
                subMatrix.mData[subi][subj] = mData[i][j];
                subj++;
            }
            subi++;
        }
        det += sign * mData[0][f] * subMatrix.Determinant();
        sign = -sign;
    }
    return det;
}


Matrix Matrix::Inverse() const {
    assert(mNumRows == mNumCols);
    assert(mNumRows > 0);
    double det = Determinant();
    assert(std::abs(det) > 1e-9);

    int n = mNumRows;
    Matrix inverse(n, n);

    if (n == 1) {
        inverse.mData[0][0] = 1.0 / mData[0][0];
        return inverse;
    }

    Matrix adj(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix subMat = GetSubMatrix(*this, i, j);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj(j + 1, i + 1) = sign * subMat.Determinant();
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse.mData[i][j] = adj.mData[i][j] / det;
        }
    }
    return inverse;
}


Matrix Matrix::PseudoInverse() const {
    Matrix At = this->Transpose();
    if (mNumRows >= mNumCols) {
        Matrix AtA = At * (*this);

        if (std::abs(AtA.Determinant()) < 1e-9 && AtA.GetNumberOfRows() > 0) {
            std::cerr << "Warning: (A^T * A) is singular in PseudoInverse. Returning empty Matrix." << std::endl;
            return Matrix(0,0);
        }

        if (AtA.GetNumberOfRows() == 0) return Matrix(0,0);
        return AtA.Inverse() * At;
    } else {
        Matrix AAt = (*this) * At;
        if (std::abs(AAt.Determinant()) < 1e-9 && AAt.GetNumberOfRows() > 0) {
            std::cerr << "Warning: (A * A^T) is singular in PseudoInverse. Returning empty Matrix." << std::endl;
            return Matrix(0,0);
        }

        if (AAt.GetNumberOfRows() == 0) return Matrix(0,0);
        return At * AAt.Inverse();
    }
}

void Matrix::SwapRows(Matrix& A, Vector& b, int r1, int r2) {
    assert(r1 >= 0 && r1 < A.GetNumberOfRows());
    assert(r2 >= 0 && r2 < A.GetNumberOfRows());
    assert(A.GetNumberOfRows() == b.GetSize());

    for (int j = 0; j < A.GetNumberOfColumns(); ++j) {
        double temp = A.mData[r1][j];
        A.mData[r1][j] = A.mData[r2][j];
        A.mData[r2][j] = temp;
    }

    double temp_b = b[r1];
    b[r1] = b[r2];
    b[r2] = temp_b;
}

void Matrix::SwapRowsMatrixOnly(Matrix& A, int r1, int r2) {
    assert(r1 >= 0 && r1 < A.GetNumberOfRows());
    assert(r2 >= 0 && r2 < A.GetNumberOfRows());
    
    for (int j = 0; j < A.GetNumberOfColumns(); ++j) {
        double temp = A.mData[r1][j];
        A.mData[r1][j] = A.mData[r2][j];
        A.mData[r2][j] = temp;
    }
}
