#pragma once

#include "Vector.h"

#include <iostream>

class Matrix {
private:
    int mNumRows;
    int mNumCols;
    double** mData;

    void AllocateMemory();
    void DeallocateMemory();

public:
    Matrix();
    Matrix(int numRows, int numCols);
    Matrix(const Matrix& otherMatrix);
    ~Matrix();

    int GetNumberOfRows() const;
    int GetNumberOfColumns() const;

    Matrix& operator=(const Matrix& otherMatrix);
    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;

    Matrix operator+() const;
    Matrix operator-() const;

    Matrix operator+(const Matrix& m1) const;
    Matrix operator-(const Matrix& m1) const;
    Matrix operator*(const Matrix& m1) const;
    Matrix operator*(double scalar) const;
    Vector operator*(const Vector& v) const;

    friend Matrix operator*(double scalar, const Matrix& m);
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

    Matrix Transpose() const;
    double Determinant() const;
    Matrix Inverse() const;
    Matrix PseudoInverse() const;

    static void SwapRows(Matrix& A, Vector& b, int r1, int r2);
    static void SwapRowsMatrixOnly(Matrix& A, int r1, int r2);
};
