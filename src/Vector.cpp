#include "Vector.h"

#include <cassert>
#include <cmath>
#include <iomanip>

Vector::Vector() : mSize(0), mData(nullptr) {}

Vector::Vector(int size) {
    assert(size >= 0);
    mSize = size;
    if (mSize > 0) {
        mData = new double[mSize];
        for (int i = 0; i < mSize; ++i) {
            mData[i] = 0.0;
        }
    } else {
        mData = nullptr;
    }
}

Vector::Vector(const Vector& otherVector) {
    mSize = otherVector.mSize;
    if (mSize > 0) {
        mData = new double[mSize];
        for (int i = 0; i < mSize; i++) {
            mData[i] = otherVector.mData[i];
        }
    } else {
        mData = nullptr;
    }
}

Vector::~Vector() {
    delete[] mData;
}

int Vector::GetSize() const {
    return mSize;
}

Vector& Vector::operator=(const Vector& otherVector) {
    if (this == &otherVector) {
        return *this;
    }
    if (mSize != otherVector.mSize && mData != nullptr) {
        delete[] mData;
        mData = nullptr;
    }
    mSize = otherVector.mSize;
    if (mSize > 0) {
        if (mData == nullptr) {
            mData = new double[mSize];
        }
        for (int i = 0; i < mSize; i++) {
            mData[i] = otherVector.mData[i];
        }
    } else {
        mData = nullptr;
    }
    return *this;
}

double& Vector::operator[](int index) {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

const double& Vector::operator[](int index) const {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

double& Vector::operator()(int index) {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

const double& Vector::operator()(int index) const {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

Vector Vector::operator+() const {
    return *this;
}

Vector Vector::operator-() const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = -mData[i];
    }
    return result;
}

Vector Vector::operator+(const Vector& v1) const {
    assert(mSize == v1.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] + v1.mData[i];
    }
    return result;
}

Vector Vector::operator-(const Vector& v1) const {
    assert(mSize == v1.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] - v1.mData[i];
    }
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] * scalar;
    }
    return result;
}

Vector operator*(double scalar, const Vector& v) {
    return v * scalar;
}

std::ostream& operator<<(std::ostream& os, const Vector& v) {
    os << "(";
    for (int i = 0; i < v.mSize; ++i) {
        os << std::fixed << std::setprecision(4) << v.mData[i];
        if (i < v.mSize - 1) {
            os << ", ";
        }
    }
    os << ")";
    return os;
}

double Vector::Norm(int p) const {
    assert(p >= 1);
    double sum = 0.0;
    
    if (p == 1) {
        for (int i = 0; i < mSize; ++i) {
            sum += std::abs(mData[i]);
        }
        return sum;
    }

    if (p == 2) {
        for (int i = 0; i < mSize; ++i) {
            sum += mData[i] * mData[i];
        }
        return std::sqrt(sum);
    }

    for (int i = 0; i < mSize; ++i) {
        sum += std::pow(std::abs(mData[i]), p);
    }

    return std::pow(sum, 1.0/p);
}
