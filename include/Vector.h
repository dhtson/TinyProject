#pragma once

#include <iostream>

class Vector {
private:
    int mSize;
    double* mData;

public:
    Vector();
    Vector(int size);
    Vector(const Vector& otherVector);
    ~Vector();

    int GetSize() const;

    Vector& operator=(const Vector& otherVector);
    double& operator[](int index);
    const double& operator[](int index) const;
    double& operator()(int index);
    const double& operator()(int index) const;

    Vector operator+() const;
    Vector operator-() const;

    Vector operator+(const Vector& v1) const;
    Vector operator-(const Vector& v1) const;
    Vector operator*(double scalar) const;

    friend Vector operator*(double scalar, const Vector& v);
    friend std::ostream& operator<<(std::ostream& os, const Vector& v);

    double Norm(int p = 2) const;
};
