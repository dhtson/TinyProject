#pragma once

#include "Matrix.h"
#include "Vector.h"

class LinearSystem {
protected:
    int mSize;
    Matrix* mpA;
    Vector* mpb;
    bool mOwnsPointers;

private:
    LinearSystem();
    LinearSystem(const LinearSystem& other);


public:
    LinearSystem(Matrix& A, Vector& b, bool copyData = true);
    virtual ~LinearSystem();

    virtual Vector Solve();

    Matrix GetMatrix() const;
    Vector GetRHSVector() const;
};
