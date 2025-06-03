#pragma once

#include "LinearSystem.h"

class PosSymLinSystem : public LinearSystem {
public:
    PosSymLinSystem(Matrix& A, Vector& b, bool copyData = true);
    ~PosSymLinSystem() override = default;

    Vector Solve() override;

private:
    bool IsSymmetric(const Matrix& A) const;
};
