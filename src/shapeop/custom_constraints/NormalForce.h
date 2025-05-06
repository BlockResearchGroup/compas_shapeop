#pragma once

#include "../Force.h"
#include "../Types.h"

namespace ShapeOp {

class NormalForce : public Force {
public:
    NormalForce(const std::vector<std::vector<int>> &faces, double magnitude);
    virtual Vector3 get(const Matrix3X &positions, int id) const override;

private:
    std::vector<std::vector<int>> faces_; // List of faces (each face is a list of vertex indices)
    double magnitude_;                    // Magnitude of the normal force
};

} // namespace ShapeOp
