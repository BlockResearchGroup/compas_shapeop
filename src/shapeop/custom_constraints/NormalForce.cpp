#include "NormalForce.h"

namespace ShapeOp {

NormalForce::NormalForce(const std::vector<std::vector<int>> &faces, double magnitude)
    : faces_(faces), magnitude_(magnitude) {}

Vector3 NormalForce::get(const Matrix3X &positions, int id) const {
    Vector3 accumulatedNormal = Vector3::Zero();

    // Iterate over all faces
    for (const auto &face : faces_) {
        if (std::find(face.begin(), face.end(), id) != face.end()) {
            // Calculate the face normal
            Vector3 v0 = positions.col(face[0]);
            Vector3 v1 = positions.col(face[1]);
            Vector3 v2 = positions.col(face[2]);
            Vector3 normal = (v1 - v0).cross(v2 - v0);

            // Weight the normal by the area of the triangle (length of the cross product)
            double area = normal.norm();
            if (area > 0) {
                normal.normalize();
                accumulatedNormal += area * normal;
            }
        }
    }

    // Normalize the accumulated normal
    if (accumulatedNormal.norm() > 0) {
        accumulatedNormal.normalize();
    }

    // Apply the force in the direction of the normal
    return magnitude_ * accumulatedNormal;
}

} // namespace ShapeOp
