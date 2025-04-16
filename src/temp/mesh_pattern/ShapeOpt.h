//
// Created by ziqwang on 2019-04-13.
//

#ifndef MESHPATTERN_SHAPEOPT_H
#define MESHPATTERN_SHAPEOPT_H

#include "iglMesh.h"
#include <libShapeOp/src/Common.h>
#include <libShapeOp/src/Solver.h>
#include <libShapeOp/src/Constraint.h>
#include <memory>
using std::make_shared;
using std::shared_ptr;
using Eigen::Vector3d;

class ShapeOpt{
public:
    ShapeOp::Matrix3X V_;
    ShapeOp::Solver solver;
    iglMesh mesh_;

public:

    void addPlanarConstraint(double weight);

    void addRegularConstraint(double weight);

    void addBoundaryConstraint(double weight);

    void runShapeOp(int iterations);

    Vector3d getCentroid(int fid);

    Vector3d getNormal(int fid);

public:

    void loadMesh(iglMesh &_mesh);

    iglMesh getMesh();
};


#endif //MESHPATTERN_SHAPEOPT_H
