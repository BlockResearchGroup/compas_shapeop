//
// Created by ziqwang on 2019-04-13.
//
// ShapeOpt.cpp - Implementation of mesh optimization algorithms using ShapeOp library
// 
// This file implements a set of architectural geometry optimization techniques that allow
// transforming general meshes into meshes with specific properties that are important for
// architectural fabrication and construction, such as planar faces, regular polygonal faces,
// and preserving boundary conditions.
//

#include "ShapeOpt.h"
#include <igl/boundary_loop.h>

/**
 * Loads a mesh from the iglMesh format and prepares the ShapeOp solver.
 * 
 * @param _mesh Input mesh in iglMesh format
 */
void ShapeOpt::loadMesh(iglMesh &_mesh)
{
    mesh_ = _mesh;
    // Convert mesh vertices to ShapeOp's expected format (transposed)
    V_ = mesh_.V_.transpose();
    // Initialize a new ShapeOp solver
    solver = ShapeOp::Solver();
    // Set the initial points for the solver
    solver.setPoints(V_);
    //std::cout << V_ << std::endl;
}

/**
 * Adds planarity constraints to all faces in the mesh.
 * 
 * Planarity is critical for architectural panels fabrication, as non-planar panels
 * are expensive to manufacture. This function enforces that all vertices of each face
 * lie on a common plane.
 * 
 * @param weight The weight of the constraint (higher values enforce stronger planarity)
 */
void ShapeOpt::addPlanarConstraint(double weight)
{
    //in plane constraint
    std::vector<int> id_vector;
    for(int id = 0; id < mesh_.faces_.size(); id++)
    {
        id_vector.clear();
        // Collect all vertex indices for this face
        for(int verID : mesh_.faces_[id]){
            id_vector.push_back(verID);
        }
        // Create and add a plane constraint for this face
        shared_ptr<ShapeOp::Constraint> constraint = make_shared<ShapeOp::PlaneConstraint>(id_vector, weight, V_);
        solver.addConstraint(constraint);
    }
}

/**
 * Adds constraints that transform mesh faces into regular polygons.
 * 
 * This is useful for architectural patterns where regular polygons provide
 * aesthetic and structural advantages. The function creates regular polygonal
 * templates aligned with each face's normal direction and uses similarity
 * constraints to make the actual geometry approach these ideal shapes.
 * 
 * @param weight The weight of the constraint (higher values enforce stronger regularity)
 */
void ShapeOpt::addRegularConstraint(double weight)
{
    //regular polygon constraint
    for(int id = 0 ;id < mesh_.faces_.size(); id++)
    {
        int num_verIDs = mesh_.faces_[id].size();

        // Calculate face center and normal
        Eigen::Vector3d center = getCentroid(id);
        Vector3d normal = getNormal(id);
        
        // Create a local coordinate system on the face
        // Find X and Y axes perpendicular to the normal
        Vector3d x_axis = Vector3d(1, 0, 0).cross(normal);
        if(x_axis.norm() < 1e-4) x_axis = Vector3d(0, 1, 0).cross(normal);
        x_axis.normalize();
        Vector3d y_axis = normal.cross(x_axis);
        y_axis.normalize();
        
        // Create a regular polygon template in the face's plane
        Eigen::MatrixXd shape(3, num_verIDs);
        vector<int> id_vector;
        for(int jd = 0; jd < num_verIDs; jd ++)
        {
            // Create evenly spaced points in a circle to form a regular polygon
            float angle = 3.1415926 * 2 / num_verIDs * jd;
            Vector3d pt = x_axis * std::cos(angle) + y_axis * std::sin(angle) + center;
            shape.col(jd) << pt[0], pt[1], pt[2];
            id_vector.push_back(mesh_.faces_[id][jd]);
        }
        
        // Create a similarity constraint that attempts to transform the actual face
        // into the regular polygon template (allowing scaling, rotation, and translation)
        shared_ptr<ShapeOp::SimilarityConstraint> constraint = make_shared<ShapeOp::SimilarityConstraint>(id_vector, weight, V_ ,true, true, true);
        std::vector<ShapeOp::Matrix3X> shapes;
        shapes.push_back(shape);
        constraint->setShapes(shapes);
        solver.addConstraint(constraint);
    }
}

/**
 * Adds constraints to fix the boundary vertices in place.
 * 
 * This is important for architectural applications where the mesh needs to
 * connect to existing structures or maintain specific boundary conditions.
 * The function detects boundary loops and adds closeness constraints to
 * keep those vertices fixed during optimization.
 * 
 * @param weight The weight of the constraint (higher values enforce stronger boundary preservation)
 */
void ShapeOpt::addBoundaryConstraint(double weight)
{
    // Detect boundary loops using libigl
    std::vector<std::vector<int> > L;
    igl::boundary_loop(mesh_.F_, L);

    // For each boundary loop
    for(int id = 0; id < L.size(); id++)
    {
        // For each vertex in the boundary
        for(int jd = 0; jd < L[id].size(); jd++)
        {
            std::vector<int> id_vector;
            int verID = L[id][jd]; 
            id_vector.push_back(verID);
            
            // Get the current position of the boundary vertex
            ShapeOp::Vector3 pt = V_.col(verID);
            
            // Create a closeness constraint to keep this vertex fixed
            shared_ptr<ShapeOp::ClosenessConstraint> constraint = make_shared<ShapeOp::ClosenessConstraint>(id_vector, weight, V_);
            constraint->setPosition(pt);
            solver.addConstraint(constraint);
        }
    }
}

/**
 * Runs the ShapeOp solver to optimize the mesh according to all added constraints.
 * 
 * This performs the actual optimization process, iteratively moving vertices to
 * satisfy the constraints as best as possible (weighted by their importance).
 * 
 * @param iterations Number of solver iterations to perform
 */
void ShapeOpt::runShapeOp(int iterations)
{
    // Initialize the solver
    solver.initialize();
    // Run the solver for the specified number of iterations
    solver.solve(iterations);
    // Get the optimized point positions
    ShapeOp::Matrix3X P = solver.getPoints();
    // Update our internal representation
    V_ = P;
    //std::cout << V_ << std::endl;
}

/**
 * Returns the optimized mesh.
 * 
 * Converts the internal optimized geometry back to the iglMesh format
 * for further processing or visualization.
 * 
 * @return The optimized mesh in iglMesh format
 */
iglMesh ShapeOpt::getMesh()
{
    // Transfer optimized vertex positions back to the mesh structure
    mesh_.V_ = V_.transpose();
    return mesh_;
}

/**
 * Calculates the centroid (geometric center) of a face.
 * 
 * Used as a helper function for various geometric operations.
 * 
 * @param fid Face index
 * @return 3D position of the face centroid
 */
Vector3d ShapeOpt::getCentroid(int fid)
{
    Vector3d centroid(0, 0, 0);
    // Sum all vertex positions
    for(int jd = 0; jd < mesh_.faces_[fid].size(); jd ++)
    {
        int vID = mesh_.faces_[fid][jd];
        Vector3d pt = V_.col(vID);
        centroid += pt;
    }
    // Divide by number of vertices to get average position
    centroid /= mesh_.faces_[fid].size();
    return centroid;
}

/**
 * Calculates the normal vector of a face using Principal Component Analysis (PCA).
 * 
 * This method is robust for non-planar faces and uses eigenanalysis of the
 * vertex covariance matrix to find the direction of least variance, which
 * corresponds to the face normal.
 * 
 * @param fid Face index
 * @return Unit normal vector of the face
 */
Vector3d ShapeOpt::getNormal(int fid) {
    // Collect all vertices of the face
    vector<vector<double>> vers;
    for(int jd = 0; jd < mesh_.faces_[fid].size(); jd ++)
    {
        int vID = mesh_.faces_[fid][jd];
        vector<double> pt;
        pt.push_back(V_(0, vID));pt.push_back(V_(1, vID));pt.push_back(V_(2, vID));
        vers.push_back(pt);
    }

    // Check if we have enough vertices for a proper face
    if(vers.size() < 3)
    {
        return Vector3d(0, 0, 0);
    }

    // Get face centroid
    Vector3d centroid = getCentroid(fid);

    // Compute the covariance matrix of vertex positions relative to centroid
    // This is used for Principal Component Analysis (PCA)
    double xx, xy, xz, yy, yz, zz;
    xx = xy = xz = yy = yz = zz = 0;
    for(vector<double> pt: vers)
    {
        Vector3d r = Vector3d(pt[0], pt[1], pt[2]) - centroid;
        xx += r[0] * r[0];
        xy += r[0] * r[1];
        xz += r[0] * r[2];
        yy += r[1] * r[1];
        yz += r[1] * r[2];
        zz += r[2] * r[2];
    }

    // Calculate determinants to find the eigenvector with smallest eigenvalue
    // This corresponds to the normal direction (direction of least variance)
    double det_x = yy*zz - yz*yz;
    double det_y = xx*zz - xz*xz;
    double det_z = xx*yy - xy*xy;

    double maxDet = std::max(det_x, std::max(det_y, det_z));
    if(maxDet <= 0){
        return Vector3d(0, 0, 0);
    }

    // Select the eigenvector corresponding to the smallest eigenvalue
    // This will be perpendicular to the best-fit plane through the vertices
    Vector3d normal(0, 0, 0);
    if(maxDet == det_x)
    {
        normal = Vector3d(det_x, xz*yz - xy*zz, xy*yz - xz*yy);
    }
    else if(maxDet == det_y)
    {
        normal = Vector3d(xz*yz - xy*zz, det_y, xy*xz - yz*xx);
    }
    else {
        normal = Vector3d(xy*yz - xz*yy, xy*xz - yz*xx, det_z);
    };
    // Normalize to unit length
    normal /= normal.norm();
    return normal;
}
