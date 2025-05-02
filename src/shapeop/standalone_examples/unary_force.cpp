#include "pch.h"
#include "Solver.h"
#include "Constraint.h"
#include "Force.h"
#include <iostream>
#include <vector>
#include <fstream>

int main() {
    // Grid size - smaller grid for faster execution
    const int rows = 14;
    const int cols = 14;
    const double spacing = 1.0;
    const double gravity_force = 0.001;

    // Helper to get index from grid coordinates
    auto index = [cols](int x, int y) { return y * cols + x; };

    // Create grid points
    ShapeOp::Matrix3X points(3, rows * cols);
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            int i = index(x, y);
            double normX = static_cast<double>(x) / (cols - 1);
            double normY = static_cast<double>(y) / (rows - 1);
            points.col(i) = ShapeOp::Vector3(
                normX * spacing * (cols - 1),
                -normY * spacing * (rows - 1),
                0.0
            );
        }
    }

    // Initialize solver with minimal configuration
    ShapeOp::Solver solver;
    solver.setPoints(points);

    // Only add necessary constraints
    // Pin corners
    {
        std::vector<int> top_left = {index(0, 0)};
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(top_left, 1e5, solver.getPoints());
        solver.addConstraint(constraint);
    }
    {
        std::vector<int> top_right = {index(cols - 1, 0)};
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(top_right, 1e5, solver.getPoints());
        solver.addConstraint(constraint);
    }
    {
        std::vector<int> bottom_left = {index(0, rows - 1)};
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(bottom_left, 1e5, solver.getPoints());
        solver.addConstraint(constraint);
    }
    {
        std::vector<int> bottom_right = {index(cols - 1, rows - 1)};
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(bottom_right, 1e5, solver.getPoints());
        solver.addConstraint(constraint);
    }

    {
        std::vector<int> center = {index(cols / 2, rows / 2)};
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(center, 1e5, solver.getPoints());
        solver.addConstraint(constraint);
    }

    // Only add essential edge constraints
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            int i = index(x, y);
            
            if (x + 1 < cols) {
                std::vector<int> edge = {i, index(x + 1, y)};
                auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(edge, 1.0, solver.getPoints());
                solver.addConstraint(constraint);
            }
            
            if (y + 1 < rows) {
                std::vector<int> edge = {i, index(x, y + 1)};
                auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(edge, 1.0, solver.getPoints());
                solver.addConstraint(constraint);
            }
        }
    }

    // Add gravity
    ShapeOp::Vector3 force(0.0, 0.0, gravity_force);
    auto gravity = std::make_shared<ShapeOp::GravityForce>(force);
    solver.addForces(gravity);

    // Static solver is faster
    solver.initialize(false);

    // Minimal iterations for speed
    std::cout << "Running simulation... ";
    const int num_iterations = 1000;
    for (int i = 0; i < num_iterations; i++) {
        solver.solve(1);
    }
    std::cout << "done." << std::endl;
    
    // Write mesh to OBJ file
    ShapeOp::Matrix3X finalPoints = solver.getPoints();
    std::ofstream meshFile("unary_force.obj");
    if (meshFile.is_open()) {
        // Write vertices
        for (int i = 0; i < finalPoints.cols(); i++) {
            meshFile << "v " 
                     << finalPoints.col(i)[0] << " " 
                     << finalPoints.col(i)[1] << " " 
                     << finalPoints.col(i)[2] << std::endl;
        }
        
        // Write faces
        for (int y = 0; y < rows - 1; ++y) {
            for (int x = 0; x < cols - 1; ++x) {
                int v1 = index(x, y) + 1;
                int v2 = index(x+1, y) + 1;
                int v3 = index(x+1, y+1) + 1;
                int v4 = index(x, y+1) + 1;
                meshFile << "f " << v1 << " " << v2 << " " << v3 << " " << v4 << std::endl;
            }
        }
        
        meshFile.close();
        std::cout << "Mesh written to unary_force.obj" << std::endl;
    }

    return 0;
}
