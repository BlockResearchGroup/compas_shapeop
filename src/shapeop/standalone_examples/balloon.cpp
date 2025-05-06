#include "pch.h"
#include "NormalForce.h"
#include "Solver.h"
#include "Constraint.h"
#include <fstream>
#include <iostream>
#include <vector>

int main() {
    const int rows = 10, cols = 10;
    double gridSize = 2.0;
    ShapeOp::Matrix3X points(3, rows * cols);

    // Helper to get the index of a grid point
    auto index = [cols](int x, int y) { return y * cols + x; };

    // Initialize grid points
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            int idx = index(x, y);
            points(0, idx) = x * gridSize / (cols - 1);
            points(1, idx) = y * gridSize / (rows - 1);
            points(2, idx) = 0.0; // Flat grid
        }
    }

    // Initialize solver
    ShapeOp::Solver solver;
    solver.setPoints(points);

    // Add constraints to pin the corners
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

    // Add edge constraints to maintain grid structure
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

    // Define faces of the grid
    std::vector<std::vector<int>> faces;
    for (int y = 0; y < rows - 1; ++y) {
        for (int x = 0; x < cols - 1; ++x) {
            int v1 = index(x, y);
            int v2 = index(x + 1, y);
            int v3 = index(x + 1, y + 1);
            int v4 = index(x, y + 1);

            // Add two triangular faces for each quad
            faces.push_back({v1, v2, v3});
            faces.push_back({v1, v3, v4});
        }
    }

    // Add normal force
    double normalForceMagnitude = 0.5;
    auto normalForce = std::make_shared<ShapeOp::NormalForce>(faces, normalForceMagnitude);
    solver.addForces(normalForce);

    // Initialize and solve
    solver.initialize(false);
    for (int i = 0; i < 1000; ++i) {
        solver.solve(1);
    }

    // Get final points
    const ShapeOp::Matrix3X &finalPoints = solver.getPoints();

    // Debug: Print final points to console
    for (int i = 0; i < finalPoints.cols(); ++i) {
        std::cout << "Point " << i << ": "
                  << finalPoints(0, i) << ", "
                  << finalPoints(1, i) << ", "
                  << finalPoints(2, i) << std::endl;
    }

    // Write mesh to OBJ file
    std::ofstream objFile("balloon_with_normal_force.obj");
    if (objFile.is_open()) {
        // Write vertices
        for (int i = 0; i < finalPoints.cols(); ++i) {
            objFile << "v " << finalPoints(0, i) << " " << finalPoints(1, i) << " " << finalPoints(2, i) << "\n";
        }

        // Write faces
        for (int y = 0; y < rows - 1; ++y) {
            for (int x = 0; x < cols - 1; ++x) {
                int v1 = index(x, y) + 1;
                int v2 = index(x + 1, y) + 1;
                int v3 = index(x + 1, y + 1) + 1;
                int v4 = index(x, y + 1) + 1;
                objFile << "f " << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";
            }
        }

        objFile.close();
        std::cout << "Mesh written to balloon_with_normal_force.obj" << std::endl;
    } else {
        std::cerr << "Failed to open OBJ file for writing." << std::endl;
    }

    return 0;
}