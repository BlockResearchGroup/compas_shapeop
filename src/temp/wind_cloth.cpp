#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>
#include "pch.h"

// Simple cloth simulation using ShapeOp
// Demonstrates cloth hanging from two corners

int main() {
    // Parameters for the cloth
    const int rows = 20, cols = 20;
    double gridSize = 1.0;
    
    // Initialize points matrix
    ShapeOp::Matrix3X points(3, rows * cols);
    
    // Convenient index to access grid vertex
    auto index = [cols](int x, int y) { return y * cols + x; };
    
    // Create flat grid
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            double posX = x * gridSize;
            double posY = 0.0;
            double posZ = y * gridSize;
            
            int idx = index(x, y);
            points(0, idx) = posX;
            points(1, idx) = posY;
            points(2, idx) = posZ;
        }
    }
    
    // Create solver
    ShapeOp::Solver solver;
    solver.setPoints(points);
    
    // Fix only two corners (diagonal corners)
    std::vector<std::pair<int, int>> fixedCorners = {
        {0, 0},           // top-left
        {cols-1, rows-1}  // bottom-right
    };
    
    double cornerWeight = 1e5; // Strong constraint
    
    for (const auto& corner : fixedCorners) {
        int x = corner.first;
        int y = corner.second;
        std::vector<int> idx = {index(x, y)};
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(
            idx, cornerWeight, solver.getPoints());
        solver.addConstraint(constraint);
        std::cout << "Fixed corner at (" << x << ", " << y << ")" << std::endl;
    }
    
    // Add edge constraints
    double edgeWeight = 10.0;
    
    // Horizontal edges
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols - 1; ++x) {
            std::vector<int> edgeIndices = {index(x, y), index(x+1, y)};
            auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(
                edgeIndices, edgeWeight, solver.getPoints(), 0.8, 1.2);
            solver.addConstraint(constraint);
        }
    }
    
    // Vertical edges
    for (int y = 0; y < rows - 1; ++y) {
        for (int x = 0; x < cols; ++x) {
            std::vector<int> edgeIndices = {index(x, y), index(x, y+1)};
            auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(
                edgeIndices, edgeWeight, solver.getPoints(), 0.8, 1.2);
            solver.addConstraint(constraint);
        }
    }
    
    // Add gravity force
    ShapeOp::Vector3 gravity(0.0, -0.1, 0.0); // Y is down
    auto gravityForce = std::make_shared<ShapeOp::GravityForce>(gravity);
    solver.addForces(gravityForce);
    
    // Initialize and solve
    solver.initialize(true);
    
    std::cout << "Solving cloth hanging from two corners..." << std::endl;
    
    // Solve for 100 iterations
    const int iterations = 100;
    for (int i = 0; i < iterations; ++i) {
        solver.solve(1);
        
        if (i % 20 == 0) {
            std::cout << "Iteration " << i << std::endl;
        }
    }
    
    std::cout << "Simulation complete." << std::endl;
    
    // Write the final result to an OBJ file
    const ShapeOp::Matrix3X& finalPoints = solver.getPoints();
    std::string filename = "hanging_cloth.obj";
    std::ofstream objFile(filename);
    
    if (!objFile.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return 1;
    }
    
    // Write header
    objFile << "# Hanging cloth mesh" << std::endl;
    objFile << "# Vertices: " << rows * cols << std::endl;
    objFile << "# Faces: " << (rows-1) * (cols-1) << std::endl;
    
    // Write vertices
    for (int i = 0; i < rows * cols; ++i) {
        objFile << "v " 
                << finalPoints(0, i) << " " 
                << finalPoints(1, i) << " " 
                << finalPoints(2, i) << std::endl;
    }
    
    // Write quad faces
    for (int y = 0; y < rows - 1; ++y) {
        for (int x = 0; x < cols - 1; ++x) {
            int i00 = index(x, y) + 1;       // +1 because OBJ indices start at 1
            int i10 = index(x+1, y) + 1;
            int i11 = index(x+1, y+1) + 1;
            int i01 = index(x, y+1) + 1;
            
            objFile << "f " << i00 << " " << i10 << " " << i11 << " " << i01 << std::endl;
        }
    }
    
    objFile.close();
    std::cout << "Wrote result to " << filename << std::endl;
    
    return 0;
}
