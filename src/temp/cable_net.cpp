#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include "pch.h"

int main() {
    // Create a simple cable net structure
    // A cable net structure is a tensile structure where cables can only resist tension
    // Cables can only resist tension, not compression
    
    // Parameters for the cable net
    double gridSize = 2.0;
    const int rows = 10, cols = 10;
    ShapeOp::Matrix3X points(3, rows * cols);
    
    // Convenient index to access grid vertex
    auto index = [cols](int x, int y) { return y * cols + x; };
    
    // Create flat grid (completely flat initially)
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            double posX = x * gridSize / (cols - 1);
            double posY = y * gridSize / (rows - 1);
            double posZ = 0.0;  // Flat grid to start
            
            int idx = index(x, y);
            points(0, idx) = posX;
            points(1, idx) = posY;
            points(2, idx) = posZ;
        }
    }
    
    // Create solver
    ShapeOp::Solver solver;
    solver.setPoints(points);
    
    // Define which corners to lift and which to keep fixed at ground level
    double liftHeight = 1.0; // Height to lift the corners
    
    // Lift two diagonal corners (top-left and bottom-right)
    std::vector<std::pair<int, int>> liftedCorners = {
        {0, 0},          // top-left corner
        {cols-1, rows-1} // bottom-right corner
    };
    
    // Keep the other two corners at ground level
    std::vector<std::pair<int, int>> groundCorners = {
        {cols-1, 0},     // top-right corner
        {0, rows-1}      // bottom-left corner
    };
    
    // Apply ClosenessConstraints to all corners
    double cornerWeight = 1e5; // Very high weight to fix corners firmly
    
    // Fix the lifted corners at their elevated positions
    for (const auto& corner : liftedCorners) {
        int x = corner.first;
        int y = corner.second;
        std::vector<int> cornerPoint = {index(x, y)};
        
        // Get initial position and modify the Z coordinate
        ShapeOp::Vector3 cornerPos = points.col(index(x, y));
        cornerPos(2) = liftHeight; // Set to the lift height
        
        // Create constraint and set the target position
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(
            cornerPoint, cornerWeight, solver.getPoints());
        constraint->setPosition(cornerPos);
        solver.addConstraint(constraint);
        
        std::cout << "Lifted corner at (" << x << ", " << y << ") to height " << liftHeight << std::endl;
    }
    
    // Fix the ground corners at z=0
    for (const auto& corner : groundCorners) {
        int x = corner.first;
        int y = corner.second;
        std::vector<int> cornerPoint = {index(x, y)};
        
        // Use the original position (which has z=0)
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(
            cornerPoint, cornerWeight, solver.getPoints());
        solver.addConstraint(constraint);
        
        std::cout << "Fixed corner at (" << x << ", " << y << ") at ground level" << std::endl;
    }
    
    // Add standard edge constraints to maintain mesh connectivity
    // But set target lengths to be 50% of original
    double shrinkFactor = 0.5; // This will shrink edges to half their length
    double edgeWeight = 100.0; // Higher weight to enforce the shrinking
    
    // Horizontal cables
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols - 1; ++x) {
            std::vector<int> edgeIndices = {index(x, y), index(x+1, y)};
            
            // Target: shrink to 50% of original length
            // Range: between 45% and 55% of original (allows small variations)
            auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(
                edgeIndices, edgeWeight, solver.getPoints(), 
                shrinkFactor - 0.05, // Min: 45% of original length
                shrinkFactor + 0.05  // Max: 55% of original length
            );
            solver.addConstraint(constraint);
        }
    }
    
    // Vertical cables
    for (int y = 0; y < rows - 1; ++y) {
        for (int x = 0; x < cols; ++x) {
            std::vector<int> edgeIndices = {index(x, y), index(x, y+1)};
            
            // Target: shrink to 50% of original length 
            // Range: between 45% and 55% of original (allows small variations)
            auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(
                edgeIndices, edgeWeight, solver.getPoints(), 
                shrinkFactor - 0.05, // Min: 45% of original length
                shrinkFactor + 0.05  // Max: 55% of original length
            );
            solver.addConstraint(constraint);
        }
    }
    
    // Initialize and solve
    solver.initialize(false);
    
    std::cout << "Optimizing cable net structure with corner lifting... ";
    const int num_iterations = 100;
    for (int i = 0; i < num_iterations; i++) {
        solver.solve(1);
        
        if (i % 10 == 0) {
            std::cout << (i * 100) / num_iterations << "% " << std::flush;
        }
    }
    std::cout << "done." << std::endl;
    
    // Get results
    const ShapeOp::Matrix3X& final_points = solver.getPoints();
    
    // Write the result to an OBJ file for visualization
    std::ofstream objFile("cable_net.obj");
    
    // Write vertices
    for (int i = 0; i < rows * cols; ++i) {
        objFile << "v " 
                << final_points(0, i) << " " 
                << final_points(1, i) << " " 
                << final_points(2, i) << std::endl;
    }
    
    // Write faces as quads (instead of triangles)
    for (int y = 0; y < rows - 1; ++y) {
        for (int x = 0; x < cols - 1; ++x) {
            // Get the four corners of the quad (in counter-clockwise order)
            int i00 = index(x, y) + 1;       // bottom-left  (+1 because OBJ indices start at 1)
            int i10 = index(x+1, y) + 1;     // bottom-right
            int i11 = index(x+1, y+1) + 1;   // top-right
            int i01 = index(x, y+1) + 1;     // top-left
            
            // Write a single quad face with 4 vertices
            objFile << "f " << i00 << " " << i10 << " " << i11 << " " << i01 << std::endl;
        }
    }
    
    objFile.close();
    std::cout << "Wrote cable_net.obj" << std::endl;
    
    return 0;
}
