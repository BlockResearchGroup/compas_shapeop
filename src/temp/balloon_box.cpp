#include "pch.h"
#include "NormalForce.h"
#include "Solver.h"
#include "Constraint.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <set>

void readOBJ(const std::string &filename, ShapeOp::Matrix3X &points, std::vector<std::vector<int>> &faces) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open OBJ file: " + filename);
    }

    std::vector<ShapeOp::Vector3> vertices;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            // Read vertex
            double x, y, z;
            iss >> x >> y >> z;
            vertices.emplace_back(x, y, z);
        } else if (prefix == "f") {
            // Read face
            std::vector<int> face;
            std::string vertex;
            while (iss >> vertex) {
                size_t pos = vertex.find('/');
                int index = std::stoi(vertex.substr(0, pos)) - 1; // OBJ indices are 1-based
                face.push_back(index);
            }
            faces.push_back(face);
        }
    }

    // Convert vertices to ShapeOp::Matrix3X
    points.resize(3, vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        points.col(i) = vertices[i];
    }
}

int main() {
    // Read vertices and faces from OBJ file
    ShapeOp::Matrix3X points;
    std::vector<std::vector<int>> faces;
    const std::string objFilePath = "data/m0.obj";

    try {
        readOBJ(objFilePath, points, faces);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // Initialize solver
    ShapeOp::Solver solver;
    solver.setPoints(points);

    // Add closeness constraints to all vertices with a small stiffness value
    double smallStiffness = 0.001; // Very small constraint value
    for (int i = 0; i < points.cols(); ++i) {
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(std::vector<int>{i}, smallStiffness, solver.getPoints());
        solver.addConstraint(constraint);
    }
    auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(std::vector<int>{0}, 1000, solver.getPoints());
    solver.addConstraint(constraint);
    constraint = std::make_shared<ShapeOp::ClosenessConstraint>(std::vector<int>{100}, 1000, solver.getPoints());
    solver.addConstraint(constraint);
    constraint = std::make_shared<ShapeOp::ClosenessConstraint>(std::vector<int>{350}, 1000, solver.getPoints());
    solver.addConstraint(constraint);

    // Add edge constraints (strings)
    std::set<std::pair<int, int>> uniqueEdges; // To avoid duplicate edges
    for (const auto &face : faces) {
        for (size_t i = 0; i < face.size(); ++i) {
            int v1 = face[i];
            int v2 = face[(i + 1) % face.size()]; // Wrap around to form a closed loop
            if (v1 > v2) std::swap(v1, v2);      // Ensure consistent ordering
            if (uniqueEdges.insert({v1, v2}).second) {
                // Add edge constraint
                std::vector<int> edge = {v1, v2};
                auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(edge, 0.1, solver.getPoints());
                solver.addConstraint(constraint);
            }
        }
    }

    // Add normal force
    double normalForceMagnitude = 0.1;
    auto normalForce = std::make_shared<ShapeOp::NormalForce>(faces, normalForceMagnitude);
    solver.addForces(normalForce);

    // Initialize and solve
    solver.initialize(false);
    for (int i = 0; i < 10; ++i) {
        solver.solve(1);
    }

    // Get final points
    const ShapeOp::Matrix3X &finalPoints = solver.getPoints();

    // Write mesh to OBJ file
    std::ofstream objFile("balloon_box_with_normal_force.obj");
    if (objFile.is_open()) {
        // Write vertices
        for (int i = 0; i < finalPoints.cols(); ++i) {
            objFile << "v " << finalPoints(0, i) << " " << finalPoints(1, i) << " " << finalPoints(2, i) << "\n";
        }

        // Write faces
        for (const auto &face : faces) {
            objFile << "f";
            for (int index : face) {
                objFile << " " << (index + 1); // OBJ indices are 1-based
            }
            objFile << "\n";
        }

        objFile.close();
        std::cout << "Mesh written to balloon_box_with_normal_force.obj" << std::endl;
    } else {
        std::cerr << "Failed to open OBJ file for writing." << std::endl;
    }

    return 0;
}