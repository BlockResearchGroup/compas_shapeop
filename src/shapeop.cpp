#include <nanobind/nanobind.h>

#include "compas.h"
#include "Solver.h" // Include the full definition of ShapeOpSolver
#include "Constraint.h"
#include "Force.h"
#include <iostream>
#include <vector>
#include <fstream>
#include "compas_api.h"



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



namespace nb = nanobind;

using namespace nb::literals;

NB_MODULE(_shapeop, m) {
    m.doc() = "This is a \"hello world\" example with nanobind";
    m.def("add", [](int a, int b) { return a + b; }, "a"_a, "b"_a);

    // Solver lifecycle
    m.def("create_solver", &shapeop_create, "Create a new ShapeOp solver");
    m.def("delete_solver", &shapeop_delete, "Delete a ShapeOp solver", "solver"_a);
    m.def("init_solver", &shapeop_init, "Initialize the solver", "solver"_a);
    m.def("solve", &shapeop_solve, "Run the solver for a given number of iterations", "solver"_a, "iterations"_a);

    // Transfer geometry
    m.def("set_points", [](ShapeOpSolver *solver, std::vector<std::vector<double>> points) {
        int nb_points = points.size();
        std::vector<ShapeOpScalar> flat_points(nb_points * 3);
        for (int i = 0; i < nb_points; ++i) {
            flat_points[i * 3 + 0] = points[i][0];
            flat_points[i * 3 + 1] = points[i][1];
            flat_points[i * 3 + 2] = points[i][2];
        }
        shapeop_setPoints(solver, flat_points.data(), nb_points);
    }, "Set the points for the solver", "solver"_a, "points"_a);
    
    m.def("get_points", [](ShapeOpSolver *solver) {
        int nb_points = 0; // You need to know the number of points beforehand
        std::vector<ShapeOpScalar> flat_points(nb_points * 3);
        shapeop_getPoints(solver, flat_points.data(), nb_points);
        std::vector<std::vector<double>> points(nb_points, std::vector<double>(3));
        for (int i = 0; i < nb_points; ++i) {
            points[i][0] = flat_points[i * 3 + 0];
            points[i][1] = flat_points[i * 3 + 1];
            points[i][2] = flat_points[i * 3 + 2];
        }
        return points;
    }, "Get the points from the solver", "solver"_a);

    // Constraints
    m.def("add_constraint", [](ShapeOpSolver *solver, const std::string &constraint_type, std::vector<int> ids, double weight) {
        return shapeop_addConstraint(solver, constraint_type.c_str(), ids.data(), ids.size(), weight);
    }, "Add a constraint to the solver", "solver"_a, "constraint_type"_a, "ids"_a, "weight"_a);
    
    m.def("edit_constraint", [](ShapeOpSolver *solver, const std::string &constraint_type, int constraint_id, std::vector<double> scalars) {
        return shapeop_editConstraint(solver, constraint_type.c_str(), constraint_id, scalars.data(), scalars.size());
    }, "Edit an existing constraint", "solver"_a, "constraint_type"_a, "constraint_id"_a, "scalars"_a);

    // Forces
    m.def("add_gravity_force", [](ShapeOpSolver *solver, std::vector<double> force) {
        if (force.size() != 3) {
            throw std::invalid_argument("Force must have exactly 3 components");
        }
        return shapeop_addGravityForce(solver, force.data());
    }, "Add a gravity force to the solver", "solver"_a, "force"_a);
    
    m.def("add_vertex_force", [](ShapeOpSolver *solver, std::vector<double> force, int id) {
        if (force.size() != 3) {
            throw std::invalid_argument("Force must have exactly 3 components");
        }
        return shapeop_addVertexForce(solver, force.data(), id);
    }, "Add a vertex force to the solver", "solver"_a, "force"_a, "id"_a);
}
