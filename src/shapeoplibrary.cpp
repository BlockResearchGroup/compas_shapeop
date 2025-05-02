#include "Solver.h"
#include "Constraint.h"
#include "Force.h"
#include "custom_constraints/normalforce.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/pair.h>
#include <nanobind/ndarray.h>
#include <nanobind/eigen/dense.h>  // Include this for Eigen support

#include <iostream>
#include <memory>
#include <vector>
#include <stdexcept>

namespace nb = nanobind;

// Wrapper class for ShapeOp::Solver with a more unique name to avoid conflicts
class DynamicSolver {
private:
    std::unique_ptr<ShapeOp::Solver> solver;
    
    // Helper method to check if solver is valid
    bool is_valid() const {
        return solver != nullptr;
    }
    
public:
    DynamicSolver() : solver(std::make_unique<ShapeOp::Solver>()) {
        if (!solver) {
            throw std::runtime_error("Failed to create ShapeOp solver");
        }
    }
    
    ~DynamicSolver() {
        // The unique_ptr will automatically release the solver
    }
    
    // Direct access to ShapeOp's internal points matrix with zero-copy
    Eigen::Ref<Eigen::MatrixXd> get_points_ref() {
        if (!is_valid()) {
            throw std::runtime_error("Invalid solver");
        }
        // Direct access to the solver's points matrix
        return solver->points_;
    }
    
    // Set points from a list of lists
    void set_points(nb::list points_list) {
        if (!is_valid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Get number of points
        size_t num_points = len(points_list);
        if (num_points == 0) {
            throw std::runtime_error("Empty points list");
        }
        
        // Create a 3×n matrix for ShapeOp
        ShapeOp::Matrix3X matrix_points(3, num_points);
        
        // Fill the matrix from the Python list
        for (size_t i = 0; i < num_points; ++i) {
            nb::list point = nb::cast<nb::list>(points_list[i]);
            if (len(point) != 3) {
                throw std::runtime_error("Each point must have 3 coordinates (x,y,z)");
            }
            
            matrix_points(0, i) = nb::cast<double>(point[0]);
            matrix_points(1, i) = nb::cast<double>(point[1]);
            matrix_points(2, i) = nb::cast<double>(point[2]);
        }
        
        // Set points in the solver
        solver->setPoints(matrix_points);
    }
    
    // Add closeness constraint
    int add_closeness_constraint(nb::list indices, double weight) {
        if (!is_valid()) {
            return -1;
        }
        
        std::vector<int> indices_vec;
        for (nb::handle h : indices) {
            indices_vec.push_back(nb::cast<int>(h));
        }
        
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(
            indices_vec, weight, solver->getPoints());
        
        return solver->addConstraint(constraint);
    }
    
    // Add closeness constraint with target position
    int add_closeness_constraint_with_position(nb::list indices, double weight, nb::list position) {
        if (!is_valid()) {
            return -1;
        }
        
        // Validate position
        if (len(position) != 3) {
            throw std::runtime_error("Position must have 3 coordinates (x,y,z)");
        }
        
        // Convert indices to vector
        std::vector<int> indices_vec;
        for (nb::handle h : indices) {
            indices_vec.push_back(nb::cast<int>(h));
        }
        
        // Create the constraint
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(
            indices_vec, weight, solver->getPoints());
        
        // Extract and set the target position
        ShapeOp::Vector3 pos;
        pos(0) = nb::cast<double>(position[0]);
        pos(1) = nb::cast<double>(position[1]);
        pos(2) = nb::cast<double>(position[2]);
        
        // Set the target position
        constraint->setPosition(pos);
        
        // Add constraint to solver
        return solver->addConstraint(constraint);
    }
    
    // Add edge strain constraint
    int add_edge_strain_constraint(nb::list indices, double weight, double min_range, double max_range) {
        if (!is_valid()) {
            return -1;
        }
        
        std::vector<int> indices_vec;
        for (nb::handle h : indices) {
            indices_vec.push_back(nb::cast<int>(h));
        }
        
        auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(
            indices_vec, weight, solver->getPoints(), min_range, max_range);
        
        return solver->addConstraint(constraint);
    }
    
    // Add shrinking edge constraint (specifically for cable nets)
    int add_shrinking_edge_constraint(nb::list indices, double weight, double shrink_factor) {
        if (!is_valid()) {
            return -1;
        }
        
        std::vector<int> indices_vec;
        for (nb::handle h : indices) {
            indices_vec.push_back(nb::cast<int>(h));
        }
        
        if (indices_vec.size() != 2) {
            throw std::runtime_error("EdgeStrain constraint requires exactly 2 indices");
        }
        
        // Calculate the min/max range based on the shrink factor
        double min_range = shrink_factor - 0.05;  // 5% below target
        double max_range = shrink_factor + 0.05;  // 5% above target
        
        auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(
            indices_vec, weight, solver->getPoints(), min_range, max_range);
        
        return solver->addConstraint(constraint);
    }
    
    // Add circle constraint (for face circularization)
    int add_circle_constraint(nb::list indices, double weight) {
        if (!is_valid()) {
            return -1;
        }
        
        std::vector<int> indices_vec;
        for (nb::handle h : indices) {
            indices_vec.push_back(nb::cast<int>(h));
        }
        
        // Circle constraint requires at least 3 vertices
        if (indices_vec.size() < 3) {
            throw std::runtime_error("Circle constraint requires at least 3 vertices");
        }
        
        auto constraint = std::make_shared<ShapeOp::CircleConstraint>(
            indices_vec, weight, solver->getPoints());
        
        return solver->addConstraint(constraint);
    }
    
    // Add plane constraint (for face planarization)
    bool add_plane_constraint(nb::list indices, double weight) {
        if (!solver) {
            throw std::runtime_error("Solver not initialized");
        }

        // Convert Python list to std::vector
        std::vector<int> ids;
        for (auto item : indices) {
            try {
                ids.push_back(nb::cast<int>(item));
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to convert list item to integer");
            }
        }

        // Validate that we have at least 3 vertices (needed for a plane)
        if (ids.size() < 3) {
            throw std::runtime_error("Plane constraint requires at least 3 vertices");
        }

        // Create the constraint
        auto constraint = std::make_shared<ShapeOp::PlaneConstraint>(ids, weight, solver->getPoints());
        auto constraint_id = solver->addConstraint(constraint);
        return constraint_id > 0;
    }

    bool add_normal_force_with_faces(nb::ndarray<int> faces_flat_array, 
                                    nb::ndarray<int> face_sizes_array, 
                                    double magnitude) {
        if (!solver) {
            throw std::runtime_error("Solver not initialized");
        }

        // Get data pointers and dimensions
        const int* faces_flat_ptr = faces_flat_array.data();
        const int* face_sizes_ptr = face_sizes_array.data();
        
        size_t faces_flat_size = faces_flat_array.size();
        size_t face_count = face_sizes_array.size();
        
        // Convert the flat array representation to vector of vectors for faces
        std::vector<std::vector<int>> faces;
        size_t idx = 0;
        
        for (size_t i = 0; i < face_count; ++i) {
            int face_size = face_sizes_ptr[i];
            std::vector<int> face;
            
            for (int j = 0; j < face_size; ++j) {
                if (idx < faces_flat_size) {
                    face.push_back(faces_flat_ptr[idx++]);
                } else {
                    throw std::runtime_error("Face index out of bounds");
                }
            }
            
            faces.push_back(face);
        }
        
        // Create the normal force
        auto force = std::make_shared<ShapeOp::NormalForce>(faces, magnitude);
        solver->addForces({force});
        return true;
    }

    // Add vertex force (for individual vertices)
    bool add_vertex_force(double force_x, double force_y, double force_z, int vertex_id) {
        if (!is_valid()) {
            return false;
        }
        
        if (vertex_id < 0 || vertex_id >= solver->getPoints().cols()) {
            throw std::runtime_error("Vertex index out of bounds");
        }
        
        // Create a force vector
        ShapeOp::Vector3 force(force_x, force_y, force_z);
        
        // Create a vertex force and add it directly to the solver using addForces
        std::shared_ptr<ShapeOp::Force> vertex_force = std::make_shared<ShapeOp::VertexForce>(force, vertex_id);
        solver->addForces(vertex_force);
        
        return true;
    }
    
    // Initialize the solver
    bool initialize() {
        if (!is_valid()) {
            return false;
        }
        
        return solver->initialize();
    }
    
    // Solve for a number of iterations
    bool solve(int iterations) {
        if (!is_valid()) {
            return false;
        }
        
        return solver->solve(iterations);
    }
};

// Define the module with a more explicit name to avoid conflicts
NB_MODULE(_shapeoplibrary, m) {
    // Give a clear docstring about this module
    m.doc() = "ShapeOp dynamic solver binding";
    
    // Define the solver class with a unique name
    nb::class_<DynamicSolver>(m, "DynamicSolver")
        .def(nb::init<>())
        .def("set_points", &DynamicSolver::set_points)
        .def("get_points_ref", &DynamicSolver::get_points_ref)
        .def("add_closeness_constraint", &DynamicSolver::add_closeness_constraint)
        .def("add_closeness_constraint_with_position", &DynamicSolver::add_closeness_constraint_with_position)
        .def("add_edge_strain_constraint", &DynamicSolver::add_edge_strain_constraint)
        .def("add_shrinking_edge_constraint", &DynamicSolver::add_shrinking_edge_constraint)
        .def("add_circle_constraint", &DynamicSolver::add_circle_constraint)
        .def("add_plane_constraint", &DynamicSolver::add_plane_constraint)
        .def("add_normal_force_with_faces", &DynamicSolver::add_normal_force_with_faces)
        .def("add_vertex_force", &DynamicSolver::add_vertex_force)
        .def("initialize", &DynamicSolver::initialize)
        .def("solve", &DynamicSolver::solve);
}
