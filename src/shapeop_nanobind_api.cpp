#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>

#include "Solver.h"

namespace nb = nanobind;

// Helper function to convert a flat Python list to an Eigen matrix
ShapeOp::Matrix3X convertFlatListToMatrix(nb::list points_list, int n_points) {
    // Create an Eigen matrix to hold the points
    ShapeOp::Matrix3X matrix(3, n_points);
    
    // Fill the matrix with data from the flat list (column-major)
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < 3; j++) {
            matrix(j, i) = nb::cast<double>(points_list[i * 3 + j]);
        }
    }
    
    return matrix;
}

// Method to set points from a flat list
void setPointsFromFlat(ShapeOp::Solver &solver, nb::list points_list, int n_points) {
    // Check input parameters
    if (points_list.size() != n_points * 3) {
        throw std::runtime_error("Flat points list should have n_points*3 elements");
    }
    
    // Convert and set the points
    ShapeOp::Matrix3X matrix = convertFlatListToMatrix(points_list, n_points);
    solver.setPoints(matrix);
}

// Method to get points as a flat list
std::tuple<nb::list, int> getPointsAsFlat(ShapeOp::Solver &solver) {
    // Get the points matrix from the solver
    const ShapeOp::Matrix3X& matrix = solver.getPoints();
    int n_points = static_cast<int>(matrix.cols());
    
    // Create a flat Python list
    nb::list result;
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < 3; j++) {
            result.append(matrix(j, i));
        }
    }
    
    return std::make_tuple(result, n_points);
}

// Define the Python module
NB_MODULE(_shapeop_direct, m) {
    // Module documentation
    m.doc() = "ShapeOp solver bindings using nanobind with efficient data transfer";

    // Expose the Solver class
    nb::class_<ShapeOp::Solver>(m, "Solver")
        // Constructor
        .def(nb::init<>(), "Create a new ShapeOp solver")
        
        //====================================
        // Core solver methods
        //====================================
        .def("initialize", &ShapeOp::Solver::initialize, "Initialize the solver", 
             nb::arg("dynamic") = false, nb::arg("masses") = 1.0, 
             nb::arg("damping") = 1.0, nb::arg("timestep") = 1.0)
        
        .def("solve", &ShapeOp::Solver::solve, "Run the solver for a given number of iterations", 
             nb::arg("iterations") = 1)
        
        //====================================
        // Point management methods
        //====================================
        .def("set_points_from_flat", [](ShapeOp::Solver &self, nb::list points_list, int n_points) {
            setPointsFromFlat(self, points_list, n_points);
        }, "Set points from a flat list of coordinates")
        
        .def("get_points_as_flat", [](ShapeOp::Solver &self) {
            return getPointsAsFlat(self);
        }, "Get points as a flat list and point count");
}