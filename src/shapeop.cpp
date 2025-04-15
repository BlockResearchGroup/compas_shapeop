#include "compas.h"

// Simple global counter for the last number of points - this approach works for simple examples
// but would need to be made solver-specific in a production environment
static int global_nb_points = 0;


// Wrap the ShapeOpSolver with a proper C++ class that handles deletion correctly
class PyShapeOpSolver {
private:
    ShapeOpSolver* solver;
    bool owned;
    
public:
    PyShapeOpSolver() : solver(shapeop_create()), owned(true) {}
    
    ~PyShapeOpSolver() {
        if (owned && solver) {
            shapeop_delete(solver);
            solver = nullptr;
        }
    }
    
    // No copying
    PyShapeOpSolver(const PyShapeOpSolver&) = delete;
    PyShapeOpSolver& operator=(const PyShapeOpSolver&) = delete;
    
    // Move constructor
    PyShapeOpSolver(PyShapeOpSolver&& other) noexcept : solver(other.solver), owned(other.owned) {
        other.solver = nullptr;
        other.owned = false;
    }
    
    // Move assignment
    PyShapeOpSolver& operator=(PyShapeOpSolver&& other) noexcept {
        if (this != &other) {
            if (owned && solver) {
                shapeop_delete(solver);
            }
            solver = other.solver;
            owned = other.owned;
            other.solver = nullptr;
            other.owned = false;
        }
        return *this;
    }
    
    // Access the underlying solver
    ShapeOpSolver* get() const { return solver; }
    
    // Check if solver is valid
    bool isValid() const { return solver != nullptr; }
};

NB_MODULE(_shapeop, m) {
    m.doc() = "Python bindings for ShapeOp library";

    // Vector binding
    nb::bind_vector<std::vector<double>>(m, "VectorDouble");
    nb::bind_vector<std::vector<int>>(m, "VectorInt");
    nb::bind_vector<std::vector<std::tuple<int, int>>>(m, "VectorTupleIntInt");
    nb::bind_vector<std::vector<std::vector<std::tuple<int, float, float, float>>>>(m, "VectorVectorTupleIntFloatFloatFloat");

    // Properly wrap our C++ class for Python
    nb::class_<PyShapeOpSolver>(m, "Solver")
        .def(nb::init<>())
        .def("is_valid", &PyShapeOpSolver::isValid);

    // Lifecycle and solver operations
    m.def("init_solver", [](PyShapeOpSolver& solver) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        return shapeop_init(solver.get());
    }, "Initialize the solver", "solver"_a);
    
    m.def("solve", [](PyShapeOpSolver& solver, unsigned int iterations) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        return shapeop_solve(solver.get(), iterations);
    }, "Run the solver for a given number of iterations", "solver"_a, "iterations"_a);

    // Transfer geometry
    m.def("set_points", [](PyShapeOpSolver& solver, nb::list points_list) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Get number of points
        int nb_points = points_list.size();
        if (nb_points == 0) {
            throw std::runtime_error("Empty points list provided");
        }
        
        // Allocate buffer for the flattened points
        std::vector<ShapeOpScalar> flat_points(nb_points * 3);
        
        // Safely extract points
        for (int i = 0; i < nb_points; ++i) {
            // Get the point (should be a list or array-like)
            auto point = points_list[i]; 
            
            // Direct indexing approach 
            try {
                flat_points[i * 3 + 0] = nb::cast<double>(point[0]);
                flat_points[i * 3 + 1] = nb::cast<double>(point[1]);
                flat_points[i * 3 + 2] = nb::cast<double>(point[2]);
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error extracting point coordinates: ") + e.what());
            }
        }
        
        // Call the C function with the flattened array
        shapeop_setPoints(solver.get(), flat_points.data(), nb_points);
        
        // Update the global point counter
        global_nb_points = nb_points;
    }, "Set the points for the solver", "solver"_a, "points"_a);

    m.def("get_points", [](PyShapeOpSolver& solver) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Get the points directly from the C++ Solver object
        const ShapeOp::Matrix3X& points_matrix = solver.get()->s->getPoints();
        int nb_points = points_matrix.cols();
        
        if (nb_points <= 0) {
            throw std::runtime_error("No points in solver");
        }
        
        // Create a Python list to return
        nb::list result;
        for (int i = 0; i < nb_points; ++i) {
            nb::list point;
            point.append(points_matrix(0, i));
            point.append(points_matrix(1, i));
            point.append(points_matrix(2, i));
            result.append(point);
        }
        
        return result;
    }, "Get the points from the solver", "solver"_a);

    // Constraints
    m.def("add_constraint", [](PyShapeOpSolver& solver, std::string constraint_type, std::vector<int> indices, double weight) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        std::vector<int> ids = indices;
        int constraint_id = shapeop_addConstraint(solver.get(), constraint_type.c_str(), ids.data(), (int)ids.size(), weight);
        return constraint_id;
    }, "Add a constraint to the solver", "solver"_a, "constraint_type"_a, "indices"_a, "weight"_a);
    
    // Get list of constraints
    m.def("get_constraints", [](PyShapeOpSolver& solver) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Since there's no direct getConstraints() method, we need to 
        // query constraints one by one until we get an invalid one
        std::vector<int> constraint_ids;
        int id = 0;
        
        // Keep trying to get constraints until we find one that's invalid
        while (true) {
            try {
                auto constraint = solver.get()->s->getConstraint(id);
                if (!constraint) {
                    break;
                }
                constraint_ids.push_back(id);
                id++;
            } catch (...) {
                // If we get an exception, we've reached the end of the constraints
                break;
            }
        }
        
        return constraint_ids;
    }, "Get list of constraint ids", "solver"_a);
    
    // Closeness constraint - set position
    m.def("set_closeness_position", [](PyShapeOpSolver& solver, int constraint_id, nb::list position) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Get the constraint
        auto constraint = solver.get()->s->getConstraint(constraint_id);
        auto closeness_constraint = std::dynamic_pointer_cast<ShapeOp::ClosenessConstraint>(constraint);
        
        if (!closeness_constraint) {
            throw std::runtime_error("Constraint is not a ClosenessConstraint");
        }
        
        // Extract the position
        if (position.size() != 3) {
            throw std::runtime_error("Position must have 3 coordinates (x,y,z)");
        }
        
        ShapeOp::Vector3 pos;
        pos(0) = nb::cast<double>(position[0]);
        pos(1) = nb::cast<double>(position[1]);
        pos(2) = nb::cast<double>(position[2]);
        
        std::cout << "Setting closeness constraint position to: [" 
                  << pos(0) << ", " << pos(1) << ", " << pos(2) << "]" << std::endl;
        
        // Set the position
        closeness_constraint->setPosition(pos);
        
        return constraint_id;
    }, "Set the target position for a closeness constraint", "solver"_a, "constraint_id"_a, "position"_a);
    
    // Add closeness constraint with specified position (all in one operation)
    m.def("add_closeness_constraint_with_position", [](PyShapeOpSolver& solver, int vertex_id, double weight, nb::list position) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Extract the position
        if (position.size() != 3) {
            throw std::runtime_error("Position must have 3 coordinates (x,y,z)");
        }
        
        // Create vertex ID array
        std::vector<int> ids = {vertex_id};
        
        // Get current points
        const ShapeOp::Matrix3X& points = solver.get()->s->getPoints();
        
        // Create the closeness constraint directly
        ShapeOp::Vector3 pos;
        pos(0) = nb::cast<double>(position[0]);
        pos(1) = nb::cast<double>(position[1]);
        pos(2) = nb::cast<double>(position[2]);
        
        // Create the constraint directly
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(ids, weight, points);
        
        // Set the target position right away
        constraint->setPosition(pos);
        
        std::cout << "Created closeness constraint with position: [" 
                  << pos(0) << ", " << pos(1) << ", " << pos(2) << "]" << std::endl;
        
        // Add to solver
        int constraint_id = solver.get()->s->addConstraint(constraint);
        
        return constraint_id;
    }, "Add a closeness constraint with a specified target position", 
       "solver"_a, "vertex_id"_a, "weight"_a, "position"_a);
    
    // EdgeStrain constraint with min/max parameters
    m.def("add_edge_strain_constraint", [](PyShapeOpSolver& solver, std::vector<int> indices, double weight, double min_range, double max_range) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        if (indices.size() != 2) {
            throw std::runtime_error("EdgeStrain constraint requires exactly 2 indices");
        }
        
        // Get the current points from the solver to calculate the edge length
        const ShapeOp::Matrix3X& points = solver.get()->s->getPoints();
        
        // Create the constraint directly with the custom ranges
        auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(
            indices, weight, points, min_range, max_range);
            
        // Add the constraint to the solver
        int constraint_id = solver.get()->s->addConstraint(c);
        return constraint_id;
    }, "Add an edge strain constraint with custom range parameters", 
       "solver"_a, "indices"_a, "weight"_a, "min_range"_a, "max_range"_a);
    
    // Add specialized constraint for shrinking edges with a shrink factor
    m.def("add_shrinking_edge_constraint", [](PyShapeOpSolver& solver, std::vector<int> indices, double weight, double shrink_factor) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        if (indices.size() != 2) {
            throw std::runtime_error("EdgeStrain constraint requires exactly 2 indices");
        }
        
        // Get the current points from the solver
        const ShapeOp::Matrix3X& points = solver.get()->s->getPoints();
        
        // Calculate the min/max range based on the shrink factor
        double min_range = shrink_factor - 0.05;  // 5% below target
        double max_range = shrink_factor + 0.05;  // 5% above target
        
        std::cout << "Creating EdgeStrainConstraint with shrink factor: " << shrink_factor 
                  << " (range: " << min_range << " to " << max_range << ")" << std::endl;
        
        // Create the constraint with explicit reference to ShapeOp::EdgeStrainConstraint
        auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(
            indices, weight, points, min_range, max_range);
            
        // Verify that the constraint is created properly
        if (!c) {
            throw std::runtime_error("Failed to create EdgeStrainConstraint");
        }
        
        // Add the constraint to the solver
        int constraint_id = solver.get()->s->addConstraint(c);
        
        // Verify constraint was added
        if (constraint_id < 0) {
            throw std::runtime_error("Failed to add constraint to solver");
        }
        
        return constraint_id;
    }, "Add a shrinking edge constraint using a shrink factor (percentage of original length)", 
       "solver"_a, "indices"_a, "weight"_a, "shrink_factor"_a);

    // Forces
    m.def("add_vertex_force", [](PyShapeOpSolver& solver, double force_x, double force_y, double force_z, int vertex_id) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        ShapeOpScalar force[3] = {force_x, force_y, force_z};
        int force_id = shapeop_addVertexForce(solver.get(), force, vertex_id);
        return force_id;
    }, "Add a vertex force to the solver", "solver"_a, "force_x"_a, "force_y"_a, "force_z"_a, "vertex_id"_a);
    
    // Add gravity force using available ShapeOp API functions
    m.def("add_gravity_force", [](PyShapeOpSolver& solver, double gx, double gy, double gz) {
        if (!solver.isValid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        ShapeOpScalar gravity[3] = {gx, gy, gz};
        
        // Create a gravity force using vertex forces on all points
        // First get the current points
        int nb_points = global_nb_points; // Use the global counter we already have
        
        // Add gravity to each point
        for (int i = 0; i < nb_points; ++i) {
            shapeop_addVertexForce(solver.get(), gravity, i);
        }
        
        return true;
    }, "Add a gravity-like force to all points in the solver", "solver"_a, "gravity_x"_a, "gravity_y"_a, "gravity_z"_a);
}
