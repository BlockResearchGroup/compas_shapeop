#include "compas.h"

// Simple global counter for the last number of points - this approach works for simple examples
// but would need to be made solver-specific in a production environment
static int global_nb_points = 0;

//------------------------------------------------------------------------------
// C++ Implementation
//------------------------------------------------------------------------------

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

// Initialize the solver
bool init_solver(PyShapeOpSolver& solver) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    return shapeop_init(solver.get());
}

// Run the solver for a given number of iterations
bool solve(PyShapeOpSolver& solver, unsigned int iterations) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    return shapeop_solve(solver.get(), iterations);
}

// Set points for the solver from a list
void set_points(PyShapeOpSolver& solver, nb::list points_list) {
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
}

// Get points from the solver
nb::list get_points(PyShapeOpSolver& solver) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    
    // Get the points directly from the C++ Solver object
    const ShapeOp::Matrix3X& points_matrix = solver.get()->s->getPoints();
    int nb_points = points_matrix.cols();
    
    if (nb_points <= 0) {
        throw std::runtime_error("No points in solver");
    }
    
    // For now, we'll use a standard list implementation that works reliably
    // Note that for high-performance applications, you might want to use the 
    // Python API to convert this list to a NumPy array immediately after returning
    
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
}

// Add a constraint to the solver
int add_constraint(PyShapeOpSolver& solver, std::string constraint_type, std::vector<int> indices, double weight) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    std::vector<int> ids = indices;
    int constraint_id = shapeop_addConstraint(solver.get(), constraint_type.c_str(), ids.data(), (int)ids.size(), weight);
    return constraint_id;
}

// Get points directly using a pre-allocated NumPy array buffer
// This function performs no copying - it directly uses the memory in the provided buffer
bool get_points_buffer(PyShapeOpSolver& solver, double* buffer_ptr, int nb_points) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    
    // Get a reference to the points matrix in the solver
    const ShapeOp::Matrix3X& points_matrix = solver.get()->s->getPoints();
    
    // Verify dimensions match
    if (points_matrix.cols() != nb_points) {
        throw std::runtime_error("Buffer size doesn't match solver points count");
    }
    
    // Create an Eigen Map to the output buffer (which is owned by Python)
    // This creates a view into the memory buffer without copying
    Eigen::Map<ShapeOp::Matrix3X> buffer_map(buffer_ptr, 3, nb_points);
    
    // Copy the data from solver to the buffer
    // This is a direct memory-to-memory copy, very efficient
    buffer_map = points_matrix;
    
    return true;
}

// Set points using a pre-allocated NumPy array buffer
// This function performs no copying - it directly uses the memory in the provided buffer
bool set_points_buffer(PyShapeOpSolver& solver, double* buffer_ptr, int nb_points) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    
    if (nb_points <= 0) {
        throw std::runtime_error("Empty points array provided");
    }
    
    // Create an Eigen Map to the input buffer (which is owned by Python)
    // This creates a view into the memory buffer without copying
    Eigen::Map<ShapeOp::Matrix3X> buffer_map(buffer_ptr, 3, nb_points);
    
    // Set the points in the solver directly from the buffer
    // The solver will copy this data into its internal storage
    solver.get()->s->setPoints(buffer_map);
    
    // Update the global point counter
    global_nb_points = nb_points;
    
    return true;
}

// Get list of constraints
std::vector<int> get_constraints(PyShapeOpSolver& solver) {
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
}

// Set closeness position for a constraint
int set_closeness_position(PyShapeOpSolver& solver, int constraint_id, nb::list position) {
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
    
    // Set the position
    closeness_constraint->setPosition(pos);
    
    return constraint_id;
}

// Add a closeness constraint with a specified target position
int add_closeness_constraint_with_position(PyShapeOpSolver& solver, int vertex_id, double weight, nb::list position) {
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

    // Add to solver
    int constraint_id = solver.get()->s->addConstraint(constraint);
    
    return constraint_id;
}

// Add an edge strain constraint with custom range parameters
int add_edge_strain_constraint(PyShapeOpSolver& solver, std::vector<int> indices, double weight, double min_range, double max_range) {
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
}

// Add a shrinking edge constraint using a shrink factor
int add_shrinking_edge_constraint(PyShapeOpSolver& solver, std::vector<int> indices, double weight, double shrink_factor) {
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
}

// Add a vertex force to the solver
int add_vertex_force(PyShapeOpSolver& solver, double force_x, double force_y, double force_z, int vertex_id) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    ShapeOpScalar force[3] = {force_x, force_y, force_z};
    int force_id = shapeop_addVertexForce(solver.get(), force, vertex_id);
    return force_id;
}

// Add a gravity-like force to all points in the solver
bool add_gravity_force(PyShapeOpSolver& solver, double gx, double gy, double gz) {
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
}

// Add a normal force (inflation) using custom face topology
bool add_normal_force_with_faces(PyShapeOpSolver& solver, const std::vector<int>& faces_flat, const std::vector<int>& face_sizes, double magnitude) {
    if (!solver.isValid()) {
        throw std::runtime_error("Invalid solver");
    }
    
    // Call our C API function
    int face_count = static_cast<int>(face_sizes.size());
    int force_id = shapeop_addNormalForce(
        solver.get(),
        const_cast<int*>(faces_flat.data()),
        const_cast<int*>(face_sizes.data()),
        face_count,
        magnitude
    );
    
    return force_id > 0;
}

//------------------------------------------------------------------------------
// Nanobind Module Definition
//------------------------------------------------------------------------------

NB_MODULE(_shapeop, m) {
    m.doc() = "Python bindings for ShapeOp library";

    // Removed explicit vector bindings since they're causing conflicts with nanobind's
    // automatic STL bindings, and they're not essential for the core functionality
    
    // Properly wrap our C++ class for Python
    nb::class_<PyShapeOpSolver>(m, "Solver")
        .def(nb::init<>())
        .def("is_valid", &PyShapeOpSolver::isValid);

    // Lifecycle and solver operations
    m.def("init_solver", &init_solver, "Initialize the solver", nb::arg("solver"));
    m.def("solve", &solve, "Run the solver for a given number of iterations", nb::arg("solver"), nb::arg("iterations"));

    // Transfer geometry
    m.def("set_points", &set_points, "Set the points for the solver", nb::arg("solver"), nb::arg("points"));
    m.def("get_points", &get_points, "Get the points from the solver", nb::arg("solver"));
    m.def("get_points_buffer", &get_points_buffer, "Get points directly using a pre-allocated NumPy array buffer", nb::arg("solver"), nb::arg("buffer_ptr"), nb::arg("nb_points"));
    m.def("set_points_buffer", &set_points_buffer, "Set points using a pre-allocated NumPy array buffer", nb::arg("solver"), nb::arg("buffer_ptr"), nb::arg("nb_points"));

    // Constraints
    m.def("add_constraint", &add_constraint, "Add a constraint to the solver", 
          nb::arg("solver"), nb::arg("constraint_type"), nb::arg("indices"), nb::arg("weight"));
    m.def("get_constraints", &get_constraints, "Get list of constraint ids", nb::arg("solver"));
    
    // Closeness constraint - set position
    m.def("set_closeness_position", &set_closeness_position, 
          "Set the target position for a closeness constraint", 
          nb::arg("solver"), nb::arg("constraint_id"), nb::arg("position"));
    
    // Add closeness constraint with specified position (all in one operation)
    m.def("add_closeness_constraint_with_position", &add_closeness_constraint_with_position, 
          "Add a closeness constraint with a specified target position", 
          nb::arg("solver"), nb::arg("vertex_id"), nb::arg("weight"), nb::arg("position"));
    
    // EdgeStrain constraint with min/max parameters
    m.def("add_edge_strain_constraint", &add_edge_strain_constraint, 
          "Add an edge strain constraint with custom range parameters", 
          nb::arg("solver"), nb::arg("indices"), nb::arg("weight"), nb::arg("min_range"), nb::arg("max_range"));
    
    // Add specialized constraint for shrinking edges with a shrink factor
    m.def("add_shrinking_edge_constraint", &add_shrinking_edge_constraint, 
          "Add a shrinking edge constraint using a shrink factor (percentage of original length)", 
          nb::arg("solver"), nb::arg("indices"), nb::arg("weight"), nb::arg("shrink_factor"));

    // Forces
    m.def("add_vertex_force", &add_vertex_force, "Add a vertex force to the solver", 
          nb::arg("solver"), nb::arg("force_x"), nb::arg("force_y"), nb::arg("force_z"), nb::arg("vertex_id"));
    
    // Add gravity force using available ShapeOp API functions
    m.def("add_gravity_force", &add_gravity_force, "Add a gravity-like force to all points in the solver", 
          nb::arg("solver"), nb::arg("gravity_x"), nb::arg("gravity_y"), nb::arg("gravity_z"));
    
    // Add NormalForce with explicit face data using a single flat array
    m.def("add_normal_force_with_faces", &add_normal_force_with_faces, 
          "Add a normal force (inflation) using custom face topology. Pass faces as a flat list along with face sizes.", 
          nb::arg("solver"), nb::arg("faces_flat"), nb::arg("face_sizes"), nb::arg("magnitude"));
}

//------------------------------------------------------------------------------
