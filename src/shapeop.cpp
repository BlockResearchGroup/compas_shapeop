#include "compas.h"

//╔═══════════════════════════════════════════════════════════════════════════╗
//║                        DYNAMIC SOLVER CLASS                               ║
//╚═══════════════════════════════════════════════════════════════════════════╝

/**
 * @brief Wrapper class for ShapeOp::Solver with dynamic memory management
 * @details Provides a C++ interface between ShapeOp library and Python with
 *          automatic memory management using std::unique_ptr
 */
class DynamicSolver {
private:
    std::unique_ptr<ShapeOp::Solver> solver;  //!< Managed ShapeOp solver instance
    
    /**
     * @brief Helper method to check if solver is valid
     * @return true if solver pointer is not null
     */
    bool is_valid() const {
        return solver != nullptr;
    }
    
public:
    /**
     * @brief Constructor for the DynamicSolver
     * @throws std::runtime_error if solver creation fails
     */
    DynamicSolver() : solver(std::make_unique<ShapeOp::Solver>()) {
        if (!solver) {
            throw std::runtime_error("Failed to create ShapeOp solver");
        }
    }
    
    /**
     * @brief Destructor for the DynamicSolver
     * @details The unique_ptr will automatically release the solver
     */
    ~DynamicSolver() {
        // The unique_ptr will automatically release the solver
    }
    
    //┌───────────────────────────────────────────────────────────────────────┐
    //│                    SOLVER CORE FUNCTIONALITY                          │
    //└───────────────────────────────────────────────────────────────────────┘
    
    /**
     * @brief Get a reference to the solver's internal points matrix
     * @details Provides zero-copy access to the solver's internal points
     * @return Reference to the solver's points matrix
     * @throws std::runtime_error if solver is invalid
     */
    Eigen::Ref<Eigen::MatrixXd> get_points_ref() {
        if (!is_valid()) {
            throw std::runtime_error("Invalid solver");
        }
        // Direct access to the solver's points matrix
        return solver->points_;
    }
    
    /**
     * @brief Set points from a list of lists
     * @param points_list List of points where each point is a list of 3 coordinates
     * @throws std::runtime_error if points list is empty or if each point does not have 3 coordinates
     */
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
    
    /**
     * @brief Initialize the solver
     * @return true if initialization is successful
     */
    bool initialize() {
        if (!is_valid()) {
            return false;
        }
        
        return solver->initialize();
    }
    
    /**
     * @brief Solve for a number of iterations
     * @param iterations Number of iterations to solve for
     * @return true if solve is successful
     */
    bool solve(int iterations) {
        if (!is_valid()) {
            return false;
        }
        
        return solver->solve(iterations);
    }
    
    //┌───────────────────────────────────────────────────────────────────────┐
    //│                    GEOMETRIC CONSTRAINTS                              │
    //└───────────────────────────────────────────────────────────────────────┘
    
    /**
     * @brief Add closeness constraint
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @return true if constraint is added successfully
     */
    bool add_closeness_constraint(nb::list indices, double weight) {
        if (!is_valid()) {
            return false;
        }
        
        std::vector<int> indices_vec;
        for (nb::handle h : indices) {
            indices_vec.push_back(nb::cast<int>(h));
        }
        
        auto constraint = std::make_shared<ShapeOp::ClosenessConstraint>(
            indices_vec, weight, solver->getPoints());
        
        return solver->addConstraint(constraint) > 0;
    }
    
    /**
     * @brief Add closeness constraint with target position
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @param position Target position
     * @return true if constraint is added successfully
     */
    bool add_closeness_constraint_with_position(nb::list indices, double weight, nb::list position) {
        if (!is_valid()) {
            return false;
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
        return solver->addConstraint(constraint) > 0;
    }
    
    /**
     * @brief Add edge strain constraint
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @param min_range Minimum range of the constraint
     * @param max_range Maximum range of the constraint
     * @return true if constraint is added successfully
     */
    bool add_edge_strain_constraint(nb::list indices, double weight, double min_range, double max_range) {
        if (!is_valid()) {
            return false;
        }
        
        std::vector<int> indices_vec;
        for (nb::handle h : indices) {
            indices_vec.push_back(nb::cast<int>(h));
        }
        
        auto constraint = std::make_shared<ShapeOp::EdgeStrainConstraint>(
            indices_vec, weight, solver->getPoints(), min_range, max_range);
        
        return solver->addConstraint(constraint) > 0;
    }
    
    /**
     * @brief Add shrinking edge constraint (specifically for cable nets)
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @param shrink_factor Shrink factor of the constraint
     * @return true if constraint is added successfully
     */
    bool add_shrinking_edge_constraint(nb::list indices, double weight, double shrink_factor) {
        if (!is_valid()) {
            return false;
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
        
        return solver->addConstraint(constraint) > 0;
    }
    
    /**
     * @brief Add circle constraint (for face circularization)
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @return true if constraint is added successfully
     */
    bool add_circle_constraint(nb::list indices, double weight) {
        if (!is_valid()) {
            return false;
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
        
        return solver->addConstraint(constraint) > 0;
    }
    
    /**
     * @brief Add plane constraint (for face planarization)
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @return true if constraint is added successfully
     */
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

    /**
     * @brief Add similarity constraint (for regular polygon formation)
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @param allow_scaling Allow scaling of the constraint
     * @param allow_rotation Allow rotation of the constraint
     * @param allow_translation Allow translation of the constraint
     * @return true if constraint is added successfully
     */
    bool add_similarity_constraint(nb::list indices, double weight, 
                                  bool allow_scaling, bool allow_rotation, bool allow_translation) {
        if (!is_valid()) {
            return false;
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
        
        // Create the constraint
        auto constraint = std::make_shared<ShapeOp::SimilarityConstraint>(
            ids, weight, solver->getPoints(), 
            allow_scaling, allow_rotation, allow_translation);
        
        return solver->addConstraint(constraint) > 0;
    }
    
    /**
     * @brief Set shape for similarity constraint - allows specifying the target shape
     * @param constraint_id ID of the constraint
     * @param points List of points defining the target shape
     * @return true if shape is set successfully
     */
    bool set_similarity_constraint_shape(int constraint_id, nb::list points) {
        if (!is_valid()) {
            return false;
        }
        
        // Validate constraint ID - ShapeOp uses 1-based constraint IDs
        if (constraint_id <= 0) {
            throw std::runtime_error("Invalid constraint ID - must be > 0");
        }
        
        // Get the constraint directly by ID
        auto constraint = solver->getConstraint(constraint_id);
        if (!constraint) {
            throw std::runtime_error("Constraint not found with the given ID");
        }
        
        // Check if it's a similarity constraint
        auto similarity_constraint = std::dynamic_pointer_cast<ShapeOp::SimilarityConstraint>(constraint);
        if (!similarity_constraint) {
            throw std::runtime_error("Constraint is not a SimilarityConstraint");
        }
        
        // Convert Python list to ShapeOp::Matrix3X
        int num_points = len(points);
        if (num_points <= 0) {
            throw std::runtime_error("Empty points list");
        }
        
        ShapeOp::Matrix3X shape(3, num_points);
        for (int i = 0; i < num_points; i++) {
            nb::list point = nb::cast<nb::list>(points[i]);
            if (len(point) != 3) {
                throw std::runtime_error("Each point must have 3 coordinates (x,y,z)");
            }
            
            shape.col(i)[0] = nb::cast<double>(point[0]);
            shape.col(i)[1] = nb::cast<double>(point[1]);
            shape.col(i)[2] = nb::cast<double>(point[2]);
        }
        
        // Set the shape
        std::vector<ShapeOp::Matrix3X> shapes;
        shapes.push_back(shape);
        similarity_constraint->setShapes(shapes);
        
        return true;
    }
    
    /**
     * @brief Add regular polygon constraint for a face (high-level convenience function)
     * @param indices List of vertex indices
     * @param weight Weight of the constraint
     * @return true if constraint is added successfully
     */
    bool add_regular_polygon_constraint(nb::list indices, double weight) {
        if (!is_valid()) {
            return false;
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
        
        // Need at least 3 vertices to define a polygon
        if (ids.size() < 3) {
            throw std::runtime_error("Regular polygon constraint requires at least 3 vertices");
        }
        
        // Calculate face centroid
        ShapeOp::Vector3 centroid = ShapeOp::Vector3::Zero();
        for (const auto& id : ids) {
            centroid += solver->getPoints().col(id);
        }
        centroid /= ids.size();
        
        // Calculate face normal using first three points
        ShapeOp::Vector3 p0 = solver->getPoints().col(ids[0]);
        ShapeOp::Vector3 p1 = solver->getPoints().col(ids[1]);
        ShapeOp::Vector3 p2 = solver->getPoints().col(ids[2]);
        ShapeOp::Vector3 normal = (p1 - p0).cross(p2 - p0).normalized();
        
        // Create a local coordinate system on the face
        ShapeOp::Vector3 x_axis = ShapeOp::Vector3(1, 0, 0);
        if (std::abs(x_axis.dot(normal)) > 0.9) {
            x_axis = ShapeOp::Vector3(0, 1, 0);
        }
        x_axis = (x_axis - normal * x_axis.dot(normal)).normalized();
        ShapeOp::Vector3 y_axis = normal.cross(x_axis).normalized();
        
        // Create a regular polygon template
        ShapeOp::Matrix3X shape(3, ids.size());
        for (size_t i = 0; i < ids.size(); i++) {
            double angle = 2.0 * M_PI * i / ids.size();
            ShapeOp::Vector3 pt = centroid + 
                                  x_axis * std::cos(angle) + 
                                  y_axis * std::sin(angle);
            shape.col(i) = pt;
        }
        
        // Create a similarity constraint
        auto constraint = std::make_shared<ShapeOp::SimilarityConstraint>(
            ids, weight, solver->getPoints(), true, true, true);
        
        // Set the regular polygon shape
        std::vector<ShapeOp::Matrix3X> shapes;
        shapes.push_back(shape);
        constraint->setShapes(shapes);
        
        return solver->addConstraint(constraint) > 0;
    }

    //┌───────────────────────────────────────────────────────────────────────┐
    //│                    FORCE CONSTRAINTS                                  │
    //└───────────────────────────────────────────────────────────────────────┘

    /**
     * @brief Add normal force with faces
     * @param faces_flat_array Flat array of face indices
     * @param face_sizes_array Array of face sizes
     * @param magnitude Magnitude of the force
     * @return true if force is added successfully
     */
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

    /**
     * @brief Add vertex force (for individual vertices)
     * @param force_x X-component of the force
     * @param force_y Y-component of the force
     * @param force_z Z-component of the force
     * @param vertex_id ID of the vertex
     * @return true if force is added successfully
     */
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
};

//╔═══════════════════════════════════════════════════════════════════════════╗
//║                        PYTHON MODULE DEFINITION                           ║
//╚═══════════════════════════════════════════════════════════════════════════╝

// Define the module with a more explicit name to avoid conflicts
NB_MODULE(_shapeop, m) {
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
        .def("add_similarity_constraint", &DynamicSolver::add_similarity_constraint)
        .def("set_similarity_constraint_shape", &DynamicSolver::set_similarity_constraint_shape)
        .def("add_regular_polygon_constraint", &DynamicSolver::add_regular_polygon_constraint)
        .def("add_normal_force_with_faces", &DynamicSolver::add_normal_force_with_faces)
        .def("add_vertex_force", &DynamicSolver::add_vertex_force)
        .def("initialize", &DynamicSolver::initialize)
        .def("solve", &DynamicSolver::solve);
}
