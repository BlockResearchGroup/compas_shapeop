from compas_shapeop import _shapeoplibrary
from compas.datastructures import Mesh

class ShapeOpSolver:
    """Optimized ShapeOp solver with direct zero-copy access to solver memory.
    
    This solver uses nanobind's Eigen integration to provide direct zero-copy 
    access to the ShapeOp solver's internal memory. This ensures maximum
    performance for dynamic simulations.
    
    The implementation maintains a direct numpy view into the C++ solver's
    Eigen matrix memory. When the solver modifies point positions, the NumPy
    array is automatically updated without any copying or data conversion.
    
    Available constraints:
    - Closeness
    - EdgeStrain
    
    Available forces:
    - VertexForce
    """
    
    def __init__(self):
        """Initialize a new optimized ShapeOpSolver.
        
        Creates a new DynamicSolver instance that handles direct
        memory sharing between C++ and Python.
        """

        self._solver = _shapeoplibrary.DynamicSolver()
        self.points = None # Direct reference to ShapeOp's points matrix
        
    
    def set_points(self, points):
        """Set the vertex positions in the solver.
        
        This method initializes the solver's internal memory with the
        provided points. After setting the points, it establishes a
        direct zero-copy view to the solver's memory.
        
        Parameters
        ----------
        points : array-like
            Array of 3D points in shape (n, 3).
        """
        self._solver.set_points(points)
    
    def get_points(self):
        """Get the current vertex positions from the solver.
        
        Returns a direct view of the solver's point data. This is a zero-copy
        operation that provides direct access to the memory used by the C++ solver.
        
        Returns
        -------
        numpy.ndarray
            Direct view of the solver's point data in shape (n, 3).
        """
        return self.points
    
    def add_closeness_constraint(self, indices, weight=1.0):
        """Add a closeness constraint to the solver.
        
        A closeness constraint tries to keep vertices close to their
        original positions. This directly adds the constraint to the C++ solver.
        
        Parameters
        ----------
        indices : list
            List of vertex indices to constrain.
        weight : float, optional
            Weight of the constraint. Higher values make the constraint stronger.
            
        Returns
        -------
        int
            ID of the added constraint.
        """
        return self._solver.add_closeness_constraint(indices, weight)
    
    def add_closeness_constraint_with_position(self, indices, weight, position):
        """Add a closeness constraint with a specified target position.
        
        A closeness constraint that keeps vertices close to a specified target position
        rather than their original positions. This is useful for fixing points in space
        or forcing points to move to specific locations.
        
        Parameters
        ----------
        indices : list
            List of vertex indices to constrain.
        weight : float
            Weight of the constraint. Higher values make the constraint stronger.
        position : list
            Target position [x, y, z] for the constraint.
            
        Returns
        -------
        int
            ID of the added constraint.
        """
        return self._solver.add_closeness_constraint_with_position(indices, weight, position)
    
    def add_edge_strain_constraint(self, indices, weight=1.0, min_range=0.9, max_range=1.1):
        """Add an edge strain constraint to the solver.
        
        An edge strain constraint tries to keep the distance between two vertices
        within a specified range relative to the original distance.
        This directly adds the constraint to the C++ solver.
        
        Parameters
        ----------
        indices : list
            List of vertex indices to constrain. Must contain exactly 2 indices.
        weight : float, optional
            Weight of the constraint. Higher values make the constraint stronger.
        min_range : float, optional
            Minimum allowed relative length (default=0.9)
        max_range : float, optional
            Maximum allowed relative length (default=1.1)
            
        Returns
        -------
        int
            ID of the added constraint.
        """
        return self._solver.add_edge_strain_constraint(indices, weight, min_range, max_range)
    
    def add_shrinking_edge_constraint(self, indices, weight=1.0, shrink_factor=0.25):
        """Add a shrinking edge constraint to the solver.
        
        A shrinking edge constraint tries to shrink the edge length by a specified factor.
        This is particularly useful for cable nets and other structures that need to
        maintain tension. The constraint creates a min/max range of ±5% around the target length.
        
        Parameters
        ----------
        indices : list
            List of vertex indices to constrain. Must contain exactly 2 indices.
        weight : float, optional
            Weight of the constraint. Higher values make the constraint stronger.
        shrink_factor : float, optional
            Target shrinking factor (default=0.25). The target length will be
            (1.0 - shrink_factor) times the original length.
            
        Returns
        -------
        int
            ID of the added constraint.
        """
        return self._solver.add_shrinking_edge_constraint(indices, weight, shrink_factor)
    
    def add_circle_constraint(self, indices, weight=1.0):
        """Add a circle constraint to the solver.
        
        A circle constraint tries to make a set of vertices lie on a circle.
        This is useful for circularizing polygonal faces in a mesh.
        
        Parameters
        ----------
        indices : list
            List of vertex indices to constrain. Must contain at least 3 vertices.
        weight : float, optional
            Weight of the constraint. Higher values make the constraint stronger.
            
        Returns
        -------
        int
            ID of the added constraint.
        """
        return self._solver.add_circle_constraint(indices, weight)
    
    def add_plane_constraint(self, indices, weight=1.0):
        """Add a plane constraint to the solver.
        
        A plane constraint tries to make a set of vertices lie on a plane.
        This is useful for planarizing polygonal faces in a mesh.
        
        Parameters
        ----------
        indices : list
            List of vertex indices to constrain. Must contain at least 3 vertices.
        weight : float, optional
            Weight of the constraint. Higher values make the constraint stronger.
            
        Returns
        -------
        int
            ID of the added constraint.
        """
        return self._solver.add_plane_constraint(indices, weight)
    
    def add_vertex_force(self, force_x, force_y, force_z, vertex_id):
        """Add a force to a specific vertex.
        
        A vertex force applies a constant force vector to a specific vertex.
        This directly adds the force to the C++ solver.
        
        Parameters
        ----------
        force_x : float
            X component of the force.
        force_y : float
            Y component of the force.
        force_z : float
            Z component of the force.
        vertex_id : int
            Index of the vertex to apply the force to.
            
        Returns
        -------
        bool
            True if the force was added successfully.
        """
        return self._solver.add_vertex_force(force_x, force_y, force_z, vertex_id)

    def _setup_direct_view(self):
        """Set up a direct view to the solver's internal memory for zero-copy operations.
        
        Gets a direct reference to the solver's points matrix and transposes it
        to match the standard (N×3) format expected by Python code. This doesn't copy 
        data, just creates a view with different strides.
        """
        # Get direct reference to the solver's points matrix and transpose it
        self.points = self._solver.get_points_ref().T
        
        # Verify the view has the expected shape
        if self.points.shape[1] != 3:
            raise ValueError(f"Expected 3 columns (x,y,z), got {self.points.shape[1]}")

    def init(self):
        """Initialize the solver.
        
        This method must be called after adding all constraints and before
        calling solve(). It builds the solver's internal structures based on
        the added constraints.
        
        After initialization, the direct view to the solver's memory is
        established to ensure zero-copy access to the point positions.
        
        Returns
        -------
        ndarray
            Direct reference to the solver's points matrix, in row-major (N×3) format.
            This is the same as accessing the `.points` attribute after initialization.
        """
        # Initialize the solver
        self._solver.initialize()
        
        # Set up the direct view after initialization when memory is fully prepared
        self._setup_direct_view()
        
        # Return the points reference for convenience
        return self.points
    
    def solve(self, iterations=10):
        """Solve the constraint problem.
        
        Runs the simulation for the specified number of iterations.
        The solver will try to satisfy all constraints while respecting
        the applied forces. All point positions are updated directly in
        shared memory.
        
        Parameters
        ----------
        iterations : int, optional
            Number of iterations to run.
        """
        self._solver.solve(iterations)
    
    # ==========================================================================
    # COMPAS Mesh Integration Methods
    # ==========================================================================
    
    @classmethod
    def from_mesh(cls, mesh):
        """Create a solver initialized with vertices from a COMPAS mesh.
        
        This is a convenience method to create a solver directly from a mesh.
        
        Parameters
        ----------
        mesh : :class:`compas.datastructures.Mesh`
            A COMPAS mesh.
            
        Returns
        -------
        :class:`ShapeOpSolver`
            A new solver instance initialized with the mesh vertices.
        """
        solver = cls()
        solver.set_points(mesh.to_vertices_and_faces()[0])
        return solver
    
    def add_mesh_edge_strain_constraint(self, mesh, weight=1.0, min_range=0.8, max_range=1.2):
        """Add edge strain constraints to all edges of a COMPAS mesh.
        
        This is a convenience method that adds edge strain constraints to all
        edges in the provided mesh.
        
        Parameters
        ----------
        mesh : :class:`compas.datastructures.Mesh`
            A COMPAS mesh.
        weight : float, optional
            The weight of the constraint.
        min_range : float, optional
            The minimum allowed edge length ratio (0.8 means 20% compression allowed).
        max_range : float, optional
            The maximum allowed edge length ratio (1.2 means 20% stretching allowed).
        """
        for u, v in mesh.edges():
            self.add_edge_strain_constraint([u, v], weight, min_range, max_range)
    
    def add_mesh_vertex_force(self, mesh, force_x, force_y, force_z, exclude_vertices=None):
        """Add a force to all vertices of a COMPAS mesh.
        
        This is a convenience method that adds the same force to all vertices
        in the provided mesh, optionally excluding specified vertices.
        
        Parameters
        ----------
        mesh : :class:`compas.datastructures.Mesh`
            A COMPAS mesh.
        force_x : float
            Force component in x direction.
        force_y : float
            Force component in y direction.
        force_z : float
            Force component in z direction.
        exclude_vertices : list, optional
            List of vertex indices to exclude from force application.
        """
        exclude_vertices = exclude_vertices or []
        for vertex in mesh.vertices():
            if vertex not in exclude_vertices:
                self.add_vertex_force(force_x, force_y, force_z, vertex)
                
    def add_mesh_closeness_constraint(self, mesh, vertices, weight=1e5):
        """Add closeness constraints to specific vertices of a COMPAS mesh.
        
        This is a convenience method that adds closeness constraints to the
        specified vertices in the provided mesh.
        
        Parameters
        ----------
        mesh : :class:`compas.datastructures.Mesh`
            A COMPAS mesh.
        vertices : list
            List of vertex indices to constrain.
        weight : float, optional
            The weight of the constraint.
        """
        for vertex in vertices:
            self.add_closeness_constraint([vertex], weight)
