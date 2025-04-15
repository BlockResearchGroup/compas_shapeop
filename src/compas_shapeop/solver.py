from compas_shapeop import _shapeop as so


class Solver:
    """A wrapper class providing an object-oriented API for the ShapeOp solver.

    Parameters
    ----------
    max_iterations : int, optional
        Maximum number of iterations to run. Default is 1000.

    Attributes
    ----------
    iterations : int
        The current number of iterations performed.
    max_iterations : int
        The maximum number of iterations the solver will run.
    """

    def __init__(self, max_iterations: int = 1000):
        """Create a new ShapeOp solver.

        Parameters
        ----------
        max_iterations : int, optional
            Maximum number of iterations to run. Default is 1000.
        """
        self._solver = so.Solver()
        self._deleted = False
        self._max_iterations = max_iterations
        self._current_iteration = 0

    def init(self) -> bool:
        """Initialize the solver.

        Returns
        -------
        bool
            True if initialization was successful.
        """
        self._current_iteration = 0
        return so.init_solver(self._solver)

    def solve(self, iterations: int = None) -> bool:
        """Run the solver for a given number of iterations.

        Parameters
        ----------
        iterations : int, optional
            Number of iterations to run. If None, will run until max_iterations is reached.

        Returns
        -------
        bool
            True if solving was successful.
        """
        if iterations is None:
            # Run until max_iterations is reached
            remaining = self._max_iterations - self._current_iteration
            if remaining <= 0:
                return 0

            iterations_to_run = remaining
        else:
            iterations_to_run = iterations

        result = so.solve(self._solver, iterations_to_run)
        self._current_iteration += iterations_to_run

        # Check if we've reached max iterations and auto-clean if needed
        if self._current_iteration >= self._max_iterations:
            # Auto-delete the solver when max iterations are reached
            self._deleted = True

        return result

    def reset_iterations(self) -> None:
        """Reset the iteration counter to zero."""
        self._current_iteration = 0

    @property
    def iterations(self) -> int:
        """Get the current number of iterations performed.

        Returns
        -------
        int
            The number of solver iterations performed.
        """
        return self._current_iteration

    @property
    def max_iterations(self) -> int:
        """Get the maximum number of iterations.

        Returns
        -------
        int
            The maximum number of iterations.
        """
        return self._max_iterations

    @max_iterations.setter
    def max_iterations(self, value: int) -> None:
        """Set the maximum number of iterations.

        Parameters
        ----------
        value : int
            The new maximum number of iterations.
        """
        self._max_iterations = value

    def set_points(self, points: list[list[float]]) -> None:
        """Set the points for the solver.

        Parameters
        ----------
        points : list[list[float]]
            List of point coordinates as [x, y, z] lists.
        """
        so.set_points(self._solver, points)

    def get_points(self) -> list[list[float]]:
        """Get the points from the solver.

        Returns
        -------
        list[list[float]]
            List of points with their current coordinates.
        """
        return so.get_points(self._solver)

    def add_constraint(self, constraint_type: str, indices: list[int], weight: float) -> int:
        """Add a constraint to the solver.

        Parameters
        ----------
        constraint_type : str
            Type of constraint to add (e.g. "Closeness", "Plane", "Circle", "EdgeStrain").
        indices : list[int]
            List of vertex indices this constraint applies to.
        weight : float
            Weight (importance) of this constraint.

        Returns
        -------
        int
            Identifier for the created constraint.
        """
        return so.add_constraint(self._solver, constraint_type, indices, weight)

    def get_constraints(self) -> list[int]:
        """Get list of constraint ids.

        Returns
        -------
        list[int]
            List of all constraint identifiers.
        """
        return so.get_constraints(self._solver)

    def set_closeness_position(self, constraint_id: int, position: list[float]) -> bool:
        """Set the target position for a closeness constraint.

        Parameters
        ----------
        constraint_id : int
            Identifier of the closeness constraint.
        position : list[float]
            Target position as [x, y, z] coordinates.

        Returns
        -------
        bool
            True if position was updated successfully.
        """
        return so.set_closeness_position(self._solver, constraint_id, position)

    def add_closeness_constraint_with_position(self, vertex_id: int, weight: float, position: list[float]) -> int:
        """Add a closeness constraint with a specified target position.

        Parameters
        ----------
        vertex_id : int
            Index of the vertex to constrain.
        weight : float
            Weight (importance) of this constraint.
        position : list[float]
            Target position as [x, y, z] coordinates.

        Returns
        -------
        int
            Identifier for the created constraint.
        """
        return so.add_closeness_constraint_with_position(self._solver, vertex_id, weight, position)

    def add_edge_strain_constraint(self, indices: list[int], weight: float, min_range: float, max_range: float) -> int:
        """Add an edge strain constraint with custom range parameters.

        Parameters
        ----------
        indices : list[int]
            List of 2 vertex indices defining the edge.
        weight : float
            Weight (importance) of this constraint.
        min_range : float
            Minimum allowed strain factor.
        max_range : float
            Maximum allowed strain factor.

        Returns
        -------
        int
            Identifier for the created constraint.
        """
        return so.add_edge_strain_constraint(self._solver, indices, weight, min_range, max_range)

    def add_shrinking_edge_constraint(self, indices: list[int], weight: float, shrink_factor: float) -> int:
        """Add a shrinking edge constraint using a shrink factor.

        Parameters
        ----------
        indices : list[int]
            List of 2 vertex indices defining the edge.
        weight : float
            Weight (importance) of this constraint.
        shrink_factor : float
            Factor by which to shrink the edge.

        Returns
        -------
        int
            Identifier for the created constraint.
        """
        return so.add_shrinking_edge_constraint(self._solver, indices, weight, shrink_factor)

    def add_vertex_force(self, force_x: float, force_y: float, force_z: float, vertex_id: int) -> bool:
        """Add a vertex force to the solver.

        Parameters
        ----------
        force_x : float
            X component of the force vector.
        force_y : float
            Y component of the force vector.
        force_z : float
            Z component of the force vector.
        vertex_id : int
            Index of the vertex to apply the force to.

        Returns
        -------
        bool
            True if the force was added successfully.
        """
        return so.add_vertex_force(self._solver, force_x, force_y, force_z, vertex_id)

    def add_gravity_force(self, gravity_x: float, gravity_y: float, gravity_z: float) -> bool:
        """Add a gravity-like force to all points in the solver.

        Parameters
        ----------
        gravity_x : float
            X component of the gravity vector.
        gravity_y : float
            Y component of the gravity vector.
        gravity_z : float
            Z component of the gravity vector.

        Returns
        -------
        bool
            True if the force was added successfully.
        """
        return so.add_gravity_force(self._solver, gravity_x, gravity_y, gravity_z)

    def add_normal_force_with_faces(self, faces_flat: list[int], face_sizes: list[int], magnitude: float) -> bool:
        """Add a normal force (inflation) using custom face topology.

        Parameters
        ----------
        faces_flat : list[int]
            Flattened list of face vertices.
        face_sizes : list[int]
            List containing the number of vertices in each face.
        magnitude : float
            Magnitude of the normal force.

        Returns
        -------
        bool
            True if the force was added successfully.
        """
        return so.add_normal_force_with_faces(self._solver, faces_flat, face_sizes, magnitude)

    def __enter__(self) -> 'Solver':
        """Support for 'with' statement (context manager).

        Returns
        -------
        Solver
            The solver instance.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Clean up resources when exiting a 'with' block."""
        self._deleted = True

    def delete(self) -> None:
        """Mark the solver as deleted to prevent further use."""
        self._deleted = True
