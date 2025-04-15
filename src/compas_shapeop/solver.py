from compas_shapeop import _shapeop as so


class Solver:
    """A wrapper class providing an object-oriented API for the ShapeOp solver.

    This class wraps the procedural API of _shapeop and provides a more Pythonic
    interface with methods directly on the solver object.
    """

    def __init__(self):
        """Create a new ShapeOp solver."""
        self._solver = so.Solver()

    def init(self):
        """Initialize the solver."""
        return so.init_solver(self._solver)

    def solve(self, iterations):
        """Run the solver for a given number of iterations."""
        return so.solve(self._solver, iterations)

    def set_points(self, points):
        """Set the points for the solver."""
        so.set_points(self._solver, points)

    def get_points(self):
        """Get the points from the solver."""
        return so.get_points(self._solver)

    def add_constraint(self, constraint_type, indices, weight):
        """Add a constraint to the solver."""
        return so.add_constraint(self._solver, constraint_type, indices, weight)

    def get_constraints(self):
        """Get list of constraint ids."""
        return so.get_constraints(self._solver)

    def set_closeness_position(self, constraint_id, position):
        """Set the target position for a closeness constraint."""
        return so.set_closeness_position(self._solver, constraint_id, position)

    def add_closeness_constraint_with_position(self, vertex_id, weight, position):
        """Add a closeness constraint with a specified target position."""
        return so.add_closeness_constraint_with_position(self._solver, vertex_id, weight, position)

    def add_edge_strain_constraint(self, indices, weight, min_range, max_range):
        """Add an edge strain constraint with custom range parameters."""
        return so.add_edge_strain_constraint(self._solver, indices, weight, min_range, max_range)

    def add_shrinking_edge_constraint(self, indices, weight, shrink_factor):
        """Add a shrinking edge constraint using a shrink factor."""
        return so.add_shrinking_edge_constraint(self._solver, indices, weight, shrink_factor)

    def add_vertex_force(self, force_x, force_y, force_z, vertex_id):
        """Add a vertex force to the solver."""
        return so.add_vertex_force(self._solver, force_x, force_y, force_z, vertex_id)

    def add_gravity_force(self, gravity_x, gravity_y, gravity_z):
        """Add a gravity-like force to all points in the solver."""
        return so.add_gravity_force(self._solver, gravity_x, gravity_y, gravity_z)

    def add_normal_force_with_faces(self, faces_flat, face_sizes, magnitude):
        """Add a normal force (inflation) using custom face topology."""
        return so.add_normal_force_with_faces(self._solver, faces_flat, face_sizes, magnitude)
