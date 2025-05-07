from .constraint import Constraint


class MeshForcesConstraint(Constraint):
    """Add a vertex force to a specific point.

    Parameters
    ----------
    index : int
        Index of the vertex to apply force to.
    force_x : float, optional
        X component of the force vector. Default is 0.0.
    force_y : float, optional
        Y component of the force vector. Default is 0.0.
    force_z : float, optional
        Z component of the force vector. Default is 0.0.

    Returns
    -------
    int
        ID of the added force.
    """

    def __init__(self, shape, vertices, force_x=0.0, force_y=0.0, force_z=0.05, weight=1e5):
        self.shape = shape
        self.points = shape.to_vertices_and_faces()[0]

        self.particle_indices = []

        self.vertices = vertices
        self.force_x = force_x
        self.force_y = force_y
        self.force_z = force_z
        self.weight = weight

        self.solver = None

    def compute(self):
        """Call the C++ Method to calculate the constraint."""
        for local_vertex in self.vertices:
            vertex = self.particle_indices[local_vertex]
            self.solver.add_vertex_force(vertex, self.force_x, self.force_y, self.force_z)

    def update(self):
        """This method is called after solve method is called to reconstruct the original shape."""
        return
