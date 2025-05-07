from .constraint import Constraint


class MeshVerticesConstraint(Constraint):
    """Add closeness constraints to multiple vertices.

    Parameters
    ----------
    vertices : list
        List of vertex indices to constrain.
    weight : float, optional
        Weight of the constraints. Default is 1e5.
    """

    def __init__(self, shape, vertices, weight=1e5):
        self.shape = shape
        self.points = shape.to_vertices_and_faces()[0]

        self.particle_indices = []

        self.weight = weight
        self.vertices = vertices

        self.solver = None

    def compute(self):
        """Call the C++ Method to calculate the constraint."""
        for local_vertex in self.vertices:
            vertex = self.particle_indices[local_vertex]
            self.solver.add_closeness_constraint(vertex, self.weight)

    def update(self):
        """This method is called after solve method is called to reconstruct the original shape."""
        return
