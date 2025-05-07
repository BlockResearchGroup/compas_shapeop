from .constraint import Constraint


class MeshEdgesConstraint(Constraint):
    """Add edge length constraints to all edges of a COMPAS mesh.

    Parameters
    ----------
    min_range : float, optional
        Minimum allowed relative length. Default is 0.9.
    max_range : float, optional
        Maximum allowed relative length. Default is 1.1.
    shrink_factor : float, optional
        Target shrinking factor (default=0). The target length will be
        (1.0 - shrink_factor) times the original length.
    exclude_edges : list, optional
        List of edges to exclude from constraints. Default is None.
    weight : float, optional
        Weight of the constraints. Default is 1.0.
    """

    def __init__(self, shape, min_range=0.99, max_range=1.01, shrink_factor=0, exclude_edges=None, weight=1e1):
        self.shape = shape
        self.points = shape.to_vertices_and_faces()[0]

        self.particle_indices = []

        self.weight = weight
        self.min_range = min_range
        self.max_range = max_range
        self.shrink_factor = shrink_factor
        self.exclude_edges = exclude_edges
        self.constraint_ids = []

        self.solver = None

    def compute(self):
        """Call the C++ Method to calculate the constraint."""
        exclude_edges = self.exclude_edges or []

        for local_u, local_v in self.shape.edges():
            u = self.particle_indices[local_u]
            v = self.particle_indices[local_v]
            if (u, v) not in exclude_edges and (v, u) not in exclude_edges:
                if self.shrink_factor > 0:
                    cid = self.solver.add_shrinking_edge_constraint(u, v, self.weight, self.shrink_factor)
                    self.constraint_ids.append(cid)
                else:
                    cid = self.solver.add_edge_strain_constraint(u, v, self.weight, self.min_range, self.max_range)
                    self.constraint_ids.append(cid)

    def update(self):
        """This method is called after solve method is called to reconstruct the original shape."""
        for i, vertex in enumerate(self.shape.vertices()):
            self.shape.vertex_attributes(vertex, "xyz", self.solver.points[self.particle_indices[i]])
