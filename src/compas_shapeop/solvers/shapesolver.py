from compas_shapeop import Solver


class ShapeSolver(Solver):
    """Create a solver from a COMPAS geometry types.

    This is a convenience method that allows you to initialize a solver
    directly from a COMPAS mesh.

    Attributes
    ----------
    shapes : :class:`compas.datastructures.Mesh`
        A COMPAS mesh.
    points : numpy.ndarray
        Direct reference to the solver's points matrix in shape (n, 3).
    """

    def __init__(self):
        super().__init__()
        self.constraints = []
        self.is_initialized = False

    def solve(self, iterations=10, t=1e-3):
        """Solve the constraint problem.

        Parameters
        ----------
        iterations : int, optional
            Number of iterations to run. Default is 10.
        t : float, optional
            Tolerance for point comparison. Default is 1e-3.
        """

        # Intialize the solver, this is needed after you add points, contraints and forces
        if not self.is_initialized:
            points = []

            def is_same_point(p1, p2, tolerance):
                return abs(p1[0] - p2[0]) < tolerance and abs(p1[1] - p2[1]) < tolerance and abs(p1[2] - p2[2]) < tolerance

            for constraint in self.constraints:
                for point in constraint.points:
                    for i, existing in enumerate(points):
                        if is_same_point(point, existing, t):
                            constraint.particle_indices.append(i)
                            break
                    else:
                        constraint.particle_indices.append(len(points))
                        points.append(point)

            self.points = points

            for constraint in self.constraints:
                constraint.compute()

            # Initialize the C++ constrains and retrieve the numpy point array
            self.init()
            self.is_initialized = True

        # Run the solver.
        self._solver.solve(iterations)

        # Update the shapes constraints.
        for constraint in self.constraints:
            constraint.update()

    def add(self, constraint):
        constraint.solver = self
        self.constraints.append(constraint)
        return constraint
