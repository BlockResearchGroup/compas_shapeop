class Constraint:
    def __init__(self, shape=None):
        self.shape = None
        self.points = []
        self.particle_indices = []
        self.weight = 1e5
        self.solver = None

    def compute(self, target=None, weight=1e5):
        """Call the C++ Method to calculate the constraint."""
        raise NotImplementedError

    def update(self):
        """This method is called after solve method is called to reconstruct the original shape."""
        raise NotImplementedError
