from .shapeop import Solver
from .solvers import MeshSolver, ShapeSolver
import os

__version__ = "0.1.0"

HERE = os.path.dirname(__file__)
HOME = os.path.abspath(os.path.join(HERE, "../../"))
DATA = os.path.join(HOME, "data")
DOCS = os.path.abspath(os.path.join(HOME, "docs"))
TEMP = os.path.abspath(os.path.join(HOME, "temp"))

__all__ = ["HOME", "DATA", "DOCS", "TEMP", "Solver", "MeshSolver", "ShapeSolver"]
