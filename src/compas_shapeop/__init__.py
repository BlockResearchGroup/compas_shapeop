from .solver import Solver
import os

# Import optimized module
try:
    from compas_shapeop import _pyshapeop
except ImportError:
    _pyshapeop = None

# Import the optimized ShapeOpSolver
try:
    from .shapeoplibrary import ShapeOpSolver
except ImportError:
    ShapeOpSolver = None

__version__ = "0.1.0"


HERE = os.path.dirname(__file__)

HOME = os.path.abspath(os.path.join(HERE, "../../"))
DATA = os.path.abspath(os.path.join(HOME, "data"))
DOCS = os.path.abspath(os.path.join(HOME, "docs"))
TEMP = os.path.abspath(os.path.join(HOME, "temp"))


__all__ = ["HOME", "DATA", "DOCS", "TEMP", "Solver", "_pyshapeop", "ShapeOpSolver"]
