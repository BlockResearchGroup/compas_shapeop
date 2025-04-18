"""
Direct ShapeOp Solver Interface
===============================

This module provides a direct interface to the ShapeOp solver using nanobind for
efficient data transfer between Python and C++.
"""

import numpy as np

try:
    from compas_shapeop import _shapeop_direct as so_direct
except ImportError:
    raise ImportError(
        "Could not import _shapeop_direct. "
        "Make sure you have built the ShapeOp nanobind API correctly."
    )


class DirectSolver:
    """
    Direct ShapeOp solver with efficient data transfer.
    
    This class provides a direct interface to the ShapeOp solver using nanobind
    for efficient data transfer between Python and C++. It uses NumPy arrays
    for point data, which allows for minimal copying between Python and C++.
    
    Examples
    --------
    >>> import numpy as np
    >>> from compas_shapeop.direct_solver import DirectSolver
    >>> 
    >>> # Create a solver with 4 points
    >>> points = np.array([
    ...     [0.0, 0.0, 0.0],
    ...     [1.0, 0.0, 0.0],
    ...     [1.0, 1.0, 0.0],
    ...     [0.0, 1.0, 0.0]
    ... ])
    >>> 
    >>> solver = DirectSolver()
    >>> solver.set_points(points)
    >>> solver.initialize()
    >>> 
    >>> # Implement constraints here...
    >>> 
    >>> # Solve
    >>> solver.solve(10)
    >>> 
    >>> # Get the updated points
    >>> updated_points = solver.get_points()
    """
    
    def __init__(self):
        """Initialize a new DirectSolver."""
        self._solver = so_direct.Solver()
        self._current_iteration = 0
        self.max_iterations = 100
    
    def set_points(self, points: np.ndarray) -> None:
        """
        Set the points in the solver.
        
        This method provides efficient data transfer between NumPy
        and the ShapeOp solver with minimal copying.
        
        Parameters
        ----------
        points : numpy.ndarray
            Array of shape (n_points, 3) containing the point positions.
        """
        if not isinstance(points, np.ndarray):
            points = np.array(points, dtype=np.float64)
        
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError(f"Points must have shape (n_points, 3), got {points.shape}")
        
        # Ensure we have contiguous, correctly-typed data
        points = np.ascontiguousarray(points, dtype=np.float64)
        
        # Flatten the array to a list for the C++ function
        flat_list = points.flatten().tolist()
        n_points = points.shape[0]
        
        # Use the flat-array method
        self._solver.set_points_from_flat(flat_list, n_points)
    
    def get_points(self) -> np.ndarray:
        """
        Get the current points from the solver.
        
        This method provides efficient data transfer between the 
        ShapeOp solver and NumPy with minimal copying.
        
        Returns
        -------
        numpy.ndarray
            Array of shape (n_points, 3) containing the current point positions.
        """
        # Get the flat points and number of points
        flat_points, n_points = self._solver.get_points_as_flat()
        
        # Convert the flat list to a NumPy array
        points_array = np.array(flat_points, dtype=np.float64)
        
        # Reshape to (n_points, 3)
        return points_array.reshape(n_points, 3)
    
    def initialize(self, dynamic: bool = False, masses: float = 1.0, 
                  damping: float = 1.0, timestep: float = 1.0) -> bool:
        """
        Initialize the solver.
        
        Parameters
        ----------
        dynamic : bool, optional
            Whether to use dynamic solving, by default False
        masses : float, optional
            Mass to use for all points, by default 1.0
        damping : float, optional
            Damping coefficient, by default 1.0
        timestep : float, optional
            Timestep for dynamic solving, by default 1.0
            
        Returns
        -------
        bool
            True if initialization was successful
        """
        return self._solver.initialize(dynamic, masses, damping, timestep)
    
    def solve(self, iterations: int = 1) -> bool:
        """
        Run the solver for a given number of iterations.
        
        Parameters
        ----------
        iterations : int, optional
            Number of iterations to run, by default 1
            
        Returns
        -------
        bool
            True if solving was successful, False otherwise.
        """
        self._current_iteration += iterations
        return self._solver.solve(iterations)
