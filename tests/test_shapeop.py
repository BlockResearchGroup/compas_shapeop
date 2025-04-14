"""
Test script for the ShapeOp Python bindings.
"""
import numpy as np
from compas_shapeop import _shapeop as so

def test_basic_solver():
    """Test the basic functionality of the ShapeOp solver."""
    print("Creating solver...")
    # Create the solver using the constructor of our Solver class
    solver = so.Solver()
    
    # Create a simple grid of points (3x3)
    print("Creating point grid...")
    points = []
    for y in range(3):
        for x in range(3):
            points.append([float(x), float(y), 0.0])
            
    print(f"Setting {len(points)} points...")
    so.set_points(solver, points)
    
    print("Adding constraints...")
    # Pin the corners
    corner_indices = [0, 2, 6, 8]  # corners of our 3x3 grid
    for idx in corner_indices:
        so.add_constraint(solver, "Closeness", [idx], 1000.0)
    
    # Add edge strain constraints to maintain edge lengths
    edges = [
        [0, 1], [1, 2],  # top row
        [3, 4], [4, 5],  # middle row
        [6, 7], [7, 8],  # bottom row
        [0, 3], [3, 6],  # left column
        [1, 4], [4, 7],  # middle column
        [2, 5], [5, 8]   # right column
    ]
    
    for edge in edges:
        so.add_constraint(solver, "EdgeStrain", edge, 1.0)
    
    # Add a downward force to the center point
    print("Adding forces...")
    so.add_vertex_force(solver, 0.0, 0.0, -0.1, 4)  # center point (index 4)
    
    print("Initializing solver...")
    so.init_solver(solver)
    
    print("Running solver...")
    iterations = 100
    so.solve(solver, iterations)
    
    print("Getting points...")
    result_points = so.get_points(solver)
    
    print("Results:")
    for i, point in enumerate(result_points):
        print(f"Point {i}: {point}")
    
    print("Done! (Memory cleanup handled automatically)")
    # No need to manually delete the solver - our C++ wrapper handles it

if __name__ == "__main__":
    test_basic_solver()
