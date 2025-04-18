"""
Direct Solver Demo
==================

This example demonstrates the basic usage of the DirectSolver class,
which provides a simple interface to the ShapeOp solver with efficient data transfer.
"""

import numpy as np
import time
from compas_shapeop.direct_solver import DirectSolver

# Create some simple test points (a cube)
vertices = np.array([
    [0.0, 0.0, 0.0],  # 0
    [1.0, 0.0, 0.0],  # 1
    [1.0, 1.0, 0.0],  # 2
    [0.0, 1.0, 0.0],  # 3
    [0.0, 0.0, 1.0],  # 4
    [1.0, 0.0, 1.0],  # 5
    [1.0, 1.0, 1.0],  # 6
    [0.0, 1.0, 1.0],  # 7
])

print(f"Created {len(vertices)} vertices")

# Create a DirectSolver
print("Creating DirectSolver...")
solver = DirectSolver()

# Set the solver parameters
solver.max_iterations = 5

# Set the points in the solver
print("Setting points...")
start_time = time.time()
solver.set_points(vertices)
print(f"Points set in {time.time() - start_time:.6f} seconds")

# Initialize the solver 
print("Initializing solver...")
solver.initialize(dynamic=False)

# Get the points back to verify
print("Getting points...")
start_time = time.time()
points = solver.get_points()
print(f"Points retrieved in {time.time() - start_time:.6f} seconds")

# Verify that the points match
print("\nVerifying points match:")
for i in range(len(vertices)):
    original = vertices[i]
    retrieved = points[i]
    print(f"Original: {original}, Retrieved: {retrieved}")
    assert np.allclose(original, retrieved), f"Point {i} doesn't match!"

# Try running the solver
print("\nRunning solver for 1 iteration...")
start_time = time.time()
solver.solve(1)
print(f"Solve completed in {time.time() - start_time:.6f} seconds")

print("\nTest completed successfully!")
