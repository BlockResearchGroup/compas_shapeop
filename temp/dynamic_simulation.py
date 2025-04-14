"""
Dynamic Grid simulation using ShapeOp with interactive visualization
Python equivalent of the C++ unary_force example with live visualization
"""
import numpy as np
import os
from compas_shapeop import _shapeop as so
from compas.datastructures import Mesh
from compas.geometry import Point
from compas.colors import Color
from compas_viewer import Viewer

# Grid parameters
rows = 14
cols = 14
spacing = 1.0
gravity_force = 0.005

# Helper to get index from grid coordinates
def index(x, y):
    return y * cols + x

# Create a mesh we'll update during the simulation
mesh = Mesh()

# Initialize viewer
viewer = Viewer()

# Create grid points
points = []
for y in range(rows):
    for x in range(cols):
        i = index(x, y)
        norm_x = float(x) / (cols - 1)
        norm_y = float(y) / (rows - 1)
        points.append([
            norm_x * spacing * (cols - 1),
            -norm_y * spacing * (rows - 1),
            0.0
        ])

# Initialize solver
print("Creating solver...")
solver = so.Solver()

# Set points to the solver
print(f"Setting {len(points)} points...")
so.set_points(solver, points)

# Add constraints
print("Adding constraints...")

# Pin corners
corner_indices = [
    index(0, 0),               # Top-left corner
    index(cols - 1, 0),        # Top-right corner
    index(0, rows - 1),        # Bottom-left corner
    index(cols - 1, rows - 1)  # Bottom-right corner
]

for idx in corner_indices:
    so.add_constraint(solver, "Closeness", [idx], 1e5)

# Pin center point
center_idx = index(cols // 2, rows // 2)
so.add_constraint(solver, "Closeness", [center_idx], 1e5)

# Add edge strain constraints
print("Adding edge constraints...")
for y in range(rows):
    for x in range(cols):
        i = index(x, y)
        
        # Horizontal edges
        if x + 1 < cols:
            edge = [i, index(x + 1, y)]
            so.add_constraint(solver, "EdgeStrain", edge, 1.0)
        
        # Vertical edges
        if y + 1 < rows:
            edge = [i, index(x, y + 1)]
            so.add_constraint(solver, "EdgeStrain", edge, 1.0)

# Add gravity force
print("Adding gravity force...")
# For a uniform gravity effect, add a small downward force to every non-pinned point
fixed_indices = corner_indices + [center_idx]
for i in range(rows * cols):
    if i not in fixed_indices:
        so.add_vertex_force(solver, 0.0, 0.0, gravity_force, i)

# Initialize solver
print("Initializing solver...")
so.init_solver(solver)

# Create initial mesh structure
# Add vertices to the mesh
for i, point in enumerate(points):
    mesh.add_vertex(i, x=point[0], y=point[1], z=point[2])

# Add faces to the mesh
for y in range(rows - 1):
    for x in range(cols - 1):
        v1 = index(x, y)
        v2 = index(x+1, y)
        v3 = index(x+1, y+1)
        v4 = index(x, y+1)
        mesh.add_face([v1, v2, v3, v4])

# Add the mesh to the viewer scene
mesh_obj = viewer.scene.add(mesh)

# Highlight fixed points with spheres
for idx in fixed_indices:
    point = Point(*points[idx])
    viewer.scene.add(point, size=0.2, color=Color.red())

# Set up the iteration counter
iteration = 0
total_iterations = 1000

@viewer.on(interval=50)  # Update every 50ms
def deform_mesh(frame):
    global iteration, mesh, solver
    
    # Exit if we've reached the maximum iterations
    if iteration >= total_iterations:
        return
    
    # Update the simulation by one step
    so.solve(solver, 1)
    
    # Get the updated points
    updated_points = so.get_points(solver)
    
    # Update the mesh vertices
    for i, point in enumerate(updated_points):
        mesh.vertex_attributes(i, 'xyz', point)
    
    # Update the mesh in the viewer
    mesh_obj.update(update_data=True)
    
    # Print progress periodically
    if iteration % 50 == 0:
        print(f"Iteration {iteration}/{total_iterations}")
    
    # When done, save the result to file
    # if iteration == total_iterations - 1:
    #     print("Simulation complete!")
    #     mesh.to_obj("python_dynamic_simulation.obj")
    
    iteration += 1

print("Starting dynamic simulation viewer...")
viewer.show()
