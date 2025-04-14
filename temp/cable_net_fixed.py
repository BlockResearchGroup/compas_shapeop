"""
Cable net simulation that directly matches the C++ implementation
Using the new add_closeness_constraint_with_position function for corner lifting
"""
import numpy as np
from compas_shapeop import _shapeop as so
from compas.datastructures import Mesh
from compas.geometry import Point
from compas.colors import Color
from compas_viewer import Viewer

# Grid parameters
rows = 10
cols = 10
grid_size = 2.0
lift_height = 1.0  # Height to lift corners

# Helper to get index from grid coordinates
def index(x, y):
    return y * cols + x

# Create flat grid (points) for the cable net
points = []
for y in range(rows):
    for x in range(cols):
        pos_x = x * grid_size / (cols - 1)
        pos_y = y * grid_size / (rows - 1)
        pos_z = 0.0  # Flat grid to start
        points.append([pos_x, pos_y, pos_z])

# Initialize solver
print("Creating solver...")
solver = so.Solver()

# Set points to the solver
print(f"Setting {len(points)} points...")
so.set_points(solver, points)

# Define which corners to lift and which to keep fixed
lifted_corners = [
    (0, 0),           # top-left corner
    (cols-1, rows-1)  # bottom-right corner
]

ground_corners = [
    (cols-1, 0),     # top-right corner
    (0, rows-1)      # bottom-left corner
]

# Apply corner constraints - DIRECTLY matching C++ implementation
corner_weight = 1e5  # Very high weight to fix corners firmly
print("Adding corner constraints...")

# Fix the lifted corners at their elevated positions
for corner in lifted_corners:
    x, y = corner
    corner_idx = index(x, y)
    
    # Get initial position
    corner_pos = points[corner_idx].copy()
    corner_pos[2] = lift_height  # Set to the lift height
    
    # Create constraint and set the target position in one operation
    so.add_closeness_constraint_with_position(solver, corner_idx, corner_weight, corner_pos)
    print(f"Lifted corner at ({x}, {y}) to height {lift_height}")

# Fix the ground corners at z=0
for corner in ground_corners:
    x, y = corner
    corner_idx = index(x, y)
    
    # Use the original position (which has z=0)
    so.add_constraint(solver, "Closeness", [corner_idx], corner_weight)
    print(f"Fixed corner at ({x}, {y}) at ground level")

# Add edge constraints with shrinking effect EXACTLY like C++ 
print("Adding edge constraints with shrinking effect...")
edge_weight = 100.0  # Higher weight to enforce the shrinking
shrink_factor = 0.2  # Shrink to 50% of original length

# Horizontal cables
for y in range(rows):
    for x in range(cols - 1):
        edge_indices = [index(x, y), index(x+1, y)]
        
        # Target: shrink to 50% of original length
        # Range: between 45% and 55% of original (allows small variations)
        so.add_shrinking_edge_constraint(solver, edge_indices, edge_weight, shrink_factor)

# Vertical cables
for y in range(rows - 1):
    for x in range(cols):
        edge_indices = [index(x, y), index(x, y+1)]
        
        # Target: shrink to 50% of original length
        # Range: between 45% and 55% of original (allows small variations)
        so.add_shrinking_edge_constraint(solver, edge_indices, edge_weight, shrink_factor)

# Initialize and solve
print("Initializing solver...")
so.init_solver(solver)

# Create mesh for visualization
mesh = Mesh()

# Create mesh vertices
for i, point in enumerate(points):
    mesh.add_vertex(i, x=point[0], y=point[1], z=point[2])

# Create mesh faces (quads)
for y in range(rows - 1):
    for x in range(cols - 1):
        # Get the four corners of the quad (in counter-clockwise order)
        v1 = index(x, y)
        v2 = index(x+1, y)
        v3 = index(x+1, y+1)
        v4 = index(x, y+1)
        mesh.add_face([v1, v2, v3, v4])

# Set up viewer
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh, show_edges=True)

# Highlight fixed and lifted corners with spheres
for corner in lifted_corners:
    x, y = corner
    idx = index(x, y)
    # Show lifted corners in red
    point = Point(*points[idx])
    viewer.scene.add(point, size=0.1, color=Color.red())

for corner in ground_corners:
    x, y = corner
    idx = index(x, y)
    # Show ground corners in blue
    point = Point(*points[idx])
    viewer.scene.add(point, size=0.1, color=Color.blue())

# Set up the iteration counter
iteration = 0
total_iterations = 1000

@viewer.on(interval=1)  # Update every 50ms
def update_cable_net(frame):
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
    if iteration % 10 == 0:
        print(f"Iteration {iteration}/{total_iterations}")
    
    # When done, save the result to file
    if iteration == total_iterations - 1:
        print("Simulation complete!")
    
    iteration += 1

print("Starting cable net simulation...")
viewer.show()
