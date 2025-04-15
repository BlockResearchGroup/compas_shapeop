"""
Balloon Simulation using ShapeOp NormalForce implementation with interactive viewer

This example demonstrates how to simulate a balloon inflation using the ShapeOp library.
It uses the newly implemented NormalForce binding to create an inflation effect with live visualization.
"""

import numpy as np
from compas_shapeop import _shapeop as so
from compas.datastructures import Mesh
from compas_viewer import Viewer

print("Starting balloon simulation...")

# Create the balloon mesh
mesh = Mesh.from_obj('data/m0.obj')
# Convert generators to lists for length calculation
vertices_list = list(mesh.vertices())
faces_list = list(mesh.faces())
print(f"Loaded mesh with {len(vertices_list)} vertices and {len(faces_list)} faces")

# Extract vertices and faces in the format needed for ShapeOp
vertices, faces = mesh.to_vertices_and_faces()
print(f"Extracted vertices and faces data: {len(vertices)} vertices")

# Create viewer instance
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)

# Create ShapeOp solver
solver = so.Solver()
print("Created solver")

# Add points to the solver
so.set_points(solver, vertices)
print("Added points to solver")

# Create constraints to maintain the shape
# Set up edge constraints to maintain the structure
weight = 1.0
edge_count = 0
for u, v in mesh.edges():
    indices = [u, v]
    so.add_edge_strain_constraint(solver, indices, weight, 0.1, 1.0)
    edge_count += 1
print(f"Added {edge_count} edge constraints")

# Fix a few vertices to prevent rigid body motion
weight = 100.0  # Strong constraint
fixed_vertices = [0]  # Fix first few vertices
for vertex in fixed_vertices:
    position = mesh.vertex_coordinates(vertex)
    so.add_constraint(solver, "Closeness", [vertex], weight)
print(f"Fixed {len(fixed_vertices)} vertices")

# Initialize the solver
result = so.init_solver(solver)
print(f"Initialized solver: {result}")

# Prepare face data for the normal force - use face indices, not coordinates!
faces_flat = []
face_sizes = []
for face in faces:
    # Just use the vertex indices directly, don't convert to coordinates
    face_size = len(face)
    face_sizes.append(face_size)
    faces_flat.extend(face)  # Add the integer indices to flat list

# Convert Python lists to numpy arrays as nanobind expects proper integer vectors
faces_flat_np = np.array(faces_flat, dtype=np.int32)
face_sizes_np = np.array(face_sizes, dtype=np.int32)
print(f"Prepared {len(face_sizes)} faces for normal force calculation")

# Simulation parameters
max_iterations = 1000
current_iteration = 0
inflation_force = 0.15  # Increased initial force magnitude
steps_per_frame = 5  # More steps per frame for faster simulation

# Apply normal force just ONCE at the beginning (instead of every frame)
print("Applying normal force once...")
normal_force_result = so.add_normal_force_with_faces(solver, faces_flat_np, face_sizes_np, inflation_force)
print(f"Normal force applied, result: {normal_force_result}")

# Define the callback for the viewer's update
@viewer.on(interval=1)
def inflate_balloon(frame):
    global current_iteration
    
    if current_iteration >= max_iterations:
        return
    
    # Run multiple simulation steps per visual update
    for _ in range(steps_per_frame):
        if current_iteration >= max_iterations:
            break
            
        # Solve one step
        so.solve(solver, 1)
        current_iteration += 1
    
    # Update the geometry in the COMPAS mesh
    updated_points = so.get_points(solver)
    
    for i, vertex in enumerate(vertices_list):
        mesh.vertex_attributes(vertex, 'xyz', updated_points[i])

    mesh_obj.update(update_data=True)



print("Starting interactive balloon simulation...")
print("Press ESC to exit")
viewer.show()

# Clean up when viewer is closed
so.delete(solver)
print("All done!")
