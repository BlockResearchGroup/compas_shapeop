"""
Mesh Face Planarization Example using COMPAS ShapeOp

This example demonstrates how to planarize the faces of a mesh using 
the Plane constraint from ShapeOp. The hypar mesh from COMPAS examples is used,
and then ShapeOp is used to optimize the face planarity.
"""

import numpy as np
import compas
import math
import random
from compas_shapeop import _shapeop as so
from compas.datastructures import Mesh
from compas.geometry import Point, Plane, Polygon
from compas.colors import Color, ColorMap
from compas_viewer import Viewer

# Load the hypar mesh from COMPAS examples
print("Loading hypar mesh...")
mesh = Mesh.from_obj(compas.get('hypar.obj'))

# Print mesh statistics
print(f"Mesh loaded with {len(list(mesh.vertices()))} vertices and {len(list(mesh.faces()))} faces")

# Initialize viewer
viewer = Viewer()

# Get mesh coordinates
vertices = list(mesh.vertices())
xyz = []
for key in vertices:
    xyz.extend(mesh.vertex_coordinates(key))

# Convert to the format needed by ShapeOp
points = []
for i in range(0, len(xyz), 3):
    if i + 2 < len(xyz):
        points.append([xyz[i], xyz[i+1], xyz[i+2]])

print(f"Prepared {len(points)} points for ShapeOp")

# Initialize solver
print("Creating solver...")
solver = so.Solver()

# Set points to the solver
print(f"Setting {len(points)} points...")
so.set_points(solver, points)

# Add constraints
print("Adding constraints...")

# Create a key to index mapping
key_index = {key: i for i, key in enumerate(vertices)}

# Find boundary vertices to fix
boundary_vertices = list(mesh.vertices_on_boundary())
print(f"Found {len(boundary_vertices)} boundary vertices")

# Find corner vertices (vertices with exactly two boundary edges)
corner_vertices = []
for key in boundary_vertices:
    nbrs = mesh.vertex_neighbors(key)
    boundary_nbrs = [nbr for nbr in nbrs if nbr in boundary_vertices]
    if len(boundary_nbrs) == 2:
        corner_vertices.append(key)

print(f"Found {len(corner_vertices)} corner vertices to fix")

# Add closeness constraints only to corners to allow more deformation
corner_weight = 1.0  # Higher weight to keep corners fixed
for key in corner_vertices:
    vertex_idx = key_index[key]
    so.add_constraint(solver, "Closeness", [vertex_idx], corner_weight)

# Add edge constraints to maintain mesh structure
edge_weight = 1.0
edge_count = 0

# Add constraints for existing edges
for u, v in mesh.edges():
    if u in key_index and v in key_index:
        i = key_index[u]
        j = key_index[v]
        so.add_constraint(solver, "EdgeStrain", [i, j], edge_weight)
        edge_count += 1

# Add diagonal constraints for quad faces
diagonal_weight = 0.5  # Lower weight for diagonals
diagonal_count = 0

for fkey in mesh.faces():
    vertices = mesh.face_vertices(fkey)
    # Check if it's a quad face
    if len(vertices) == 4:
        # Add diagonal constraints (0-2 and 1-3)
        if all(v in key_index for v in vertices):
            # First diagonal
            i = key_index[vertices[0]]
            j = key_index[vertices[2]]
            so.add_constraint(solver, "EdgeStrain", [i, j], diagonal_weight)
            
            # Second diagonal
            i = key_index[vertices[1]]
            j = key_index[vertices[3]]
            so.add_constraint(solver, "EdgeStrain", [i, j], diagonal_weight)
            
            diagonal_count += 2

print(f"Added {edge_count} edge constraints and {diagonal_count} diagonal constraints")

# Add Plane constraint to each face
plane_weight = 10000.0  # Increased weight for stronger planarization effect
plane_count = 0

# Collect faces and apply plane constraints
faces = []
for fkey in mesh.faces():
    face_vertices = [key_index[key] for key in mesh.face_vertices(fkey)]
    # Add the face for later planarity measurement
    faces.append(face_vertices)
    # Apply Plane constraint
    so.add_constraint(solver, "Plane", face_vertices, plane_weight)
    plane_count += 1

print(f"Added Plane constraint to {plane_count} faces")

# Initialize the solver
print("Initializing solver...")
init_result = so.init_solver(solver)
print(f"Solver initialized: {init_result}")

# Add the mesh to the scene
mesh_obj = viewer.scene.add(mesh)

# Color the mesh: boundary vertices in red, other points in blue
for key in mesh.vertices():
    if key in boundary_vertices:
        mesh.vertex_attribute(key, 'color', Color.from_hex("#e74c3c"))  # Red for boundary
    else:
        mesh.vertex_attribute(key, 'color', Color.from_hex("#3498db"))  # Blue for other points

# Create a colormap for visualizing planarity
colormap = ColorMap.from_mpl('viridis')

# Function to measure face planarity
def measure_face_planarity(points, faces):
    """Measure the maximum distance from any vertex to the best-fit plane for each face."""
    if not faces:
        return {}
    
    planarity = {}
    
    for i, face in enumerate(faces):
        # Get face vertices
        face_points = [points[idx] for idx in face]
        if len(face_points) < 3:
            planarity[i] = 0.0
            continue
        
        # Calculate face center
        center = [sum(p[j] for p in face_points) / len(face_points) for j in range(3)]
        
        # Calculate normal using cross product of two edges
        v0 = [face_points[1][j] - face_points[0][j] for j in range(3)]
        v1 = [face_points[2][j] - face_points[0][j] for j in range(3)]
        normal = [
            v0[1] * v1[2] - v0[2] * v1[1],
            v0[2] * v1[0] - v0[0] * v1[2],
            v0[0] * v1[1] - v0[1] * v1[0]
        ]
        
        # Normalize normal
        length = math.sqrt(sum(n*n for n in normal))
        if length > 0:
            normal = [n/length for n in normal]
        else:
            planarity[i] = 0.0
            continue
        
        # Find max distance to plane
        max_dist = 0.0
        for p in face_points:
            # Vector from center to point
            vec = [p[j] - center[j] for j in range(3)]
            # Distance to plane is dot product with normal
            dist = abs(sum(vec[j] * normal[j] for j in range(3)))
            max_dist = max(max_dist, dist)
        
        planarity[i] = max_dist
    
    return planarity

# Function to calculate average planarity
def calculate_average_planarity(planarity_dict):
    """Calculate the average planarity of all faces"""
    if not planarity_dict:
        return 0.0
    return sum(planarity_dict.values()) / len(planarity_dict)

# Measure initial planarity
initial_planarity = measure_face_planarity(points, faces)
avg_initial_planarity = calculate_average_planarity(initial_planarity)
print(f"Initial average face planarity: {avg_initial_planarity:.6f} (0.0 = perfectly planar)")

# Color the faces by planarity
max_planarity_for_color = 0.2  # Scale for coloring
for i, face in enumerate(faces):
    if i in initial_planarity:
        planar_value = min(1.0, initial_planarity[i] / max_planarity_for_color)
        color = colormap(1.0 - planar_value)  # Invert so red = high deviation
        # Get the actual face key from mesh
        face_key = list(mesh.faces())[i] if i < len(list(mesh.faces())) else None
        if face_key is not None:
            mesh.face_attribute(face_key, 'color', color)

# Simulation parameters
max_iterations = 1000
current_iteration = 0
steps_per_frame = 1

@viewer.on(interval=1)
def update(frame):
    global current_iteration
    
    if current_iteration >= max_iterations:
        return
    
    # Run multiple simulation steps per frame
    for _ in range(steps_per_frame):
        if current_iteration >= max_iterations:
            break
        
        # Solve one step
        so.solve(solver, 1)
        current_iteration += 1
        
        # Print progress every 10 iterations
        if current_iteration % 1 == 0:
            updated_points = so.get_points(solver)
            planarity = measure_face_planarity(updated_points, faces)
            avg_planarity = calculate_average_planarity(planarity)
            
            improvement = (avg_initial_planarity - avg_planarity) / avg_initial_planarity * 100.0
            
            print(f"Iteration {current_iteration}/{max_iterations} - " 
                  f"Avg Face Planarity: {avg_planarity:.6f} " 
                  f"(Improvement: {improvement:.1f}%)")
    
    # Update mesh with new point positions
    updated_points = so.get_points(solver)
    
    # Update vertex positions
    for key, idx in key_index.items():
        if idx < len(updated_points):
            mesh.vertex_attributes(key, 'xyz', updated_points[idx])
    
    # Update face colors based on new planarity
    planarity = measure_face_planarity(updated_points, faces)
    for i, face in enumerate(faces):
        if i in planarity:
            planar_value = min(1.0, planarity[i] / max_planarity_for_color)
            color = colormap(1.0 - planar_value)  # Invert so red = high deviation
            # Get the actual face key from mesh
            face_key = list(mesh.faces())[i] if i < len(list(mesh.faces())) else None
            if face_key is not None:
                mesh.face_attribute(face_key, 'color', color)
    
    # Update the mesh in the viewer
    mesh_obj.update(update_data=True)

print("Starting mesh face planarization...")
print("Press ESC to exit viewer")
viewer.show()

# Measure final planarity
final_points = so.get_points(solver)
final_planarity = measure_face_planarity(final_points, faces)
avg_final_planarity = calculate_average_planarity(final_planarity)

# Calculate improvement
improvement = (avg_initial_planarity - avg_final_planarity) / avg_initial_planarity * 100.0

print("\nFinal Results:")
print(f"Initial average face planarity: {avg_initial_planarity:.6f}")
print(f"Final average face planarity: {avg_final_planarity:.6f}")
print(f"Planarity improvement: {improvement:.1f}%")

print("Simulation complete!")
