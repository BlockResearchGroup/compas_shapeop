"""
Grid Circularization Example using COMPAS ShapeOp

This example demonstrates the Circle constraint from ShapeOp applied to 
every face of a grid mesh with random height variations.
A small vertical force is applied to create a dynamic effect.
"""

import numpy as np
import math
import random
from compas_shapeop import _shapeop as so
from compas.datastructures import Mesh
from compas.geometry import Point, Circle, Polygon
from compas.colors import Color
from compas_viewer import Viewer

# Grid parameters
rows = 8
cols = 8
spacing = 1.0
grid_scale = 0.8  # Scale to make the circle constraint more visible
height_noise = 0.5  # Maximum random height variation
gravity_force = 0.1  # Vertical force magnitude (smaller than dynamic_simulation to maintain stability)

# Set random seed for reproducibility
random.seed(42)

# Helper to get index from grid coordinates
def index(x, y):
    return y * cols + x

# Create a mesh we'll update during the simulation
mesh = Mesh()

# Initialize viewer
viewer = Viewer()

# Create grid points with random height variations
points = []
for y in range(rows):
    for x in range(cols):
        i = index(x, y)
        norm_x = float(x) / (cols - 1)
        norm_y = float(y) / (rows - 1)
        
        # Add random height variation (z coordinate)
        # Keep perimeter points flat for better Circle constraint application
        is_perimeter = (x == 0 or y == 0 or x == cols-1 or y == rows-1)
        random_height = 0.0 if is_perimeter else random.uniform(-height_noise, height_noise)
        
        points.append([
            (norm_x - 0.5) * spacing * (cols - 1) * grid_scale,  # Center the grid
            (norm_y - 0.5) * spacing * (rows - 1) * grid_scale,
            random_height*0
        ])

print(f"Created {len(points)} grid points with random height variations")

# Initialize solver
print("Creating solver...")
solver = so.Solver()

# Set points to the solver
print(f"Setting {len(points)} points...")
so.set_points(solver, points)

# Add constraints
print("Adding constraints...")

# 1. Pin the four corners to maintain structure
corner_indices = [
    index(0, 0),               # Top-left corner
    index(cols - 1, 0),        # Top-right corner
    index(0, rows - 1),        # Bottom-left corner
    # index(cols - 1, rows - 1)  # Bottom-right corner
]

corner_weight = 1000000.0
for corner_idx in corner_indices:
    so.add_constraint(solver, "Closeness", [corner_idx], corner_weight)
print(f"Fixed {len(corner_indices)} corner points")

# 2. Add edge constraints for grid structure (optional but helps stability)
edge_weight = 5.0
edge_count = 0

# Horizontal edges
for y in range(rows):
    for x in range(cols - 1):
        i = index(x, y)
        j = index(x + 1, y)
        so.add_constraint(solver, "EdgeStrain", [i, j], edge_weight)
        edge_count += 1

# Vertical edges
for x in range(cols):
    for y in range(rows - 1):
        i = index(x, y)
        j = index(x, y + 1)
        so.add_constraint(solver, "EdgeStrain", [i, j], edge_weight)
        edge_count += 1

print(f"Added {edge_count} edge constraints")

# 3. Add Circle constraint to each face of the mesh (the key change)
circle_weight = 500.0  # Higher weight makes the circle constraint stronger
circle_count = 0

# Create the face structure of the mesh while applying circle constraints
faces = []
for y in range(rows - 1):
    for x in range(cols - 1):
        v1 = index(x, y)
        v2 = index(x + 1, y)
        v3 = index(x + 1, y + 1)
        v4 = index(x, y + 1)
        
        # Add face to mesh structure
        face = [v1, v2, v3, v4]
        faces.append(face)
        
        # Apply Circle constraint to each face
        so.add_constraint(solver, "Circle", face, circle_weight)
        circle_count += 1

print(f"Added Circle constraint to {circle_count} faces")

# 4. Add vertical force (gravity) to all non-fixed points
print("Adding vertical force...")
for i in range(rows * cols):
    if i not in corner_indices:
        # Apply a small downward force (positive Z in this case)
        so.add_vertex_force(solver, 0.0, 0.0, gravity_force, i)

# We still track the perimeter for measurement purposes
perimeter_indices = []
# Top edge
for x in range(cols):
    perimeter_indices.append(index(x, 0))
# Right edge
for y in range(1, rows):
    perimeter_indices.append(index(cols - 1, y))
# Bottom edge
for x in range(cols - 2, -1, -1):
    perimeter_indices.append(index(x, rows - 1))
# Left edge
for y in range(rows - 2, 0, -1):
    perimeter_indices.append(index(0, y))

# Initialize the solver
print("Initializing solver...")
init_result = so.init_solver(solver)
print(f"Solver initialized: {init_result}")

# Add vertices to mesh
for i, p in enumerate(points):
    mesh.add_vertex(i, x=p[0], y=p[1], z=p[2])

# Add faces to create a grid mesh
for face in faces:
    mesh.add_face(face)

# Set up visualization with colored vertices
mesh_obj = viewer.scene.add(mesh, show_points=True)

# Color the mesh: corners in red, other points in blue
for i in range(len(points)):
    if i in corner_indices:
        mesh.vertex_attribute(i, 'color', Color.from_hex("#e74c3c"))  # Red for corners
    else:
        mesh.vertex_attribute(i, 'color', Color.from_hex("#3498db"))  # Blue for other points

# Simulation parameters
max_iterations = 1000
current_iteration = 0
steps_per_frame = 5

# Function to measure how circular the faces are
def measure_face_circularity(points, faces):
    """Measure the average circularity of all faces in the mesh (1.0 = perfect circle)"""
    if not faces:
        return 0.0
    
    total_circularity = 0.0
    
    for face in faces:
        # For each face, calculate its circularity
        face_points = [points[idx] for idx in face]
        
        # Calculate centroid of the face
        center_x = sum(p[0] for p in face_points) / len(face_points)
        center_y = sum(p[1] for p in face_points) / len(face_points)
        
        # Calculate radii (ignore z-component for circularity measurement)
        radii = []
        for p in face_points:
            dx = p[0] - center_x
            dy = p[1] - center_y
            radii.append(math.sqrt(dx*dx + dy*dy))
        
        avg_radius = sum(radii) / len(radii)
        std_dev = np.std(radii) if radii else 0.0
        
        # Calculate circularity for this face
        face_circularity = max(0.0, 1.0 - min(1.0, (std_dev / (avg_radius + 1e-10))))
        total_circularity += face_circularity
    
    # Return average circularity across all faces
    return total_circularity / len(faces)

# Measure initial circularity
init_circularity = measure_face_circularity(points, faces)
print(f"Initial average face circularity score: {init_circularity:.2f} (1.0 = perfect circles)")

# Create initial circle meshes
circles = []
circle_objs = []

# Create initial circle meshes once
for face_idx, face in enumerate(faces):
    face_points = [Point(*points[idx]) for idx in face]
    circle = Circle.from_points(face_points)
    circle_poly = circle.to_polygon(50)
    circle_mesh = Mesh.from_vertices_and_faces(*circle_poly.to_vertices_and_faces())
    circles.append(circle_mesh)
    circle_objs.append(viewer.scene.add(circle_mesh, show_faces=False, linecolor=Color.red()))


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
        
        # Print progress every 25 iterations
        if current_iteration % 25 == 0:
            updated_points = so.get_points(solver)
            cur_circularity = measure_face_circularity(updated_points, faces)
            
            improvement = (cur_circularity - init_circularity) / max(0.01, 1.0 - init_circularity) * 100.0
            
            print(f"Iteration {current_iteration}/{max_iterations} - " 
                  f"Avg Face Circularity: {cur_circularity:.2f} " 
                  f"(Improvement: {improvement:.1f}%)")
    
    # Update mesh with new point positions
    updated_points = so.get_points(solver)
    
    for i, vertex in enumerate(mesh.vertices()):
        if i < len(updated_points):
            mesh.vertex_attributes(vertex, 'xyz', updated_points[i])
    
    mesh_obj.update(update_data=True)

    # # Update circle meshes
    # for idx, face in enumerate(faces):
    #     face_points = [Point(*updated_points[i]) for i in face]
        
    #     circle = Circle.from_points(face_points)
    #     circle_poly = circle.to_polygon(50)
        
    #     vertices, poly_faces = circle_poly.to_vertices_and_faces()
    #     circles[idx].clear()
    #     for i, coords in enumerate(vertices):
    #         circles[idx].add_vertex(i, x=coords[0], y=coords[1], z=coords[2])
            
    #     for f in poly_faces:
    #         circles[idx].add_face(f)

    #     circle_objs[idx].update(update_data=True)


print("Starting grid faces to circles transformation with vertical force...")
print("Press ESC to exit viewer")
viewer.show()

# Measure final circularity
final_points = so.get_points(solver)
final_circularity = measure_face_circularity(final_points, faces)

# Calculate improvement
improvement = (final_circularity - init_circularity) / max(0.01, 1.0 - init_circularity) * 100.0

print("\nFinal Results:")
print(f"Final average face circularity score: {final_circularity:.2f} (1.0 = perfect circles)")
print(f"Circularity improvement: {improvement:.1f}%")

print("Simulation complete!")
