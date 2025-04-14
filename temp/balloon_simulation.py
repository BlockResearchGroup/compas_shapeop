"""
Balloon simulation using ShapeOp with interactive visualization
Based on balloon_box.cpp example demonstrating normal force inflation
"""
import os
import numpy as np
from compas_shapeop import _shapeop as so
from compas.datastructures import Mesh
from compas.colors import Color
from compas_viewer import Viewer
import time

# Check if we have a sample mesh file or create a simple box
MESH_FILE = "data/m0.obj"


# Initialize the mesh and viewer
mesh = Mesh()
viewer = Viewer()

def read_obj(filename):
    """Read vertices and faces from an OBJ file"""
    vertices = []
    faces = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('v '):
                # Read vertex coordinates (split line and convert to float)
                parts = line.strip().split()
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                vertices.append([x, y, z])
            elif line.startswith('f '):
                # Read face indices (convert to 0-based indexing)
                parts = line.strip().split()
                # Handle different OBJ face formats (v, v/vt, v/vt/vn)
                face = []
                for p in parts[1:]:
                    # Get just the vertex index (ignore texture and normal indices)
                    idx = int(p.split('/')[0]) - 1  # OBJ uses 1-based indexing
                    face.append(idx)
                faces.append(face)
    
    return vertices, faces

# Load the mesh
vertices, faces = read_obj(MESH_FILE)
print(f"Loaded mesh with {len(vertices)} vertices and {len(faces)} faces")

# Add vertices to the COMPAS mesh
for i, (x, y, z) in enumerate(vertices):
    mesh.add_vertex(i, x=x, y=y, z=z)

# Add faces to the COMPAS mesh
for i, face in enumerate(faces):
    mesh.add_face(face)

# Initialize ShapeOp solver
print("Creating solver...")
solver = so.Solver()

# Set points to the solver
print("Setting points...")
so.set_points(solver, vertices)

# Remove solver parameter settings since they're not available in the current API
# so.set_solver_parameter(solver, "Iterations", 5)
# so.set_solver_parameter(solver, "Omega", 1.9)

# Add weak closeness constraints to all vertices to maintain stability
print("Adding closeness constraints...")
small_stiffness = 0.00001
for i in range(len(vertices)):
    so.add_constraint(solver, "Closeness", [i], small_stiffness)

# Add strong closeness constraints to keep some vertices fixed
# For our simple box, let's fix a few corners
# Fixing more vertices for a more controlled inflation
fixed_vertices = [0, 150,350, 75]

print(f"Fixing {len(fixed_vertices)} vertices")
for idx in fixed_vertices:
    so.add_constraint(solver, "Closeness", [idx], 1000.0)
    
# Add edge strain constraints
print("Adding edge constraints...")
# Get the edges from the mesh to avoid duplicates
processed_edges = set()
for face in faces:
    for i in range(len(face)):
        v1 = face[i]
        v2 = face[(i + 1) % len(face)]
        
        # Ensure consistent edge ordering
        edge = tuple(sorted([v1, v2]))
        
        if edge not in processed_edges:
            processed_edges.add(edge)
            so.add_constraint(solver, "EdgeStrain", list(edge), 0.09)

# Initialize solver
print("Initializing solver...")
so.init_solver(solver)

# Add normal forces (inflation forces)
print("Adding normal forces...")
# We'll compute normals per vertex for each iteration

# Calculate face normals
def calculate_face_normals(vertices, faces):
    """Calculate normal vectors for each face"""
    face_normals = []
    
    for face in faces:
        if len(face) >= 3:  # Ensure face has at least 3 vertices
            v0 = np.array(vertices[face[0]])
            v1 = np.array(vertices[face[1]])
            v2 = np.array(vertices[face[2]])
            
            # Calculate normal using cross product
            normal = np.cross(v1 - v0, v2 - v0)
            
            # Normalize
            norm = np.linalg.norm(normal)
            if norm > 0:
                normal = normal / norm
            
            face_normals.append(normal)
        else:
            # Degenerate face, add zero normal
            face_normals.append(np.array([0, 0, 0]))
    
    return face_normals

# Associate faces with vertices for efficient normal force calculation
vertex_to_faces = {}
for i, face in enumerate(faces):
    for v_idx in face:
        if v_idx not in vertex_to_faces:
            vertex_to_faces[v_idx] = []
        vertex_to_faces[v_idx].append(i)

# Add the mesh to the viewer
mesh_obj = viewer.scene.add(mesh)

# Highlight fixed points with spheres
# for idx in fixed_vertices:
#     point = Point(*vertices[idx])
#     viewer.scene.add(point, size=0.1, color=Color.red())

# Setup simulation parameters
inflation_force = 0.05 # Increased for more visible effect
max_iterations = 10000
current_iteration = 0
steps_per_frame = 1  # Reduced to slow down simulation

# Get current point positions
current_points = so.get_points(solver)

# Calculate face normals
face_normals = calculate_face_normals(current_points, faces)

# Apply normal forces to each vertex (inflation)
for vertex_idx in range(len(current_points)):
    if vertex_idx in vertex_to_faces:
        # Get all faces this vertex belongs to
        vertex_faces = vertex_to_faces[vertex_idx]
        
        # Calculate average normal for this vertex
        if vertex_faces:
            avg_normal = np.zeros(3)
            for face_idx in vertex_faces:
                if face_idx < len(face_normals):
                    avg_normal += face_normals[face_idx]
            
            # Normalize
            norm = np.linalg.norm(avg_normal)
            if norm > 0:
                avg_normal = avg_normal / norm
                
                # Apply normal force
                force = inflation_force * avg_normal
                so.add_vertex_force(solver, force[0], force[1], force[2], vertex_idx)

@viewer.on(interval=1)  # Increased interval to slow down updates
def inflate_balloon(frame):
    global current_iteration, vertices, mesh
    
    if current_iteration >= max_iterations:
        return
    
    # Overall timing
    start_time = time.time()
    
    # Run multiple simulation steps per visual update
    solve_start = time.time()
    for _ in range(steps_per_frame):
        if current_iteration >= max_iterations:
            break
        
        # Step the simulation
        so.solve(solver, 1)
        current_iteration += 1
    solve_time = time.time() - solve_start
    
    # Update the geometry in the COMPAS mesh
    get_points_start = time.time()
    updated_points = so.get_points(solver)
    get_points_time = time.time() - get_points_start
    
    vertex_update_start = time.time()
    for i, point in enumerate(updated_points):
        mesh.vertex_attributes(i, 'xyz', point)
    vertex_update_time = time.time() - vertex_update_start
    
    # Update the visualization
    vis_update_start = time.time()
    mesh_obj.update(update_data=True)
    vis_update_time = time.time() - vis_update_start
    
    # Total frame time
    total_time = time.time() - start_time
    
    # Print progress and timing info
    if current_iteration % 20 != 0 or current_iteration <= max_iterations:
        print(f"Iteration {current_iteration}/{max_iterations}")
        print(f"Timing (seconds):")
        print(f"  Solve:         {solve_time:.6f}")
        print(f"  Get points:    {get_points_time:.6f}")
        print(f"  Update mesh:   {vertex_update_time:.6f}")
        print(f"  Update view:   {vis_update_time:.6f}")
        print(f"  Total frame:   {total_time:.6f}")
        print(f"  FPS:           {1.0/total_time:.2f}")
    
print("Starting balloon simulation...")
viewer.show()
