import numpy as np
from compas_shapeop import Solver
from compas.datastructures import Mesh
from compas_viewer import Viewer

###############################################################################################
# Create balloon mesh and prepare for solver
###############################################################################################
mesh = Mesh.from_obj('data/circular_inflation.obj')
vertices_list = list(mesh.vertices())
faces_list = list(mesh.faces())
vertices, faces = mesh.to_vertices_and_faces()

###############################################################################################
# Initialize solver and set points
###############################################################################################
solver = Solver(1000)
solver.set_points(vertices)

###############################################################################################
# Add constraints to maintain mesh structure
###############################################################################################
edge_weight = 100.0
edge_count = 0
for u, v in mesh.edges():
    solver.add_edge_strain_constraint([u, v], edge_weight, 0.95, 1.05)
    edge_count += 1

anchor_weight = 0.1
for vertex in vertices_list:
    solver.add_constraint("Closeness", [vertex], anchor_weight)

###############################################################################################
# Initialize solver
###############################################################################################
solver.init()

###############################################################################################
# Prepare normal force data
###############################################################################################
faces_flat = []
face_sizes = []
for face in faces:
    face_size = len(face)
    face_sizes.append(face_size)
    faces_flat.extend(face)

faces_flat_np = np.array(faces_flat, dtype=np.int32)
face_sizes_np = np.array(face_sizes, dtype=np.int32)

inflation_force = 100
solver.add_normal_force_with_faces(faces_flat_np, face_sizes_np, inflation_force)

###############################################################################################
# Setup viewer and animation
###############################################################################################
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)

current_iteration = 0

@viewer.on(interval=1)
def inflate_balloon(frame):
    global current_iteration
    
    if current_iteration >= solver.max_iterations:
        return

    solver.solve(100)
    current_iteration += 100
    
    updated_points = solver.get_points()
    
    for i, vertex in enumerate(vertices_list):
        mesh.vertex_attributes(vertex, 'xyz', updated_points[i])

    mesh_obj.update(update_data=True)

viewer.show()