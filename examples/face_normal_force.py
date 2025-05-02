from compas.datastructures import Mesh
from compas_viewer import Viewer
from compas_shapeop.shapeop import Solver

# ==========================================================================
# Create balloon mesh and prepare for solver
# ==========================================================================

# Load the mesh
mesh = Mesh.from_obj('data/m0.obj')
vertices_list = list(mesh.vertices())
faces_list = list(mesh.faces())

# ==========================================================================
# Initialize solver directly from mesh
# ==========================================================================

solver = Solver.from_mesh(mesh)

# ==========================================================================
# Add constraints to maintain mesh structure
# ==========================================================================

# Add edge strain constraints to maintain mesh structure
edge_weight = 100.0
for u, v in mesh.edges():
    solver.add_edge_strain_constraint([u, v], edge_weight, 0.95, 1.05)

# Add weak closeness constraints to all vertices to maintain stability
anchor_weight = 10
for vertex in vertices_list:
    solver.add_closeness_constraint([vertex], anchor_weight)

# ==========================================================================
# Initialize solver
# ==========================================================================

points_ref = solver.init()

# ==========================================================================
# Prepare normal force data
# ==========================================================================

faces_flat = []
face_sizes = []
for face in faces_list:
    face_vertices = mesh.face_vertices(face)
    face_size = len(face_vertices)
    face_sizes.append(face_size)
    faces_flat.extend(face_vertices)

inflation_force = 100
solver.add_normal_force_with_faces(faces_flat, face_sizes, inflation_force)

# ==========================================================================
# Setup viewer and animation
# ==========================================================================

viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)

@viewer.on(interval=1)
def update(frame):
    solver.solve(10)
    
    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, 'xyz', points_ref[i])

    mesh_obj.update(update_data=True)

viewer.show()
