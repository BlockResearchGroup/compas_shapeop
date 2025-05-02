import compas
from compas.datastructures import Mesh
from compas_viewer import Viewer
from compas_shapeop import Solver

# ==========================================================================
# Create mesh
# ==========================================================================

mesh = compas.json_load("data/hex_mesh.json")
mesh = Mesh.from_polygons(mesh.to_polygons())
mesh.unify_cycles()
mesh.weld()
mesh.scale(10)
mesh.translate([-20, -20, -1])

# ==========================================================================
# Initialize solver with mesh vertices
# ==========================================================================

solver = Solver.from_mesh(mesh)

# ==========================================================================
# Add constraints
# ==========================================================================

# solver.add_mesh_closeness_constraint(mesh, mesh.vertices(), weight=0.01)
solver.add_mesh_edge_strain_constraint(mesh, weight=0.1, min_range=0.7, max_range=1.2)


# Add plane constraints to all faces for planarization
face_weight = 10
for fkey in mesh.faces():
    face_vertices = mesh.face_vertices(fkey)
    if len(face_vertices) > 3:  # Triangles are already planar
        solver.add_plane_constraint(face_vertices, face_weight)

# Add regular polygon constraints to regularize face shapes
polygon_weight = 0.5
for fkey in mesh.faces():
    face_vertices = mesh.face_vertices(fkey)
    if len(face_vertices) > 3:  # Only apply to non-triangular faces
        solver.add_regular_polygon_constraint(face_vertices, polygon_weight)

# Fix boundary vertices with closeness constraints
boundary_weight = 0.01
boundary_vertices = set(mesh.vertices_on_boundary())
for idx in boundary_vertices:
    solver.add_closeness_constraint([idx], boundary_weight)


# ==========================================================================
# Initialize the solver
# ==========================================================================

points_ref = solver.init()

# ==========================================================================
# Visualization
# ==========================================================================
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)

for line in mesh.to_lines():
    viewer.scene.add(compas.geometry.Line(*line))

for v in boundary_vertices:
    viewer.scene.add(compas.geometry.Point(*mesh.vertex_coordinates(v)))

@viewer.on(interval=1)  
def deform_mesh(frame):

    solver.solve(10)
    
    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, 'xyz', points_ref[i])
    
    mesh_obj.update(update_data=True)

viewer.show()
