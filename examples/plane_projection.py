from compas.datastructures import Mesh
from compas_viewer import Viewer
from compas_shapeop.shapeop import Solver
import compas

# ==========================================================================
# Load mesh and prepare for solver
# ==========================================================================

mesh = Mesh.from_obj(compas.get("hypar.obj"))
mesh.translate([-2.5, -2.5, -1])
mesh.scale(2)

# ==========================================================================
# Initialize solver directly from mesh
# ==========================================================================

solver = Solver.from_mesh(mesh)

# ==========================================================================
# Find boundary vertices for constraints
# ==========================================================================

# Find boundary vertices - they will be kept fixed
boundary_vertices = set(mesh.vertices_on_boundary())

# Define corner vertices (a subset of boundary vertices at the corners)
corner_vertices = []
for key in boundary_vertices:
    nbrs = mesh.vertex_neighbors(key)
    boundary_nbrs = [nbr for nbr in nbrs if nbr in boundary_vertices]
    if len(boundary_nbrs) == 2:
        corner_vertices.append(key)

# ==========================================================================
# Add constraints
# ==========================================================================

# Add closeness constraints to corners
corner_weight = 1.0
for idx in corner_vertices:
    solver.add_closeness_constraint(idx, corner_weight)

# Add edge strain constraints to maintain mesh topology
edge_weight = 1.0
for u, v in mesh.edges():
    solver.add_edge_strain_constraint(u, v, edge_weight)

# Add diagonal constraints to quads to prevent shearing
diagonal_weight = 0.2
for fkey in mesh.faces():
    vertices = mesh.face_vertices(fkey)
    if len(vertices) == 4:  # Only apply to quads
        solver.add_edge_strain_constraint(vertices[0], vertices[2], diagonal_weight)
        solver.add_edge_strain_constraint(vertices[1], vertices[3], diagonal_weight)

# Add Plane constraint to each face
plane_weight = 10000.0
faces = []
for fkey in mesh.faces():
    face_vertices = list(mesh.face_vertices(fkey))
    faces.append(face_vertices)
    solver.add_plane_constraint(face_vertices, plane_weight)

# ==========================================================================
# Initialize the solver
# ==========================================================================

solver.init()
points_ref = solver.points

# ==========================================================================
# Visualization
# ==========================================================================

viewer = Viewer()

mesh_obj = viewer.scene.add(mesh)


@viewer.on(interval=1)
def update(frame):
    # Run solver iteration
    solver.solve(1)

    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, "xyz", points_ref[i])

    mesh_obj.update(update_data=True)


viewer.show()
