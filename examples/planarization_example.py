import compas
from compas_shapeop import Solver
from compas.datastructures import Mesh
from compas.colors import Color
from compas_viewer import Viewer

###############################################################################################
# Load mesh and prepare for solver
###############################################################################################
mesh = Mesh.from_obj(compas.get('hypar.obj'))
points = mesh.to_vertices_and_faces()[0]

###############################################################################################
# Initialize solver and set points
###############################################################################################
solver = Solver(1000)
solver.set_points(points)

###############################################################################################
# Add constraints
###############################################################################################
# Find boundary vertices to fix
boundary_vertices = list(mesh.vertices_on_boundary())

# Find corner vertices (vertices with exactly two boundary edges)
corner_vertices = []
for key in boundary_vertices:
    nbrs = mesh.vertex_neighbors(key)
    boundary_nbrs = [nbr for nbr in nbrs if nbr in boundary_vertices]
    if len(boundary_nbrs) == 2:
        corner_vertices.append(key)

# Add closeness constraints to corners
corner_weight = 1.0
for idx in corner_vertices:
    solver.add_constraint("Closeness", [idx], corner_weight)

# Add edge constraints to maintain mesh structure
edge_weight = 1.0
for u, v in mesh.edges():
    solver.add_constraint("EdgeStrain", [u, v], edge_weight)

# Add diagonal constraints for quad faces
diagonal_weight = 0.5
for fkey in mesh.faces():
    vertices = mesh.face_vertices(fkey)
    if len(vertices) == 4:
        # First diagonal
        solver.add_constraint("EdgeStrain", [vertices[0], vertices[2]], diagonal_weight)
        # Second diagonal
        solver.add_constraint("EdgeStrain", [vertices[1], vertices[3]], diagonal_weight)

# Add Plane constraint to each face
plane_weight = 10000.0
for fkey in mesh.faces():
    face_vertices = list(mesh.face_vertices(fkey))
    solver.add_constraint("Plane", face_vertices, plane_weight)

###############################################################################################
# Initialize solver
###############################################################################################
solver.init()

###############################################################################################
# Setup visualization
###############################################################################################
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)

# Color the mesh: boundary vertices in red, other points in blue
for key in mesh.vertices():
    color = Color.from_hex("#e74c3c") if key in boundary_vertices else Color.from_hex("#3498db")
    mesh.vertex_attribute(key, 'color', color)

###############################################################################################
# Animation
###############################################################################################
max_iterations = 1000
current_iteration = 0

@viewer.on(interval=1)
def update(frame):
    global current_iteration
    
    if current_iteration >= solver.max_iterations:
        return

    solver.solve(1)
    current_iteration += 1
    
    updated_points = solver.get_points()

    for i, vertex in enumerate(mesh.vertices()):
        if i < len(updated_points):
            mesh.vertex_attributes(vertex, 'xyz', updated_points[i])
    
    mesh_obj.update(update_data=True)

viewer.show()