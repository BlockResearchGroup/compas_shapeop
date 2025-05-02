from compas.datastructures import Mesh
from compas.colors import Color
from compas.geometry import Circle, Plane, Point
from compas_viewer import Viewer
from compas_shapeop.shapeop import Solver

# ==========================================================================
# Create mesh grid and prepare for solver
# ==========================================================================

mesh = Mesh.from_meshgrid(10.0, 7, 10.0, 7)

# ==========================================================================
# Initialize solver directly from mesh
# ==========================================================================

solver = Solver.from_mesh(mesh)

# ==========================================================================
# Add constraints
# ==========================================================================

# Find corner vertices
corners = []
for v in mesh.vertices():
    if len(mesh.vertex_neighbors(v)) == 2:
        corners.append(v)

# Pin the corners with high weight
corner_weight = 1000.0
for idx in corners:
    solver.add_closeness_constraint([idx], corner_weight)

# Add edge constraints with medium weight
edge_weight = 5.0
for u, v in mesh.edges():
    solver.add_edge_strain_constraint([u, v], edge_weight)

# Add Circle constraint to each face with high weight
circle_weight = 100.0
faces = []
for face in mesh.faces():
    face_vertices = list(mesh.face_vertices(face))
    faces.append(face_vertices)
    solver.add_circle_constraint(face_vertices, circle_weight)

# Add gravity force
gravity_force = 0.5
solver.add_mesh_vertex_force(mesh, 0.0, 0.0, gravity_force)

# ==========================================================================
# Initialize the solver
# ==========================================================================

points_ref = solver.init()

# ==========================================================================
# Visualization
# ==========================================================================

viewer = Viewer()
mesh_obj = viewer.scene.add(mesh, show_points=True)

# Create circle visualizations
circles = []
circle_objs = []

for face_vertices in faces:
    face_points = [Point(*mesh.vertex_coordinates(idx)) for idx in face_vertices]
    circle = Circle.from_points(face_points)
    circles.append(circle)
    circle_objs.append(viewer.scene.add(circle, linecolor=Color.red(), u=32))

iteration = 0

@viewer.on(interval=1)
def deform_mesh(frame):
    global iteration
    if iteration >= 100:
        return

    solver.solve(5)
    
    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, 'xyz', points_ref[i])

    mesh_obj.update(update_data=True)
    
    # Update the circles
    for idx, face_vertices in enumerate(faces):
        face_points = [Point(*mesh.vertex_coordinates(vertex)) for vertex in face_vertices]
        circle = Circle.from_three_points(face_points[0], face_points[1], face_points[2])
        circles[idx].frame = circle.frame
        circles[idx].radius = circle.radius
        circle_objs[idx].update(update_data=True)

    
    iteration += 1

viewer.show()
