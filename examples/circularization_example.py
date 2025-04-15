from compas_shapeop import Solver
from compas.datastructures import Mesh
from compas.geometry import Point, Circle
from compas.colors import Color
from compas_viewer import Viewer

###############################################################################################
# Create mesh grid and prepare for solver
###############################################################################################
mesh = Mesh.from_meshgrid(10.0, 7, 10.0, 7)

###############################################################################################
# Initialize solver and set points
###############################################################################################
points = mesh.to_vertices_and_faces()[0]
solver = Solver()
solver.set_points(points)

###############################################################################################
# Add constraints
###############################################################################################
# Find corner vertices
corners = []
for v in mesh.vertices():
    if len(mesh.vertex_neighbors(v)) == 2:
        corners.append(v)

# Pin the corners
corner_weight = 1000.0
for idx in corners:
    solver.add_constraint("Closeness", [idx], corner_weight)

# Add edge constraints
edge_weight = 5.0
for u, v in mesh.edges():
    solver.add_constraint("EdgeStrain", [u, v], edge_weight)

# Add Circle constraint to each face
circle_weight = 100.0
faces = []
for face in mesh.faces():
    face_vertices = [key for key in mesh.face_vertices(face)]
    faces.append(face_vertices)
    solver.add_constraint("Circle", face_vertices, circle_weight)

# Add gravity force
gravity_force = 0.1
for i in mesh.vertices():
    solver.add_vertex_force(0.0, 0.0, gravity_force, i)

###############################################################################################
# Initialize solver
###############################################################################################
solver.init()

###############################################################################################
# Setup visualization
###############################################################################################
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh, show_points=True)

# Color vertices
for i in mesh.vertices():
    color = Color.from_hex("#e74c3c") if i in corners else Color.from_hex("#3498db")
    mesh.vertex_attribute(i, 'color', color)

# Create circle visualizations
circles = []
circle_objs = []

for face in faces:
    face_points = [Point(*points[idx]) for idx in face]
    circle = Circle.from_points(face_points)
    circle_poly = circle.to_polygon(20)
    circle_mesh = Mesh.from_vertices_and_faces(*circle_poly.to_vertices_and_faces())
    circles.append(circle_mesh)
    circle_objs.append(viewer.scene.add(circle_mesh, show_faces=False, linecolor=Color.red()))

###############################################################################################
# Animation
###############################################################################################
max_iterations = 1000
current_iteration = 0
steps_per_frame = 5

@viewer.on(interval=1)
def update(frame):
    global current_iteration
    
    if current_iteration >= max_iterations:
        return
    
    for _ in range(steps_per_frame):
        if current_iteration >= max_iterations:
            break
        
        solver.solve(1)
        current_iteration += 1
    
    updated_points = solver.get_points()
    
    for i, vertex in enumerate(mesh.vertices()):
        if i < len(updated_points):
            mesh.vertex_attributes(vertex, 'xyz', updated_points[i])
    
    mesh_obj.update(update_data=True)
    
    # Update circles
    for idx, face in enumerate(faces):
        face_points = [Point(*updated_points[i]) for i in face]
        circle = Circle.from_points(face_points)
        circle_poly = circle.to_polygon(20)
        
        vertices, poly_faces = circle_poly.to_vertices_and_faces()
        circles[idx].clear()
        for i, coords in enumerate(vertices):
            circles[idx].add_vertex(i, x=coords[0], y=coords[1], z=coords[2])
            
        for f in poly_faces:
            circles[idx].add_face(f)

        circle_objs[idx].update(update_data=True)

viewer.show()

###############################################################################################
# Cleanup
###############################################################################################
solver.delete()
