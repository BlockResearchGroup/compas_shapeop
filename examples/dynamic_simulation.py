from compas_shapeop import Solver
from compas.datastructures import Mesh
from compas.geometry import Point
from compas.colors import Color
from compas_viewer import Viewer

###############################################################################################
# Create mesh grid and prepare for solver
###############################################################################################
rows, cols = 14, 14
spacing = 1.0
gravity_force = 0.005

# Create a mesh grid directly, similar to other examples
mesh = Mesh.from_meshgrid(spacing * (cols - 1), cols-1, spacing * (rows - 1), rows-1)

# Get the mesh points for the solver
points = mesh.to_vertices_and_faces()[0]

###############################################################################################
# Initialize solver and set points
###############################################################################################
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

# Pin corners
corner_weight = 1e5
for idx in corners:
    solver.add_constraint("Closeness", [idx], corner_weight)

# Find center vertex (approximate)
center_vertex = None
center_x, center_y = (cols-1)/2, (rows-1)/2
min_distance = float('inf')

for v in mesh.vertices():
    pos = mesh.vertex_coordinates(v)
    dist = ((pos[0] - center_x * spacing)**2 + (pos[1] - center_y * spacing)**2)**0.5
    if dist < min_distance:
        min_distance = dist
        center_vertex = v

# Pin center point
solver.add_constraint("Closeness", [center_vertex], corner_weight)

# Add edge strain constraints
for u, v in mesh.edges():
    solver.add_constraint("EdgeStrain", [u, v], 1.0)

# Add gravity force to non-fixed vertices
fixed_indices = corners + [center_vertex]
for v in mesh.vertices():
    if v not in fixed_indices:
        solver.add_vertex_force(0.0, 0.0, gravity_force, v)

###############################################################################################
# Initialize solver
###############################################################################################
solver.init()

###############################################################################################
# Setup visualization
###############################################################################################
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)

# Highlight fixed points with spheres
for idx in fixed_indices:
    pos = mesh.vertex_coordinates(idx)
    point = Point(*pos)
    viewer.scene.add(point, size=0.2, color=Color.red())

###############################################################################################
# Animation
###############################################################################################
iteration = 0

@viewer.on(interval=1)  
def deform_mesh(frame):
    global iteration, mesh, solver
    
    if iteration >= solver.max_iterations:
        return
    
    solver.solve(1)
    updated_points = solver.get_points()
    
    for i, vertex in enumerate(mesh.vertices()):
        if i < len(updated_points):
            mesh.vertex_attributes(vertex, 'xyz', updated_points[i])
    
    mesh_obj.update(update_data=True)
    
    iteration += 1

viewer.show()