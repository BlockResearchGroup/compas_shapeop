from compas_shapeop import Solver
from compas.datastructures import Mesh
from compas_viewer import Viewer

###############################################################################################
# Create mesh grid and prepare for solver
###############################################################################################
mesh = Mesh.from_meshgrid(10.0, 10, 10.0, 10)

###############################################################################################
# Initialize solver and set points
###############################################################################################
solver = Solver()
points = mesh.to_vertices_and_faces()[0]
solver.set_points(points)

###############################################################################################
# Corners
# First and fourth corners lifted
# Second and third corners fixed
# Add edge constraints with shrinking
###############################################################################################
corners = []
for v in mesh.vertices():
    if len(mesh.vertex_neighbors(v)) == 2:
        corners.append(v)

for c in [0, 3]:
    idx = corners[c]
    pos = list(mesh.vertex_point(idx))
    pos[2] = 5.0
    solver.add_closeness_constraint_with_position(idx, 1000.0, pos)


for c in [1, 2]:
    idx = corners[c]
    solver.add_constraint("Closeness", [idx], 1000.0)


for u, v in mesh.edges():
    solver.add_shrinking_edge_constraint([u, v], 10.0, 0.25)

###############################################################################################
# Initialize solver
###############################################################################################
solver.init()

###############################################################################################
# Setup viewer and animation
###############################################################################################
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh, show_edges=True)
iteration = 0


@viewer.on(interval=0)
def update(frame):
    global iteration
    if iteration >= solver.max_iterations:
        return
        
    solver.solve(1)
    updated_points = solver.get_points()
    
    for key, idx in enumerate(mesh.vertices()):
        if idx < len(updated_points):
            mesh.vertex_attributes(key, 'xyz', updated_points[idx])
    
    mesh_obj.update(update_data=True)
    iteration += 1

viewer.show()