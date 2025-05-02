from compas.datastructures import Mesh
from compas_viewer import Viewer
from compas_shapeop.shapeoplibrary import ShapeOpSolver

# ==========================================================================
# Create mesh
# ==========================================================================

mesh = Mesh.from_meshgrid(10.0, 10, 10.0, 10)

# ==========================================================================
# Initialize solver from mesh
# ==========================================================================

solver = ShapeOpSolver.from_mesh(mesh)

# ==========================================================================
# Find corner vertices
# ==========================================================================

corners = []
for v in mesh.vertices():
    if len(mesh.vertex_neighbors(v)) == 2:
        corners.append(v)

# ==========================================================================
# Add constraints with target positions
# ==========================================================================

# Fix two corners at their original position
for c in [corners[1], corners[2]]:
    pos = list(mesh.vertex_point(c))
    solver.add_closeness_constraint_with_position([c], 1000.0, pos)

# Lift two corners by setting a target position with z=5.0
for c in [corners[0], corners[3]]:
    pos = list(mesh.vertex_point(c))
    pos[2] = 5.0  # Raise corners to z=5.0
    solver.add_closeness_constraint_with_position([c], 1000.0, pos)

# Add shrinking edge constraints for tension effect
for u, v in mesh.edges():
    solver.add_edge_strain_constraint([u, v], 10.0, 0.2, 0.3)  # Targets ~25% shrinkage

# ==========================================================================
# Initialize the solver
# ==========================================================================

points_ref = solver.init()

# ==========================================================================
# Visualization
# ==========================================================================

viewer = Viewer()
mesh_obj = viewer.scene.add(mesh, show_edges=True)
iteration = 0

@viewer.on(interval=0)
def update(frame):
    global iteration
    if iteration >= 100:
        return
        
    # Run solver iteration
    solver.solve(1)
    
    # Update mesh vertices from the solver's points
    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, 'xyz', points_ref[i])
    
    # Update the viewer
    mesh_obj.update(update_data=True)
    iteration += 1

viewer.show()
