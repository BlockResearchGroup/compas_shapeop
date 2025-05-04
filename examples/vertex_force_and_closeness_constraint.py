from compas.datastructures import Mesh
from compas_viewer import Viewer
from compas_shapeop.shapeop import Solver

# ==========================================================================
# Create mesh
# ==========================================================================

rows, cols = 14, 14
mesh = Mesh.from_meshgrid(nx=cols - 1, ny=rows - 1, dx=1.0 * (cols - 1), dy=1.0 * (rows - 1))
mesh.translate([-5, -5, 0])

# ==========================================================================
# Initialize solver directly from mesh
# ==========================================================================

solver = Solver.from_mesh(mesh)

# ==========================================================================
# Add constraints
# ==========================================================================

corner_vertices = [v for v in mesh.vertices() if len(mesh.vertex_neighbors(v)) == 2]
solver.add_closeness_constraints(corner_vertices, weight=1e5)
solver.add_mesh_edge_strain_constraint(mesh, weight=1.0, min_range=0.8, max_range=1.1)
solver.add_gravity_force(fz=0.001)

# ==========================================================================
# Initialize the solver
# ==========================================================================

solver.init(dynamic=True, masses=1.0, damping=1.0, timestep=7)
points_ref = solver.points

# ==========================================================================
# Visualization
# ==========================================================================
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)


@viewer.on(interval=1)
def deform_mesh(frame):
    # Run one solver iteration
    solver.solve(20)

    # Update mesh vertices directly from the points_ref array
    # No need to manually copy values - the shared memory is automatically updated
    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, "xyz", points_ref[i])

    # Update the viewer
    mesh_obj.update(update_data=True)


viewer.show()
