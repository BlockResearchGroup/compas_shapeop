from compas.datastructures import Mesh
from compas_viewer import Viewer
from compas_shapeop.shapeop import Solver

# ==========================================================================
# Create triangulated mesh
# ==========================================================================

# First create a basic quad mesh
rows, cols = 14, 14
mesh = Mesh.from_meshgrid(nx=cols - 1, ny=rows - 1, dx=1.0 * (cols - 1), dy=1.0 * (rows - 1))
mesh.translate([-5, -5, 0])
mesh.quads_to_triangles()

# ==========================================================================
# Initialize solver directly from mesh
# ==========================================================================

solver = Solver.from_mesh(mesh)

# ==========================================================================
# Add constraints
# ==========================================================================

# Fix the corners in place
corner_vertices = [v for v in mesh.vertices() if len(mesh.vertex_neighbors(v)) < 4]
corner_vertices.pop()

solver.add_closeness_constraints(corner_vertices, weight=1e5)

# Add edge strain to maintain edge lengths
solver.add_mesh_edge_strain_constraint(mesh, weight=1.0, min_range=0.8, max_range=1.2)

# Add bending constraints to control the bending between adjacent triangle faces
solver.add_mesh_bending_constraints(mesh, weight=0.03, min_range=0.9, max_range=1.1)

# Add a gentle gravity force
solver.add_gravity_force(fz=0.0005)

# ==========================================================================
# Initialize the solver
# ==========================================================================

solver.init(dynamic=True, masses=1.0, damping=0.95, timestep=5)
points_ref = solver.points

# ==========================================================================
# Visualization
# ==========================================================================
viewer = Viewer()
mesh_obj = viewer.scene.add(mesh)


@viewer.on(interval=1)
def deform_mesh(frame):
    # Run solver iterations
    solver.solve(15)

    # Update mesh vertices directly from the points_ref array
    # No need to manually copy values - the shared memory is automatically updated
    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, "xyz", points_ref[i])

    # Update the viewer
    mesh_obj.update(update_data=True)


viewer.show()
