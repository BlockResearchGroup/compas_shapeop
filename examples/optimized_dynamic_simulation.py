from compas.datastructures import Mesh
from compas_viewer import Viewer
from compas_shapeop.shapeoplibrary import ShapeOpSolver

# ==========================================================================
# Create mesh
# ==========================================================================

rows, cols = 14, 14
mesh = Mesh.from_meshgrid(
    nx=cols-1, 
    ny=rows-1, 
    dx=1.0 * (cols - 1), 
    dy=1.0 * (rows - 1)
)

# ==========================================================================
# Initialize solver directly from mesh
# ==========================================================================

solver = ShapeOpSolver.from_mesh(mesh)

# ==========================================================================
# Add constraints
# ==========================================================================

corner_vertices = [v for v in mesh.vertices() if len(mesh.vertex_neighbors(v)) == 2]
solver.add_mesh_closeness_constraint(mesh, corner_vertices, weight=1e5)
solver.add_mesh_edge_strain_constraint(mesh, weight=1.0, min_range=0.8, max_range=1.1)
solver.add_mesh_vertex_force(mesh, 0, 0, 0.001)

# ==========================================================================
# Initialize the solver
# ==========================================================================

points_ref = solver.init()

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
        mesh.vertex_attributes(vertex, 'xyz', points_ref[i])
    
    # Update the viewer
    mesh_obj.update(update_data=True)

viewer.show()
