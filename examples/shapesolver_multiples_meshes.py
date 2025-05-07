from compas.datastructures import Mesh
from compas_shapeop.solvers import ShapeSolver
from compas_shapeop.constraints import MeshVerticesConstraint, MeshForcesConstraint, MeshEdgesConstraint
from compas_viewer import Viewer

mesh0 = Mesh.from_meshgrid(10, 8, 10, 8)
mesh0.translate([-10, -5, 0])

mesh1 = Mesh.from_meshgrid(10, 8, 10, 8)
mesh1.translate([0, -5, 0])

s = ShapeSolver()

mesh_edges_constraint0 = MeshEdgesConstraint(mesh0)
mesh_vertices_constraint0 = MeshVerticesConstraint(mesh0, mesh0.vertices_where({"vertex_degree": 2}))
mesh_forces_constraint0 = MeshForcesConstraint(mesh0, mesh0.vertices(), force_z=0.05)
s.add(mesh_edges_constraint0)
s.add(mesh_vertices_constraint0)
s.add(mesh_forces_constraint0)

mesh_edges_constraint1 = MeshEdgesConstraint(mesh1)
mesh_vertices_constraint1 = MeshVerticesConstraint(mesh1, mesh1.vertices_where({"vertex_degree": 2}))
mesh_forces_constraint1 = MeshForcesConstraint(mesh1, mesh1.vertices(), force_z=-0.1)
s.add(mesh_edges_constraint1)
s.add(mesh_vertices_constraint1)
s.add(mesh_forces_constraint1)


viewer = Viewer()
mesh_obj0 = viewer.scene.add(mesh_edges_constraint0.shape)
mesh_obj1 = viewer.scene.add(mesh_edges_constraint1.shape)

@viewer.on(interval=1)
def update(frame):
    s.solve(1)
    mesh_obj0.update(update_data=True)
    mesh_obj1.update(update_data=True)

viewer.show()
