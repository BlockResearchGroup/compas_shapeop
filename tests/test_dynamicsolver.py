import numpy as np
from compas.datastructures import Mesh
from compas_shapeop.shapeop import Solver


def test_solver_basic_functionality():
    """Test that the SolverWrapper can perform basic solving operations without C++ errors."""

    rows, cols = 14, 14
    mesh = Mesh.from_meshgrid(nx=cols - 1, ny=rows - 1, dx=1.0 * (cols - 1), dy=1.0 * (rows - 1))
    mesh.translate([-5, -5, 0])
    points = mesh.to_vertices_and_faces()[0]
    
    solver = Solver()
    solver.points = points
    corner_vertices = [v for v in mesh.vertices() if len(mesh.vertex_neighbors(v)) == 2]
    solver.add_closeness_constraints(corner_vertices, weight=1e5)
    for edge in mesh.edges():
        solver.add_edge_strain_constraint(edge[0], edge[1], weight=1.0, min_range=0.8, max_range=1.1)
    solver.add_gravity_force(0, 0, 0.001)
    solver.init()
    original_positions = solver.points.copy()
    solver.solve(20)

    # Check that points have changed (solver should update points_ref in-place)
    assert not np.array_equal(original_positions, solver.points), "Points should change after solving"

    # Specifically check that at least some z-coordinates have changed due to the z-force
    z_changes = np.abs(solver.points[2, :] - original_positions[2, :])
    assert np.max(z_changes) > 0, "At least some z-coordinates should change due to the force"
