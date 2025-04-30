import numpy as np
from compas_shapeop.shapeop_optimized import ShapeOpSolver

def test_basic_grid():
    """Test a simple grid with closeness and edge strain constraints"""
    # Create a simple 3x3 grid of points
    x = np.linspace(0, 2, 3)
    y = np.linspace(0, 2, 3)
    xx, yy = np.meshgrid(x, y)
    zz = np.zeros_like(xx)
    
    # Reshape into a list of 3D points
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    
    print("Points shape:", points.shape)
    print("Initial points F-contiguous:", np.isfortran(points))
    
    # Create solver
    solver = ShapeOpSolver()
    
    # Set points - this will convert to F-order
    solver.set_points(points)
    
    # Add constraints - fix the corners
    for i in [0, 2, 6, 8]:  # corner indices in a 3x3 grid
        solver.add_closeness_constraint([i], 1000.0)
    
    # Add edge constraints to maintain the structure
    edges = [
        [0, 1], [1, 2],  # top row
        [3, 4], [4, 5],  # middle row
        [6, 7], [7, 8],  # bottom row
        [0, 3], [3, 6],  # left column
        [1, 4], [4, 7],  # middle column
        [2, 5], [5, 8]   # right column
    ]
    
    for edge in edges:
        solver.add_edge_strain_constraint(edge, 1.0, 0.9, 1.1)
    
    # Initialize and solve
    solver.init()
    solver.solve(100)
    
    # Get results
    result_points = solver.get_points()
    print("Result points shape:", result_points.shape)
    print("Result points F-contiguous:", np.isfortran(result_points))
    
    # Check that corners are fixed
    corner_indices = [0, 2, 6, 8]
    for i in corner_indices:
        np.testing.assert_array_almost_equal(
            result_points[i], 
            points[i], 
            decimal=5, 
            err_msg=f"Corner point {i} moved"
        )
    
    print("All tests passed!")
    return result_points

def test_cloth_simulation():
    """Test a cloth-like grid with gravity and fixed top corners"""
    rows, cols = 10, 10
    x = np.linspace(0, 1, cols)
    y = np.linspace(0, 1, rows)
    xx, yy = np.meshgrid(x, y)
    zz = np.zeros_like(xx)
    
    # Reshape into a list of 3D points
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    
    # Create solver
    solver = ShapeOpSolver()
    
    # Set points
    solver.set_points(points)
    
    # Add constraints - fix the top corners
    top_left = 0
    top_right = cols - 1
    solver.add_closeness_constraint([top_left], 1000.0)
    solver.add_closeness_constraint([top_right], 1000.0)
    
    # Add edge constraints to maintain the cloth structure
    # Horizontal and vertical connections
    edges = []
    for i in range(rows):
        for j in range(cols):
            idx = i * cols + j
            
            # Right connection
            if j < cols - 1:
                edges.append([idx, idx + 1])
            
            # Down connection
            if i < rows - 1:
                edges.append([idx, idx + cols])
    
    # Add edge strain constraints with some elasticity
    for edge in edges:
        solver.add_edge_strain_constraint(edge, 1.0, 0.9, 1.1)
    
    # Initialize solver
    solver.init()
    
    # Solve with "gravity" by applying a downward force
    for _ in range(20):
        # Get current positions
        current_positions = solver.get_points()
        
        # Modify Z coordinates to simulate gravity (except for fixed points)
        for i in range(len(current_positions)):
            if i != top_left and i != top_right:
                current_positions[i, 2] -= 0.01
        
        # Update positions
        solver.set_points(current_positions)
        
        # Solve step
        solver.solve(5)
    
    # Get final result
    result_points = solver.get_points()
    
    # Check that top corners remained fixed
    np.testing.assert_array_almost_equal(
        result_points[top_left, :2], 
        points[top_left, :2], 
        decimal=5, 
        err_msg="Top left corner moved in X-Y plane"
    )
    
    np.testing.assert_array_almost_equal(
        result_points[top_right, :2], 
        points[top_right, :2], 
        decimal=5, 
        err_msg="Top right corner moved in X-Y plane"
    )
    
    # Check that points have moved down due to "gravity"
    assert np.any(result_points[:, 2] < -0.1), "Points didn't move down enough"
    
    print("Cloth simulation test passed!")
    return result_points
    
if __name__ == "__main__":
    print("\n=== Running basic grid test ===")
    grid_points = test_basic_grid()
    
    print("\n=== Running cloth simulation test ===")
    cloth_points = test_cloth_simulation()
    
    print("\nAll tests completed successfully!")
