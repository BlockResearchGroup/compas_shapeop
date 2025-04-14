import numpy as np
from compas_shapeop import _shapeop

# Grid size - smaller grid for faster execution
rows = 14
cols = 14
spacing = 1.0
gravity_force = 0.001

# Helper to get index from grid coordinates
def index(x, y):
    return y * cols + x

# Create grid points using NumPy
points_array = np.zeros((rows * cols, 3), dtype=np.float64)
for y in range(rows):
    for x in range(cols):
        i = index(x, y)
        normX = float(x) / (cols - 1)
        normY = float(y) / (rows - 1)
        points_array[i] = [
            normX * spacing * (cols - 1),
            -normY * spacing * (rows - 1),
            0.0
        ]

# Convert to a regular Python list (no VectorDouble needed now)
points = points_array.tolist()

# Create a solver
solver = _shapeop.ShapeOpSolver()

# Now the C++ binding directly handles Python lists
_shapeop.set_points(solver, points)

# Add constraints to pin the corners
corner_indices = [
    0,  # Top-left corner
    cols - 1,  # Top-right corner
    (rows - 1) * cols,  # Bottom-left corner
    rows * cols - 1  # Bottom-right corner
]
for idx in corner_indices:
    _shapeop.add_constraint(solver, "Closeness", [idx], weight=1e5)

# Add edge constraints for structural integrity
for y in range(rows):
    for x in range(cols):
        i = index(x, y)
        if x + 1 < cols:  # Horizontal edge
            _shapeop.add_constraint(solver, "EdgeStrain", [i, i + 1], weight=1.0)
        if y + 1 < rows:  # Vertical edge
            _shapeop.add_constraint(solver, "EdgeStrain", [i, i + cols], weight=1.0)

# Add a gravity force
_shapeop.add_gravity_force(solver, [0.0, 0.0, gravity_force])

# Initialize the solver
_shapeop.init_solver(solver)

# Run the solver for 1000 iterations
_shapeop.solve(solver, 1000)

# Get the final points
final_points = _shapeop.get_points(solver)

# Write the points to an OBJ file
with open("output.obj", "w") as f:
    # Write vertices
    for point in final_points:
        f.write(f"v {point[0]} {point[1]} {point[2]}\n")
    
    # Write faces (quads)
    for y in range(rows - 1):
        for x in range(cols - 1):
            i = index(x, y) + 1  # OBJ indices are 1-based
            j = index(x + 1, y) + 1
            k = index(x + 1, y + 1) + 1
            l = index(x, y + 1) + 1
            f.write(f"f {i} {j} {k} {l}\n")

print(f"Results written to output.obj")

# Clean up
_shapeop.delete_solver(solver)