import compas_shapeop
from compas_shapeop import solver

print("Successfully imported compas_shapeop")
print("Available modules:", dir(compas_shapeop))
print("Solver module contents:", dir(solver))

# Try to create a solver instance
try:
    s = solver.Solver()
    print("Successfully created a Solver instance")
except Exception as e:
    print(f"Error creating Solver: {e}")

print("Build and import successful!")
