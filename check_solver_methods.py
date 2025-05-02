"""Check available methods in the DynamicSolver class"""

# Try importing directly from the module
try:
    from compas_shapeop._shapeoplibrary import DynamicSolver
    print("Successfully imported DynamicSolver directly from _shapeoplibrary")
except ImportError:
    print("Could not import DynamicSolver directly")

# Try importing through the wrapper
try:
    from compas_shapeop.shapeop_optimized import ShapeOpSolver
    print("\nSuccessfully imported ShapeOpSolver wrapper")
    
    # Create an instance
    solver = ShapeOpSolver()
    
    # Print internal solver attributes
    print("\nShapeOpSolver._solver attributes:")
    print(dir(solver._solver))
except ImportError:
    print("Could not import ShapeOpSolver wrapper")
