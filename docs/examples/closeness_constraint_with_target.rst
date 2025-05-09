********************************************************************************
Closeness Constraint with Target
********************************************************************************

This example demonstrates how to use COMPAS ShapeOp to create a tensioned cable net structure with control over specific vertex positions. The example is similar to "vertex_force_and_closeness_constraint" example, but with the addition of closeness constraints with target positions and shrinking edge constraints to create tension in the structure.

.. figure:: /_images/closeness_constraint_with_target.png
    :figclass: figure
    :class: figure-img img-fluid

.. literalinclude:: ../../examples/meshsolver_closeness_constraint_with_target.py
    :language: python
