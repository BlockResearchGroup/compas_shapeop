********************************************************************************
Vertex Force and Closeness Constraint
********************************************************************************

This example demonstrates a dynamic cloth simulation using vertex force, edge strain constraint, and closeness constraints.
First we create a grid mesh and initialize the solver from it. Then we add closeness constraints to pin the corners of the mesh in place. Next we add edge strain constraints to make the mesh more regular. Finally we add vertex force to make the mesh more dynamic. We run the solver for 100 iterations and update the mesh vertices in real-time. The simulation is interactive and can be stopped at any time.

.. figure:: /_images/vertex_force_and_closeness_constraint.gif
    :figclass: figure
    :class: figure-img img-fluid

.. literalinclude:: ../../examples/meshsolver_vertex_force_and_closeness_constraint.py
    :language: python
