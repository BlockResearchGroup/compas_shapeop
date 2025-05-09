********************************************************************************
Circularization
********************************************************************************

This example demonstrates how to use COMPAS ShapeOp to apply circle constraints to meshes. Because faces with vertices on a circle have inscribed circles, this constraint is often used alongside planarization constraints for architectural panelization.

.. figure:: /_images/circularization.gif
    :figclass: figure
    :class: figure-img img-fluid

.. literalinclude:: ../../examples/meshsolver_circularization.py
    :language: python