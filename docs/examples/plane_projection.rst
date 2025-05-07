********************************************************************************
Plane Projection
********************************************************************************

This example demonstrates how to project mesh faces vertices onto a plane. Edges of the mesh use edge strain constraints of a) mesh edges and b) diagonals of each face to preserve the overall shape and reduce distortion while flattening the mesh.

.. figure:: /_images/plane_projection.gif
    :figclass: figure
    :class: figure-img img-fluid

.. literalinclude:: ../../examples/meshsolver_plane_projection.py
    :language: python
