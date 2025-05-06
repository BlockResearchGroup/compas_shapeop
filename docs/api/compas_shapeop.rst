********************************************************************************
compas_shapeop
********************************************************************************

.. currentmodule:: compas_shapeop

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Solver

Functions and Methods
=====================

Mesh Integration
----------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Solver.from_mesh
    Solver.add_mesh_edge_strain_constraint
    Solver.add_mesh_vertex_force

Constraints
-----------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Solver.add_closeness_constraint
    Solver.add_closeness_constraint_with_position
    Solver.add_edge_strain_constraint
    Solver.add_shrinking_edge_constraint
    Solver.add_circle_constraint
    Solver.add_plane_constraint
    Solver.add_similarity_constraint
    Solver.add_regular_polygon_constraint
    Solver.add_bending_constraint
    Solver.add_shape_constraint

Forces
------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Solver.add_vertex_force
    Solver.add_normal_force_with_faces
    Solver.add_gravity_force

Core Methods
------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Solver.points
    Solver.init
    Solver.solve

.. toctree::
    :maxdepth: 1
