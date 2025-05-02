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
    Solver.add_mesh_closeness_constraint

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
    Solver.set_similarity_constraint_shape
    Solver.add_regular_polygon_constraint

Forces
------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Solver.add_vertex_force
    Solver.add_normal_force_with_faces

Core Methods
------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Solver.set_points
    Solver.get_points
    Solver.init
    Solver.solve

.. toctree::
    :maxdepth: 1
