********************************************************************************
compas_shapeop
********************************************************************************

.. currentmodule:: compas_shapeop.solvers

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshSolver

Functions and Methods
=====================

Mesh Integration
----------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshSolver.from_obj
    MeshSolver.from_grid

Constraints
-----------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshSolver.constrain_edge_lengths
    MeshSolver.constrain_face_diagonals
    MeshSolver.constrain_face_planarity
    MeshSolver.constrain_face_regularization
    MeshSolver.constrain_triface_bending
    MeshSolver.fix_vertex
    MeshSolver.fix_vertices

Forces
------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshSolver.add_gravity
    MeshSolver.inflate

Core Methods
------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshSolver.solve

.. toctree::
    :maxdepth: 1
