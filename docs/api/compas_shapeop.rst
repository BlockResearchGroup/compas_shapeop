********************************************************************************
compas_shapeop
********************************************************************************

.. currentmodule:: compas_shapeop

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    meshsolver.MeshSolver

Functions and Methods
=====================

Mesh Integration
----------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    meshsolver.MeshSolver.from_obj
    meshsolver.MeshSolver.from_grid

Constraints
-----------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    meshsolver.MeshSolver.constrain_edge_lengths
    meshsolver.MeshSolver.constrain_face_diagonals
    meshsolver.MeshSolver.constrain_face_planarity
    meshsolver.MeshSolver.constrain_face_regularization
    meshsolver.MeshSolver.constrain_triface_bending
    meshsolver.MeshSolver.fix_vertex
    meshsolver.MeshSolver.fix_vertices

Forces
------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    meshsolver.MeshSolver.add_gravity
    meshsolver.MeshSolver.inflate

Core Methods
------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    meshsolver.MeshSolver.solve

.. toctree::
    :maxdepth: 1
