********************************************************************************
compas_shapeop
********************************************************************************

.. figure:: /_images/vertex_force_and_closeness_constraint.gif
    :figclass: figure
    :class: figure-img img-fluid

.. rst-class:: lead

COMPAS ShapeOp provides Python bindings for `ShapeOp <https://shapeop.org/>`_, a robust and efficient C++ implementation of a physics solver for geometry processing.
The binding is generated with `Nanobind <https://nanobind.readthedocs.io/>`_ and offers zero-copy memory sharing between Python and C++ for high-performance simulations.

Key features include:

* Mesh planarization and regularization
* Physically-based cloth and cable simulations
* Geometry optimization with constraints
* Zero-copy integration with NumPy and Eigen
* Seamless integration with COMPAS meshes
* Interactive mesh manipulation in real-time


Table of Contents
=================

.. toctree::
   :maxdepth: 3
   :titlesonly:

   Introduction <self>
   installation
   devguide
   tutorial
   examples
   api
   license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
