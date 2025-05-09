********************************************************************************
Conda Environment
********************************************************************************

There are two ways to set up the development environment:

Using environment.yml (recommended)
-----------------------------------

.. code-block:: bash

    conda env create -f environment.yml
    conda activate compas_shapeop

Manual setup
------------

.. code-block:: bash

    conda create -n compas_shapeop -c conda-forge python=3.9 compas -y
    pip install -r requirements-dev.txt
    pip install --no-build-isolation -ve . -Ceditable.rebuild=true

Both methods will create and configure the same development environment.

Building Wheels Locally
-----------------------

To test package building locally before PyPI release on all recent Python versions:

.. code-block:: bash

    pip install cibuildwheel
    cibuildwheel --output-dir wheelhouse .