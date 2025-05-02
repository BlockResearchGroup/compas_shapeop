# COMPAS ShapeOp

COMPAS bindings for ShapeOp

## Installation

Stable releases can be installed from PyPI.

```bash
pip install compas_shapeop
```

To install the latest version for development, do:

```bash
git clone https://github.com/blockresearchgroup/compas_shapeop.git
cd compas_shapeop
pip install --no-build-isolation -ve .
```

## Documentation

For further "getting started" instructions, a tutorial, examples, and an API reference,
please check out the online documentation here: [COMPAS ShapeOp docs](https://blockresearchgroup.github.io/compas_shapeop)

## Issue Tracker

If you find a bug or if you have a problem with running the code, please file an issue on the [Issue Tracker](https://github.com/blockresearchgroup/compas_shapeop/issues).


## Installation with editable rebuild

```bash
rm -rf build/
/home/pv/anaconda3/envs/compas_shapeop/bin/pip install --no-build-isolation -ve . -Ceditable.rebuild=true --no-cache-dir
python examples/optimized_dynamic_simulation.py
```