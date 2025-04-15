// pch.h - Precompiled Header for ShapeOp
#pragma once

// We want to use ShapeOp as a compiled library, not header-only because it takes 30 seconds to compile...
#undef SHAPEOP_HEADER_ONLY

// Standard library includes
#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <cmath>
#include <unordered_set>

// Eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// ShapeOp includes - use relative paths as resolved by CMake
#include "Solver.h"
#include "Constraint.h"
#include "Force.h"
#include "Types.h"
#include "Common.h"
#include "API.h"
#include "normalforce.h"  // Custom class for inflation force

// Nanobind includes
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

// C++ includes
#include <vector>
#include <memory>
#include <iostream>
#include <stdexcept>

namespace nb = nanobind;
using namespace nb::literals;

// Include ShapeOp API
#include "shapeop_api.h"
