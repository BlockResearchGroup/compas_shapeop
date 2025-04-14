// pch.h - Precompiled Header for ShapeOp
#pragma once

// This is a key change - we want to use ShapeOp as a compiled library, not header-only
#undef SHAPEOP_HEADER_ONLY

// Standard library includes
#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <cmath>

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