cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 11)
add_library(MeshPatternLib STATIC iglMesh.cpp ShapeOpt.cpp)
target_link_libraries(MeshPatternLib PRIVATE igl::core ShapeOp)

pybind11_add_module(MeshPattern MODULE MeshPattern.cpp)
target_link_libraries(MeshPattern PRIVATE MeshPatternLib)