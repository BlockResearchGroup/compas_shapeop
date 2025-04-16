#include <vector>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "iglMesh.h"
#include "ShapeOpt.h"

namespace py = pybind11;
using namespace pybind11::literals;

int add(int i, int j){
    return i + j;
}

PYBIND11_MODULE(MeshPattern, m)
{
    m.doc() = "Mesh Pattern";
    m.def("add", &add, "A function which adds two numbers");
    m.def("createMesh", [](){
        return iglMesh();
    });
    
    py::class_<iglMesh>(m, "iglMesh")
            .def(py::init())
            .def("loadMesh", &iglMesh::loadMesh)
            .def("parametrization_simple", &iglMesh::parametrization_simple)
            .def("parametrization_lscm", &iglMesh::parametrization_lscm)
            .def("mapMesh3D_simple", &iglMesh::mapMesh3D_simple)
            .def("mapMesh3D_AABB", &iglMesh::mapMesh3D_AABB)
            .def("saveWireFrame", &iglMesh::saveWireFrame)
            .def("getUVs", &iglMesh::getUVs)
            .def("getVertices",  &iglMesh::getVertices)
            .def("getFaces",  &iglMesh::getFaces);


    py::class_<ShapeOpt>(m, "ShapeOpt")
            .def(py::init())
            .def("loadMesh", &ShapeOpt::loadMesh)
            .def("getMesh", &ShapeOpt::getMesh)
            .def("addPlanarConstraint", &ShapeOpt::addPlanarConstraint)
            .def("addRegularConstraint", &ShapeOpt::addRegularConstraint)
            .def("addBoundaryConstraint", &ShapeOpt::addBoundaryConstraint)
            .def("runShapeOp", &ShapeOpt::runShapeOp);

}
