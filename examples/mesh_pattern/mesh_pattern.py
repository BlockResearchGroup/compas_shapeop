from tessagon.types.hex_tessagon import HexTessagon
from tessagon.adaptors.list_adaptor import ListAdaptor

import pathlib

from compas.datastructures import Mesh
from mesh_pattern import cppMesh


def create2DPattern():
    options = {
        "function": lambda u, v: [u, v, 0],
        "u_range": [0, 1.0],
        "v_range": [0, 1.0],
        "u_num": 16,
        "v_num": 10,
        "u_cyclic": False,
        "v_cyclic": False,
        "adaptor_class": ListAdaptor,
    }
    tessagon = HexTessagon(**options)
    vertices = tessagon.create_mesh()["vert_list"]
    faces = tessagon.create_mesh()["face_list"]
    return [vertices, faces]


def processMesh(vertices, faces, pv, pf):
    # load iglMesh
    mesh = MeshPattern.iglMesh()
    mesh.loadMesh(vertices, faces)

    # do parametrization
    mesh.parametrization_simple()
    vertices = mesh.getUVs()
    faces = mesh.getFaces()

    # do lift
    pattern2D = MeshPattern.iglMesh()
    pattern2D.loadMesh(pv, pf)
    mesh.mapMesh3D_AABB(pattern2D)
    vertices = pattern2D.getVertices()
    faces = pattern2D.getFaces()

    # do shapeop
    solver = MeshPattern.ShapeOpt()
    solver.loadMesh(pattern2D)
    solver.addPlanarConstraint(1.0)
    solver.addRegularConstraint(1.0)
    solver.addBoundaryConstraint(1.0)
    solver.runShapeOp(10)
    pattern3D = solver.getMesh()
    vertices = pattern3D.getVertices()
    faces = pattern3D.getFaces()

    # output wireframe
    wire = pattern3D.saveWireFrame(0.01, 10)
    vertices = wire.getVertices()
    faces = wire.getFaces()

    return [vertices, faces]


ROOT = pathlib.Path(__file__).parent.parent
DATA = ROOT / "data"
FILE = DATA / "minimal_surface.obj"

surface = Mesh.from_obj(FILE)

# getting the vertices and faces list from compas mesh
vertices, faces = surface.to_vertices_and_faces()

pv, pf = hex_tessagon.create2DPattern()
vertices, faces = cppMesh.processMesh(vertices, faces, pv, pf)

mesh = Mesh.from_vertices_and_faces(vertices, faces)