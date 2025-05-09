from typing import List
from typing import Optional
from typing import Tuple

import numpy as np

from compas.datastructures import Mesh
from compas.geometry import Transformation
from compas_shapeop import Solver


class MeshSolver(Solver):
    """Create a solver from a COMPAS mesh.

    This is a convenience method that allows you to initialize a solver
    directly from a COMPAS mesh.

    Attributes
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        A COMPAS mesh.
    points : numpy.ndarray
        Direct reference to the solver's points matrix in shape (n, 3).
    """

    def __init__(self, mesh: Mesh) -> None:
        super().__init__()
        self.mesh: Mesh = mesh
        self.is_initialized: bool = False
        # This calls the setter in Solver class
        self.points: np.ndarray = mesh.to_vertices_and_faces()[0]

    def solve(self, iterations: int = 10) -> None:
        """Solve the constraint problem.

        Parameters
        ----------
        iterations : int, optional
            Number of iterations to run. Default is 10.
        """

        # Intialize the solver, this is needed after you add points, contraints and forces
        if not self.is_initialized:
            self.init()
            self.is_initialized = True

        # Run the solver.
        self._solver.solve(iterations)

        # Update the mesh.
        for i, vertex in enumerate(self.mesh.vertices()):
            self.mesh.vertex_attributes(vertex, "xyz", self.points[i])

    @classmethod
    def from_grid(cls, dx: float, nx: int, dy: float, ny: int) -> "MeshSolver":
        """Create a MeshSolver from a rectangular grid mesh.

        Parameters
        ----------
        dx : float
            Size of a grid cell in the x direction.
        nx : int
            Number of grid cells in the x direction.
        dy : float
            Size of a grid cell in the y direction.
        ny : int
            Number of grid cells in the y direction.

        Returns
        -------
        MeshSolver
            A new mesh solver instance.
        """
        mesh: Mesh = Mesh.from_meshgrid(dx, nx, dy, ny).translated([dy * -0.5, dx * -0.5, 0])
        ms: "MeshSolver" = cls(mesh)
        return ms

    @classmethod
    def from_obj(cls, obj: str, transformation: Optional[Transformation] = None) -> "MeshSolver":
        """Create a MeshSolver from an OBJ file.

        Parameters
        ----------
        obj : str
            Path to the OBJ file.
        transformations : list of :class:`compas.geometry.Transformation`, optional
            Transformations to apply to the mesh after loading.

        Returns
        -------
        MeshSolver
            A new mesh solver instance.
        """
        mesh: Mesh = Mesh.from_obj(obj)
        if transformation:
            mesh.transform(transformation)
        ms: "MeshSolver" = cls(mesh)
        return ms

    def constrain_edge_lengths(
        self, min_range: float = 0.99, max_range: float = 1.01, shrink_factor: float = 0, exclude_edges: Optional[List[Tuple[int, int]]] = None, weight: float = 1e1
    ) -> List[int]:
        """Add edge length constraints to all edges of a COMPAS mesh.

        Parameters
        ----------
        min_range : float, optional
            Minimum allowed relative length. Default is 0.9.
        max_range : float, optional
            Maximum allowed relative length. Default is 1.1.
        shrink_factor : float, optional
            Target shrinking factor (default=0). The target length will be
            (1.0 - shrink_factor) times the original length.
        exclude_edges : list, optional
            List of edges to exclude from constraints. Default is None.
        weight : float, optional
            Weight of the constraints. Default is 1.0.

        Returns
        -------
        list
            IDs of all added constraints.
        """
        exclude_edges = exclude_edges or []
        constraint_ids: List[int] = []

        for u, v in self.mesh.edges():
            if (u, v) not in exclude_edges and (v, u) not in exclude_edges:
                if shrink_factor > 0:
                    cid: int = self.add_shrinking_edge_constraint(u, v, weight, shrink_factor)
                    constraint_ids.append(cid)
                else:
                    cid: int = self.add_edge_strain_constraint(u, v, weight, min_range, max_range)
                    constraint_ids.append(cid)

        return constraint_ids

    def constrain_face_diagonals(self, min_range: float = 0.99, max_range: float = 1.01, weight: float = 1e1) -> List[int]:
        """Add diagonal constraints to quads to prevent shearing.

        Parameters
        ----------
        min_range : float, optional
            Minimum allowed relative length. Default is 0.9.
        max_range : float, optional
            Maximum allowed relative length. Default is 1.1.
        weight : float, optional
            Weight of the constraints. Default is 1.0.

        Returns
        -------
        list[int]
            IDs of all added constraints.
        """
        constraint_ids: List[int] = []
        for fkey in self.mesh.faces():
            vertices: List[int] = self.mesh.face_vertices(fkey)
            if len(vertices) == 4:  # Only apply to quads
                cid: int = self.add_edge_strain_constraint(vertices[0], vertices[2], 0.2)
                constraint_ids.append(cid)
                cid = self.add_edge_strain_constraint(vertices[1], vertices[3], 0.2)
                constraint_ids.append(cid)
        return constraint_ids

    def constrain_face_planarity(self, type: str = "plane", weight: float = 1e5) -> List[int]:
        """Add plane constraints to all faces of a COMPAS mesh.

        Parameters
        ----------
        type : str, optional
            "plane" or "circle". Default is "plane".
        weight : float, optional
            Weight of the constraints. Default is 1e5.

        Returns
        -------
        list[int]
            IDs of all added constraints.
        """

        constraint_ids: List[int] = []
        faces: List[List[int]] = []
        for fkey in self.mesh.faces():
            face_vertices: List[int] = list(self.mesh.face_vertices(fkey))
            faces.append(face_vertices)
            if type == "plane":
                cid: int = self.add_plane_constraint(face_vertices, weight)
            elif type == "circle":
                cid: int = self.add_circle_constraint(face_vertices, weight)
            constraint_ids.append(cid)
        return constraint_ids

    def constrain_face_regularization(self, weight: float = 1e3) -> List[int]:
        """Add equalize face constraints to all faces of a COMPAS mesh.

        Parameters
        ----------
        weight : float, optional
            Weight of the constraints. Default is 1e5.

        Returns
        -------
        list[int]
            IDs of all added constraints.
        """
        constraint_ids: List[int] = []
        for fkey in self.mesh.faces():
            face_vertices: List[int] = self.mesh.face_vertices(fkey)
            if len(face_vertices) > 3:  # Only apply to non-triangular faces
                cid: int = self.add_regular_polygon_constraint(face_vertices, weight)
                constraint_ids.append(cid)
        return constraint_ids

    def constrain_triface_bending(self, weight: float = 1.0, min_range: float = 1.0, max_range: float = 1.0) -> List[int]:
        """Add bending constraints to all pairs of adjacent triangular faces in a COMPAS mesh.

        Parameters
        ----------
        weight : float, optional
            Weight of the constraints. Default is 1.0.
        min_range : float, optional
            Minimum bend factor. Default is 1.0.
        max_range : float, optional
            Maximum bend factor. Default is 1.0.

        Returns
        -------
        list
            IDs of all added constraints.
        """
        constraint_ids: List[int] = []
        for edge in self.mesh.edges():
            u: int = edge[0]
            v: int = edge[1]

            # connected faces
            faces: List[Optional[int]] = list(self.mesh.edge_faces(edge))
            if faces[0] is None or faces[1] is None:
                continue

            # Get vertices for both faces
            face1_vertices: List[int] = self.mesh.face_vertices(faces[0])
            face2_vertices: List[int] = self.mesh.face_vertices(faces[1])

            # Find vertices that are not part of the shared edge
            id2: Optional[int] = None
            for vertex in face1_vertices:
                if vertex != u and vertex != v:
                    id2 = vertex
                    break

            id3: Optional[int] = None
            for vertex in face2_vertices:
                if vertex != u and vertex != v:
                    id3 = vertex
                    break

            # Now we need to order them as [id2, id0, id1, id3]
            # where id0-id1is the shared edge
            if id2 is not None and id3 is not None:
                constraint_ids.append(self.add_bending_constraint([u, v, id2, id3], weight, min_range, max_range))

        return constraint_ids

    def fix_vertices(self, vertices: Optional[List[int]] = None, targets: Optional[List[List[float]]] = None, weight: float = 1e5) -> List[int]:
        """Add closeness constraints to vertices of a COMPAS mesh.

        Parameters
        ----------
        vertices : list[int], optional
            List of vertices to constrain. If None, all vertices will be constrained.
        targets : list[[float, float, float]], optional
            Target positions for the vertices. Must match length of vertices if provided.
        weight : float, optional
            Weight of the constraints. Default is 1e5.

        Returns
        -------
        list[int]
            IDs of all added constraints.
        """
        constraint_ids: List[int] = []

        if vertices is None:
            vertices = self.mesh.vertices()

        # Convert to list if targets are provided to ensure length matching
        vertices = list(vertices) if targets is not None else vertices

        for idx, v in enumerate(vertices):
            target: Optional[List[float]] = targets[idx] if isinstance(targets, list) and len(targets) == len(vertices) else None
            cid: Optional[int] = self.fix_vertex(v, target, weight)
            if cid is not None:
                constraint_ids.append(cid)

        return constraint_ids

    def fix_vertex(self, vertex: Optional[int] = None, target: Optional[List[float]] = None, weight: float = 1e5) -> Optional[int]:
        """Add closeness constraints to all vertices of a COMPAS mesh to the boundary edges.
        Valence of the vertices allows to select the corners.

        Parameters
        ----------
        vertex : int, optional
            Vertex to constrain. Default is None.
        target : [float, float, float], optional
            Target position for the vertex. Default is None.
        weight : float, optional
            Weight of the constraints. Default is 1e5.

        Returns
        -------
        int
            ID of the added constraint.
        """

        if isinstance(vertex, int):
            if target is not None:
                cid = self.add_closeness_constraint_with_position(vertex, weight, target[0], target[1], target[2])
                return cid
            else:
                cid = self.add_closeness_constraint(vertex, weight)
                return cid

        return None

    def add_gravity(self, fz: float = 0.01) -> int:
        """Add gravity to the mesh.

        Parameters
        ----------
        fz : float, optional
            Z component of the gravity force vector. Default is -0.001.

        Returns
        -------
        int
            ID of the added force.
        """
        return self._solver.add_gravity_force(0, 0, fz)

    def inflate(self, weight: float = 1e1) -> None:
        """Add inflation force to the mesh.

        Parameters
        ----------
        weight : float, optional
            Weight of the force. Default is 1e5.

        Returns
        -------
        int
            ID of the added force.
        """
        faces_flat: List[int] = []
        face_sizes: List[int] = []
        for face in self.mesh.faces():
            face_vertices: List[int] = self.mesh.face_vertices(face)
            face_size: int = len(face_vertices)
            face_sizes.append(face_size)
            faces_flat.extend(face_vertices)

        self.add_normal_force_with_faces(faces_flat, face_sizes, weight)
