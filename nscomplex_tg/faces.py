from sage.all import vector, QQ, Polyhedron, PowerSeriesRing
from . import regina_util, surfaces, enumerate_surfaces
import regina
import networkx as nx
import json

def zero_set_quad(surface):
    """
    The indices of the quads of a normal surface with weight 0.
    """
    vector = getattr(surface, 'quad_vector', surface)
    return frozenset(i for i, v in enumerate(vector) if v == 0)

def zero_set_standard(surface):
    """
    The indices of the quads and tri of a normal surface with weight 0.
    """
    vector = getattr(surface, 'tri_quad_vector', surface)
    return frozenset(i for i, v in enumerate(vector) if v == 0)

def zero_set_quad_oct(surface):
    """
    The indices of the quads and octs of an almost normal surface with weight 0.
    """
    vector = getattr(surface, 'quad_oct_vector', surface)
    return frozenset(i for i, v in enumerate(vector) if v == 0)

def is_admissible_quad(zero_set, num_tet):
    """
    A zero set is admissible when at least two quad weights are
    zero in every tetrahedron.
    """
    if len(zero_set) == 0:
        return False
    for i in range(0, num_tet):
        if len(zero_set.intersection({3*i, 3*i + 1, 3*i + 2})) < 2:
            return False
    return True

def is_admissible_standard(zero_set, num_tet):
    """
    A zero set is admissible when at least two quad weights are
    zero in every tetrahedron.
    """
    if len(zero_set) == 0:
        return False
    for i in range(0, num_tet):
        if len(zero_set.intersection({7*i + 4,  7*i + 5,  7*i + 6})) < 2:
            return False
    return True


def is_admissible_quad_oct(zero_set, num_tet):
    """
    A zero set is admissible when at least five quad and oct weights
    are zero in every tetrahedron, and all but oct weights is zero.

    Assumes we're looking at quad-oct coordinates (i.e. no triangles).
    """
    if len(zero_set) == 0:
        return False

    non_zero_octs = 0
    for i in range(0, num_tet):
        if len(zero_set.intersection({6*i + k for k in range(6)})) < 5:
            return False
        oct_pos = {6*i + k for k in range(3, 6)}
        non_zero_octs += len(oct_pos - zero_set)
        if non_zero_octs > 1:
            return False
    return True

def create_admissibility_test_excluding_obvious_compressions(triangulation):
    """
    In standard coordinates, we can use the strengthened admissibility
    criterion that there is no chain of quads forming an annulus
    around a thin edge.
    """
    quads_around_edges = regina_util.quad_types_around_thin_edges(triangulation)

    def is_admissible(zero_set, num_tet):
        if not is_admissible_standard(zero_set, num_tet):
            return False
        for quads_around in quads_around_edges:
            indices = {7*tet + 4 + quad for tet, quad in quads_around}
            if indices.isdisjoint(zero_set):
                return False
        return True

    return is_admissible

class AdmissibleFace(object):
    """
    An admissible face of the normal surface solution space.

    WARNING: Because we are only interested in faces that to not carry
    the vertex link, the associated cone is stored in quad coordinates
    rather than triangle-quad coordinates to keep the dimension down.
    """
    def __init__(self, dim, zero_set, vertex_surfaces):
        self.dim = dim  # dimension of the face, not the cone
        self.zero_set = zero_set
        self.F0 = F0 = vertex_surfaces[0]
        self.num_tets = F0.triangulation.countTetrahedra()
        self.vertex_surfaces = vertex_surfaces
        self.euler_one_vertices = []
        for F in self.vertex_surfaces:
            if F.euler < 0:
                vec = vector(QQ, F.quad_vector)/F.euler
            else:
                vec = None  # Vertex link
            self.euler_one_vertices.append(vec)

    def __repr__(self):
        return 'AFace(d=%d, %s)' % (self.dim, self.vertex_surfaces)

    def default_surface_in_interior(self):
        F = self.vertex_surfaces[0]
        for G in self.vertex_surfaces[1:]:
            F = F.add(G, 'standard')
        return F

    def contains_vertex_link(self):
        F = self.default_surface_in_interior()
        return regina_util.contains_vertex_link(F.surface)

    def euler_slice_polyhedron(self, euler):
        verts = [euler*v for v in self.euler_one_vertices]
        P = Polyhedron(vertices=verts, backend='normaliz', base_ring=QQ)
        return P

    def max_euler_surface_in_interior(self):
        euler = -2
        interior = []
        while len(interior) == 0:
            P = self.euler_slice_polyhedron(euler)
            interior = [p for p in P.integral_points() if P.relative_interior_contains(p)]
            euler += -2
        T = self.F0.triangulation
        S = regina_util.normal_surface_from_quads(T, interior[0])
        return surfaces.NormalSurface(S, -1)

    def cone(self):
        vertex_rays = [-2*v for v in self.euler_one_vertices]
        P = Polyhedron(rays=vertex_rays, backend='normaliz', base_ring=QQ)
        return P

    def ehrhart_series(self):
        """
        The generating function for b(-2n) for this face.
        """
        verts = [-2*v for v in self.euler_one_vertices]
        P = Polyhedron(vertices=verts, backend='normaliz', base_ring=QQ)
        return P.ehrhart_series(variable='x')

    def ehrhart_series_of_interior(self):
        """
        The generating function for b(-2n) for the interior of this face.
        """
        verts = [-2*v for v in self.euler_one_vertices]
        P = Polyhedron(vertices=verts, backend='normaliz', base_ring=QQ)
        # Using Ehrhart-Macdonald reciprocity, see e.g. Thm 4.6.26 of Stanley.
        p = P.ehrhart_series(variable='x')
        x = p.parent().gen()
        return ((-1)**(self.dim + 1))*p(1/x)

    def _quad_vectors_in_interior(self, euler):
        P = self.euler_slice_polyhedron(euler)
        interior = [p for p in P.integral_points()
                    if P.relative_interior_contains(p)]
        return interior

    def num_of_genus_in_interior(self, genus):
        euler = 2 - 2*genus
        T = self.F0.triangulation
        ans = []
        for quad_vec in self._quad_vectors_in_interior(euler):
            S = regina_util.normal_surface_from_quads(T, quad_vec)
            assert regina_util.to_int(S.eulerChar()) == euler
            if S.isConnected():
                ans.append(S)
        return len(ans)

    def surfaces_genus_in_interior(self, genus):
        euler = 2 - 2*genus
        T = self.F0.triangulation
        ans = []
        for quad_vec in self._quad_vectors_in_interior(euler):
            S = regina_util.normal_surface_from_quads(T, quad_vec)
            assert regina_util.to_int(S.eulerChar()) == euler
            if S.isConnected():
                ans.append(S)
        return ans

    def surfaces_euler_in_interior(self, euler):
        T = self.F0.triangulation
        ans = []
        for quad_vec in self._quad_vectors_in_interior(euler):
            S = regina_util.normal_surface_from_quads(T, quad_vec)
            assert regina_util.to_int(S.eulerChar()) == euler
            if S.isConnected():
                ans.append(S)
        return ans




def admissible_faces(surfaces,
                     coordinates='standard',
                     num_tet=None,
                     zero_set_fn=None,
                     admissible_face_fn=None,
                     return_type=AdmissibleFace):

    """
    Returns the admissible faces of the (almost) normal surface
    polytope in the specified coordinates. It is assumed that surfaces
    are all vertex surfaces in these coordinates, and that surfaces
    includes all such that satisfy the admissibility condition.

    Closely follows Algorithm 3.2 of Burton (2014):

    https://doi.org/10.1137/1.9781611973198.11

    >>> T = regina_util.as_regina('K10n10')
    >>> raw_surfaces = enumerate_surfaces.vertex_surfaces(T)
    >>> verts = [surfaces.NormalSurface(S, i) for i, S in enumerate(raw_surfaces)]
    >>> faces = admissible_faces(verts, 'quad')
    >>> sorted(face.dim for face in faces)
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
    >>> len(faces[-1].vertex_surfaces)
    2
    """
    if coordinates == 'quad':
        zero_set, admissible = zero_set_quad, is_admissible_quad
    elif coordinates == 'standard':
        zero_set, admissible = zero_set_standard, is_admissible_standard
    elif coordinates == 'quad_oct':
        zero_set, admissible = zero_set_quad_oct, is_admissible_quad_oct
    else:
        assert coordinates == 'custom'
        zero_set, admissible = zero_set_fn, admissible_face_fn

    if len(surfaces) == 0:
        return []
    if num_tet is None:
        num_tet = surfaces[0].triangulation.countTetrahedra()

    surfaces = [F for F in surfaces if admissible(zero_set(F), num_tet)]
    S = {1:{zero_set(F) for F in surfaces}}
    zeros_to_surface = {zero_set(F):F for F in surfaces}
    assert len(S[1]) == len(surfaces)

    maximal = list()
    k = 1
    while len(S[k]) > 0:
        k += 1
        S[k] = set()
        for z in S[k - 1]:
            success = False
            for v in S[1]:
                v_cap_z = v.intersection(z)
                if (not z.issubset(v)) and admissible(v_cap_z, num_tet):
                    S[k].add(v_cap_z)
                    success = True
            if not success:
                maximal.append(z)
        # Remove nonmaximal elements
        for z in S[k].copy():
            for w in S[k]:
                if z.issubset(w) and w != z:
                    S[k].discard(z)
                    break

    all_faces = []
    for k in sorted(S.keys()):
        for zeros in S[k]:
            verts = [v for v in S[1] if zeros.issubset(v)]
            surfaces = [zeros_to_surface[v] for v in verts]
            surfaces.sort(key=lambda x:x.index)
            face = return_type(k - 1, zeros, surfaces)
            face.maximal = zeros in maximal
            all_faces.append(face)

    return all_faces

class FacesComplex(object):
    def __init__(self, faces):
        self.faces = faces
        self.maximal = [C for C in faces if C.maximal]
        if len(faces) > 0:
            vertex_surfaces = set.union(*[set(face.vertex_surfaces) for face in self.maximal])
            self.vertex_surfaces = sorted(vertex_surfaces, key=lambda S:S.index)
            self.dim = max(face.dim for face in self.maximal)
            self.euler = sum([(-1)**face.dim for face in faces])
        else:
            self.vertex_surfaces = []
            self.dim = -1
            self.euler = 0

    def ehrhart_series(self):
        return sum(C.ehrhart_series_of_interior() for C in self.faces)

    def num_of_genus(self, genus, upper_bound=None):
        if not upper_bound is None:
            return [self.num_of_genus(g) for g in range(genus, upper_bound)]
        return sum(C.num_of_genus_in_interior(genus) for C in self.faces)

    def check_ehrhart_series(self, num_to_check):
        R = PowerSeriesRing(QQ, 'x', default_prec=num_to_check)
        f = R(self.ehrhart_series())
        for n in range(1, num_to_check + 1):
            euler = -2*n
            vecs = set()
            for C in self.faces:
                for vec in C._quad_vectors_in_interior(euler):
                    vecs.add(tuple(vec))
            assert len(vecs) == f[n]

            alt_vecs = set()
            for C in self.maximal:
                P = C.euler_slice_polyhedron(euler)
                for vec in P.integral_points():
                    alt_vecs.add(tuple(vec))
            assert vecs == alt_vecs


    def __len__(self):
        return len(self.faces)

    def __getitem__(self, i):
        return self.faces[i]

    def __repr__(self):
        info = (self.dim, len(self.maximal), len(self.vertex_surfaces), self.euler)
        return 'FacesComplex(dim=%d, num max=%d, num vert=%d, euler=%d)' % info

if __name__ == '__main__':
    import doctest
    doctest.testmod()
