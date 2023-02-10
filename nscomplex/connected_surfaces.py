import snappy
import regina
import sys, os
import networkx as nx
import collections
from .surfaces import (connected_surfaces_to_euler_char,
                       is_disjoint, vertex_location,
                       link_of_all_vertices)
from .enumerate_surfaces import (min_euler_vert_std_sans_obvious_compression,
                                 num_vert_std_sans_obvious_compression)
from . import regina_util, faces, skeleta

class MoreSurfacesNeedToFindFaces(Exception):
    pass

class ConnectedSurfaces(object):
    """
    Given an ideal triangulation of a 1-cusped hyperbolic 3-manifold M
    with H_2(M; F_2) = 0, applies the algorithm in Section 6 of [DGR] to
    enumerate all normal and almost normal surfaces to the specified
    lower bound on Euler characteristic and determine exactly which of
    the normal ones are essential.

    >>> CS = ConnectedSurfaces('t12071', -6)
    >>> len(CS.normal), len(CS.almost_normal), len(CS.tubed)
    (44, 28, 175)
    >>> len(CS.incompressible), len(CS.least_weight)
    (30, 30)
    >>> LW = CS.essential_faces_of_normal_polytope(); LW
    FacesComplex(dim=2, num max=6, num vert=10, euler=1)
    >>> LW.ehrhart_series()
    (-x^3 + x^2 - 10*x)/(x^3 - 3*x^2 + 3*x - 1)
    >>> LW.num_of_genus(2, 4)
    [10, 8]

    If you don't specify the euler_bound, it defaults to the minimum
    over vertex normal surfaces in standard coordinates without
    obvious compressions about thin edges.
    """
    def __init__(self, triangulation, euler_bound=None):
        regina_tri = regina_util.as_regina(triangulation)
        name = regina_tri.isoSig()

        # Needed to ensure that all closed surfaces are seperating and
        # orientable.
        cusps = regina_tri.countVertices()
        assert cusps == 1  # Needed by Regina for NS_QUAD_CLOSED coor
        assert regina_tri.homologyH2Z2() == cusps - 1

        isosig = regina_tri.isoSig()
        if name == isosig:
            self.name = isosig
        else:
            self.name = name + '_' + isosig
        self.triangulation = regina_tri

        threshold = min_euler_vert_std_sans_obvious_compression(regina_tri)
        if euler_bound is None:
            euler_bound = threshold
        else:
            assert euler_bound <= threshold
        self.euler_bound = euler_bound
        self._setup_surfaces()
        self._build_graphs()

    def _setup_surfaces(self):
        T, e = self.triangulation, self.euler_bound
        surfaces = connected_surfaces_to_euler_char(T, e)
        self.surfaces = sum(surfaces, [])
        self.normal, self.almost_normal, self.tubed = surfaces
        self.surfaces_by_name = {S.name:S for S in self.surfaces}
        self.normal_by_weights = {N.edge_weights:N for N in self.normal}
        self.almost_normal_by_weights = {N.edge_weights:N for N in self.almost_normal}
        n = T.countEdges()
        self.normal_by_weights[n*(0,)] = 'empty'
        self.normal_by_weights[n*(2,)] = 'link'
        if e is None:
            self.euler_bound = min([0] + [S.euler for S in self.surfaces])
        assert self.euler_bound % 2 == 0
        self.max_genus = (2 - self.euler_bound) // 2

    def _build_graphs(self):
        self.pre = G = nx.DiGraph()
        # For recording non-normal isotopies of normal surfaces:
        self.NG = NG = nx.Graph()
        NG.add_nodes_from(self.normal)
        normal_that_compress = []

        G.add_nodes_from(self.surfaces)
        for i, A in enumerate(self.surfaces):
            for B in self.surfaces[i + 1:]:
                if not (A.kind == B.kind == 'tubed') and is_disjoint(A, B):
                    location = vertex_location(A, B)
                    if location <= 0:
                        G.add_edge(B, A)
                    if location >= 0:
                        G.add_edge(A, B)

        self.simple_edges = simple = []
        self.non_simple_edges = non_simple = []
        edges = G.edges()
        for A, B in edges:
            if (B, A) in edges:
                non_simple.append((A,B))
            else:
                simple.append((A, B))

        # Now add edges and decorations according to the
        # normalizations of almost normal surfaces.
        by_weights = self.normal_by_weights
        for S in self.almost_normal + self.tubed:
            wt_a, wt_b = S.normalizations()
            if wt_a in by_weights:
                A = by_weights[wt_a]
                if A == 'empty':
                    compresses_away = True
                else:
                    assert A != 'link'
                    attrs = {'color':'#5a59d7', 'penwidth':3}
                    compresses_away = A.genus != S.genus
                    if not compresses_away:
                        attrs['dir'] = 'both'
                    G.add_edge(A, S, **attrs)

            else:
                compresses_away = True

            if wt_b in by_weights:
                B = by_weights[wt_b]
                if B == 'link':
                    compresses_toward = True
                else:
                    assert B != 'empty'
                    attrs = {'color':'#5a59d7', 'penwidth':3}
                    compresses_toward = B.genus != S.genus
                    if not compresses_toward:
                        attrs['dir'] = 'both'
                    G.add_edge(S, B, **attrs)
            else:
                compresses_toward = True

            shapes = {(False, False): 'rectangle',
                       (True, False): 'trapezium',
                       (False, True): 'invtrapezium',
                        (True, True): 'hexagon'}
            S.shape = shapes[compresses_away, compresses_toward]

            if not compresses_away and not compresses_toward:
                NG.add_edge(A, B, label=S.name)
            if compresses_away and not compresses_toward:
                normal_that_compress.append(B)
            if not compresses_away and compresses_toward:
                normal_that_compress.append(A)

        # Now finally work out which normal surfaces are compressible.
        isotopy_classes = list(nx.connected_components(NG))
        for iso_class in isotopy_classes:
            compresses = any(N in normal_that_compress for N in iso_class)
            for N in iso_class:
                N.compresses = compresses
                if not compresses:
                    N.shape = 'rectangle'
        self.incompressible = [N for N in self.normal if N.compresses==False]
        self.incomp_graph = NG.subgraph(self.incompressible)
        self._build_least_weight_graph()

    def _build_least_weight_graph(self):
        IG = self.incomp_graph
        all_lw_surfaces = set()
        for component in nx.connected_components(IG):
            min_wt = min(N.total_weight for N in component)
            lw_surfaces = {N for N in component if N.total_weight == min_wt}
            assert nx.is_connected(IG.subgraph(lw_surfaces))
            all_lw_surfaces.update(lw_surfaces)
        self.lw_graph = self.LG = IG.subgraph(all_lw_surfaces)
        self.least_weight = sorted(all_lw_surfaces, key=lambda S:S.index)

    def edge_weight_info(self, out=sys.stdout):
        for kind in ['normal', 'octogon', 'tubed']:
            for S in [S for S in self.surfaces if S.kind == kind]:
                out.write(S.name + ': ' + repr(S.edge_weights))
                if kind == 'tubed':
                    v, w = S.tube_ends
                    out.write(' %s >--< %s' % (v, w))
                out.write(' g=%d' % S.genus)
                if kind == 'normal':
                    if S.compresses:
                        out.write(' compressible')
                    else:
                        out.write(' incompressible')
                out.write('\n')
            if kind != 'tubed':
                out.write('\n')

    def save_info(self, root_dir='data_new', force=False, open_pdfs=False):
        import draw_graph
        root_dir = os.path.expanduser(root_dir)
        dir = root_dir + '/' + self.name + '/'
        if os.path.exists(dir):
            if force:
                os.system('rm -rf ' + dir)
            else:
                raise ValueError('The directory ' + dir + ' already exists')
        os.mkdir(dir)
        draw_graph.save_graph(self.pre, dir + 'pre.pdf', self.name + ' pregraph')
        if len(self.incompressible):
            draw_graph.save_graph(self.incomp_graph, dir + 'incomp.pdf',
                                  self.name + ' incompressible graph')
        file = open(dir + 'surfaces.txt', 'w')
        self.edge_weight_info(file)
        file.close()
        if open_pdfs: # macOS only
            os.system('open -a Preview ' + dir + '*.pdf')

    def check_essential_surfaces_against_regina(self):
        essential = self.incompressible
        if len(essential) == 0:
            return regina_util.is_small(self.triangulation)
        else:
            for N in self.normal:
                N._determine_compressibility_via_regina()
                if N.compresses and N.regina_compressible == 'incompressible':
                    return False
                if not N.compresses and N.regina_compressible != 'incompressible':
                    return False
        return True

    def num_incompressible_by_genus(self, concise=False):
        counts = collections.Counter()
        for component in nx.connected_components(self.incomp_graph):
            N = next(iter(component))
            counts.update([N.genus])
        genera = range(2, self.max_genus + 1)
        if concise:
            ans = [counts[g] for g in genera]
        else:
            ans = collections.OrderedDict([(g, counts[g]) for g in genera])
        return ans

    def lookup_surface(self, surface):
        if surface.kind == 'normal':
            d = self.normal_by_weights
        elif surface.kind == 'octogon':
            d = self.almost_normal_by_weights
        else:
            raise ValueError('Lookup for tubed surfaces not implemented')
        return d.get(surface.edge_weights, None)


    def no_weight_changing_isotopies(self):
        for surfaces in nx.connected_components(self.incomp_graph):
            if len({N.total_weight for N in surfaces}) > 1:
                return False
        return True

    def has_isotopy_of_incompressible(self):
        for surfaces in nx.connected_components(self.incomp_graph):
            if len(surfaces) > 1:
                return True
        return False

    def has_isotopy_of_least_weight(self):
        for surfaces in nx.connected_components(self.lw_graph):
            if len(surfaces) > 1:
                return True
        return False


    def _connected_components_hack(self, surface):
        skel = skeleta.OneSkeletonOfSurface(self.triangulation, surface.edge_weights)
        ans = set()
        for edge_weights in skel.connected_components()[0]:
            if not edge_weights in self.normal_by_weights:
                raise MoreSurfacesNeedToFindFaces
            N =  self.normal_by_weights[edge_weights]
            ans.add(N)
        return ans

    def essential_faces_of_normal_polytope(self, include_all_faces=False):
        T = self.triangulation
        lw_surfaces = set(self.least_weight)

        # We consider vertex surfaces with respect to the standard
        # coordinates that lack obvious compressions.
        #
        # Only change needed to handle more vertices would be to add
        # the vertex-linking surface for each vertex.

        assert T.countVertices() == 1
        L = link_of_all_vertices(T)
        vertex_surfaces = [L] + [S for S in self.normal
                                 if (S.is_vertex_standard and
                                     not S.has_obvious_compression())]
        assert len(vertex_surfaces) == num_vert_std_sans_obvious_compression(T)

        # Find all faces of P_T spanned by this class of vertex surfaces.

        is_admissible = faces.create_admissibility_test_excluding_obvious_compressions(T)
        all_faces = faces.admissible_faces(vertex_surfaces,
                                           zero_set_fn=faces.zero_set_standard,
                                           admissible_face_fn=is_admissible)

        # For each face whose vertices are all essential lw-surfaces,
        # check a surface in its interior to determine whether it is
        # an essential lw-face.

        essential_faces = []
        for C in all_faces:
            if lw_surfaces.issuperset(C.vertex_surfaces):
                F = C.default_surface_in_interior()
                try:
                    pieces = self._connected_components_hack(F)
                except MoreSurfacesNeedToFindFaces:
                    F = C.max_euler_surface_in_interior()
                    pieces = self._connected_components_hack(F)
                C.determinative_surfaces = sorted(pieces, key=lambda x:x.index)
                if lw_surfaces.issuperset(C.determinative_surfaces):
                    essential_faces.append(C)

        # Mark the maximal faces as such and return the answer.

        def is_maximal(C):
            for D in essential_faces:
                u, v = C.zero_set, D.zero_set
                if u != v and u.issuperset(v):
                    return False
            return True

        for C in essential_faces:
            C.maximal = is_maximal(C)

        if include_all_faces:
            return essential_faces, all_faces
        return faces.FacesComplex(essential_faces)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
