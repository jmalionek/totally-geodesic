"""
>>> M = snappy.Manifold('K7_20')
>>> T = regina.Triangulation3(M._to_string())
>>> normal, almost, tubed = connected_surfaces_to_euler_char(T, -6)
>>> len(normal), len(almost), len(tubed)
(10, 3, 44)
>>> normal[0], almost[0], tubed[0]
(N0, A0, T0)
>>> normal = [S for S in normal if S.connected]
>>> surfaces = normal + almost + tubed
>>> all(S.connected for S in surfaces)
True
>>> Counter(S0.is_compatible(S1) for S0 in normal for S1 in almost)
Counter({False: 24, True: 6})


>>> surfaces_copy = sum(connected_surfaces_to_euler_char(T, -6), [])
>>> Counter(S0==S1 for S0 in surfaces for S1 in surfaces_copy)
Counter({False: 3192, True: 57})
>>> pairs = [(S0, S1) for S0 in surfaces for S1 in surfaces if S0 != S1]
>>> good_pairs = [(S0, S1) for (S0, S1) in pairs if S0.kind != 'tubed' or S1.kind != 'tubed']
>>> len(good_pairs)
1300
>>> dis = [(S0, S1) for S0, S1 in good_pairs if is_disjoint(S0, S1)]
>>> len(dis)
142
>>> [(S0, S1) for S0, S1 in dis if vertex_location(S0, S1) != -vertex_location(S1, S0)]
[]
>>> Counter(vertex_location(S0, S1) for S0, S1 in dis) == Counter({1:71, -1:71})
True

>>> A0, A1, A2 = almost
>>> A0.normalizations()
[(0, 0, 0, 0, 0, 0, 0), (2, 2, 2, 2, 2, 2, 0)]
>>> A1.normalizations()
[(0, 0, 0, 0, 0, 0, 0), (2, 2, 2, 2, 2, 2, 0)]
>>> A2.normalizations()
[(0, 0, 0, 0, 0, 0, 0), (2, 2, 2, 2, 2, 2, 2)]

>>> norm_tubed = [T.normalizations() for T in tubed]

Testing haken sum of surfaces in standard vs quad coors

>>> hsum = regina_util.haken_sum
>>> edge_wts = regina_util.edge_weights
>>> remove_links = regina_util.remove_vertex_links
>>> N0, N4 = normal[0].surface, normal[4].surface
>>> edge_wts(N0), edge_wts(N4)
((2, 2, 2, 2, 2, 2, 0), (0, 2, 2, 2, 2, 2, 2))
>>> sum_std = hsum(N0, N4, 'standard'); edge_wts(sum_std)
(2, 4, 4, 4, 4, 4, 2)
>>> sum_quad = hsum(N0, N4, 'quad'); edge_wts(sum_quad)
(0, 2, 2, 2, 2, 2, 0)
>>> remove_links(sum_std) == sum_quad
True
>>> (normal[0] + normal[4]).edge_weights == edge_wts(sum_quad)
True
"""

import snappy
import regina
from collections import Counter
from . import regina_util, skeleta, enumerate_surfaces

class Surface(object):
    def __repr__(self):
        return self.name


class NormalSurface(Surface):
    def __init__(self, regina_surface, index):
        self.surface = S = regina_surface
        self.triangulation = S.triangulation()
        self.index = index
        self.connected = S.isConnected()
        self.orientable = S.isOrientable()
        self.edge_weights = regina_util.edge_weights(S)
        self.quad_vector = regina_util.quad_vector(S)
        self.tri_vector = regina_util.tri_vector(S)
        self.oct_vector = regina_util.oct_vector(S)
        self.tri_quad_vector = regina_util.tri_quad_vector(S)
        self.quad_oct_vector = regina_util.extract_vector(S, coordinates='quad')
        self.full_vector = regina_util.extract_vector(S, coordinates='standard')
        self.euler = regina_util.to_int(S.eulerChar())
        self.total_weight = sum(self.edge_weights)
        if self.connected:
            self.genus = (2 - self.euler)//2
        self.kind = 'normal'
        self._format_name()
        self.is_fundamental = None  # means unknown
        self.is_vertex = None  # means unknown
        self.is_vertex_standard = None  # means unknown
        self._hash = hash((self.kind, tuple(self.full_vector)))
        
        # For display as a node in a graph
        self.fillcolor = '#efe1ff'  # light purple
        self.style = 'filled'
        self.shape = 'ellipse'

    def _format_name(self):
        index, kind = self.index, self.kind
        prefixes = {'normal':'N', 'octogon':'A', 'tubed':'T'}
        if isinstance(index, tuple):
            pieces = [prefixes[kind] + repr(i) for kind, i in index]
            self.name = '(' + ' + '.join(pieces) + ')'
        else:
            self.name = prefixes[kind] + repr(index)

    def _index_as_tuple(self):
        index = self.index
        if isinstance(index, tuple):
            return index
        else:
            return ((self.kind, self.index),)

    def __eq__(self, other):
        if not isinstance(other, Surface):
            return False
        if self.kind != other.kind:
            return False
        else:
            return self.edge_weights == other.edge_weights

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._hash

    def add(self, other, coordinates):
        kinds = {self.kind, other.kind}
        untubed = {'normal', 'octogon'}
        if not kinds.issubset(untubed):
            raise ValueError('Surfaces not of correct types')
        if kinds == {'octogon'}:
            raise ValueError('Sum would have two octogons')
        if not self.is_compatible(other):
            raise ValueError('Surfaces are not compatible')
        the_sum = regina_util.haken_sum(self.surface, other.surface, coordinates)
        index = self._index_as_tuple() + other._index_as_tuple()
        if 'octogon' in kinds:
            return AlmostNormalSurface(the_sum, index)
        return NormalSurface(the_sum, index)

    def __add__(self, other):
        return self.add(other, 'quad')

    def __rmul__(self, scalar):
        assert scalar >= 1
        if scalar == 1:
            return self
        return self + self.__rmul__(scalar - 1)

    def is_compatible(self, other):
        untubed = ['normal', 'octogon']
        if self.kind in untubed and other.kind in untubed:
            return self.surface.locallyCompatible(other.surface)
        raise ValueError('Surfaces not of correct types')

    def has_obvious_compression(self):
        return regina_util.has_annulus_of_quads_around_a_thin_edge(self.surface)

    def _determine_compressibility_via_regina(self):
        assert self.connected
        X = self.surface.cutAlong()
        X.idealToFinite()   # Important!          
        X.simplifyToLocalMinimum()
        assert len(X.components()) == 2
        X.splitIntoComponents()
        A, B = regina_util.children(X)
        if len(A.boundaryComponents()) != 1:
            A, B = B, A
        assert len(A.boundaryComponents()) == 1
        assert len(B.boundaryComponents()) == 2
        a_comp = A.hasCompressingDisc()
        b_comp = B.hasCompressingDisc()
        if not a_comp and not b_comp:
            self.regina_compressible = 'incompressible'
            #self.shape = 'rectangle'
        elif a_comp and not b_comp:
            self.regina_compressible = 'away from vertex'
            #self.shape = 'trapezium'
        elif not a_comp and b_comp:
            self.regina_compressible = 'towards vertex'
            #self.shape = 'invtrapezium'
        else:
            self.regina_compressible = 'both sides'
            #self.shape = 'hexagon'

def link_of_all_vertices(triangulation):
    S = regina_util.vertex_link(triangulation)
    L = NormalSurface(S, index=-1)
    L.__repr__ = lambda self:'N_oo'
    return L

class AlmostNormalSurface(NormalSurface):
    def __init__(self, regina_surface, index):
        NormalSurface.__init__(self, regina_surface, index)

        # For display as a node in a graph
        self.kind = 'octogon'
        self._format_name()
        self.shape = 'ellipse'
        self.fillcolor = 'white'

    def normalizations(self):
        """
        Return the two normalizations of self as edge weight vectors.
        When self is connected and orientable, the normalizations are
        returned in the order

        (neg tightening (away from cusp), pos tightening (toward cusp))
        """
        saddles = regina_util.saddles_of_almost_normal(self.surface)
        if self.connected and self.orientable:
            def away_from_cusp(saddle):
                (i, x), (j, y) = saddle
                assert x + 1 == y
                return x % 2 == 0
            
            if not away_from_cusp(saddles[0]):
                saddles.reverse()
            assert away_from_cusp(saddles[0]) and not away_from_cusp(saddles[1])

        ans = []
        for saddle in saddles:
            skel = skeleta.OneSkeletonOfSurface(self.triangulation, self.edge_weights)
            skel.merge_vertices(*saddle)
            skel.normalize()
            ans.append(tuple(skel.edge_weights))
        return ans


    def _determine_compressibility_via_regina(self):
        raise ValueError('Not implemented')

    def has_obvious_compression(self):
        raise ValueError('Not implemented')
    
class TubedSurface(NormalSurface):
    def __init__(self, normal_surface, surface_index, tube_ends):
        if isinstance(normal_surface, NormalSurface):
            normal_surface = normal_surface.surface
        NormalSurface.__init__(self, normal_surface, surface_index)
        self.euler += -2
        self.kind = 'tubed'
        self._format_name()

        # Make sure that tube spec is valid.
        v, w = tube_ends
        if v > w:
            v, w = w, v
        assert v[0] == w[0] and v[1] + 1 == w[1]
        assert w[1] < self.edge_weights[v[0]]
        self.tube_ends = (v, w)

        if not self.connected:
            skel = skeleta.OneSkeletonOfSurface(self.triangulation, self.edge_weights)
            pieces, edge_intersections = skel.connected_components()
            if len(pieces) == 2:
                intersections = edge_intersections[v[0]]
                if intersections[v[1]] != intersections[w[1]]:
                    self.connected = True

        if self.connected:
            self.genus = (2 - self.euler)//2

        self._hash = hash((self.kind, self.tube_ends, tuple(self.full_vector)))
        # For display as a node in a graph
        self.shape = 'diamond'
        self.style = 'filled'
        self.fillcolor = '#dfffd7'  # lime green

    def normalizations(self):
        """
        Return the two normalizations of self as edge weight vectors.
        When self is connected and orientable, the normalizations are
        returned in the order

        (neg tightening (away from cusp), pos tightening (toward cusp))
        """
        trivial = regina_util.edge_weights(self.surface)
        skel = skeleta.OneSkeletonOfSurface(self.triangulation, self.edge_weights)
        skel.merge_vertices(*self.tube_ends)
        skel.normalize()
        hard = tuple(skel.edge_weights)
        ans = [trivial, hard]
        if self.connected and self.orientable:
            if self.tube_ends[0][1] % 2 == 0:
                ans.reverse()
        return ans

    def __rmul__(self, scalar):
        raise ValueError('Not implemented')

    def __eq__(self, other):
        if not isinstance(other, Surface):
            return False
        if self.kind != other.kind:
            return False
        else:
            weights_match = self.edge_weights == other.edge_weights
            tubes_match = self.tube_ends == other.tube_ends
            return weights_match and tubes_match

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._hash
    
    def _determine_compressibility_via_regina(self):
        raise ValueError('Not implemented')

    def has_obvious_compression(self):
        raise ValueError('Not implemented')



def connected_surfaces_to_euler_char(regina_triangulation, euler_char):
    """
    All connected orientable (almost) normal surfaces in the given
    1-vertex 1-efficient triangulation with euler characteristic
    bounded below by "euler_char", with the exception of the vertex
    link.  Includes both the octagonal and tubed varieties of almost
    normal surfaces.

    When euler_char is not given it defaults to the minimum euler
    characteristic of any fundamental normal or almost surface with
    octagon in the triangulation.
    """
    T = regina_triangulation
    L = regina_util.vertex_link(T)
    wt = regina_util.edge_weights    
    raw_surfaces = enumerate_surfaces.surfaces_to_euler_char(T, euler_char)
    if euler_char is None:
        euler_char = raw_surfaces.min_euler_of_orientable_connected()
        euler_char_standard = enumerate_surfaces.min_euler_vertex_surface_standard(T)
        if euler_char_standard < euler_char:
            euler_char = euler_char_standard
            raw_surfaces = enumerate_surfaces.surfaces_to_euler_char(T, euler_char)
    raw_normal = raw_surfaces.normal(orientable=True)
    tubed = []

    edge_weights = {wt(N) for N in raw_normal}
    def new_surface(S):
        return wt(S) not in edge_weights
    extra_for_tubes = [L] + [regina_util.haken_sum(S, L, 'standard') for S in raw_normal]
    extra_for_tubes = [N for N in extra_for_tubes if new_surface(N)]
    for N in extra_for_tubes + raw_normal:
        N_euler = regina_util.to_int(N.eulerChar())
        if N_euler - 2 >= euler_char:
            skel = skeleta.OneSkeletonOfSurface(T, wt(N))
            for tube_ends in skel.normal_isotopy_classes_of_interstitial_arcs():
                S = TubedSurface(N, len(tubed), tube_ends)
                if S.connected:
                    tubed.append(S)

    raw_normal = [N for N in raw_normal if N.isConnected()]
    raw_almost = [A for A in raw_surfaces.almost_normal(orientable=True)
                  if A.isConnected()]
    normal = [NormalSurface(S, i) for i, S in enumerate(raw_normal)]
    almost = [AlmostNormalSurface(S, i) for i, S in enumerate(raw_almost)]

    # Mark the fundamental surfaces for later reference
    weights_of_fundamental = {wt(S) for S in raw_surfaces.initial_surfaces}
    for S in normal + almost:
        S.is_fundamental = S.edge_weights in weights_of_fundamental

    # Same for vertex normal surfaces:
    weights_of_vertex = {wt(S) for S in enumerate_surfaces.vertex_surfaces(T)}
    for S in normal:
        S.is_vertex = S.edge_weights in weights_of_vertex

    # Same for vertex normal surfaces in standard coordinates:
    weights_of_vertex_standard = {wt(S) for S in enumerate_surfaces.vertex_surfaces_standard(T)}
    for S in normal:
        S.is_vertex_standard = S.edge_weights in weights_of_vertex_standard
        
    return normal, almost, tubed


def standardize_surface_pair(A, B):
    """
    Return the pair of surfaces in the order (C, D) so that C has the
    minimum type in [normal, almost normal, tubed].
    """
    return sorted([A, B], key=lambda X:X.kind)

def is_disjoint(A, B):
    assert isinstance(A, Surface) and isinstance(B, Surface)
    A, B = standardize_surface_pair(A, B)
    if not A.surface.disjoint(B.surface):
        return False
    if not isinstance(B, TubedSurface):
        return True
    
    # Ok we have at least one tube to deal with.
    if not (A.connected and B.connected):
        raise ValueError('Sorry, disconnected case not fully implemented')
    if isinstance(A, TubedSurface):
        raise ValueError('Sorry, dual tube case not implemented')

    # Corner case: one surface is a tubed version of the other
    w_a, w_b = A.edge_weights, B.edge_weights
    if w_a == w_b:
        assert isinstance(A, NormalSurface)
        return True

    # Typical case when B is tubed but A is not.
    T = A.triangulation
    edge_inter = skeleta.edge_intersections(T, w_a, w_b)
    for s, S in enumerate([A, B]):
        if isinstance(S, TubedSurface):
            v, w = S.tube_ends
            S_pos = [i for i, x in enumerate(edge_inter[v[0]]) if x == s]
            if abs(S_pos[v[1]] - S_pos[w[1]]) > 1:
                return False

    return True
    

def vertex_location(A, B):
    """
    Given two separating connected disjoint (almost) normal surfaces
    in a one-vertex triangulation T, there are three possible
    locations for said vertex in X = T - (A U B) which has three
    connected components:

    a) In the component of X whose boundary consists entirely of A.

    b) In the component of X whose boundary contains both A and B.

    c) In the component of X whose boundary consists entirely of B.

    This function returns [-1, 0, 1] corresponding to the three
    options above.
    """
    assert isinstance(A, Surface) and isinstance(B, Surface)
    assert A.connected and B.connected
    assert A.orientable and B.orientable
    assert is_disjoint(A, B) and A != B
    T = A.triangulation
    assert T.vertices().size() == 1
    assert T.homology(2).isTrivial()
    w_a, w_b = A.edge_weights, B.edge_weights
    edge_intersections = skeleta.edge_intersections(T, w_a, w_b)

    # A tricky case is when you are comparing a tubed surface to its
    # underlying normal surface.
    if w_a == w_b:
        if isinstance(A, TubedSurface) and isinstance(B, NormalSurface):
            v, w = A.tube_ends
            # Check if the tube on inside common normal surface
            return 1 if v[1] % 2 == 0 else -1
        else:
            assert isinstance(A, NormalSurface) and isinstance(B, TubedSurface)
            v, w = B.tube_ends
            # Check if the tube on inside common normal surface
            return -1 if v[1] % 2 == 0 else 1

    # One has to think this through, but basically for parity reasons
    # it's OK that we're ignoring tubes given that the two almost
    # normal surfaces have unique positions relative to each other.
    for intersects in edge_intersections:
        if 0 in intersects and 1 in intersects:
            if intersects[0] != intersects[-1]:
                return 0
            s0 = intersects[0]
            if intersects.index(1 - s0) % 2 == 0:
                return 0
            return 1 if s0 == 1 else -1
    
    # In this case, we can reach the unique vertex from each surface
    # along a segment of an edge disjoint from the other surface.
    return 0
        
    


if __name__ == '__main__':
    import doctest
    doctest.testmod()
