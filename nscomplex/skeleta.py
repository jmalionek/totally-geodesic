import snappy
import regina
import networkx as nx
import collections
from . import regina_util

# In order to build the 1-skeleton of a Regina normal surface, we need
# to extract the information about the 2-skeleton of the
# triangulation of the ambient 3-manifold.

class Face(object):
    """
    An oriented triangle in a Regina triangulation
    """
    def __init__(self, index, oriented_edges):
        assert len(oriented_edges) == 3
        self.oriented_edges = oriented_edges
        self.edges = [e for e, o in oriented_edges]
        self.relative_orientations = [o for e, o in oriented_edges]
        
    def __getitem__(self, i):
        return self.oriented_edges[i]
    
class TwoSkeletonOfTriangulation(object):
    """
    >>> T = regina.Triangulation3('hLALPkcbbegffghhroxhsj')
    >>> X = TwoSkeletonOfTriangulation(T)
    >>> X.thick_edges()
    [0, 2, 3, 5]
    >>> X.homology()
    []
    """
    
    def __init__(self, regina_triangulation):
        self.triangulation = T = regina_triangulation
        triangles = T.triangles()
        self.num_triangles = triangles.size()
        self.num_edges = T.edges().size()
        self.faces = []
        for t, tri in enumerate(triangles):
            oriented_edges = []
            for i in range(3):
                edge = tri.edge(i)
                mapping = tri.edgeMapping(i)
                if (mapping[0], mapping[1]) in [(0, 1), (1, 2), (2, 0)]:
                    orientation = 1
                else:
                    orientation = -1
                oriented_edges.append((int(edge.index()), orientation))
            self.faces.append(Face(t, oriented_edges))

    def homology(self):
        """
        For a 1-vertex triangulation, return H_1 as a list of elementary
        divisors.  
        """
        import sage.all
        assert self.triangulation.vertices().size() == 1
        A = sage.all.matrix(sage.all.ZZ, self.num_edges, self.num_triangles)
        for i, face in enumerate(self.faces):
            for edge, sign in face:
                A[edge, i] += sign

        d = [x for x in A.elementary_divisors() if x != 1]
        return d

    def thick_edges(self):
        thick = set()
        for face in self.faces:
            edges = [e[0] for e in face]
            for e in edges:
                if edges.count(e) > 1:
                    thick.add(e)

        return sorted(thick)
    

def test_two_skeleton(num_samples=None):
    """
    A check that we are extracting the information on
    triangulation correctly by computing homology.

    >>> test_two_skeleton(20)
    """
    census = snappy.OrientableClosedCensus
    if num_samples is not None:
        sample = [census.random() for i in range(num_samples)]
    else:
        sample = census
    
    for M in sample:
        T = regina.Triangulation3(M.filled_triangulation()._to_string())
        T.intelligentSimplify()
        skel = TwoSkeletonOfTriangulation(T)
        assert M.homology().elementary_divisors() == skel.homology()

# Now on to building the 1-skeleton of the a given Regina normal
# surface.
    
def weights_to_flow(a, b, c):
    """
    Given integer weights along the three successive edges of a
    triangle that satisfy the triangular inequality, compute the
    associated pattern of arcs.
    """
    assert a <= b + c and b <= a + c and c <= a + b
    assert (a + b + c) % 2 == 0
    u = (a + b - c)//2
    v = (b + c - a)//2
    w = (a + c - b)//2
    assert u + w == a and u + v == b and v + w == c
    return {(0,1):u, (1,2):v, (2,0):w}

class ArcsInTriangle(object):
    """
    This class describes an embedded family of proper arcs in a
    triangle.

    The sides of the triangle are labelled 0, 1, 2 anticlockwise. The
    ends of the arcs are numbered with respect to the same
    orientation.  Arcs are stored in an oriented fashion so that their
    vertices are increasing.

    >>> A = ArcsInTriangle([5, 3, 4])
    >>> sorted(A.arcs)   # doctest: +NORMALIZE_WHITESPACE
    [((0, 3), (1, 1)), ((0, 4), (1, 0)), ((1, 2), (2, 0)), 
     ((2, 1), (0, 2)), ((2, 2), (0, 1)), ((2, 3), (0, 0))]
    >>> A.normal()
    True
    """
    def __init__(self, weights, arcs=None):
        self.weights = weights
        self.arc_vertices = []
        for side, weight in enumerate(weights):
            self.arc_vertices += [(side, v) for v in range(weight)]
        if arcs is None:  # Assume arcs are normal, i.e. no bigons.
            arcs = []
            for (i, j), flow in weights_to_flow(*weights).items():
                v_i = range(weights[i])[-flow:]
                v_j = reversed(range(weights[j])[:flow])
                arcs += [((i, x), (j, y)) for x, y in zip(v_i, v_j)]
        self.arcs = arcs
        # self._check_consistency()

    @staticmethod
    def _is_normalized_arc(arc_info):
        (i, x), (j, y) = arc_info
        if i == j:
            return x < y
        return (i, j) in [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)]

    @staticmethod
    def _normalized_arc(arc):
        v, w = arc
        if ArcsInTriangle._is_normalized_arc((v, w)):
            return (v, w)
        else:
            return (w, v)

    def _check_consistency(self):
        assert self.arc_vertices == [(s, v) for s in range(3)
                                     for v in range(self.weights[s])]
        assert len(self.arc_vertices) == 2*len(self.arcs)
        assert all(self._is_normalized_arc(arc) for arc in self.arcs)
        seen = set()
        for arc in self.arcs:
            seen.update(arc)
        assert seen == set(self.arc_vertices)
        
    def an_innermost_bigon(self):
        """
        >>> ArcsInTriangle([3, 4, 5]).an_innermost_bigon()
        >>> A = ArcsInTriangle([4, 0, 0], [((0,0), (0,3)), ((0, 1), (0, 2))])
        >>> A.normal()
        False
        >>> A.an_innermost_bigon()
        ((0, 1), (0, 2))
        """
        for arc in self.arcs:
            ((i, x), (j, y)) = arc
            if i == j and x + 1 == y:
                return arc
            
    def along_arc(self, vertex):
        """
        Returns the other end point of the arc starting at the given
        vertex.

        >>> A = ArcsInTriangle([5, 3, 4])
        >>> A.along_arc((0, 2))
        (2, 1)
        >>> A.along_arc((2, 1))
        (0, 2)
        """
        for arc in self.arcs:
            if vertex in arc:
                return [v for v in arc if v != vertex][0]

    def merge_vertices(self, v, w):
        """
        Join the arcs attached to the given vertices together and push
        into the interior of the triangle.

        >>> A = ArcsInTriangle([2, 2, 2])
        >>> A.merge_vertices((2, 1), (2, 0))
        >>> A.merge_vertices((0, 0), (0, 1))
        >>> A.arcs
        [((1, 0), (1, 1))]
        >>> A.merge_vertices((1, 0), (1, 1))
        >>> A.weights
        [0, 0, 0]
        """
        if v > w:
            v, w = w, v
        if not (v[0] == w[0] and v[1] + 1 == w[1]):
            raise ValueError('Vertices are not adjacent')
        a, b = self.along_arc(v), self.along_arc(w)

        new_arcs = [arc for arc in self.arcs if not v in arc and not w in arc]
        if a == w and b == v:
            # We will get a loop in the interior of the triangle which
            # will be promptly forgotten.
            pass
        elif a[0] == v[0] or b[0] == v[0]:
            # In this situation, merging the vertices would create a
            # "cusp".  This can't happen in our application to
            # normalizing surfaces, since "all folds are to one side".
            raise ValueError('Merging vertices would create a cusp')
        else:
            new_arcs.append(self._normalized_arc((a, b)))
            

        reindexed = {c:c for c in self.arc_vertices}
        i, x = v
        reindexed.pop(v)
        reindexed.pop(w)
        for y in range(x + 2, self.weights[i]):
            reindexed[(i, y)] = (i, y - 2)

        self.weights[i] += -2
        self.arc_vertices = sorted(reindexed.values())
        self.arcs = [(reindexed[c], reindexed[d]) for c, d in new_arcs]
        #self._check_consistency()

    def normal(self):
        return all(v[0] != w[0] for v, w in self.arcs)

    def normal_pairing(self):
        """
        Each pair of adjacent parallel normal arcs gives a normal isotopy
        of the corresponding segments of the sides of the triangle.

        >>> A = ArcsInTriangle([4, 3, 5])
        >>> A.normal_pairing()   # doctest: +NORMALIZE_WHITESPACE
        [(((0, 0), (0, 1)), ((2, 4), (2, 3))), 
         (((0, 1), (0, 2)), ((2, 3), (2, 2))),
         (((2, 0), (2, 1)), ((1, 2), (1, 1)))]
        """
        assert self.normal()
        ans = []
        for i, j in [(0, 2), (1, 0), (2, 1)]:
            starts = [v for v in self.arc_vertices
                    if v[0] == i and self.along_arc(v)[0] == j]
            starts.sort()
            assert [x for _,x in starts] == list(range(len(starts)))
            for k in range(len(starts) - 1):
                v, w = starts[k], starts[k+1]
                ans.append(((v, w), (self.along_arc(v), self.along_arc(w))))
        return sorted(ans)
            


class FaceWithArcs(object):
    """
    This is primarily a helper class for going between the "global"
    labeling of vertices via:
    
    (edge index in 3-manifold triangulation, vertex number along edge)

    and the "local" labeling of vertices with respect to the edges of
    the triangle.
    """
    def __init__(self, face, arcs):
        self.face, self.arcs = face, arcs

    def local_to_global(self, local_edge_and_vertex):
        local_edge, local_vertex = local_edge_and_vertex
        global_edge, orient = self.face.oriented_edges[local_edge]
        if orient == 1:
            return (global_edge, local_vertex)
        else:
            weight = self.arcs.weights[local_edge]
            return (global_edge, weight - local_vertex - 1)

    def global_to_local(self, global_edge_and_vertex, local_edge):
        global_edge, global_vertex = global_edge_and_vertex
        edge, orient = self.face.oriented_edges[local_edge]
        assert edge == global_edge
        if orient == 1:
            return (local_edge, global_vertex)
        else:
            weight = self.arcs.weights[local_edge]
            return (local_edge, weight - global_vertex - 1)
                   
    def global_arcs(self):
        to_global = self.local_to_global
        return [(to_global(v), to_global(w)) for v, w in self.arcs.arcs]

    def merge_vertices(self, v, w):
        # First, we must figure out which local vertices are to be
        # merged. As multiple sides of this face can be identified
        # with the same global edge, there can be several steps.
        face, arcs = self.face, self.arcs
        to_local = self.global_to_local
        e = v[0]
        for i in range(3):
            if face.edges[i] == e:
                a, b = to_local(v, i), to_local(w, i)
                arcs.merge_vertices(a, b)

    def an_innermost_bigon(self):
        ans = self.arcs.an_innermost_bigon()
        if ans is not None:
            v, w = ans
            return tuple(sorted((self.local_to_global(v), self.local_to_global(w))))


class OneSkeletonOfSurface(object):
    """
    The 1-skeleton of a normal or almost normal surface given by its
    edge weights. 

    >>> M = snappy.Manifold('K5_1')
    >>> T = regina.Triangulation3(M._to_string())
    >>> N0 = OneSkeletonOfSurface(T, (2, 4, 2, 0, 2))
    >>> len(N0.connected_components()[0])
    1
    >>> A0 =  OneSkeletonOfSurface(T, (2, 4, 2, 2, 2))
    >>> len(A0.connected_components()[0])
    1
    >>> N0A0 = OneSkeletonOfSurface(T, (4, 8, 4, 2, 4))
    >>> surfaces, edge_intersections = N0A0.connected_components()
    >>> sorted(surfaces) == [N0.edge_weights, A0.edge_weights]
    True
    >>> edge_intersections[1]
    [0, 1, 1, 0, 0, 1, 1, 0]
    """
    def __init__(self, tri3, edge_weights):
        self.triangulation = T = tri3
        self.twoskeleton_of_tri = X = TwoSkeletonOfTriangulation(T)
        self._setup_vertices(edge_weights)
        self.faces_with_arcs = faces_with_arcs = []
        for face in X.faces:
            arcs = ArcsInTriangle([edge_weights[e] for e in face.edges])
            faces_with_arcs.append(FaceWithArcs(face, arcs))            
        self._build_graph()

    def _setup_vertices(self, edge_weights):
        self.edge_weights = tuple(edge_weights)
        self.verts_by_edge = verts_by_edge = dict()
        for e, weight in enumerate(self.edge_weights):
            verts_by_edge[e] = [(e, i) for i in range(weight)]
        self.vertices = sum(verts_by_edge.values(), [])
        
    def _build_graph(self):
        self.graph = G = nx.MultiGraph()
        G.add_nodes_from(self.vertices)
        for face_arcs in self.faces_with_arcs:
            G.add_edges_from(face_arcs.global_arcs())

    def merge_vertices(self, v, w):
        """
        The basic opperation in the normalization process is to eliminate
        two adjacent vertices along an edge that part of an innermost
        bigon in some face.  In terms of the one skeleton, this is the
        same as adding a small tube to the surface and removing any
        circle intersections with the interiors of the faces.

        In this example, we create the one skeleton of a tubing of the
        vertex link along edge 2.

        >>> M = snappy.Manifold('K8_291')
        >>> T = regina.Triangulation3(M._to_string())
        >>> L = OneSkeletonOfSurface(T, 8*(2,))
        >>> L.merge_vertices((2,0), (2, 1))
        >>> L.edge_weights
        (2, 2, 0, 2, 2, 2, 2, 2)
        >>> face_arcs = L.faces_with_arcs[1]
        >>> face_arcs.face.edges
        [3, 2, 0]
        >>> face_arcs.arcs.arcs
        [((2, 1), (0, 0)), ((2, 0), (0, 1))]

        This surface is in fact normal.

        >>> N = OneSkeletonOfSurface(T, (2, 2, 0, 2, 2, 2, 2, 2))
        >>> L_arcs = [sorted(fa.arcs.arcs) for fa in L.faces_with_arcs]
        >>> N_arcs = [sorted(fa.arcs.arcs) for fa in N.faces_with_arcs]
        >>> L_arcs == N_arcs
        True
        """
        if v > w:
            v, w = w, v
        assert v[0] == w[0] and v[1] + 1 == w[1]
        edge_weights = list(self.edge_weights)
        edge_weights[v[0]] += -2
        self._setup_vertices(edge_weights)
        for face_arcs in self.faces_with_arcs:
            face_arcs.merge_vertices(v, w)

    def an_innermost_bigon(self):
        for face_arcs in self.faces_with_arcs:
            bigon = face_arcs.an_innermost_bigon()
            if bigon is not None:
                return bigon

    def normalize(self):
        """
        Do merges until the one skeleton is a normal family of arcs in the
        two skeleton of the ambient triangulation.  Here is an example
        where we tube the vertex link along an edge for the initial
        surface. 

        >>> M = snappy.Manifold('K8_291')
        >>> T = regina.Triangulation3(M._to_string())
        >>> L = OneSkeletonOfSurface(T, 8*(2,))
        >>> L.merge_vertices((1,0), (1, 1))
        >>> L.edge_weights
        (2, 0, 2, 2, 2, 2, 2, 2)
        >>> L.an_innermost_bigon()
        ((2, 0), (2, 1))
        >>> L.normalize()

        The result is a certain normal surface

        >>> N = OneSkeletonOfSurface(T, (2, 0, 0, 2, 2, 2, 2, 2))
        >>> L_arcs = [sorted(fa.arcs.arcs) for fa in L.faces_with_arcs]
        >>> N_arcs = [sorted(fa.arcs.arcs) for fa in N.faces_with_arcs]
        >>> L_arcs == N_arcs
        True

        Now for an example coming from an almost normal surface with
        an octagon. There are two different ways we can push the
        octagonal saddle and then normalize. The first one falls of
        the 1-skeleton completely.

        >>> weights = (0, 2, 2, 2, 4, 2, 2, 4)
        >>> N = OneSkeletonOfSurface(T, weights)
        >>> N.merge_vertices((1, 0), (1, 1))
        >>> N.normalize()
        >>> N.edge_weights
        (0, 0, 0, 0, 0, 0, 0, 0)

        The other normalization becomes a certain normal surface
        
        >>> N = OneSkeletonOfSurface(T, weights)
        >>> N.merge_vertices((4, 1), (4, 2))
        >>> N.normalize()
        >>> N.edge_weights
        (0, 2, 2, 2, 2, 2, 2, 4)
        """
        bigon = self.an_innermost_bigon()
        while bigon is not None:
            self.merge_vertices(*bigon)
            bigon = self.an_innermost_bigon()
        
    def connected_components(self):
        """
        Returns the connected components of self as two lists:

        a. The first gives the surfaces themselves, where each
        surface is given by its vector of edge weights.

        b. The second describes how the surfaces meet each oriented
        edge of the ambient triangulation T.  Specifically, there is
        one sublist for each edge of T gives the indices of the
        surfaces hit by the edge in order.
        """

        components = list(nx.connected_components(self.graph))
        vertex_to_surface_index = dict()
        for i, component in enumerate(components):
            for vertex in component:
                vertex_to_surface_index[vertex] = i

        n = len(self.edge_weights)
        edges = range(n)
        surface_edge_weights = [list() for component in components] 
        edge_intersections = []
        for e in edges:
            intersections =  [vertex_to_surface_index[v]
                              for v in self.verts_by_edge[e]]
            edge_intersections.append(intersections)
            for i in range(len(components)):
                surface_edge_weights[i].append(intersections.count(i))

        surfaces = [tuple(weights) for weights in surface_edge_weights]
        return (surfaces, edge_intersections)

    def normal_isotopy_classes_of_interstitial_arcs(self):
        """
        Intersecting the edges of the 3-manifold triangulation with our
        surface gives a family of "interstitial arcs".  Often, pairs
        of these arcs are properly normally isotopic by sliding them
        along a face.  When considering normal surfaces with tubes, we
        can restrict to only tubing along one representative per
        "normal isotopy class", hence this function.
        
        >>> M = snappy.Manifold('K5_1')
        >>> T = regina.Triangulation3(M._to_string())
        >>> L = OneSkeletonOfSurface(T, (2, 2, 2, 2, 2))
        >>> len(L.normal_isotopy_classes_of_interstitial_arcs())
        5
        >>> N0 = OneSkeletonOfSurface(T, (2, 4, 2, 0, 2))
        >>> len(N0.normal_isotopy_classes_of_interstitial_arcs())
        2
        """
        interstitial_arcs = []
        for vertices in self.verts_by_edge.values():
            for i in range(len(vertices) - 1):
                interstitial_arcs.append((vertices[i], vertices[i+1]))

        G = nx.Graph()
        G.add_nodes_from(interstitial_arcs)
        for face_arcs in self.faces_with_arcs:
            def to_global(pairing):
                a, b = [face_arcs.local_to_global(x) for x in pairing]
                if a > b:
                    a, b = b, a
                return a, b
            
            for arc0, arc1 in face_arcs.arcs.normal_pairing():
                G.add_edge(to_global(arc0), to_global(arc1))

        return sorted(sorted(C)[0] for C in nx.connected_components(G))
        


def edge_intersections(T, w_0, w_1):
    """
    Given weight vectors of two connected (almost) normal surfaces,
    return a list of the pattern of intersections of the two surfaces
    along the edges of the 3-manifold triangulation T.

    >>> M = snappy.Manifold('K5_1')
    >>> T = regina.Triangulation3(M._to_string())
    >>> ans = edge_intersections(T, (2, 4, 2, 0, 2), (2, 4, 2, 2, 2))
    >>> ans[1]
    [1, 0, 0, 1, 1, 0, 0, 1]
    >>> ans = edge_intersections(T, (2, 4, 2, 2, 2), (2, 4, 2, 0, 2))
    >>> ans[1]
    [0, 1, 1, 0, 0, 1, 1, 0]
    """
    w_sum = regina_util.add_weights(w_0, w_1)
    C = OneSkeletonOfSurface(T, w_sum)
    surfaces, edge_intersections = C.connected_components()
    if len(surfaces) != 2:
        raise ValueError('Original surfaces were not connected')
    if [w_0, w_1] != surfaces:
        assert [w_1, w_0] == surfaces
        edge_intersections = [[1 - x for x in xs]
                              for xs in edge_intersections]
    return edge_intersections


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

    >>> M = snappy.Manifold('K5_1')
    >>> T = regina.Triangulation3(M._to_string())
    >>> normal, almost = regina_util.fundamental_genus2_surfaces(T)
    >>> N0, N1, N2 = normal
    >>> A0 = almost[0]
    >>> vertex_location(N0, N1)
    1
    >>> vertex_location(N1, N0)
    -1
    >>> vertex_location(N0, A0)
    1
    """

    T = A.triangulation()
    assert T.vertices().size() == 1
    assert T.homology(2).isTrivial()
    assert A.isConnected() and B.isConnected() and A.disjoint(B)

    w_a = regina_util.edge_weights(A)
    w_b = regina_util.edge_weights(B)
    
    # To figure out the location of the vertex, we look for an edge of
    # T that meets both A and B.
    for edge in edge_intersections(T, w_a, w_b):
        if edge.count(0) > 0 and edge.count(1) > 0:
            if edge[0] == 0:
                assert edge[-1] == 0
                return -1
            else:
                assert edge[0] == edge[-1] == 1
                return 1

    # In this case, we can reach the unique vertex from each surface
    # along a segment of an edge disjoint from the other surface.
    return 0


if __name__ == '__main__':
    import doctest
    doctest.testmod()
