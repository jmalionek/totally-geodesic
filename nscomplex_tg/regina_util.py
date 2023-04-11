"""
>>> T = regina.Triangulation3('iLLALQcbcdegghhhtsfxjgpxt')
>>> normal, almost = fundamental_genus2_surfaces(T)
>>> N0, N1 = normal[:2]
>>> A0, A1 = almost[:2]
>>> edge_weights(N1)
(0, 2, 0, 2, 2, 2, 2, 2)
>>> num_octagons(N0), num_octagons(A1)
(0, 1)
>>> N0.locallyCompatible(A0), N0.locallyCompatible(A1)
(True, False)
>>> S = haken_sum(N0, A0); num_octagons(S)
1
"""

import snappy, regina

def as_regina(triangulation):
    """
    >>> as_regina('m004').isoSig()
    'cPcbbbiht'
    """
    if isinstance(triangulation, str):
        triangulation = snappy.Triangulation(triangulation)
    if hasattr(triangulation, '_to_string'):
        # SnapPy Triangulation
        regina_tri = regina.Triangulation3(triangulation._to_string())
    else:
        regina_tri = triangulation
    return regina_tri

def to_int(regina_integer):
    return int(regina_integer.stringValue())

def children(packet):
    child = packet.firstChild()
    while child:
        yield child
        child = child.nextSibling()

def edge_weights(regina_surface):
    S = regina_surface
    edges = S.triangulation().edges()
    return tuple(to_int(S.edgeWeight(i)) for i, e in enumerate(edges))

def num_octagons(surface):
    ans = 0
    T = surface.triangulation()
    for i in range(T.countTetrahedra()):
        for j in range(3):
            ans += to_int(surface.octs(i, j))
    return ans

def is_normal(surface):
    return surface.octPosition().type == -1

def add_weights(V, W):
    """
    Add two equal-length tuples together, coordinate by coordinate.
    """
    assert len(V) == len(W)
    return tuple(v + w for v, w in zip(V, W))

def is_link_of_thin_edge(weights):
    ans = weights.count(0) == 1 and weights.count(2) == len(weights) - 1
    edge = weights.index(0) if ans else None
    return ans, edge

def contains_vertex_link(surface):
    S = surface
    T = S.triangulation()
    assert T.countVertices() == 1
    return all(S.triangles(i, j) > 0
               for i in range(T.countTetrahedra()) for j in range(4))

def quad_types_around_thin_edges(triangulation):
    """
    >>> T = as_regina('s783')
    >>> quad_types_around_thin_edges(T)
    [[(0, 1), (3, 2), (5, 1), (1, 2)], [(2, 1), (5, 1), (3, 2), (4, 0)]]
    """
    ans = []
    for edge in triangulation.edges():
        # If the edge is not thing, we skip it.
        tets_around = {emb.tetrahedron().index() for emb in edge.embeddings()}
        if edge.degree() != len(tets_around):
            continue

        quads = []
        for edge_embedding in edge.embeddings():
            # I couldn't figure out the correct way to get at this info.
            tet, edge = edge_embedding.str().split()
            tet = int(tet)
            assert len(edge) == 4 and edge[0] == '(' and edge[3] == ')'
            a, b = int(edge[1]), int(edge[2])
            quad_type = int(regina.quadSeparating[a][b])
            quads.append((tet, quad_type))

        ans.append(quads)

    return ans


def has_annulus_of_quads_around_a_thin_edge(surface):
    S = surface
    for quads_around in quad_types_around_thin_edges(S.triangulation()):
        if all([to_int(S.quads(tet, quad)) > 0 for tet, quad in quads_around]):
            return True

    return False

regina_coor_names = {'standard':regina.NS_AN_STANDARD,
                     'quad':regina.NS_AN_QUAD_OCT_CLOSED}

def extract_vector(surface, coordinates='standard'):
    """
    Extract the raw vector of the (almost) normal surface in either Regina's
    NS_AN_STANDARD or NS_AN_QUAD_OCT_CLOSED coordinate system.
    """
    S = surface
    T = S.triangulation()
    n = T.countTetrahedra()
    ans = []
    for i in range(n):
        if coordinates == 'standard':
            for j in range(4):
                ans.append(S.triangles(i, j))
        for j in range(3):
            ans.append(S.quads(i, j))
        for j in range(3):
            ans.append(S.octs(i, j))

    return [to_int(a) for a in ans]

def in_standard_coordinates(surface):
    """
    Return a copy of surface as a regina normal surface in standard
    coordinates.
    """
    T = surface.triangulation()
    v = extract_vector(surface, 'standard')
    return regina.NormalSurface(T, regina.NS_AN_STANDARD, v)

def haken_sum(S1, S2, coordinates='standard'):
    """
    Returns the surface corresponding to the sum of the vectors
    encoding the inputs in the requested coordinate system.  In
    standard coordinates, this is the geometric opperation of Haken
    sum.  In quad coordinates, it is the slightly less geometric
    operation of Haken sum *followed by removal of all vertex linking
    components*.

    Regardless of the coordinates used for the sum, the returned
    surface is encoded in regina.NS_AN_STANDARD coordinates for
    consistency.
    """
    T = S1.triangulation()
    assert S1.locallyCompatible(S2)
    assert is_normal(S1) or is_normal(S2)
    v1 = extract_vector(S1, 'standard')
    v2 = extract_vector(S2, 'standard')
    sum_vec = [x1 + x2 for x1, x2 in zip(v1, v2)]
    A = regina.NormalSurface(T, regina.NS_AN_STANDARD, sum_vec)
    if coordinates=='quad':
        A = remove_vertex_links(A)
    assert S1.locallyCompatible(A) and S2.locallyCompatible(A)
    assert S1.eulerChar() + S2.eulerChar() == A.eulerChar()
    return A

def quad_vector(surface):
    S = surface
    T = S.triangulation()
    n = T.countTetrahedra()
    return [to_int(S.quads(i, j)) for i in range(n) for j in range(3)]

def tri_vector(surface):
    S = surface
    T = S.triangulation()
    n = T.countTetrahedra()
    return [to_int(S.triangles(i, j)) for i in range(n) for j in range(4)]

def oct_vector(surface):
    S = surface
    T = S.triangulation()
    n = T.countTetrahedra()
    return [to_int(S.octs(i, j)) for i in range(n) for j in range(3)]

def quad_oct_vector(surface):
    S = surface
    T = S.triangulation()
    n = T.countTetrahedra()
    ans = []
    for i in range(n):
        ans += [to_int(S.quads(i, j))  for j in range(3)] + [to_int(S.octs(i, j))  for j in range(3)]
    return ans

def tri_quad_vector(surface):
    S = surface
    T = S.triangulation()
    n = T.countTetrahedra()
    ans = []
    for i in range(n):
        ans += [to_int(S.triangles(i, j))  for j in range(4)]
        ans += [to_int(S.quads(i, j))  for j in range(3)]
    return ans

def vertex_link(triangulation):
    """
    >>> T = regina.Triangulation3('iLALPMcbcbfffghhxxjxhqxhw')
    >>> L = vertex_link(T)
    >>> L.isConnected()
    True
    """
    n = triangulation.countTetrahedra()
    vector = n*[1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    L = regina.NormalSurface(triangulation, regina.NS_AN_STANDARD, vector)
    assert L.eulerChar() == 0
    return L

def contains_vertex_link(surface):
    vector = [v for v in extract_vector(surface)]
    T = surface.triangulation()
    n = T.countTetrahedra()
    tri_indices = [10*i + j for j in range(4) for i in range(n)]
    min_tri_weight = min(vector[i] for i in tri_indices)
    return min_tri_weight > 0

def remove_vertex_links(surface):
    """
    Return a copy of the surface with all trivial vertex linking
    components removed.
    """
    vector = [v for v in extract_vector(surface)]
    T = surface.triangulation()
    n = T.countTetrahedra()
    tri_indices = [10*i + j for j in range(4) for i in range(n)]
    min_tri_weight = min(vector[i] for i in tri_indices)
    for i in tri_indices:
        vector[i] += -min_tri_weight
    return regina.NormalSurface(T, regina.NS_AN_STANDARD, vector)

def normal_surface_from_quads(triangulation, quad_vector):
    quad_vector = [int(v) for v in quad_vector]
    return regina.NormalSurface(triangulation, regina.NS_QUAD_CLOSED, quad_vector)

def normal_surface_from_quads_and_octs(triangulation, quad_oct_vector):
    quad_oct_vector = [int(v) for v in quad_oct_vector]
    return regina.NormalSurface(triangulation, regina.NS_AN_QUAD_OCT_CLOSED, quad_oct_vector)

def saddles_of_almost_normal(surface):
    """
    Given an almost normal surface with one octagon, returns the two
    edge locations where we can push the octagonal saddle through to
    start the normalization process.

    >>> M = snappy.Manifold('K3_1')
    >>> T = regina.Triangulation3(M._to_string())
    >>> vector = [0,1,1,2,0,0,1,0,0,0,0,1,2,1,0,1,0,0,0,0,0,0,1,1,0,0,0,1,0,0]
    >>> A = regina.NormalSurface(T, regina.NS_AN_STANDARD, vector)
    >>> edge_weights(A)
    (2, 2, 4)
    >>> saddles_of_almost_normal(A)
    [((1, 0), (1, 1)), ((2, 1), (2, 2))]
    """
    S = surface
    T = S.triangulation()
    assert num_octagons(S) == 1
    pos = S.octPosition()
    tet_index, oct_type = pos.tetIndex, pos.type
    tet = T.tetrahedron(tet_index)

    # As per engine/triangulation/facenumbering.h, edges in Regina are
    # ordered lexicographically: 01, 02, 03, 12, 13, 23
    #
    # Additionally, quad type i separates edge[i] from edge[-i] and
    # octagon type i separates intersects edge[i] and edge[-i] twice.

    edges_per_regina = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    saddles = []
    for i in [oct_type, 5-oct_type]:
        v, w = edges_per_regina[i]
        edge_index = int(tet.edge(i).index())
        mapping = tet.edgeMapping(i)
        if (mapping[0], mapping[1]) == (v, w):
            start_vert = v
        else:
            assert (mapping[0], mapping[1]) == (w, v)
            start_vert = w
        w = to_int(S.triangles(tet_index, start_vert))
        saddles.append(((edge_index, w), (edge_index, w+1)))

    return sorted(saddles)


def empty_tetrahedra(surface):
    quads = quad_vector(surface)
    octs = oct_vector(surface)
    n = len(quads)/3
    positions = [i for i in range(n) if quads[3*i:3*i+3]==[0,0,0]
                 and octs[3*i:3*i+3]==[0,0,0]]
    return positions

def in_threes(listlike):
    n = len(listlike)
    assert n % 3 == 0
    return [listlike[i:i + 3] for i in range(0, n, 3)]

def is_mildly_singular(surface):
    assert is_normal(surface)
    return all(0 in q for q in in_threes(quad_vector(surface)))

def is_non_singular(surface):
    assert is_normal(surface)
    return all(q.count(0) >= 2 for q in in_threes(quad_vector(surface)))

def is_small(triangulation):
    """
    Decides if an ideal triangulation is small.

    >>> T = as_regina('m004')
    >>> is_small(T)
    True

    >>> T = as_regina('m137')
    >>> is_small(T)            # doctest: +SKIP
    False
    """
    R = triangulation
    surfaces = regina.NormalSurfaces(R, regina.NS_QUAD_CLOSED)
    if surfaces is None:
        return True
    
    for i in range(surfaces.size()):
        S = surfaces.surface(i)
        if S.isVertexLinking() or S.isThinEdgeLink() != (None, None):
            continue
        X = S.cutAlong()
        X.idealToFinite()
        X.simplifyToLocalMinimum()
        if not X.hasCompressingDisc():
            return False
    return True

def least_genus_essential_surface(triangulation):
    """
    Returns a least genus essential surface in the given ideal
    triangulation, if any.

    >> T = regina.Triangulation3('eLPkbcdddlwnnv')
    >> S = least_genus_essential_surface(T)
    >> S.eulerChar()
    -2
    """
    R = triangulation
    regina_surfaces = regina.NormalSurfaces(R, regina.NS_STANDARD)
    surfaces = [regina_surfaces.surface(i) for i in range(regina_surfaces.size())]
    surfaces.sort(key=lambda S:-S.eulerChar())
    for S in surfaces:
        if S.isVertexLinking() or S.isThinEdgeLink() != (None, None):
            continue
        X = S.cutAlong()
        X.idealToFinite()
        X.simplifyToLocalMinimum()
        if not X.hasCompressingDisc():
            return S

def fundamental_genus2_surfaces(regina_triangulation):
    T = regina_triangulation
    surfaces = regina.NormalSurfaces(T,
                        regina.NS_AN_STANDARD,
                        regina.NS_FUNDAMENTAL | regina.NS_EMBEDDED_ONLY)
    normal, almost_normal = [], []
    for i in range(surfaces.size()):
        S = surfaces.surface(i)
        if S.isOrientable() and S.eulerChar() == -2:
            if is_normal(S):
                normal.append(S)
            elif num_octagons(S) == 1:
                almost_normal.append(S)

    return normal, almost_normal

if __name__ == '__main__':
    import doctest
    doctest.testmod()
