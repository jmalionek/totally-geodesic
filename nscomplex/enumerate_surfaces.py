"""
Enumerating closed surfaces in a 1-efficient ideal triangulation.
"""

import snappy, regina
from .regina_util import (to_int, is_normal, haken_sum,
                         extract_vector, num_octagons,
                         remove_vertex_links,
                         in_standard_coordinates,
                         has_annulus_of_quads_around_a_thin_edge)

class SurfacesByEuler(object):
    """
    A bag for holding (almost) normal surfaces, sorted by Euler
    characteristic.  The surfaces are not allowed to have components
    that are vertex linking tori. At this stage, only almost normal
    surfaces with octogons are allowed, not those with tubes.
    """
    def __init__(self, surface_list):
        self.surfaces = dict()
        for S in surface_list:
            self.add(S)
        self.initial_surfaces = surface_list

    def __getitem__(self, euler):
        return self.surfaces.get(euler, list())

    def add(self, surface):
        S = remove_vertex_links(surface)
        e = to_int(S.eulerChar())
        current = self[e]
        if not S.isEmpty() and not any(S == C for C in current):
            self.surfaces[e] = current + [S]
            return True
        return False

    def __len__(self):
        return sun(len(surfaces) for surfaces in self.surfaces.values())

    def saturate_to_euler_bound(self, euler):
        initial = sum(self.surfaces.values(), list())
        progress = True
        while progress:
            progress = False
            current = sum(self.surfaces.values(), list())
            for I in initial:
                for C in current:
                    if to_int(I.eulerChar() + C.eulerChar()) >= euler:
                       if I.locallyCompatible(C):
                           if is_normal(I) or is_normal(C):
                               # The geometric Haken sum could have
                               # vertex linking tori, but these will
                               # be eliminated by the "add" method,
                               # making the result equivalent to taking
                               # the haken_sum with coordinates='quad'
                               S = haken_sum(I, C, coordinates='standard')
                               if self.add(S):
                                   progress = True

    def genus(self, genus):
        ans = []
        for S in self[2 - 2*genus]:
            if S.isOrientable() and S.isConnected():
                ans.append(S)
        return ans

    def normal(self, orientable=False):
        ans = []
        for euler in sorted(self.surfaces.keys(), reverse=True):
            for S in self[euler]:
                if is_normal(S):
                    if not orientable or S.isOrientable():
                        ans.append(S)
        return ans

    def almost_normal(self, orientable=False):
        ans = []
        for euler in sorted(self.surfaces.keys(), reverse=True):
            for S in self[euler]:
                if not is_normal(S):
                    if not orientable or S.isOrientable():
                        ans.append(S)
        return ans

    def connected_normal(self, euler_bound):
        return [S for S in self.normal()
                if S.isConnected() and to_int(S.eulerChar()) >= euler_bound]

    def connected_almost_normal(self, euler_bound):
        return [S for S in self.almost_normal()
                if S.isConnected() and to_int(S.eulerChar()) >= euler_bound]

    def min_euler_of_orientable_connected(self):
        surfaces = self.normal(orientable=True)
        surfaces += self.almost_normal(orientable=True)
        surfaces = [S for S in surfaces if S.isConnected()]
        if len(surfaces) == 0:
            return 0
        else:
            return min(to_int(S.eulerChar()) for S in surfaces)

def regina_enumerate_surfaces(triangulation):
    assert triangulation.isOriented()
    surfaces = regina.NormalSurfaces(triangulation,
                               regina.NS_AN_QUAD_OCT_CLOSED,
                               regina.NS_FUNDAMENTAL | regina.NS_EMBEDDED_ONLY)
    return [in_standard_coordinates(surfaces.surface(i))
            for i in range(surfaces.size())]

def min_euler_vertex_surface(triangulation):
    assert triangulation.isOriented()
    surfaces = regina.NormalSurfaces(triangulation,
                                    regina.NS_QUAD_CLOSED)
    surfaces = [surfaces.surface(i) for i in range(surfaces.size())]
    if len(surfaces) == 0:
        return 0
    else:
        return min(to_int(S.eulerChar()) for S in surfaces)

def min_euler_vertex_surface_standard(triangulation):
    assert triangulation.isOriented()
    surfaces = regina.NormalSurfaces(triangulation,
                                    regina.NS_STANDARD)
    surfaces = [surfaces.surface(i) for i in range(surfaces.size())]
    if len(surfaces) == 0:
        return 0
    else:
        return min(to_int(S.eulerChar()) for S in surfaces)

def min_euler_vert_std_sans_obvious_compression(triangulation):
    assert triangulation.isOriented()
    surfaces = regina.NormalSurfaces(triangulation,
                                    regina.NS_STANDARD)
    surfaces = [surfaces.surface(i) for i in range(surfaces.size())]
    if len(surfaces) == 0:
        return 0
    else:
        return min(to_int(S.eulerChar())
                   for S in surfaces
                   if not has_annulus_of_quads_around_a_thin_edge(S))

def num_vert_std_sans_obvious_compression(triangulation):
    assert triangulation.isOriented()
    surfaces = regina.NormalSurfaces(triangulation,
                                    regina.NS_STANDARD)
    surfaces = [surfaces.surface(i) for i in range(surfaces.size())]
    return len([S for S in surfaces
                if not has_annulus_of_quads_around_a_thin_edge(S)])

def vertex_surfaces(triangulation):
    assert triangulation.isOriented()
    surfaces = regina.NormalSurfaces(triangulation,
                                    regina.NS_QUAD_CLOSED)
    return [in_standard_coordinates(surfaces.surface(i))
            for i in range(surfaces.size())]

def vertex_surfaces_standard(triangulation):
    assert triangulation.isOriented()
    surfaces = regina.NormalSurfaces(triangulation,
                                    regina.NS_STANDARD)
    return [in_standard_coordinates(surfaces.surface(i))
            for i in range(surfaces.size())]

def surfaces_to_euler_char(regina_triangulation, euler_bound=None):
    """
    Returns all normal and almost normal surfaces with

    euler_char >= euler_bound

    in a triangulation where every normal surface of euler
    characteristic 0 is the vertex link.

    When euler_bound is not given it defaults to the minimum euler
    characteristic of any fundamental normal or almost surface with
    octogon in the triangulation.
    """
    T = regina_triangulation
    n = T.countTetrahedra()
    assert T.countVertices() == 1
    surfaces = regina_enumerate_surfaces(T)
    surfaces = [S for S in surfaces if num_octagons(S) <= 1]
    if len(surfaces) == 0:
        return SurfacesByEuler([])

    if euler_bound == None:
        euler_bound = min(to_int(S.eulerChar()) for S in surfaces)

    surfaces = [S for S in surfaces if to_int(S.eulerChar()) >= euler_bound]
    surfaces = SurfacesByEuler(surfaces)
    assert len(surfaces[0]) == 0  # Down with vertex links!

    # combine the surfaces in all possible ways
    surfaces.saturate_to_euler_bound(euler_bound)
    return surfaces


def surfaces_of_genus(regina_triangulation, genus):
    """
    >>> M = snappy.Manifold('m372')
    >>> N = regina.Triangulation3(M._to_string())
    >>> len(surfaces_of_genus(N, 2))
    7
    >>> len(surfaces_of_genus(N, 3))
    0

    >>> M = snappy.Manifold('m004')
    >>> N = regina.Triangulation3(M._to_string())
    >>> len(surfaces_of_genus(N, 2))
    0
    """
    surfaces = surfaces_to_euler_char(regina_triangulation, 2 - 2*genus)
    return surfaces.genus(genus)

def surfaces_of_euler_char(regina_triangulation, euler_char):
    """
    Return all *orientable* normal and almost normal surfaces in the
    given triangulation of the specified euler_characteristic.  The
    surfaces need *not* be connected.

    >>> M = snappy.Manifold('v1234')
    >>> N = regina.Triangulation3(M._to_string())
    >>> len(surfaces_of_euler_char(N, -4))
    11
    >>> len(surfaces_of_euler_char(N, -5))
    0
    >>> sorted({S.isConnected() for S in surfaces_of_euler_char(N, -6)})
    [False, True]
    """
    surfaces = surfaces_to_euler_char(regina_triangulation, euler_char)
    return [S for S in surfaces[euler_char] if S.isOrientable()]

def surfaces_of_genus_normal(triangulation, genus):
    """
    sage: [surfaces_of_genus_normal('K6_32',g) for g in range(2,10)]
    [(4, 0), (1, 1), (0, 1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]
    """
    M = snappy.Manifold(triangulation)
    N = regina.Triangulation3(M._to_string())
    surfaces = surfaces_of_genus(N, genus)
    normal = [S for S in surfaces if is_normal(S)]
    return len(normal), len(surfaces)-len(normal)

def surfaces_of_euler_char_normal(triangulation, euler_char):
    """
    sage: surfaces_of_euler_char_normal('v1234', -4)
    (15, 2)
    """
    M = snappy.Manifold(triangulation)
    N = regina.Triangulation3(M._to_string())
    surfaces = surfaces_of_euler_char(N, euler_char)
    normal = [S for S in surfaces if is_normal(S)]
    return len(normal), len(surfaces)-len(normal)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
