import pickle
import snappy, regina
import math
import time

import nscomplex_tg.regina_util
from normal_surfaces import *
from itertools import combinations
from nscomplex_tg import regina_util, surfaces, enumerate_surfaces
from sage.all import block_matrix, matrix, vector, CC, FreeGroup

def obviously_compressible(surface):
    """
    Given a Regina NormalSurface object determines whether the surface has an obvious compressing disc.
    """
    S = surfaces.NormalSurface(surface, 0)
    return S.has_obvious_compression()

def detect_totally_geodesic(manifold, filename=None):
    """
    Given a Snappy Manifold object, finds all totally geodesic surfaces and possibly totally geodesic surfaces (surfaces
    that reach Algorithm2 step g) then pickles relevant information.
    The name for the resulting pickle file can be given with the filename argument; by default the file will be saved as
    "info_(name of manifold provided by Snappy)".
    """
    M = manifold
    if 'degenerate' in M.solution_type():
        return

    euler_bd = -math.floor(M.volume()/ (4 * 0.29156 * math.pi))

    # use vertex surfaces to enumerate surfaces

    tik_total = time.perf_counter()

    vertex_surfaces = regina.NormalSurfaces(regina.Triangulation3(M), regina.NS_QUAD_CLOSED)  # algHints=regina.NS_VERTEX_DD not needed
    vertex_surfaces_ess = [S for S in vertex_surfaces if not obviously_compressible(S)]
    vertex_surfaces_vec = [nscomplex_tg.regina_util.extract_vector(S) for S in vertex_surfaces]

    tok = time.perf_counter()

    find_vertex_sfces = tok - tik_total

    tik = time.perf_counter()

    all_surfaces = []
    for e in range(-1, euler_bd - 1, -1):
        all_surfaces.extend(find_surfaces(vertex_surfaces_ess, e))

    incomp_surfaces = [S for S in all_surfaces if not obviously_compressible(S)]
    incomp_our_surfaces = [from_regina_normal_surface(S, M) for S in incomp_surfaces]
    incomp_vec = [S.get_vector() for S in incomp_our_surfaces]

    tok_total = time.perf_counter()

    enumerate_from_vertex_sfces = tok_total - tik

    # finding totally geodesic surfaces

    tik = time.perf_counter()

    tot_geo_surfaces = []
    potential_tot_geo_surfaces = []  # these are surfaces that are orientable and whose traces are all real,
                                     # need to apply steps g, h of Algorithm 2 using nscomplex to check for totally geodesic

    for surface in incomp_our_surfaces:
        orientable = surface.surface.isOrientable()
        all_real = is_Fuchsian(M, surface)
        if orientable and all_real:
            potential_tot_geo_surfaces.append(surface.get_vector())
        elif not orientable and all_real:
            tot_geo_surfaces.append(surface.get_vector())

    tok = time.perf_counter()

    tot_geo_time = tok - tik

    result = {'manifold': regina.Triangulation3(M).tightEncoding(),
              'manifold_name': M.name(),
              'runtime_vertex_surfaces': find_vertex_sfces,
              'runtime_enumerate_surfaces': enumerate_from_vertex_sfces,
              'runtime_tot_geo': tot_geo_time,
              'vertex_surfaces_vec': vertex_surfaces_vec,
              'all_surfaces_vec': incomp_vec,
              'tot_geo': tot_geo_surfaces,
              'potential_tot_geo': potential_tot_geo_surfaces}

    if filename == None:
        filename = f'info_{M.name()}'

    with open(filename, 'wb') as file:
        pickle.dump(result, file)


def is_Fuchsian(M, surface):
    """
    Determines whether the given normal_surfaces.NormalSurface object in the given Snappy manifold is Fuchsian.
    If the given surface is non-orientable, this checks whether its double is Fuchsian.
    """
    G = M.fundamental_group(simplify_presentation=False)
    orientable = surface.surface.isOrientable()
    all_real = True
    if orientable:
        gens = surface.simplified_generators(False)
    else:
        vec = surface.get_vector()
        double_vec = tuple(2 * x for x in vec)
        double_surface = vec_to_NormalSurface(double_vec, M)
        gens = double_surface.simplified_generators(False)
    gens_matrix = [Tietze_to_matrix(gen, G) for gen in gens]
    comb = list(combinations(list(range(len(gens))), 1)) \
           + list(combinations(list(range(len(gens))), 2)) \
           + list(combinations(list(range(len(gens))), 3))
    for c in comb:
        if not all_real:
            break
        gen = matrix.identity(CC, 2)
        for n in c:
            gen *= gens_matrix[n]
        if gen.trace().imag().abs() > 0.001:
            all_real = False
            break
    return all_real
