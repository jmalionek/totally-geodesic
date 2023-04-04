#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --array=1-10
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SPATCH --output=/data/keeling/a/chaeryn2/totally_geodesic/htlinkexterior%j
#SPATCH --error=/data/keeling/a/chaeryn2/totally_geodesic_error/htlinkexterior%j

import os
import pickle

import snappy, regina
import math
import normal_surfaces
import time
from itertools import combinations
from nscomplex import regina_util, surfaces
from sage.all import block_matrix, matrix, vector, CC, FreeGroup
def has_annulus_of_quads_around_a_thin_edge(surface):
    S = surfaces.NormalSurface(surface)
    for quads_around in regina_util.quad_types_around_thin_edges(S.triangulation()):
        if all([regina_util.to_int(S.quads(tet, quad)) > 0 for tet, quad in quads_around]):
            return True
    return False

if __name__ == '__main__':
    # index = int(os.environ['SLURM_ARRAY_TASK_ID'])
    index = 10
    with open('htlinkexterior.txt', 'r') as data:
        i = 0
        while i < index:
            data.readline()
        name = data.readline()
    M = snappy.Manifold(name)
    genus_bd = math.floor(M.volume() / (4 * math.pi) + 1)

    tik = time.perf_counter()

    vertex_surfaces = regina.NormalSurfaces(regina.Triangulation3(M), regina.NS_QUAD_CLOSED, algHints=regina.NS_VERTEX_DD)
    all_surfaces = []
    for g in range(2, genus_bd + 1):
        all_surfaces.extend(normal_surfaces.find_surfaces(vertex_surfaces))

    tok = time.perf_counter()

    incomp_surfaces = [S for S in all_surfaces if not has_annulus_of_quads_around_a_thin_edge(S)]
    incomp_our_surfaces = [from_regina_normal_surface(S, M) for S in incomp_surfaces]
    incomp_vec = [S.get_vector() for S in incomp_our_surfaces]

    find_surface_time = tok - tik

    tik = time.perf_counter()

    tot_geo_surfaces = []
    for surface in incomp_our_surfaces:
        all_real = True
        gens = surface.fundamental_group_embedding()
        gens_matrix = [Tietze_to_matrix(gen, M) for gen in gens]
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
        if all_real:
            tot_geo_surfaces.append(surface.get_vector())

    tok = time.perf_counter()

    surface_fun_gp_time = tok - tik

    result = {'manifold': M.isometry_signature(),
              'runtime_surfaces': find_surface_time,
              'runtime_gp': surface_fun_gp_time,
              'all_surfaces': incomp_vec,
              'tot_geo': tot_geo_surfaces}
    print(result)

    # directory = '/data/keeling/a/chaeryn2/totally_geodesic/'
    # filename = 'totally_geodesic_info_manifold%i' % index
    # with open(directory+filename, 'wb') as file:
    #     pickle.dump(result, file)

