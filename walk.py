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
import time
from normal_surfaces import *
from itertools import combinations
from nscomplex import regina_util, surfaces, enumerate_surfaces
from sage.all import block_matrix, matrix, vector, CC, FreeGroup

def obviously_compressible(surface):
    '''
    surface: regina normal surface
    '''
    S = surfaces.NormalSurface(surface, 0)
    return S.has_obvious_compression()

def detect_totally_geodesic(k=0):
    # index = int(os.environ['SLURM_ARRAY_TASK_ID'])
    index = k
    with open('htlinkexterior.txt', 'r') as data:
        i = 0
        while i < index:
            data.readline()
            i += 1
        name = data.readline()
    M = snappy.Manifold(name)
    genus_bd = math.floor(M.volume() / (4 * math.pi) + 1)

    print(M.num_tetrahedra())

    # use vertex surfaces to enumerate surfaces

    tik = time.perf_counter()

    vertex_surfaces = regina.NormalSurfaces(regina.Triangulation3(M), regina.NS_QUAD_CLOSED, algHints=regina.NS_VERTEX_DD)
    vertex_surfaces_ess = [S for S in vertex_surfaces if not obviously_compressible(S)]

    all_surfaces = []
    for g in range(2, genus_bd + 1):
        all_surfaces.extend(find_surfaces(vertex_surfaces_ess, g))

    tok = time.perf_counter()

    incomp_surfaces = [S for S in all_surfaces if not obviously_compressible(S)]
    incomp_our_surfaces = [from_regina_normal_surface(S, M) for S in incomp_surfaces]
    incomp_vec = [S.get_vector() for S in incomp_our_surfaces]

    find_surface_time = tok - tik

    print('time taken for vertex surfaces', find_surface_time)

    # use fundamental surfaces to enumerate surfaces

    tik = time.perf_counter()

    all_surfaces_from_fundamental_dict = enumerate_surfaces.surfaces_to_euler_char(regina.Triangulation3(M), 2 - 2* genus_bd)
    all_surfaces_from_fundamental = all_surfaces_from_fundamental_dict.connected_normal(2 - 2 * genus_bd)

    incomp_surfaces_from_fundamental = [S.removeOcts() for S in all_surfaces_from_fundamental if not obviously_compressible(S)]  # list of regina normal surfaces
    incomp_our_surfaces_from_fundamental = [from_regina_normal_surface(S, M) for S in incomp_surfaces_from_fundamental]
    incomp_vec_from_fundamental = [S.get_vector() for S in incomp_our_surfaces_from_fundamental]

    vec_set = set(incomp_vec_from_fundamental)
    incomp_surf_no_duplicate = []
    for vec in vec_set:
        incomp_surf_no_duplicate.append(vec_to_NormalSurface(vec, M))

    tok = time.perf_counter()
    find_fundamental_surface_time = tok - tik

    print('time taken for fundamental surfaces', find_fundamental_surface_time)

    print('num surfaces from vertex surfaces', len(incomp_vec))
    print('num surfaces from fundamental surfaces', len(incomp_vec_from_fundamental))
    print('duplicates?', len(vec_set.union(set(incomp_vec))))

    # finding totally geodesic surfaces (for now using vertex surface enumerations)
    #
    # tik = time.perf_counter()
    #
    # tot_geo_surfaces = []
    # for surface in incomp_our_surfaces:
    #     all_real = True
    #     gens = surface.fundamental_group_embedding()
    #     gens_matrix = [Tietze_to_matrix(gen, M) for gen in gens]
    #     comb = list(combinations(list(range(len(gens))), 1)) \
    #            + list(combinations(list(range(len(gens))), 2)) \
    #            + list(combinations(list(range(len(gens))), 3))
    #     for c in comb:
    #         if not all_real:
    #             break
    #         gen = matrix.identity(CC, 2)
    #         for n in c:
    #             gen *= gens_matrix[n]
    #         if gen.trace().imag().abs() > 0.001:
    #             all_real = False
    #             break
    #     if all_real:
    #         tot_geo_surfaces.append(surface.get_vector())
    #
    # tok = time.perf_counter()
    #
    # surface_fun_gp_time = tok - tik
    #
    # result = {'manifold': M.isometry_signature(),
    #           'runtime_surfaces': find_surface_time,
    #           'runtime_gp': surface_fun_gp_time,
    #           'all_surfaces': incomp_vec,
    #           'tot_geo': tot_geo_surfaces}
    # print(result)

    # directory = '/data/keeling/a/chaeryn2/totally_geodesic/'
    # filename = 'totally_geodesic_info_manifold%i' % index
    # with open(directory+filename, 'wb') as file:
    #     pickle.dump(result, file)

if __name__ == '__main__':
    for k in range(20, 30):
        print('at', k)
        detect_totally_geodesic(k)