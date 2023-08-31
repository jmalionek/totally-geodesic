#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --array=0-199
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=/data/keeling/a/chaeryn2/SLURM_print/right_htlinkexterior%A_%a
#SBATCH --error=/data/keeling/a/chaeryn2/SLURM_error/right_htlinkexterior%A_%a

import os
import pickle
import snappy, regina
import math
import time
import multiprocessing
import gc
from normal_surfaces import *
from itertools import combinations
from nscomplex_tg import regina_util, surfaces, enumerate_surfaces
from sage.all import block_matrix, matrix, vector, CC, FreeGroup

def obviously_compressible(surface):
    '''
    surface: regina normal surface
    '''
    S = surfaces.NormalSurface(surface, 0)
    return S.has_obvious_compression()

def detect_totally_geodesic(manifold, index):
    M = manifold
    # genus_bd_wrong = math.floor(M.volume() / (4 * math.pi) + 1)
    genus_bd_right = math.floor(M.volume()/ (8 * 0.29156 * math.pi) + 1)

    # use vertex surfaces to enumerate surfaces

    tik_total = time.perf_counter()

    vertex_surfaces = regina.NormalSurfaces(regina.Triangulation3(M), regina.NS_QUAD_CLOSED)  # algHints=regina.NS_VERTEX_DD not needed
    vertex_surfaces_ess = [S for S in vertex_surfaces if not obviously_compressible(S)]

    tok = time.perf_counter()

    find_vertex_sfces = tok - tik_total

    tik = time.perf_counter()

    all_surfaces = []
    for g in range(2, genus_bd_right + 1):
        all_surfaces.extend(find_surfaces(vertex_surfaces_ess, g))

    tok_total = time.perf_counter()

    enumerate_from_vertex_sfces = tok_total - tik

    incomp_surfaces = [S for S in all_surfaces if not obviously_compressible(S)]
    incomp_our_surfaces = [from_regina_normal_surface(S, M) for S in incomp_surfaces]
    incomp_vec = [S.get_vector() for S in incomp_our_surfaces]

    find_surface_time = tok_total - tik_total

    # (might need som time in the future?)
    # use fundamental surfaces to enumerate surfaces
    #
    # tik = time.perf_counter()
    #
    # all_surfaces_from_fundamental_dict = enumerate_surfaces.surfaces_to_euler_char(regina.Triangulation3(M), 2 - 2* genus_bd)
    # all_surfaces_from_fundamental = all_surfaces_from_fundamental_dict.connected_normal(2 - 2 * genus_bd)
    #
    # incomp_surfaces_from_fundamental = [S.removeOcts() for S in all_surfaces_from_fundamental if not obviously_compressible(S)]  # list of regina normal surfaces
    # incomp_our_surfaces_from_fundamental = [from_regina_normal_surface(S, M) for S in incomp_surfaces_from_fundamental]
    # incomp_vec_from_fundamental = [S.get_vector() for S in incomp_our_surfaces_from_fundamental]
    #
    # vec_set = set(incomp_vec_from_fundamental)
    # incomp_surf_no_duplicate = []
    # for vec in vec_set:
    #     incomp_surf_no_duplicate.append(vec_to_NormalSurface(vec, M))
    #
    # tok = time.perf_counter()
    # find_fundamental_surface_time = tok - tik

    # finding totally geodesic surfaces (using vertex surface enumerations)

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

    result = {'manifold': regina.Triangulation3(M).tightEncoding(),
              'runtime_surfaces': find_surface_time,
              'runtime_gp': surface_fun_gp_time,
              'all_surfaces': incomp_vec,
              'tot_geo': tot_geo_surfaces,
              'runtime_vertex_surfaces': find_vertex_sfces,
              'runtime_enumerating_surfaces_from_vertex_surfaces': enumerate_from_vertex_sfces}

    directory = '/data/keeling/a/chaeryn2/tg_computation_outputs/'
    filename = 'totally_geodesic_info_link%i' % index
    with open(directory+filename, 'wb') as file:
        pickle.dump(result, file)

if __name__ == '__main__':
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])

    mfld_list = []
    with open(f'taking_too_long_manifolds_{task}.txt', 'r') as f:
        too_long_list = f.read().split()
        too_long_list = [int(num) for num in too_long_list]

    for i in range(task, 279649, 200):
        if 'totally_geodesic_info_link%i'%i not in os.listdir('/data/keeling/a/chaeryn2/tg_computation_outputs/') and i not in too_long_list:
            M = snappy.HTLinkExteriors(alternating=False)[7.2:][i]
            mfld_list.append([i, M, M.num_tetrahedra()])

    mfld_list_num_tet = sorted(mfld_list, key=lambda manifold:manifold[2])
    for manifold_info in mfld_list_num_tet:
        gc.collect()
        p = multiprocessing.Process(target=detect_totally_geodesic, args=(manifold_info[1], manifold_info[0]))
        p.start()
        p.join(5000)
        if p.is_alive():
            with open(f'taking_too_long_manifolds_{task}.txt', 'a') as file:
                file.write(str(manifold_info[0]) + "\n")
            p.terminate()
            continue