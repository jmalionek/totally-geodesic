#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --array=0-49
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=/data/keeling/a/chaeryn2/SLURM_print/long_htlinkexterior%A_%a
#SBATCH --error=/data/keeling/a/chaeryn2/SLURM_error/long_htlinkexterior%A_%a

import os
import pickle
import snappy, regina
import math
import time
import multiprocessing
import gc

import nscomplex_tg.regina_util
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
    G = M.fundamental_group(simplify_presentation=False)
    for surface in incomp_our_surfaces:
        orientable = surface.surface.isOrientable()
        all_real = True
        if orientable:
            gens = surface.simplified_generators(False)
        else:
            vec = S.get_vector()
            double_vec = tuple(2*x for x in vec)
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

    directory = '/data/keeling/a/chaeryn2/computation_outputs/'
    filename = f'link_info_{i}_{M.name()}'
    with open(directory+filename, 'wb') as file:
        pickle.dump(result, file)


if __name__ == '__main__':
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])

    with open('taking_too_long_manifolds/taking_way_too_long_manifolds.txt', 'r') as f:
        way_too_long = f.read().split()
        way_too_long = [int(num) for num in way_too_long]

    long_mfld_list = []
    for i in range(200):
        with open(f'taking_too_long_manifolds/taking_too_long_manifolds_{i}.txt', 'r') as f:
            too_long_list = f.read().split()
            too_long_list = [int(num) for num in too_long_list]
            long_mfld_list.extend(too_long_list)

    wrong_saved = []
    for num in long_mfld_list:
        if num > 279648:
            wrong_saved.append(num)
    with open('taking_too_long_manifolds/wrong_saved.txt', 'a') as f:
        for num in wrong_saved:
            f.write('\n' + str(num) + '\n')

    long_mfld_list = [num for num in long_mfld_list if num < 279649]
    mfld_list = []

    for index in range(task, len(long_mfld_list), 50):
        i = long_mfld_list[index]
        if 'totally_geodesic_info_link%i'%i not in os.listdir('/data/keeling/a/chaeryn2/tg_computation_outputs/'):
            if i not in way_too_long:
                M = snappy.HTLinkExteriors(alternating=False)[7.2:][i]
                mfld_list.append([i, M, M.num_tetrahedra()])

    mfld_list_num_tet = sorted(mfld_list, key=lambda manifold:manifold[2])
    for manifold_info in mfld_list_num_tet:
        gc.collect()
        p = multiprocessing.Process(target=detect_totally_geodesic, args=(manifold_info[1], manifold_info[0]))
        p.start()
        p.join(10000)
        if p.is_alive():
            flag = True
            while flag:
                if 'dummy_file' not in os.listdir():
                    with open('dummy_file', 'w') as f:
                        pass
                    with open(f'taking_too_long_manifolds/taking_way_too_long_manifolds.txt', 'a') as file:
                        file.write("\n" + str(manifold_info[0]) + "\n")
                    os.remove('dummy_file')
                    p.terminate()
                    flag = False
                else:
                    continue
            continue