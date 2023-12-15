#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --array=0-99
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=/data/keeling/a/chaeryn2/SLURM_print/left_OCCcover%A_%a
#SBATCH --error=/data/keeling/a/chaeryn2/SLURM_error/left_OCCcover%A_%a

import os, sys
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
            if 'simplify' in sys.argv:
                gens = surface.simplified_generators()
            else:
                gens = surface.fundamental_group_embedding()
        else:
            vec = S.get_vector()
            double_vec = tuple(2*x for x in vec)
            double_surface = vec_to_NormalSurface(double_vec, M)
            if 'simplify' in sys.argv:
                gens = double_surface.simplified_generators()
            else:
                gens = double_surface.fundamental_group_embedding()

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
    if 'simplify' in sys.argv:
        filename = f'cover_info_{index[0]}_{index[1]}_{index[2]}_{M.name()}_simplify'
    else:
        filename = f'cover_info_{index[0]}_{index[1]}_{index[2]}_{M.name()}'
    with open(directory+filename, 'wb') as file:
        pickle.dump(result, file)

if __name__ == '__main__':
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])

    census = snappy.OrientableCuspedCensus()

    for M_index, M in enumerate(census):
        if M_index % 100 == task:
            vol = M.volume()
            low_bd = math.ceil(3.6 / vol)
            up_bd = math.floor(20 / vol)
            indices = range(max(low_bd, 2), up_bd + 1)
            for index in indices:
                for cover_index, cover in enumerate(M.covers(index)):
                    if f'cover_info_{M_index}_{index}_{cover_index}_' not in os.listdir('/data/keeling/a/chaeryn2/computation_outputs/'):
                        gc.collect()
                        p = multiprocessing.Process(target=detect_totally_geodesic,
                                                    args=(cover, (M_index, index, cover_index)))
                        p.start()
                        p.join(5000)
                        if p.is_alive():
                            p.terminate()
                            continue