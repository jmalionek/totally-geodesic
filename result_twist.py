#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=/data/keeling/a/chaeryn2/totally-geodesic/result_S21twist
#SBATCH --error=/data/keeling/a/chaeryn2/totally_geodesic_error/result_S21twist

import snappy, regina
import os
import pickle
import pandas as pd

def sort_results():
    manifold_data = pd.DataFrame({'index':pd.Series(dtype='int'),
                                  'manifold':pd.Series(dtype='str'),
                                  'runtime_surfaces_total':pd.Series(dtype='float'),
                                  'runtime_gp':pd.Series(dtype='float'),
                                  'num_tot_geo':pd.Series(dtype='int'),
                                  'runtime_vertex_surfaces':pd.Series(dtype='float'),
                                  'runtime_enumerating_surfaces_from_vertex_surfaces':pd.Series(dtype='float')}, index=[0])

    # manifold_data = pd.read_csv('result_S21twist.csv', delimiter=',')
    file_list = [name for name in os.listdir() if 'S21bundle' in name]
    for file in file_list:
    #    if i not in manifold_data['index']:
    #        if file in os.listdir('/data/keeling/a/chaeryn2/totally_geodesic/'):
        with open('/data/keeling/a/chaeryn2/totally_geodesic/' + file, 'rb') as f:
            info = pickle.load(f)
        if 'runtime_vertex_surfaces' in info.keys():
            manifold_data = pd.concat([manifold_data, pd.DataFrame({'index': i,
                                      'manifold': info['manifold'],
                                      'runtime_surfaces_total': info['runtime_surfaces'],
                                      'runtime_gp': info['runtime_gp'],
                                      'num_tot_geo': len(info['tot_geo']),
                                      'runtime_vertex_surfaces': info['runtime_vertex_surfaces'],
                                      'runtime_enumerating_surfaces_from_vertex_surfaces': info['runtime_enumerating_surfaces_from_vertex_surfaces']}, index=[0])])
        else:
            manifold_data = pd.concat([manifold_data, pd.DataFrame({'index': i,
                                      'manifold': info['manifold'],
                                      'runtime_surfaces_total': info['runtime_surfaces'],
                                      'runtime_gp': info['runtime_gp'],
                                      'num_tot_geo': len(info['tot_geo'])}, index=[0])])

    manifold_data.to_csv('/data/keeling/a/chaeryn2/totally-geodesic/result_S21twist.csv', index=False)

if __name__ == '__main__':
    sort_results()