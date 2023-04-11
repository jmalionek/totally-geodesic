import snappy, regina
import os
import pickle
import pandas as pd

def sort_results():
    manifold_data = pd.DataFrame({'index':pd.Series(dtype='int'),
                                  'manifold':pd.Series(dtype='str'),
                                  'runtime_surfaces':pd.Series(dtype='float'),
                                  'runtime_gp':pd.Series(dtype='float'),
                                  'num_tot_geo':pd.Series(dtype='int')})
    for i in range(271636):
        if 'totally_geodesic_info_manifold%i'%i in os.listdir('/data/keeling/a/chaeryn2/totally_geodesic/'):
            with open('/data/keeling/a/chaeryn2/totally_geodesic/totally_geodesic_info_manifold%i'%i, 'rb') as f:
                info = pickle.load(f)
            pd.concat([manifold_data, {'index':i,
                                      'manifold':info['manifold'],
                                      'runtime_surfaces':info['runtime_surfaces'],
                                      'runtime_gp':info['runtime_gp'],
                                      'num_tot_geo':len(info['tot_geo'])}])
    manifold_data.to_csv('/data/keeling/a/chaeryn2/totally-geodesic/result.csv', index=False)

if __name__ == '__main__':
    sort_results()