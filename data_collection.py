#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=/data/keeling/a/chaeryn2/SLURM_print/data_collection%A
#SBATCH --error=/data/keeling/a/chaeryn2/SLURM_error/data_collection%A


import os
import regina
import snappy
import matplotlib.pyplot as plt
import normal_surfaces
import pickle
import numpy as np
import shutil


def get_all_results(location, filename_word = None, get_manifolds = False):
	dir = location
	if filename_word is None:
		files = os.listdir(dir)
	else:
		files = [file for file in os.listdir(dir) if filename_word in file]

	keys = ['runtime_surfaces', 'runtime_gp', 'num_surfaces', 'manifold']


	all_results = {key:[] for key in keys}
	for file in files:
		with open(dir + file, 'rb') as f:
			results = pickle.load(f)
		for key in keys:
			if key not in results.keys():
				if key == 'runtime_surfaces':
					all_results[key].append(results['runtime_vertex_surfaces'] + results['runtime_enumerate_surfaces'])
				elif key == 'runtime_gp':
					all_results[key].append(results['runtime_tot_geo'])
				elif key == 'num_surfaces':
					if 'all_surfaces_vec' in results.keys():
						all_results[key].append(len(results['all_surfaces_vec']))
					else:
						all_results[key].append(len(results['all_surfaces']))
				else:
					print('BAD THINGS HAVE HAPPENED')
					assert False
			else:
				all_results[key].append(results[key])

	if get_manifolds:
		manifolds = []
		for code in all_results['manifold']:
			Mr = regina.Triangulation3.tightDecoding(code)
			manifolds.append(snappy.Manifold(Mr))
		return all_results, manifolds

	return all_results



def make_runtime_histogram(attribute, num_bins = 30, filename_word = None, image_name = None, **kwargs):
	all_results = get_all_results(filename_word)

	if image_name is None:
		image_name = 'histogram.png'

	fig, ax = plt.subplots()

	ax.histogram(all_results[attribute], all_num_bins = num_bins, **kwargs)

	fig.save(image_name)


def make_scatterplot(x_name, y_name, filename_word = None, image_name = None, **kwargs):
	all_results = get_all_results(filename_word)

	if image_name is None:
		image_name = 'scatter.png'

	fig, ax = plt.subplots()

	ax.histogram(all_results[x_name], all_results[y_name], num_bins=num_bins, **kwargs)

	fig.save(image_name)


def volume_scatterplot(y_name, filename_word = None, image_name = None, **kwargs):
	all_results, manifolds = get_all_results(filename_word, True)

	if image_name is None:
		image_name = 'scatter.png'

	fig, ax = plt.subplots()

	volumes = [M.volume() for M in manifolds]

	ax.scatter(volumes, all_results[y_name], num_bins=num_bins, **kwargs)

	fig.save(image_name)

def weekly_check():
	keywords = ['link', 'cover']

	for keyword in keywords:
		print(keyword)
		results, manifolds = get_all_results(keyword, get_manifolds = True)
		print('number done:', len(results['tot_geo']))
		print('"number with totally geodesic surfaces":', len([l for l in results['tot_geo'] if len(l) > 0]))
		actual = [i for i in range(len(manifolds)) if len(results['tot_geo'][i]) > 0 and manifolds[i].solution_type(enum=True) <= 2]
		print('Number of manifolds which ACTUALLY have tot geo', len(actual))
		print()


def number_unfinished():
	num_done = len([name for name in os.listdir('/data/keeling/a/chaeryn2/tg_computation_outputs/') if 'link' in name])
	num_total = len(snappy.HTLinkExteriors(alternating=False)[7.2:])
	return num_total - num_done

def unfinished_list():
	unfinished = set(range(len(snappy.HTLinkExteriors(alternating=False)[7.2:])))
	for name in os.listdir('/data/keeling/a/chaeryn2/results_links_names_fixed/'):
		if 'link' in name:
			index = int(name.split('_')[2])
			unfinished.discard(index)
	return unfinished

def discard_nonhyperbolic():
	uf = unfinished_list()
	hyperbolic = set()
	HTlinkexterior = list(snappy.HTLinkExteriors(alternating=False)[7.2:])
	for i in uf:
		M = HTlinkexterior[i]
		if M.solution_type(enum=True) == 1 or M.solution_type(enum=True) == 2:
			hyperbolic.add(i)
	return hyperbolic

def sort_file_names():
	htlinkexterior = snappy.HTLinkExteriors(alternating=False)[7.2:]
	htlinkexterior_names = [M.name() for M in htlinkexterior]

	dirname = '/data/keeling/a/chaeryn2/results_links_names_fixed/'
	for file in os.listdir('/data/keeling/a/chaeryn2/tg_computation_outputs/'):
		if 'link' in file:
			# print('tg_computation_outputs/', file)
			with open('/data/keeling/a/chaeryn2/tg_computation_outputs/' + file, 'rb') as f:
				data = pickle.load(f)
			index = int(file.split('_')[3][4:])

			MR = regina.Triangulation3.tightDecoding(data['manifold'])
			N = snappy.Manifold(MR)
			M = htlinkexterior[index]

			filename = f'link_info_{index}_{M.name()}'

			with open(dirname + filename, 'wb') as f:
				pickle.dump(data, f)
			# print('changed to', filename)

	for file in os.listdir('/data/keeling/a/chaeryn2/computation_outputs/'):
		if 'link' in file:
			# print('computation_outputs/', file)
			manifold_name = file.split('_')[-1]
			true_index = htlinkexterior_names.index(manifold_name)
			true_name = f'link_info_{true_index}_{manifold_name}'
			shutil.copy('/data/keeling/a/chaeryn2/computation_outputs/' + file, dirname + true_name)
			# print('changed to', true_name)

def sort_file_names_check_done():
	htlinkexterior = snappy.HTLinkExteriors(alternating=False)[7.2:]
	htlinkexterior_names = [M.name() for M in htlinkexterior]

	dirname = '/data/keeling/a/chaeryn2/results_links_names_fixed/'
	for file in os.listdir('/data/keeling/a/chaeryn2/tg_computation_outputs/'):
		if 'link' in file:
			print(file)
			index = int(file.split('_')[3][4:])
			M = htlinkexterior[index]
			filename = f'link_info_{index}_{M.name()}'

			if filename not in os.listdir('/data/keeling/a/chaeryn2/results_links_names_fixed/'):
				print('changed to ', filename)
				shutil.copy('/data/keeling/a/chaeryn2/tg_computation_outputs/' + file, dirname + filename)

	for file in os.listdir('/data/keeling/a/chaeryn2/computation_outputs/'):
		if 'link' in file:
			print(file)
			manifold_name = file.split('_')[-1]
			true_index = htlinkexterior_names.index(manifold_name)
			true_name = f'link_info_{true_index}_{manifold_name}'
			if true_name not in os.listdir('/data/keeling/a/chaeryn2/results_links_names_fixed/'):
				print('changed to ', true_name)
				shutil.copy('/data/keeling/a/chaeryn2/computation_outputs/' + file, dirname + true_name)

def save_plots():
	dir = '/data/keeling/a/chaeryn2/totally-geodesic/data_plots/'
	link_data, link_manifolds = get_all_results('/data/keeling/a/chaeryn2/results_links_names_fixed/',
	                                            get_manifolds=True)
	cover_data, cover_manifolds = get_all_results('/data/keeling/a/chaeryn2/computation_outputs/', filename_word='cover', get_manifolds=True)
	manifolds = np.array(link_manifolds + cover_manifolds)

	with open('data_collected', 'wb') as f:
		pickle.dump([link_data, link_manifolds, cover_data, cover_manifolds], f)


	# in case of repeated use
	# with open('data_collected', 'rb') as file:
	# 	f = pickle.load(file)
	# 	link_data, link_manifolds, cover_data, cover_manifolds = f

	times_enum = np.array(link_data['runtime_surfaces'] + cover_data['runtime_surfaces'])
	times_tg = np.array(link_data['runtime_gp'] + cover_data['runtime_gp'])
	times = times_tg + times_enum
	num_sfces = np.array(link_data['num_surfaces'] + cover_data['num_surfaces'])

	# We exclude outliers where runtime was too large
	good_times = (times < 5000)
	times_enum = times_enum[good_times]
	times_tg = times_tg[good_times]
	times = times[good_times]
	manifolds = manifolds[good_times]

	print('average runtime: ', np.average(times))

	fig, ax = plt.subplots()
	ax.hist(times_enum, bins=30, edgecolor='black')
	ax.set_xlabel('Runtime of enumerating surfaces in seconds')
	fig.savefig(dir + 'runtime_enum_histogram.png')
	plt.close(fig)

	# When considering runtimes of Algorithm2 we exclude any manifolds where there were no surfaces
	# (i.e. manifolds where runtime for Algorithm2 is 0)
	good_times = (num_sfces > 0)
	times_enum = times_enum[good_times]
	times_tg = times_tg[good_times]
	times = times[good_times]
	manifolds = manifolds[good_times]

	fig, ax = plt.subplots()
	ax.hist(times_tg, bins=30, edgecolor='black')
	ax.set_xlabel('Runtime of Algorithm 2 in seconds')
	fig.savefig(dir + 'runtime_tot_geo_histogram.png')
	plt.close(fig)

	fig, ax = plt.subplots()
	volumes = []
	tet_dict = {}
	tetrahedra = []
	for i, M in enumerate(manifolds):
		if 'solution' not in M.solution_type():
			volumes.append(M.volume())
			tet = M.num_tetrahedra()
			tetrahedra.append(tet)
			if tet in tet_dict.keys():
				tet_dict[tet].append(times[i])
			else:
				tet_dict[tet] = [times[i]]
		else:
			M.randomize()
			volumes.append(M.volume())
			tet = M.num_tetrahedra()
			tetrahedra.append(tet)
			if tet in tet_dict.keys():
				tet_dict[tet].append(times[i])
			else:
				tet_dict[tet] = [times[i]]

	indices = np.random.choice(np.arange(len(volumes)), size=5000)
	volumes = np.array(volumes)
	tetrahedra = np.array(tetrahedra)
	ax.scatter(volumes[indices], times[indices], s=5, alpha=.2)
	ax.set_yscale('log')
	ax.set_xlabel('Manifold volume')
	ax.set_ylabel('Log of algorithm runtime in seconds')
	fig.savefig(dir + 'volume_runtime_scatter.png')
	plt.close(fig)

	# fig, ax = plt.subplots()
	# indices = np.random.choice(np.arange(len(volumes)), size=5000)
	# volumes = np.array(volumes)
	# ax.scatter(volumes[indices], times_enum[indices], s=5, alpha=.2, c='red', marker="o", label='Enumerating surfaces')
	# ax.scatter(volumes[indices], times_tg[indices], s=5, alpha=.2, c='blue', marker="X", label='Algorithm 2')
	# # ax.set_yscale('log')
	# ax.set_xlabel('Manifold volume')
	# ax.set_ylabel('Runtime in seconds')
	# ax.legend()
	# fig.savefig(dir + 'volume_runtime_scatter_both.png')
	#
	# fig, ax = plt.subplots()
	# num_tets_list = sorted(list(tet_dict.keys()))
	# data = [tet_dict[i] for i in num_tets_list]
	# ax.boxplot(data, labels = num_tets_list)
	# ax.set_xlabel('Number of tetrahedra')
	# ax.set_ylabel('Algorithm runtime in seconds')
	# ax.set_yscale('log')
	# fig.savefig(dir + 'box_plot_num_tets_runtime.png')

	fig, ax = plt.subplots()
	indices = np.random.choice(np.arange(len(volumes)), size=5000)
	ax.scatter(tetrahedra[indices], times[indices], s=5)
	ax.set_xticks(np.unique(tetrahedra[indices]))
	ax.set_yscale('log')
	ax.set_xlabel('Number of tetrahedra')
	ax.set_ylabel('Log of algorithm runtime in seconds')
	fig.savefig(dir + 'tetrahedra_runtime_scatter.png')
	plt.close(fig)

	print('Average of enumeration runtime ratios', np.average(times_enum/times))
	print('Average of algorithm2 runtime ratios', np.average(times_tg/times))

	fig, ax = plt.subplots()
	indices = np.random.choice(np.arange(len(times_tg)), size=5000)
	ax.scatter(times_enum[indices], times_tg[indices], s=5, alpha=.2)
	ax.set_aspect('equal')
	ax.set_xlabel('Runtime of enumerating surfaces in seconds')
	ax.set_ylabel('Runtime of Algorithm 2 in seconds')
	fig.savefig(dir + 'runtime_enumvstg_scatter.png')
	plt.close(fig)


if __name__ == '__main__':
	save_plots()
	# downloading files from keeling: scp chaeryn2@keeling.earth.illinois.edu:/path_to_file /path_to_file_on_your_laptop