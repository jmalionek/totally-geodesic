import os
import regina
import snappy
import matplotlib.pyplot as plt
import normal_surfaces
import pickle
import numpy as np


def get_all_results(filename_word = None, get_manifolds = False):
	dir = '/data/keeling/a/chaeryn2/tg_computation_outputs/'
	if filename_word is None:
		files = os.listdir(dir)
	else:
		files = [file for file in os.listdir(dir) if filename_word in file]

	all_results = {}
	for file in files:
		with open(dir + file, 'rb') as f:
			results = pickle.load(f)
		for key in results.keys():
			if key not in all_results.keys():
				all_results[key] = []
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


def main():

	dir = '/data/keeling/a/chaeryn2/totally-geodesic/'
	# data, manifolds = get_all_results(get_manifolds=True)
	# with open('data_collected', 'wb') as f:
	# 	pickle.dump([data, manifolds], f)

	# in case of repeated use
	with open('data_collected', 'rb') as file:
		f = pickle.load(file)
		data = f[0]
		manifolds = f[1]

	times_enum = np.array(data['runtime_surfaces'])
	times_tg = np.array(data['runtime_gp'])
	times = np.array(data['runtime_surfaces']) + np.array(data['runtime_gp'])

	print('average runtime: ', np.average(times))

	fig, ax = plt.subplots()
	ax.hist(times_enum, bins = 30, edgecolor = 'black')
	ax.set_xlabel('Runtime of enumerating surfaces in seconds')
	fig.savefig(dir + 'runtime_enum_histogram.png')

	fig, ax = plt.subplots()
	ax.hist(times_tg, bins = 30, edgecolor = 'black')
	ax.set_xlabel('Runtime of Algorithm 2 in seconds')
	fig.savefig(dir + 'runtime_tot_geo_histogram.png')

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
	ax.scatter(volumes[indices], times[indices], s = 5, alpha = .2)
	ax.set_yscale('log')
	ax.set_xlabel('Manifold volume')
	ax.set_ylabel('Log of algorithm runtime in seconds')
	fig.savefig(dir + 'volume_runtime_scatter.png')

	fig, ax = plt.subplots()
	indices = np.random.choice(np.arange(len(volumes)), size=5000)
	volumes = np.array(volumes)
	ax.scatter(volumes[indices], times_enum[indices], s=5, alpha=.2, c='red', marker="o", label='Enumerating surfaces')
	ax.scatter(volumes[indices], times_tg[indices], s=5, alpha=.2, c='blue', marker="X", label='Algorithm 2')
	# ax.set_yscale('log')
	ax.set_xlabel('Manifold volume')
	ax.set_ylabel('Log of algorithm runtime in seconds')
	ax.legend()
	fig.savefig(dir + 'volume_runtime_scatter_both.png')

	fig, ax = plt.subplots()
	num_tets_list = sorted(list(tet_dict.keys()))
	data = [tet_dict[i] for i in num_tets_list]
	ax.boxplot(data, labels = num_tets_list)
	ax.set_xlabel('Number of tetrahedra')
	ax.set_ylabel('Algorithm runtime in seconds')
	ax.set_yscale('log')
	fig.savefig(dir + 'box_plot_num_tets_runtime.png')

	fig, ax = plt.subplots()
	indices = np.random.choice(np.arange(len(volumes)), size=5000)
	volumes = np.array(volumes)
	ax.scatter(volumes[indices], times[indices], s=5, alpha=.2)
	# ax.set_yscale('log')
	ax.set_xlabel('Manifold volume')
	ax.set_ylabel('Algorithm runtime in seconds')
	ax.legend()
	fig.savefig(dir + 'tetrahedra_runtime_scatter.png')

	print('Average of runtime ratios', np.average(times_tg/times_enum))

if __name__ == '__main__':
	main()