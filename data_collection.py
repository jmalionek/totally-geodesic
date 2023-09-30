import os
import regina
import snappy
import matplotlib.pyplot as plt
import normal_surfaces
import pickle


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

	ax.histogram(volumes, all_results[y_name], num_bins=num_bins, **kwargs)

	fig.save(image_name)

def weekly_check():
	keywords = ['link', 'cover']

	for keyword in keywords:
		print(keyword)
		results, manifolds = get_all_results(keyword, get_manifolds = True)
		print('number done:', len(results[0]))
		print('"number with totally geodesic surfaces":', len([l for l in results['tot_geo'] if len(l) > 0]))
		print('number which ACTUALLY have tot geo', len([i for i in range(len(manifolds)) if len(results['tot_geo'][i]) > 0 and manifolds[i].solution_type(enum=True) <= 2]))
		print()

