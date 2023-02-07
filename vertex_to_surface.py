import snappy
import regina
import itertools
from sage.all import vector



def surface_to_vec(surface):
	return vector([int(i.stringValue()) for i in surface.vector()], immutable=True)



def find_surfaces_bad(M):
	"""
	Takes in a snappy manifold and gives normal surfaces in M
	"""
	T = regina.Triangulation3(M)
	NS = regina.NormalSurfaces(T, regina.NS_STANDARD, regina.NS_VERTEX)
	print(f'Length of NS: {len(NS)}')
	faces = []
	N = len(NS)
	for bits in itertools.product(*([[0, 1]] * len(NS))):
		surfaces_list = [NS[i] for i in range(N) if bits[i] == 1]
		surfaces = {surface_to_vec(NS[i]) for i in range(N) if bits[i] == 1}
		compatible = True
		for a, b in itertools.combinations(surfaces_list, 2):
			if not a.locallyCompatible(b):
				compatible = False
				break
		if compatible:
			faces.append(surfaces)
	faces_original = faces.copy()
	# TODO: The list ends up empty. WHY?
	print(faces_original)
	N = len(faces)
	for i in range(N):
		for j in range(N):
			if faces_original[i].issubset(faces_original[j]):
				faces[i] = None
	return faces




# def find_surfaces_helper(surfaces, indices_not_used):
	# print(f'indices not used:{indices_not_used}')
	# N = len(surfaces)
	# # Check all pairs of surfaces to see if this is a face
	# if N - len(indices_not_used) <= 1:
	# 	for i in range(N):
	# 		if i not in indices_not_used:
	# 			return []
	# compatible = True
	# for i in range(len(surfaces)):
	# 	if i not in indices_not_used:
	# 		for j in range(i+1, N):
	# 			if j not in indices_not_used:
	# 				if not surfaces[i].locallyCompatible(surfaces[j]):
	# 					compatible = False
	# 			if not compatible:
	# 				break
	# 	if not compatible:
	# 		break
	# if compatible:
	# 	return [[surfaces[i] for i in range(N) if i not in indices_not_used]]
	# else:
	# 	# If it's not a face, then remove the indices one at a time and see if these are subfaces
	# 	faces = []
	# 	for i in range(N):
	# 		if i not in indices_not_used:
	# 			faces.extend(find_surfaces_helper(surfaces, indices_not_used + [i]))
	# 	return faces


def main():
	M = snappy.Manifold('K10n10')
	faces = find_surfaces_bad(M)
	print(faces)



if __name__ == '__main__':
	main()