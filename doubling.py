import snappy
import regina
from snappy.snap.t3mlite import Mcomplex


def read_triangulation_info(filename):
	with open(filename, 'r') as file:
		words = file.read().split()
		ograph_data = words
	# print(ograph_data)
	gluing_data = []
	for i in range(0, len(ograph_data), 8):
		tet_gluings = [int(t) for t in ograph_data[i:i+4]]
		perms = [[int(v) for v in perm] for perm in ograph_data[i+4:i+8]]
		gluing_data.append((tet_gluings, perms))
	T = Mcomplex(gluing_data)
	return T


def double_regina_triangulation(T):
	"""
	T is a regina triangulation with exactly one non-torus boundary component.
	"""
	n = T.size()
	T.insertTriangulation(T)
	boundary_component = None
	for component in T.boundaryComponents():
		if component.eulerChar() < 0:
			boundary_component = component
			print(boundary_component)
			break
	print(T.boundaryComponents())
	tet_indices = []
	face_indices = []
	for tri in boundary_component.triangles():
		# print('get inside?')
		tet_indices.append(tri.embedding(0).simplex().index())
		face_indices.append(tri.embedding(0).face())
		# print(tet, face)
		# print(T.tetrahedron(tet) in T.tetrahedra())
		# T.tetrahedron(tet).join(face, T.tetrahedron(tet + n), regina.Perm4(0, 1, 2, 3))
	for i in range(len(tet_indices)):
		# print('get inside?')
		# tet = tri.embedding(0).simplex().index()
		# face = tri.embedding(0).face()
		# print(tet, face)
		# print(T.tetrahedron(tet) in T.tetrahedra())
		tet = tet_indices[i]
		face = face_indices[i]
		T.tetrahedron(tet).join(face, T.tetrahedron(tet + n), regina.Perm4(0, 1, 2, 3))
		# print('get here too?')
	# print(T.detail())
	return T



def double_snappy_triangulation(T):
	# TODO: Change from snappy to all regina using regina.Triangulation3.fromGluings
	Tr = T.regina_triangulation()
	Tr.idealToFinite()
	Tr.intelligentSimplify()
	print(Tr.size())
	gluing_info = []
	num_tets = Tr.size()
	for tet in Tr.tetrahedra():
		adj_tets = []
		adj_perms = []
		for i in range(4):
			if tet.adjacentSimplex(i) is None:
				# THIS LINE is the doubling
				adj_tets.append(tet.index() + num_tets)
			else:
				adj_tets.append(tet.adjacentSimplex(i).index())
		for i in range(4):
			if adj_tets[i] < 0:
				# THIS LINE is the doubling
				adj_perms.append([0, 1, 2, 3])
			else:
				perm = tet.adjacentGluing(i)
				perm_list = []
				for char in perm.str():
					perm_list.append(int(char))
				adj_perms.append(perm_list)
		gluing_info.append((adj_tets, adj_perms))
	gluing_info2 = []
	for i in range(num_tets):
		adj_tets = []
		for face in range(4):
			# THESE LINES are doubling
			if gluing_info[i][0][face] >= num_tets:
				adj_tets.append(gluing_info[i][0][face] - num_tets)
			else:
				adj_tets.append(gluing_info[i][0][face] + num_tets)
		gluing_info2.append((adj_tets, gluing_info[i][1]))
	gluing_info = gluing_info + gluing_info2
	gluing_info[0][0][0] = None
	gluing_info[0][1][0] = None
	Tdouble = Mcomplex(gluing_info)
	Trdouble = Tdouble.regina_triangulation()
	print(list(Trdouble.boundaryComponents()))
	return Trdouble


def main():
	T = read_triangulation_info('example_manifold_2nd_last.txt')
	Tr = T.regina_triangulation()
	Tr.idealToFinite()
	Tr.intelligentSimplify()
	Tr = double_regina_triangulation(Tr)
	print(Tr.detail())
	Tr.finiteToIdeal()
	print(Tr.detail())
	# print(Tr)
	# print(Tr.snapPea())
	# M = snappy.Manifold(Tr.snapPea())
	# print(M)
	# print(M.volume())

	# print(Tr.detail())

	# Tr = T.regina_triangulation()
	# print(Tr.detail())
	# print(Tr.isOrientable())
	# print(Tr.fundamentalGroup())


if __name__ == '__main__':
	main()

