import nscomplex
import regina
import snappy
import networkx as nx
import matplotlib.pyplot as plt
from sage.all import block_matrix
from complex_reps import preserves_hermitian_form


class NormalSurface:
	def __init__(self, surface, manifold):
		self.polygons_list = []
		# A list of lists of lists, where the innermost list contains all of the discs which are
		# of the specific tetrahedron and type. i.e. polygons[2][3] is a list of all the discs in
		# the normal surface which are in tetrahedron 2 and of type 3.
		self.polygons = [[[] for i in range(7)] for j in range(manifold.num_tetrahedra())]
		self.surface = surface
		self.manifold = manifold

	def add_disc(self, disc):
		self.polygons_list.append(disc)
		self.polygons[disc.tetrahedron][disc.disc_type].append(disc)

	def get_polygon(self, tet_number, disc_type, disc_index):
		return self.polygons[tet_number][disc_type][disc_index]

	def dual_graph(self):
		"""
		Given a NormalSurface object from this package, returns the dual graph
		"""
		G = nx.MultiDiGraph()
		for polygon in self.polygons_list:
			G.add_node(polygon.get_id_numbers())
		for polygon in self.polygons_list:
			for index, adj_disc in enumerate(polygon.adjacent_discs):
				if not G.has_edge(adj_disc.get_id_numbers(), polygon.get_id_numbers(), key=polygon.faces[index].index()):
					G.add_edge(polygon.get_id_numbers(), adj_disc.get_id_numbers(), key=polygon.faces[index].index())
		return G

	def fundamental_group_generators(self, return_tree=False):
		undirected_dual_graph = self.dual_graph().to_undirected()
		tree_graph = nx.minimum_spanning_tree(undirected_dual_graph)
		# print(f'size of the spanning tree is {len(tree_graph.edges())}')
		diff = nx.difference(undirected_dual_graph, tree_graph)
		# print(f'the number of edges in the difference {len(diff.edges())}')
		nx.draw(tree_graph)
		if return_tree:
			return list(diff.edges(keys=True)), tree_graph
		return diff.edges

	def fundamental_group_embedding(self):
		self.manifold._choose_generators(False,False)
		gen_info = self.manifold._choose_generators_info()
		basepoint_tet = -1
		for entry in gen_info:
			if entry['generator_path'] == -1:
				basepoint_tet = entry['index']
				break
		basepoint_disc = self.polygons_list[0]
		for disc in self.polygons_list:
			if disc.tetrahedron == basepoint_tet:
				basepoint_disc = disc
				break
		cycles = []
		node_paths = []
		generators, tree = self.fundamental_group_generators(True)
		for edge in generators:
			path_to_edge = nx.shortest_path(tree, basepoint_disc.get_id_numbers(), edge[0])
			path_from_edge = nx.shortest_path(tree, edge[1], basepoint_disc.get_id_numbers())
			edge_path = []
			for i in range(len(path_to_edge) - 1):
				edge_path.append((path_to_edge[i], path_to_edge[i+1], list(tree.get_edge_data(path_to_edge[i], path_to_edge[i+1]).keys())[0]))
			edge_path.append(edge)
			for i in range(len(path_from_edge) - 1):
				edge_path.append((path_from_edge[i], path_from_edge[i+1], list(tree.get_edge_data(path_from_edge[i], path_from_edge[i+1]).keys())[0]))
			cycles.append(edge_path)
			node_paths.append(path_to_edge + path_from_edge)
		gens_in_M = []
		for i, cycle in enumerate(cycles):
			gen_in_M = []
			for j, this_edge in enumerate(cycle):
				if this_edge[0] == node_paths[i][j]:
					gen = self.find_manifold_generator(this_edge)
					if gen != 0:
						gen_in_M.append(gen)
				elif this_edge[1] == node_paths[i][j]:
					gen = -self.find_manifold_generator(this_edge)
					if gen != 0:
						gen_in_M.append(gen)
			gens_in_M.append(gen_in_M)
		return gens_in_M

	def find_manifold_generator(self, edge):
		T = self.surface.triangulation()
		starting_tet = self.get_polygon(*edge[0]).tetrahedron

		if T.triangle(edge[2]).embedding(0).tetrahedron().index() == starting_tet:
			initial_face_embedding = T.triangle(edge[2]).embedding(0)
		else:
			initial_face_embedding = T.triangle(edge[2]).embedding(1)
		face_index = initial_face_embedding.face()
		self.manifold._choose_generators(True, True)
		info = self.manifold._choose_generators_info()
		gen = info[starting_tet]['generators'][face_index]
		return gen

	def intersecting_edges(self, normal_disc, return_edge_embedding=False):
		# normal_disc should be of a Polygon class
		T = self.surface.triangulation
		tet_num = normal_disc.tetrahedron
		disc_type = normal_disc.disc_type
		edge_list = []
		if normal_disc.is_triangle():
			for n in range(0, 4):
				if disc_type != n:
					edge_list.append((disc_type, n))
		if normal_disc.is_quad():
			if disc_type == 4:
				edge_list = [(0, 2), (1, 2), (1, 3), (0,3)]
			if disc_type == 5:
				edge_list = [(0, 1), (1, 2), (2, 3), (0, 3)]
			if disc_type == 6:
				edge_list = [(0, 1), (0, 2), (2, 3), (1, 3)]
		if return_edge_embedding:
			return edge_list
		edge_indices = []
		for e in edge_list:
			edge_indices.append(T.tetrahedron(tet_num).edge(*e))
		return edge_indices

	def relations(self):
		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
		for normal_disc in self.polygons_list:
			for i, edge_index in enumerate(intersecting_edges(normal_disc)):
				edge_embedding = intersecting_edges(normal_disc, True)[i]
				edge = T.Edges[edge_index]
				arrow = edge.get_arrow()
				while True:
					binary_list = [0, 0, 0, 0]
					for n in edge_embedding:
						binary_list[n] = 1
					if arrow.Edge == snappy.snap.t3mlite.simplex.bitmap(binary_list) and arrow.Tetrahedron.Index == normal_disc.tetrahedron:
						break
					else:
						arrow = arrow.next()
		# To-do:
		# 1. for each arrow(=face of triangulation) we find the disk adjacent along that arrow
		# 2. continue edge_valence times to get a cycle of discs
		# 3. check if start = finish so that it is indeed a relator
		pass

class Polygon:
	def __init__(self, manifold):
		# Stored in edge order based on regina's weird permuation system
		self.adjacent_discs = []
		# The tetrahedron that is adjacent to this one across the edge in edges
		self.adjacent_tets = []
		# The face that the edge lies on in the given triangulation
		self.faces = []
		# The tetrahedron that this normal disk sits inside
		self.tetrahedron = None
		# The manifold that this normal disc lies in
		self.manifold = manifold
		# The disc type (a number from 0-6 inclusive, no octs or tubes!)
		self.disc_type = None
		# The index of the disc (there can be multiple discs which have the same type and tet number)
		self.disc_index = None

	def get_id_numbers(self):
		return self.tetrahedron, self.disc_type, self.disc_index

	def is_triangle(self):
		return False

	def is_quad(self):
		return False

	def __repr__(self):
		d_type = 'no_d_type'
		if self.is_quad():
			d_type = 'quad'
		if self.is_triangle():
			d_type = 'tri'
		return f'{d_type}:{self.tetrahedron}.{self.disc_type}.{self.disc_index}'


class Triangle (Polygon):
	def __init__(self, manifold):
		super().__init__(manifold)
		self.adjacent_discs = [None] * 3
		self.adjacent_tets = [None] * 3
		self.faces = [None] * 3

	def is_triangle(self):
		return True


class Quad (Polygon):
	def __init__(self, manifold):
		super().__init__(manifold)
		self.adjacent_discs = [None] * 4
		self.adjacent_tets = [None] * 4
		self.faces = [None] * 4

	def is_quad(self):
		return True


def from_regina_normal_surface(surface, manifold):
	"""
	Given a regina normal surface, and a snappy manifold, returns a NormalSurface object
	"""
	DSS = regina.DiscSetSurface(surface)
	our_surface = NormalSurface(surface, manifold)
	T = surface.triangulation()
	face_list = list(T.faces(2))
	face_gluings = []
	for i in range(len(face_list)):
		face_gluings.append(list(face_list[i].embeddings()))
	for tet_index in range(manifold.num_tetrahedra()):
		for i in range(4):
			num_tris = surface.triangles(tet_index, i)
			for tri_index in range(num_tris.longValue()):
				triangle = Triangle(manifold)
				triangle.disc_type = i
				triangle.tetrahedron = tet_index
				triangle.disc_index = tri_index
				our_surface.add_disc(triangle)

		for i in range(3):
			num_quads = surface.quads(tet_index, i)
			for quad_index in range(num_quads.longValue()):
				quad = Quad(manifold)
				quad.disc_type = 4 + i
				quad.tetrahedron = tet_index
				quad.disc_index = quad_index
				our_surface.add_disc(quad)

	for disc in our_surface.polygons_list:
		if disc.is_triangle():
			discSpec = regina.DiscSpec(disc.tetrahedron, disc.disc_type, disc.disc_index)
			for edge_index, perm in enumerate(regina.triDiscArcs[disc.disc_type]):
				adj_disc, other_perm = DSS.adjacentDisc(discSpec, perm)
				face_index = perm[3]
				for i, embeddings in enumerate(face_gluings):
					for embedding in embeddings:
						if embedding.face() == face_index and embedding.simplex().index() == disc.tetrahedron:
							disc.faces[edge_index] = face_list[i]
				# WHY BEN?!
				disc.adjacent_discs[edge_index] = our_surface.get_polygon(adj_disc.tetIndex, adj_disc.type, adj_disc.number)
				disc.adjacent_tets[edge_index] = adj_disc.tetIndex
		elif disc.is_quad():
			discSpec = regina.DiscSpec(disc.tetrahedron, disc.disc_type, disc.disc_index)
			for edge_index, perm in enumerate(regina.quadDiscArcs[disc.disc_type - 4]):
				adj_disc, other_perm = DSS.adjacentDisc(discSpec, perm)
				face_index = perm[3]
				for i, embeddings in enumerate(face_gluings):
					for embedding in embeddings:
						if embedding.face() == face_index and embedding.simplex().index() == disc.tetrahedron:
							disc.faces[edge_index] = face_list[i]
				disc.adjacent_discs[edge_index] = our_surface.get_polygon(adj_disc.tetIndex, adj_disc.type, adj_disc.number)
				disc.adjacent_tets[edge_index] = adj_disc.tetIndex
	return our_surface


def Tietze_to_string(word):
	alphabet = 'abcdefghijklmnopqrstuvwxyz'
	our_alphabet = '?' + alphabet + alphabet.upper()[::-1]
	return ''.join([our_alphabet[index] for index in word])


def Tietze_to_matrix(word, M):
	return M.fundamental_group(simplify_presentation=False).SL2C(Tietze_to_string(word))


def surface_generators(surface, manifold, SL2C=True):
	our_surface = from_regina_normal_surface(surface, manifold)
	generators = our_surface.fundamental_group_embedding()
	if SL2C:
		return [Tietze_to_matrix(word, manifold) for word in generators]
	else:
		return [Tietze_to_string(word) for word in generators]

def surface_group_in_PSL2R(surface, manifold):
	gens = surface_generators(surface, manifold)
	return preserves_hermitian_form(gens)[0]


def run_on_knots():
	# number of alternating knots starting at crossing number 4
	num_alternating = [1, 2, 3, 7, 18, 41, 123, 367, 1288, 4878, 19536, 85263, 379799, 1769979, 8400285]
	for crossing, num in zip(range(4, 20), num_alternating):
		knot_strs = [f'K{crossing}a{i+1}' for i in range(num)]
		for string in knot_strs:
			M = snappy.Manifold(string)
			print(M)
			if M.volume() < 2:
				print('Manifold not hyperbolic')
				continue
			print(f'num tetrahedra: {M.num_tetrahedra()}')
			euler = -2
			while(True):
				try:
					cs = nscomplex.ConnectedSurfaces(M, euler)
					break
				except AssertionError:
					euler = euler-2
			LW = cs.essential_faces_of_normal_polytope()
			surfaces = []
			for face in LW.maximal:
				surfaces = surfaces + face.vertex_surfaces
			if len(surfaces) == 0:
				print('Manifold has no closed essential surfaces which are not boundary parallel')
			for surface in surfaces:
				result = surface_group_in_PSL2R(surface.surface, M)
				print(result, M)

def boundary_to_surface(M):
	# M: regina triangulation, make sure to run it only on knot complements
	num_tet = M.size()
	surface_vector = [1, 1, 1, 1, 0, 0, 0] * num_tet
	return regina.NormalSurface(M, regina.NS_STANDARD, surface_vector)


def main2():
	import doubling
	T = doubling.read_triangulation_info('example_manifold_2nd_last.txt')
	Tr = T.regina_triangulation()
	Tr.idealToFinite()
	Tr.intelligentSimplify()
	Tr = doubling.double_regina_triangulation(Tr)
	Tr.finiteToIdeal()
	Tr.intelligentSimplify()
	M = snappy.Manifold(Tr.snapPea())
	print('M:', M.fundamental_group(simplify_presentation=False))
	# print(M)
	# print(M.num_tetrahedra())
	# M.randomize()
	# M.simplify()
	# print(M.num_tetrahedra())
	surfaces = regina.NormalSurfaces(regina.Triangulation3(M), regina.NS_QUAD_CLOSED, algHints=regina.NS_VERTEX_DD)
	surfaces_gen2 = [S for S in surfaces if S.eulerChar() == -2]
	for surface in surfaces_gen2:
		print(surface_group_in_PSL2R(surface, M))
		our_surface = from_regina_normal_surface(surface, M)
		gens = our_surface.fundamental_group_embedding()
		gens_matrix = [Tietze_to_matrix(gen, M).trace().imag().abs() for gen in gens]
		if max(gens_matrix) < 0.1:
		   print(surface)
		   for gen in gens:
			   print(Tietze_to_matrix(gen, M))
		# print(surface_group_in_PSL2R(surface, M))
	# # print(M.num_cusps())
	# N = snappy.Manifold('4_1')
	# print(N)
	# N.dehn_fill((1,0), 0)
	# print(N.fundamental_group())
	# print(N)
	# print(N.num_cusps())
	# M.dehn_fill((92341, 700), 1)
	# M.randomize()
	# D = M.dirichlet_domain()
	# M = snappy.Manifold(D)
	# print(M.num_tetrahedra())
	# print(M.num_cusps())
	# print(M.identify())
	# Tr.normalSurfaces()
	# cs = nscomplex.ConnectedSurfaces(M, -2)
	# for surface in S:
	#	print(surface_group_in_PSL2R(surface, M))


def main():
	import matplotlib.pyplot as plt
	import doubling
	M = snappy.Manifold('10_135')
	print(f'num tetrahedra: {M.num_tetrahedra()}')
	surfaces = nscomplex.ConnectedSurfaces(M, -8)
	print(regina.NormalSurfaces(regina.Triangulation3(M), regina.NS_AN_QUAD_OCT_CLOSED))
	print(surfaces.essential_faces_of_normal_polytope())
	print(f'number incompressible: {len(surfaces.incompressible)}')
	regina_surface = surfaces.normal[0].surface
	print(regina_surface.vector())
	our_surface = from_regina_normal_surface(regina_surface, M)
	surface_group_in_PSL2R()
	# G = our_surface.dual_graph()
	# nx.draw(G)
	# print(len(G.edges()))
	# for edge in G.edges():
	#	print(edge)
	# plt.savefig('/home/joseph/Desktop/graph example.png')
	# print(type(list(our_surface.fundamental_group_generators())[0]))
	# print(surface_group_in_PSL2R(regina_surface, M))


def main3():
	M = snappy.Manifold('4_1')
	G = M.fundamental_group(simplify_presentation=False)
	print('meridian:', G.meridian())
	print('longitude:', G.longitude())
	boundary_surface = from_regina_normal_surface(boundary_to_surface(regina.Triangulation3(M)), M)
	gens = boundary_surface.fundamental_group_embedding()
	gens_matrix = [Tietze_to_matrix(gen, M) for gen in gens]
	gens_string = [Tietze_to_string(gen) for gen in gens]
	print(gens_string)
	for i in range(len(gens)):
		for j in range(len(gens)):
			if i == j:
				break
			else:
				print(gens_matrix[i] * gens_matrix[j] * gens_matrix[i].inverse() * gens_matrix[j].inverse())


if __name__ == '__main__':
	main2()
