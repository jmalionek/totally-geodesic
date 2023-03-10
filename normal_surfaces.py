import nscomplex
import regina
import snappy
import pickle
import random
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from sage.all import block_matrix, matrix, vector, CC, FreeGroup
from complex_reps import preserves_hermitian_form
from nscomplex import faces, regina_util, surfaces
from itertools import combinations

class NormalSurface:
	def __init__(self, surface, manifold):
		self.polygons_list = []
		# A list of lists of lists, where the innermost list contains all of the discs which are
		# of the specific tetrahedron and type. i.e. polygons[2][3] is a list of all the discs in
		# the normal surface which are in tetrahedron 2 and of type 3.
		self.polygons = [[[] for i in range(7)] for j in range(manifold.num_tetrahedra())]
		self.surface = surface
		self.manifold = manifold
		self.basepoint = None

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
		self.basepoint = basepoint_disc
		cycles = []
		node_paths = []
		generators, tree = self.fundamental_group_generators(True)
		# add paths from/to basepoint to given edge (generator corresponds to a single edge outside of the tree)
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
			gens_in_M.append(self.find_manifold_generator_in_tree(node_paths[i], cycle))
		return gens_in_M

	def find_manifold_generator(self, edge):
		T = self.surface.triangulation()
		self.manifold._choose_generators(True, True)
		info = self.manifold._choose_generators_info()

		starting_tet = self.get_polygon(*edge[0]).tetrahedron
		ending_tet = self.get_polygon(*edge[1]).tetrahedron
		starting_disc = self.get_polygon(*edge[0])
		ending_disc = self.get_polygon(*edge[1])
		index = starting_disc.adjacent_discs.index(ending_disc)

		if starting_disc.is_triangle():
			disc_type = starting_disc.disc_type
			face_index = regina.triDiscArcs[disc_type][index][3]  # last digit of permutation gives the index of the face inside the tetrahedron (0 ~ 3)
		else:
			disc_type = starting_disc.disc_type - 4
			face_index = regina.quadDiscArcs[disc_type][index][3]  # last digit of permutation gives the index of the face inside the tetrahedron (0 ~ 3)
		gen = info[starting_tet]['generators'][face_index]
		return gen

	def find_manifold_generator_in_tree(self, node_path, edge_path):
		# given a path in a graph retrieve the list of edge and list of nodes, this function compares these two lists to find the right generator
		generators, tree = self.fundamental_group_generators(True)
		list_gens = []
		for i, edge in enumerate(edge_path):
			if edge[0] == node_path[i]:
				gen = self.find_manifold_generator(edge)
				if gen != 0:
					list_gens.append(gen)
			elif edge[1] == node_path[i]:
				gen = -self.find_manifold_generator(edge)
				if gen != 0:
					list_gens.append(gen)
		return list_gens

	def intersecting_edges(self, normal_disc, return_edge_embedding=False):
		# normal_disc should be of a Polygon class
		T = self.surface.triangulation()
		tet_num = normal_disc.tetrahedron
		disc_type = normal_disc.disc_type
		edge_list = []
		if normal_disc.is_triangle():
			for n in range(0, 4):
				if disc_type != n:
					edge_list.append((disc_type, n))
		if normal_disc.is_quad():
			if disc_type == 4:
				edge_list = [(0, 2), (1, 2), (1, 3), (0, 3)]
			if disc_type == 5:
				edge_list = [(0, 1), (1, 2), (2, 3), (0, 3)]
			if disc_type == 6:
				edge_list = [(0, 1), (0, 2), (2, 3), (1, 3)]
		if return_edge_embedding:
			return edge_list
		edge_indices = []
		for e in edge_list:
			edge_indices.append(T.tetrahedron(tet_num).edge(*e).index())
		return edge_indices

	def relations(self):
		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
		T_regina = regina.Triangulation3(self.manifold)
		relators = []

		self.manifold._choose_generators(True, True)
		info = self.manifold._choose_generators_info()

		for normal_disc in self.polygons_list:
			for i, edge_index in enumerate(self.intersecting_edges(normal_disc)):
				relator = []  # contains indices of faces that correspond to a single relation
				# the edges on the tetrahedron that intersect our normal disc as embeddings
				edge_embedding = self.intersecting_edges(normal_disc, True)[i]
				# actual snappy edge
				edge = T.Edges[edge_index]
				arrow = edge.get_arrow()
				while True:
					if arrow.Edge == snappy.snap.t3mlite.simplex.bitmap(edge_embedding) and arrow.Tetrahedron.Index == normal_disc.tetrahedron:
						break
					else:
						arrow = arrow.next()
				current_disc = normal_disc
				current_arrow = arrow

				for n in range(edge.valence()):
					# we look at which face we are gluing our normal disc across
					# we take the index of this face with respect to the tetrahedron that it lies in
					# corresponds to information stored in 'arrow'
					face_index_in_tetrahedron = snappy.snap.t3mlite.simplex.FaceIndex[arrow.Face]  # is a decimal, not binary!
					face = T_regina.tetrahedron(arrow.Tetrahedron.Index).triangle(face_index_in_tetrahedron)
					face_index = face.index()
					# find index of tetrahedron on which our normal disc lies in
					tetrahedron_index = arrow.Tetrahedron.Index
					# find generator
					gen = info[tetrahedron_index]['generators'][face_index_in_tetrahedron]
					# stuff needed to find the next disc
					current_disc_face_index = [face.index() for face in current_disc.faces]
					next_index = current_disc_face_index.index(face_index)
					current_disc = current_disc.adjacent_discs[next_index]
					next_disc = current_disc.adjacent_discs[next_index]
					current_arrow = current_arrow.next()
					if gen != 0:
						relator.append(gen)
					if n == edge.valence() - 1:
						assert current_disc == normal_disc  # the last disc does not glue back to the first disc  # the last disc does not glue back to the first disc

				# find paths from/to basepoint and append to relators
				tree = self.fundamental_group_generators(return_tree=True)[1]
				node_path = nx.shortest_path(tree, self.basepoint.get_id_numbers(), normal_disc.get_id_numbers())  #list of nodes in tree
				edge_path = []
				for i in range(len(node_path) - 1):
					edge_path.append((node_path[i], node_path[i + 1],
									  list(tree.get_edge_data(node_path[i], node_path[i + 1]).keys())[0]))

				start_path_relator = self.find_manifold_generator_in_tree(node_path, edge_path)
				end_path_relator = [-n for n in start_path_relator[::-1]]
				relator = start_path_relator + relator + end_path_relator
				relators.append(relator)
		return relators

	def surface_relations(self):
		"""
		Get the relations which would be in a presentation of the fundamental group of the surface (when used with the
		generators from the fundamental_group_generators function
		"""
		generators, tree = self.fundamental_group_generators(return_tree=True)

		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
		T_regina = regina.Triangulation3(self.manifold)
		relators = []

		for normal_disc in self.polygons_list:
			for i, edge_index in enumerate(self.intersecting_edges(normal_disc)):
				relator = []  # contains indices of faces that correspond to a single relation
				# the edges on the tetrahedron that intersect our normal disc as embeddings
				edge_embedding = self.intersecting_edges(normal_disc, True)[i]
				# actual snappy edge
				edge = T.Edges[edge_index]
				arrow = edge.get_arrow()
				while True:
					if arrow.Edge == snappy.snap.t3mlite.simplex.bitmap(edge_embedding) and arrow.Tetrahedron.Index == normal_disc.tetrahedron:
						break
					else:
						arrow = arrow.next()
				current_disc = normal_disc
				current_arrow = arrow

				for n in range(edge.valence()):
					# we look at which face we are gluing our normal disc across
					# we take the index of this face with respect to the tetrahedron that it lies in
					# corresponds to information stored in 'arrow'
					face_index_in_tetrahedron = snappy.snap.t3mlite.simplex.FaceIndex[arrow.Face]  # is a decimal, not binary!
					face = T_regina.tetrahedron(arrow.Tetrahedron.Index).triangle(face_index_in_tetrahedron)
					face_index = face.index()
					# find index of tetrahedron on which our normal disc lies in
					tetrahedron_index = arrow.Tetrahedron.Index

					# stuff needed to find the next disc
					current_disc_face_index = [face.index() for face in current_disc.faces]
					next_index = current_disc_face_index.index(face_index)
					next_disc = current_disc.adjacent_discs[next_index]

					# find generator
					gen = self.surface_generator_of_edge(current_disc.get_id_numbers(), next_disc.get_id_numbers(), face_index)

					if gen != 0:
						relator.append(gen)
					if n == edge.valence() - 1:
						assert next_disc == normal_disc  # the last disc does not glue back to the first disc
					current_arrow = current_arrow.next()
					current_disc = next_disc
				relators.append(relator)
		return relators

	def surface_generator_of_edge(self, initial, end, face):
		if isinstance(initial, Polygon):
			initial = initial.get_id_numbers()
		if isinstance(end, Polygon):
			end = end.get_id_numbers()
		generators, tree = self.fundamental_group_generators(True)
		for index, gen in enumerate(generators):
			if initial == gen[0] and end == gen[1] and face == gen[2]:
				return index + 1
			elif initial == gen[1] and end == gen[0] and face == gen[2]:
				return -(index + 1)
		for gen in tree.edges(keys=True):
			if initial == gen[0] and end == gen[1] and face == gen[2]:
				return 0
			elif initial == gen[1] and end == gen[0] and face == gen[2]:
				return 0
		raise RuntimeError("An edge was given which was not in the surface")

	def simplify_representation(self):
		'''
		TO-DO(?): write a method to simplify representation of surface given its generators and relations
		(may have to use magma?)
		'''
		gen = self.fundamental_group_embedding()
		relation = self.surface_relations()
		F = FreeGroup(len(gen))
		rel_in_G = [F(rel) for rel in relation]
		G = F/rel_in_G
		iso = G.simplification_isomorphism()

		#goal: find simplified_group = []
		pass

	def plot_limit_set(self, name=None):
		gens = self.fundamental_group_embedding()
		gens_matrix = [Tietze_to_matrix(gen, self.manifold) for gen in gens]
		gens_excludeI = []
		I = matrix.identity(CC, 2)
		for gen in gens_matrix:
			if (gen - I).norm() > 0.01:
				gens_excludeI.append(gen)
		points_real = []
		points_complex = []
		for i in range(10000):
			pt = vector(CC, [1, 0])
			n = random.randint(1000, 2000)
			for k in range(n):
				mat = random.choice(gens_excludeI)
				pt = mat * pt
			points_real.append((pt[0] / pt[1]).real())
			points_complex.append((pt[0] / pt[1]).imag())

		fig, ax = plt.subplots()
		ax.plot(points_real, points_complex, 'bo', markersize=0.5)
		ax.set_xticks(np.linspace(min(points_real), max(points_real), 10))
		ax.set_yticks(np.linspace(min(points_complex), max(points_complex), 10))
		if name is None:
			fig.savefig('limit_set')
		elif isinstance(name, str):
			fig.savefig(name)
		else:
			raise TypeError('Name must be a string')



class Polygon:
	def __init__(self, manifold):
		# Stored in edge order based on regina's weird permuation system
		self.adjacent_discs = []  # are our Polygon classes
		# The indices of the tetrahedron that is adjacent to this one across the edge in edges
		self.adjacent_tets = []
		# The actual faces that the edges lie on in the given triangulation
		self.faces = []
		# The index of the tetrahedron that this normal disk sits inside
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

	def __eq__(self, other):
		if isinstance(other, Polygon):
			if self.get_id_numbers() == other.get_id_numbers():
				return True
			else:
				return False
		else:
			try:
				return other.__eq__(self)
			except TypeError as r:
				raise TypeError(f"Comparison between Polygon and {type(other)} not allowed")


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
	gens_traces = [gen.trace().imag().abs() for gen in gens]
	if max(gens_traces) < 0.1:
		good_gens = []
		identity = matrix.identity(gens[0].base_ring(), 2)
		for mat in gens:
			if (mat - identity).norm('frob') >= 10**(-6):
				good_gens.append(mat)
		return preserves_hermitian_form(good_gens)[0], len(good_gens)
	else:
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

def find_surfaces(vertex_surfaces, g):
    '''
    vertex_surfaces: list of regina normal surfaces
    '''
    NS = vertex_surfaces
    NS_nscomplex = [surfaces.NormalSurface(S, i) for i, S in enumerate(NS) if S.eulerChar() < 0]
    all_faces = faces.admissible_faces(NS_nscomplex)
    print(len(all_faces))
    all_surfaces = []
    for AF in all_faces:
        all_surfaces = all_surfaces + AF.surfaces_genus_in_interior(g)
    return all_surfaces

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


def main2():
	# import doubling
	# T = doubling.read_triangulation_info('example_manifold_2nd_last.txt')
	# Tr = T.regina_triangulation()
	# Tr.idealToFinite()
	# Tr.intelligentSimplify()
	# Tr = doubling.double_regina_triangulation(Tr)
	# Tr.finiteToIdeal()
	# Tr.intelligentSimplify()
	# M = snappy.Manifold(Tr.snapPea())
	M = snappy.Manifold(b"pickle: \x16\x02\x00\t\x01\x01\x01l\xd8\xb4l\x00\x01\x01\x01\x03\xcc\xff\x01\x01\xff\x01\xff\x00\x00\x03Z\xff\x01\xff\x01\x01\xff\x00\x00\x0f\x00\x00\x00ll\xd8\xb4\x00\x01\x01\x010\xca\xff\x01\x01\xff\xff\x01\x00\x00\x00\xc6\xff\x01\x01\xff\x00\x00\x05\x08\x0c\x0el\x93\xd8'\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x04\x0c\x04\x089\xe1\xb4\xe1\x01\x01\x01\x01\x03\x00\x01\xff\x00\x00\x00\n\x01\xff\x00\x00\x07\x03\x05\x03l\x93\xe1\xb4\x01\x01\x01\x010\x00\x01\xff\x00\x00\x00\x00\x00\x00\x02\x06\x04\nl\xd8\xe1\x1e\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x06\x06\x05\x0e\x8dr\xd8\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x04\r\x0f\x15l\xe1\x93\x8d\x01\x01\x01\x01\x00\x90\xff\x01\x00\x00\x00\x00\x00\x00\x02\x12\n\x039\x93\x1e\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\xc0\xff\x01\x00\x00\x00\n\r\x0bl9l'\x01\x01\x01\x01\t\x00\xff\x01\x00\x003\x00\xff\x01\x01\xff\x00\x00\x05\x08\t\x14KK\x93\x1e\x01\x01\x01\x01\x00\x00\x00\x00`\x0c\x01\xff\x01\xff\x00\x00\t\x0c\x12\x11'l\xb49\x01\x01\x01\x01\x03\xc0\xff\x01\xff\x01\x00\x00\x00\x00\x00\x00\x03\x02\r\x0b\xe1\xd8\xd8l\x01\x01\x01\x01\t\x00\x01\xff\x00\x00\x00P\xff\x01\x00\x00\x07\x0c\t\x13\xe1\xd8ll\x01\x01\x01\x01\x00\x00\x00\x00\n\x00\xff\x01\x00\x00\x02\x0f\x10\x06'l\xd8\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x01\x07\x11\x0el9\xd8l\x01\x01\x01\x01\x00P\x01\xff\x00\x00\x00\x00\x00\x00\x11\x0e\x11\x15r\xd8\xb4\x1e\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x0b\x0f\x10\x10\x93\xd8\x8d\xb4\x01\x01\x01\x01\x03\x00\xff\x01\x00\x00\x00\x00\x00\x00\x08\x13\x15\x0b9l\x1e\xb4\x01\x01\x01\x01\x00\xc0\xff\x01\x00\x00\x00\x00\x00\x00\x14\r\x14\x12'l\x93l\x01\x01\x01\x01\x00\x00\x00\x00\x03\x00\x01\xff\x00\x00\n\x13\x15\x13K9\xb4'\x01\x01\x01\x01\x00\x00\x00\x00\t\x00\x01\xff\x00\x00\x10\x12\x07\x14KKr\xb4\x01\x01\x01\x01`\x00\x01\xff\x00\x00\x00\x00\x00\x00Regina_Triangulation")
	print('M:', M.fundamental_group(simplify_presentation=False))
	print(M.num_tetrahedra())
	# M.randomize()
	# M.simplify()
	surfaces = regina.NormalSurfaces(regina.Triangulation3(M), regina.NS_QUAD_CLOSED, algHints=regina.NS_VERTEX_DD)
	print('num surfaces', len(surfaces))
	surfaces_vec = []
	for S in surfaces:
		V = S.vector()
		V_int = [int(n.stringValue()) for n in V]
		surfaces_vec.append(V_int)
	with open("surfaces_vertex.pickle", "wb") as file:
		pickle.dump(surfaces_vec, file)
	surfaces_gen2 = find_surfaces(surfaces, 2)
	print('num genus 2 surfaces', len(surfaces_gen2))
	surfaces_gen2_vec = []
	for S in surfaces_gen2:
		V = S.vector()
		V_int = [int(n.stringValue()) for n in V]
		surfaces_vec.append(V_int)
	with open("surfaces_all.pickle", "wb") as file:
		pickle.dump(surfaces_gen2_vec, file)

		# print(surface_group_in_PSL2R(surface, M))
		# gens_matrix = [Tietze_to_matrix(gen, M).trace().imag().abs() for gen in gens]
		# if max(gens_matrix) < 0.1:
		#    print(surface)
		#    for gen in gens:
		# 	   print(Tietze_to_matrix(gen, M))
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


def main3():
	M = snappy.Manifold('4_1')
	G = M.fundamental_group(simplify_presentation=False)
	Tr = regina.Triangulation3(M)

	# print('meridian:', G.meridian())
	# print('longitude:', G.longitude())
	boundary_surface = from_regina_normal_surface(boundary_to_surface(regina.Triangulation3(M)), M)

	gens = boundary_surface.fundamental_group_embedding()
	relations = boundary_surface.relations()
	gens_matrix = [Tietze_to_matrix(gen, M) for gen in gens]
	gens_string = [Tietze_to_string(gen) for gen in gens]
	# print(len(boundary_surface.polygons_list))

	# for gen in boundary_surface.fundamental_group_generators():
	# 	print(boundary_surface.surface_generator_of_edge(*gen))

	# print(boundary_surface.polygons_list)

	for rel in relations:
		print(Tietze_to_string(rel))
		print(Tietze_to_matrix(rel, M))
	# for i in range(len(gens)):
	# 	for j in range(len(gens)):
	# 		if i == j:
	# 			break
	# 		else:
	# 			print(gens_matrix[i] * gens_matrix[j] * gens_matrix[i].inverse() * gens_matrix[j].inverse())

def main4():
	M = snappy.Manifold('4_1')
	G = M.fundamental_group(simplify_presentation=False)
	Tr = regina.Triangulation3(M)

	# print('meridian:', G.meridian())
	# print('longitude:', G.longitude())
	boundary_surface = from_regina_normal_surface(boundary_to_surface(regina.Triangulation3(M)), M)
	DSS = regina.DiscSetSurface(boundary_surface.surface)

	gens, tree = boundary_surface.fundamental_group_generators(True)
	# print(gens)
	# for i, gen in enumerate(gens):
	# 	print(i+1)
	# 	print(gen)
	# 	print(Tr.faces(2)[gen[2]])
	# 	print()
	#
	# print(tree.edges(keys=True))
	# print(Tr.faces(2))
	#
	# for polygon in boundary_surface.polygons_list:
	# 	print(polygon)
	# 	print(regina.triDiscArcs[polygon.disc_type])
	# 	for perm in regina.triDiscArcs[polygon.disc_type]:
	# 		print(DSS.adjacentDisc(regina.DiscSpec(*polygon.get_id_numbers()), perm), end=', ')
	# 	print('\n')
	print(gens)
	print(boundary_surface.surface_relations())


def main5():
	# M = snappy.Manifold(b"pickle: \x16\x02\x00\t\x01\x01\x01l\xd8\xb4l\x00\x01\x01\x01\x03\xcc\xff\x01\x01\xff\x01\xff\x00\x00\x03Z\xff\x01\xff\x01\x01\xff\x00\x00\x0f\x00\x00\x00ll\xd8\xb4\x00\x01\x01\x010\xca\xff\x01\x01\xff\xff\x01\x00\x00\x00\xc6\xff\x01\x01\xff\x00\x00\x05\x08\x0c\x0el\x93\xd8'\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x04\x0c\x04\x089\xe1\xb4\xe1\x01\x01\x01\x01\x03\x00\x01\xff\x00\x00\x00\n\x01\xff\x00\x00\x07\x03\x05\x03l\x93\xe1\xb4\x01\x01\x01\x010\x00\x01\xff\x00\x00\x00\x00\x00\x00\x02\x06\x04\nl\xd8\xe1\x1e\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x06\x06\x05\x0e\x8dr\xd8\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x04\r\x0f\x15l\xe1\x93\x8d\x01\x01\x01\x01\x00\x90\xff\x01\x00\x00\x00\x00\x00\x00\x02\x12\n\x039\x93\x1e\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\xc0\xff\x01\x00\x00\x00\n\r\x0bl9l'\x01\x01\x01\x01\t\x00\xff\x01\x00\x003\x00\xff\x01\x01\xff\x00\x00\x05\x08\t\x14KK\x93\x1e\x01\x01\x01\x01\x00\x00\x00\x00`\x0c\x01\xff\x01\xff\x00\x00\t\x0c\x12\x11'l\xb49\x01\x01\x01\x01\x03\xc0\xff\x01\xff\x01\x00\x00\x00\x00\x00\x00\x03\x02\r\x0b\xe1\xd8\xd8l\x01\x01\x01\x01\t\x00\x01\xff\x00\x00\x00P\xff\x01\x00\x00\x07\x0c\t\x13\xe1\xd8ll\x01\x01\x01\x01\x00\x00\x00\x00\n\x00\xff\x01\x00\x00\x02\x0f\x10\x06'l\xd8\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x01\x07\x11\x0el9\xd8l\x01\x01\x01\x01\x00P\x01\xff\x00\x00\x00\x00\x00\x00\x11\x0e\x11\x15r\xd8\xb4\x1e\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x0b\x0f\x10\x10\x93\xd8\x8d\xb4\x01\x01\x01\x01\x03\x00\xff\x01\x00\x00\x00\x00\x00\x00\x08\x13\x15\x0b9l\x1e\xb4\x01\x01\x01\x01\x00\xc0\xff\x01\x00\x00\x00\x00\x00\x00\x14\r\x14\x12'l\x93l\x01\x01\x01\x01\x00\x00\x00\x00\x03\x00\x01\xff\x00\x00\n\x13\x15\x13K9\xb4'\x01\x01\x01\x01\x00\x00\x00\x00\t\x00\x01\xff\x00\x00\x10\x12\x07\x14KKr\xb4\x01\x01\x01\x01`\x00\x01\xff\x00\x00\x00\x00\x00\x00Regina_Triangulation")
	# MR = regina.Triangulation3(M)
	# with open("surfaces_vertex.pickle", "rb") as file:
	# 	surface_list = pickle.load(file)
	# print(len(surface_list))
	# regina_list = []
	# for S in surface_list:
	# 	SR = regina.NormalSurface(MR, regina.NS_STANDARD, regina.VectorLarge(S))
	# 	regina_list.append(SR)
	M = snappy.Manifold('K12a1')
	MR = regina.Triangulation3(M)
	print('num tet', MR.size())
	regina_list = regina.NormalSurfaces(MR, regina.NS_QUAD_CLOSED, regina.NS_FUNDAMENTAL)
	print(len(regina_list))
	surfaces_real = []
	surfaces_complex = []
	for surface in regina_list:
		all_real = True
		our_surface = from_regina_normal_surface(surface, M)
		gens = our_surface.fundamental_group_embedding()
		gens_matrix = [Tietze_to_matrix(gen, M) for gen in gens]
		comb = combinations(list(range(len(gens))), 3)
		for c in list(comb):
			gen = gens_matrix[c[0]] * gens_matrix[c[1]] * gens_matrix[c[2]]
			if gen.trace().imag().abs() > 0.001:
				surfaces_complex.append(surface)
				all_real = False
				print('complex', surface)
				break
		if all_real:
			surfaces_real.append(surface)
			break
	print("real trace:", surfaces_real)
	print("at least one complex trace:", surfaces_complex)


def main6():
	M = snappy.Manifold(b"pickle: \x16\x02\x00\t\x01\x01\x01l\xd8\xb4l\x00\x01\x01\x01\x03\xcc\xff\x01\x01\xff\x01\xff\x00\x00\x03Z\xff\x01\xff\x01\x01\xff\x00\x00\x0f\x00\x00\x00ll\xd8\xb4\x00\x01\x01\x010\xca\xff\x01\x01\xff\xff\x01\x00\x00\x00\xc6\xff\x01\x01\xff\x00\x00\x05\x08\x0c\x0el\x93\xd8'\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x04\x0c\x04\x089\xe1\xb4\xe1\x01\x01\x01\x01\x03\x00\x01\xff\x00\x00\x00\n\x01\xff\x00\x00\x07\x03\x05\x03l\x93\xe1\xb4\x01\x01\x01\x010\x00\x01\xff\x00\x00\x00\x00\x00\x00\x02\x06\x04\nl\xd8\xe1\x1e\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x06\x06\x05\x0e\x8dr\xd8\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x04\r\x0f\x15l\xe1\x93\x8d\x01\x01\x01\x01\x00\x90\xff\x01\x00\x00\x00\x00\x00\x00\x02\x12\n\x039\x93\x1e\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\xc0\xff\x01\x00\x00\x00\n\r\x0bl9l'\x01\x01\x01\x01\t\x00\xff\x01\x00\x003\x00\xff\x01\x01\xff\x00\x00\x05\x08\t\x14KK\x93\x1e\x01\x01\x01\x01\x00\x00\x00\x00`\x0c\x01\xff\x01\xff\x00\x00\t\x0c\x12\x11'l\xb49\x01\x01\x01\x01\x03\xc0\xff\x01\xff\x01\x00\x00\x00\x00\x00\x00\x03\x02\r\x0b\xe1\xd8\xd8l\x01\x01\x01\x01\t\x00\x01\xff\x00\x00\x00P\xff\x01\x00\x00\x07\x0c\t\x13\xe1\xd8ll\x01\x01\x01\x01\x00\x00\x00\x00\n\x00\xff\x01\x00\x00\x02\x0f\x10\x06'l\xd8\xe1\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x01\x07\x11\x0el9\xd8l\x01\x01\x01\x01\x00P\x01\xff\x00\x00\x00\x00\x00\x00\x11\x0e\x11\x15r\xd8\xb4\x1e\x01\x01\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x0b\x0f\x10\x10\x93\xd8\x8d\xb4\x01\x01\x01\x01\x03\x00\xff\x01\x00\x00\x00\x00\x00\x00\x08\x13\x15\x0b9l\x1e\xb4\x01\x01\x01\x01\x00\xc0\xff\x01\x00\x00\x00\x00\x00\x00\x14\r\x14\x12'l\x93l\x01\x01\x01\x01\x00\x00\x00\x00\x03\x00\x01\xff\x00\x00\n\x13\x15\x13K9\xb4'\x01\x01\x01\x01\x00\x00\x00\x00\t\x00\x01\xff\x00\x00\x10\x12\x07\x14KKr\xb4\x01\x01\x01\x01`\x00\x01\xff\x00\x00\x00\x00\x00\x00Regina_Triangulation")
	MR = regina.Triangulation3(M)
	with open("surfaces_vertex.pickle", "rb") as file:
		surface_list = pickle.load(file)
	n = 28
	SR = regina.NormalSurface(MR, regina.NS_STANDARD, regina.VectorLarge(surface_list[n]))
	S = from_regina_normal_surface(SR, M)
	name = 'limit set-surface' + str(n)
	S.plot_limit_set(name)
	print('surface', n, SR.eulerChar())

		# print('surface', n)
		# print('euler char:', SR.eulerChar())
		# print('compressible:', not SR.isIncompressible())
		# print('connected:', SR.isConnected())
		# print('orientable:', SR.isOrientable())
		# print('embedded:', SR.embedded())
	    # index of sfces with euler char < -2: [15, 17, 21, 22, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 85, 86, 87, 88, 89, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265]

if __name__ == '__main__':
	main6()
