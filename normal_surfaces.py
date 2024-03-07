import nscomplex_tg
import regina
import snappy
import pickle
import random
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from sage.all import block_matrix, matrix, vector, CC, FreeGroup, MatrixSpace, PermutationGroup
from complex_reps import preserves_hermitian_form
from nscomplex_tg import faces, regina_util, surfaces
from itertools import combinations

class NormalSurface:
	"""
	A class used to contain information and do computations for a normal surface

	Construct the boundary torus of the exterior of the figure 8 knot as a NormalSurface
	>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))

	The fundamental group of S
	>>> S.sage_group()
	Finitely presented group < x0, x3 | x3^-1*x0*x3*x0^-1 >

	The vector corresponding to S
	>>> S.get_vector()
	(1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
	"""
	def __init__(self, surface, manifold):
		"""
		Creates a normal surface from a Regina normal surface and Snappy manifold
		"""
		self.polygons_list = []
		# A list of lists of lists, where the innermost list contains all of the discs which are
		# of the specific tetrahedron and type. i.e. polygons[2][3] is a list of all the discs in
		# the normal surface which are in tetrahedron 2 and of type 3.
		self.polygons = [[[] for i in range(7)] for j in range(manifold.num_tetrahedra())]
		self.surface = surface
		self.manifold = manifold
		self.basepoint = None
		self.edges = None
		self.tree = None

	def get_vector(self):
		"""
		Returns the normal coordinate of the normal surface as a tuple
		"""
		V = self.surface.vector()
		V_int = [int(n.stringValue()) for n in V]
		surface_vec = tuple(V_int)
		return surface_vec

	def add_disc(self, disc):
		"""
		Adds the specified disc to this normal surface
		"""
		self.polygons_list.append(disc)
		self.polygons[disc.tetrahedron][disc.disc_type].append(disc)

	def get_polygon(self, tet_number, disc_type, disc_index):
		"""
		Gets the specified normal disc from the tetrahedron number, disc type, and the index of the specific disc
		(this follows Regina's normal disc indexing conventions)
		"""
		return self.polygons[tet_number][disc_type][disc_index]

	def dual_graph(self):
		"""
		Given a NormalSurface object from this package, returns the dual graph as a networkx MultiDiGraph object
		The vertices correspond to normal discs in the normal surface and are indexed by the id numbers (see Polygon.get_id_numbers for more information)
		The edges correspond to edges of the normal surface triangulation and are indexed by the face of the tetrahedron
		the edge passes through as well as the index of the triangle it passes through in the triangulation of the
		3-manifold

		This additionally chooses a basepoint for the fundamental group of the normal surface and ensures that the basepoint
		has a directed path from the basepoint to every other normal disc.

		Construct the dual graph of the boundary torus of the exterior of the figure 8 knot
		>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
		>>> G = S.dual_graph()
		>>> G.nodes()
		NodeView(((0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 3, 0), (1, 0, 0), (1, 1, 0), (1, 2, 0), (1, 3, 0)))
		>>> G.edges(keys=True)
		OutMultiEdgeView([((0, 0, 0), (1, 2, 0), (3, 0)), ((0, 1, 0), (1, 1, 0), (0, 3)), ((0, 3, 0), (1, 3, 0), (2, 1)), ((1, 0, 0), (0, 2, 0), (3, 2)), ((1, 0, 0), (0, 1, 0), (1, 0)), ((1, 0, 0), (0, 1, 0), (2, 1)), ((1, 1, 0), (0, 0, 0), (2, 1)), ((1, 1, 0), (0, 0, 0), (3, 2)), ((1, 2, 0), (0, 3, 0), (3, 2)), ((1, 2, 0), (0, 3, 0), (0, 3)), ((1, 3, 0), (0, 2, 0), (0, 3)), ((1, 3, 0), (0, 2, 0), (1, 0))])

		Shortest path from the basepoint to (0, 3, 0)
		>>> nx.shortest_path(G, (1, 0 ,0), (0, 3, 0))
		[(1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 0), (1, 2, 0), (0, 3, 0)]
		"""
		G = nx.MultiDiGraph()
		DSS = regina.DiscSetSurface(self.surface)
		T = regina.Triangulation3(self.manifold)

		# Choose a basepoint so that the basepoint of the normal surface lies in the same tetrahedron as the basepoint
		# Snappy has chosen for the manifold (if possible).
		# If there are no discs inside this tetrahedron just choose the first normal disc in the list of discs of the normal surface.
		self.manifold._choose_generators(False, False)
		gen_info = self.manifold._choose_generators_info()
		basepoint_tet = -1
		if self.basepoint is None:
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

		# Use a breadth first search starting from the basepoint disc to construct the graph by filling in normal discs
		# step-by-step until all the edges between normal discs have been added
		for polygon in self.polygons_list:
			G.add_node(polygon.get_id_numbers())
		finished_nodes = []
		unvisited_nodes = [self.basepoint.get_id_numbers()]
		while len(unvisited_nodes) != 0:
			polygon = self.get_polygon(*unvisited_nodes.pop())
			finished_nodes.append(polygon.get_id_numbers())
			if polygon.is_triangle():
				arc_list = regina.triDiscArcs[polygon.disc_type]
			else:
				arc_list = regina.quadDiscArcs[polygon.disc_type - 4]

			for arc in arc_list:
				tet_face_number = arc[3]
				adjacent_disc, adjacent_perm = DSS.adjacentDisc(regina.DiscSpec(*polygon.get_id_numbers()), arc)

				# The information for the current edge we are trying to add
				tail = polygon.get_id_numbers()
				head = id_numbers_from_DiscSpec(adjacent_disc)
				key = tet_face_number, T.tetrahedron(polygon.tetrahedron).triangle(tet_face_number).index()

				# The information for the same edge but oriented in the opposite direction
				opp_tail = id_numbers_from_DiscSpec(adjacent_disc)
				opp_head = polygon.get_id_numbers()
				opp_key = adjacent_perm[3], T.tetrahedron(adjacent_disc.tetIndex).triangle(adjacent_perm[3]).index()

				# Checking if the opposite edge is in the graph already so we don't have duplicates
				if not G.has_edge(opp_tail, opp_head, key=opp_key):
					G.add_edge(tail, head, key=key)
				if head not in unvisited_nodes and head not in finished_nodes:
					unvisited_nodes.append(head)
		return G

	def fundamental_group_generators(self, return_tree=False):
		"""
		Retrieves the edges that lie outside a maximal tree of the dual graph of the normal surface.
		These edges correspond to a set of generators of the fundamental group of the normal surface.
		This function optionally returns the maximal tree if return_tree is True.

		Generators of the boundary torus of the exterior of the figure 8 knot
		>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
		>>> G = S.dual_graph()
		>>> S.fundamental_group_generators()
		[((1, 0, 0), (0, 1, 0), (2, 1)), ((1, 1, 0), (0, 0, 0), (3, 2)), ((1, 2, 0), (0, 3, 0), (0, 3)), ((1, 3, 0), (0, 2, 0), (0, 3)), ((1, 3, 0), (0, 2, 0), (1, 0))]

		The edges corresponding to the generators and the maximal tree make up the entire dual graph
		>>> generators, T = S.fundamental_group_generators(return_tree=True)
		>>> all_edges = set(G.edges(keys=True))
		>>> all_edges == set(generators).union(set(T.edges(keys=True)))
		True
		"""
		if self.edges is None:
			digraph = self.dual_graph()
			tree_graph = nx.minimum_spanning_arborescence(digraph, preserve_attrs=True)
			edges = list(tree_graph.edges())
			for a, b in edges:
				tree_graph.remove_edge(a,b)
				tree_graph.add_edge(a, b, list(dict(digraph[a][b]).keys())[0])
			diff = nx.difference(digraph, tree_graph)
			self.edges = list(diff.edges(keys=True))
			self.tree = tree_graph

		if return_tree:
			return self.edges, self.tree
		return self.edges

	def fundamental_group_embedding(self):
		"""
		Finds a set of generators of the fundamental group of the normal surface written in terms of the generators of the fundamental group
		of the manifold it lies in. These generators are given as numbers that come from the unreduced presentation of the fundamental group
		of the manifold computed by Snappy.

		Generators of the boundary torus of the exterior of the figure 8 knot written in terms of the generators of the manifold
		>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
		>>> S.fundamental_group_embedding()
		[[2, -1], [1, 3, -2, -1], [1, 2, -1, -3, 1, -2, -1], [1, 2, -1, 3, -2, -3], [1, 2, -1, 3, -2, 1, -3]]
		"""
		self.manifold._choose_generators(False, False)
		gen_info = self.manifold._choose_generators_info()

		generators, tree = self.fundamental_group_generators(True)
		# Add paths from/to the basepoint to the given edge (note that every generator corresponds to a single edge outside the tree)
		# and creates a loop in the fundamental group corresponding to that edge
		gens_in_M = []
		for edge in generators:

			# Create the path from the basepoint to the tail of the given edge
			path_to_tail = nx.shortest_path(tree, self.basepoint.get_id_numbers(), edge[0])
			path_to_head = nx.shortest_path(tree, self.basepoint.get_id_numbers(), edge[1])
			gen_path = []
			first_half_path = list(path_to_tail)
			for i, tail in enumerate(first_half_path[:-1]):
				head = first_half_path[i + 1]
				key = list(tree.get_edge_data(tail, head).keys())[0][0]
				tetrahedron = tail[0]
				gen = gen_info[tetrahedron]['generators'][key]
				if gen != 0:
					gen_path.append(gen)

			# Add the given edge
			edge_gen = gen_info[edge[0][0]]['generators'][edge[2][0]]
			if edge_gen != 0:
				gen_path.append(edge_gen)

			# Add the path back to the basepoint from the head of the given edge
			gen_path_to_tail = []
			for i, tail in enumerate(path_to_head[:-1]):
				head = path_to_head[i + 1]
				key = list(tree.get_edge_data(tail, head).keys())[0][0]
				tetrahedron = tail[0]
				# We're gonna reverse this later so that it actually becomes the path back
				gen = -gen_info[tetrahedron]['generators'][key]
				if gen != 0:
					gen_path_to_tail.append(gen)
			gen_path = gen_path + gen_path_to_tail[::-1]

			gens_in_M.append(gen_path)
		return gens_in_M


	def intersecting_edges(self, normal_disc, return_vertex_pairs=False):
		"""
		Given a normal disc (Polygon object) returns the list of the indices of the edges inside the Snappy triangulation that the normal disc intersects.
		Note that the indices for the edges in the Regina and Snappy triangulations do NOT match.
		If return_vertex_pairs is set to True then it returns the edge on the tetrahedron as a list of pairs of vertices (0, 1, 2, 3) that are
		its endpoints.
		Note further that the vertex pairs of the edges in Snappy and Regina DO match up.
		"""
		Tr = self.surface.triangulation()
		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
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
		if return_vertex_pairs:
			return edge_list
		edge_indices = []
		for e in edge_list:
			# We've designed this code so that it returns the indices of the edges that Snappy uses.
			# However, we have included the corresponding line for the Regina edge indices as a comment below.
			# edge_indices.append(Tr.tetrahedron(tet_num).edge(*e).index()) # REGINA ONLY
			edge_indices.append(T.Tetrahedra[tet_num].Class[snappy.snap.t3mlite.bitmap(e)].Index)
		return edge_indices

	def relations(self, surface_relations = True):
		"""
		Returns the relations of the presentation of the fundamental group of this surface.
		If surface_relations is True, this returns the relations in terms of the generators of the surface corresponding
		to the edges outside the dual spanning tree.
		If surface_relations is False, this returns the relations written in terms of the generators of the fundamental
		group of the surrounding Snappy manifold.

		Relations of the boundary torus of the exterior of the figure 8 knot
		>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
		>>> S.relations()
		[[-3, -2], [-1, -5, 4], [1, 2], [-4, 3, 5]]

		Relations written in terms of the generators of the manifold
		>>> S.relations(surface_relations=False)
		[[1, 2, -1, 3, 1, -3, 2, -2, -1], [1, -2, 3, -1, -3, 1, -1], [1, 2, -2, -1, 2, 3, -2, -1], [3, 2, -3, -2, 1, -3]]

		The abelianization of the fundamental group of this surface
		>>> FG = FreeGroup(len(S.fundamental_group_generators()))
		>>> G = FG / [FG(rel) for rel in S.relations()]
		>>> G.abelian_invariants()
		(0, 0)
		"""
		all_relations = []  # list of all relations that will be returned

		start_discs = []

		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
		Tr = regina.Triangulation3(self.manifold)
		DSS = regina.DiscSetSurface(self.surface)

		self.manifold._choose_generators(True, True)
		info = self.manifold._choose_generators_info()

		self.fundamental_group_generators(return_tree=True)

		num_edges = len(Tr.edges())

		swap_last = regina.Perm4(0, 1, 3, 2)  # switches last two digits of a permutation when multiplied on the right
		swap_both = regina.Perm4(1, 0, 3, 2)  # switches first two digits and last two digits of a permutation when multiplied on the right,
											  # e.g. swaps (1, 3, 2, 0) to (3, 1, 0, 2)

		# edge_disc_list is a list of lists which maps an edge index in the triangulation to pairs (disc, disc corner)
		# (and by disc corner we mean the corner on the disc). The disc corner is recorded as the unique edge of the
		# tetrahedron (given as a pair of endpoints) which touches that specific corner.
		# From these endpoints (and information about the normal disc), we can get the arcs on the normal disc which
		# would be glued across as you go along getting the relation.
		edge_disc_list = [[] for i in range(num_edges)]

		for disc in self.polygons_list:
			tet = disc.tetrahedron
			edge_pairs = self.intersecting_edges(disc, return_vertex_pairs=True)
			edge_indices = [Tr.tetrahedron(tet).edge(*edge).index() for edge in edge_pairs]
			for i, edge_index in enumerate(edge_indices):
				edge_disc_list[edge_index].append((disc, set(edge_pairs[i])))

		# Find the relations
		for edge in range(num_edges):
			edge_valence = Tr.edge(edge).degree()
			disc_list = edge_disc_list[edge]
			while len(disc_list) > 0:
				disc, corner = disc_list[0]
				start_discs.append(disc)
				# Find the two arcs on normal disc that intersect the given edge of the triangulation,
				# start_arc will be the one where we find the next disc.
				# At the end of finding a cycle of discs we check whether the last disc glues back to the end_arc
				start_arc = None
				end_arc = None
				if disc.is_triangle():
					arc_list = regina.triDiscArcs[disc.disc_type]
					for arc in arc_list:
						if {arc[0], arc[1]} == corner:
							start_arc = arc
							end_arc = start_arc * swap_last
							break
				elif disc.is_quad():
					arc_list = regina.quadDiscArcs[disc.disc_type - 4]
					for arc in arc_list:
						if {arc[0], arc[1]} == corner:
							start_arc = arc
							end_arc = start_arc * swap_both
							break

				current_arc = start_arc
				current_disc = disc
				relation = []
				for i in range(edge_valence):
					next_disc, next_arc = DSS.adjacentDisc(regina.DiscSpec(*current_disc.get_id_numbers()), current_arc)
					next_our_disc = self.polygons[next_disc.tetIndex][next_disc.type][next_disc.number]  # our Polygon class instead of a Regina one
					next_corner = {next_arc[0], next_arc[1]}

					if surface_relations:
						# get relations in terms of the fundamental group of the surface
						regina_tet = Tr.tetrahedron(next_disc.tetIndex)
						gen = self.surface_generator_of_edge(current_disc, next_our_disc, (current_arc[3], next_arc[3]))
						if gen != 0:
							relation.append(gen)
					else:
						# get relations in terms of the fundamental group of the manifold
						gen = info[current_disc.tetrahedron]['generators'][current_arc[3]]
						if gen != 0:
							relation.append(gen)

					disc_list.remove((next_our_disc, next_corner))

					current_disc = next_our_disc
					if current_disc.is_triangle():
						current_arc = next_arc * swap_last
					elif current_disc.is_quad():
						current_arc = next_arc * swap_both
				all_relations.append(relation)
				assert current_disc == disc  # make sure we come back to the disc, arc that we started with
				assert current_arc == start_arc

		# If we are getting the relations in terms of the generators of the fundamental group of the manifold,
		# we need to add a path from/to the basepoint of the manifold to/from the disc where we start the relation.
		if not surface_relations:
			all_unbased_relations = all_relations
			all_relations = []
			gens, tree = self.fundamental_group_generators(return_tree=True)
			for i, relation in enumerate(all_unbased_relations):
				start_disc = start_discs[i]
				node_path = nx.shortest_path(tree, self.basepoint.get_id_numbers(), start_disc.get_id_numbers())
				gen_path = []
				for i in range(len(node_path) - 1):
					tail = node_path[i]
					head = node_path[i+1]
					key = list(tree.get_edge_data(tail, head).keys())[0][0]
					tetrahedron = tail[0]
					gen = info[tetrahedron]['generators'][key]
					if gen != 0:
						gen_path.append(gen)
				relation = gen_path + relation + [-num for num in gen_path[::-1]]
				all_relations.append(relation)
		return all_relations

	def simplified_generators(self, surface_generators = True):
		"""
		Finds a simplified set of generators of the fundamental group of the normal surface written in terms of the generators of the fundamental group
		of the manifold it lies in. These generators are given as numbers that come from the unreduced presentation of the fundamental group
		of the manifold computed by Snappy.

		Simplified version of the generators of the fundamental group of boundary torus of the exterior of the figure 8 knot
		>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
		>>> S.simplified_generators()
		[[4], [5]]

		Generators above written as a loop in terms of the fundamental group of the manifold
		>>> S.simplified_generators(surface_generators=False)
		[[1, 2, -1, 3, -2, -3], [1, 2, -1, 3, -2, 1, -3]]
		"""
		simpG = self.regina_group()
		iso = simpG.intelligentSimplify()  # from the unsimplified group to the simplified one
		inv_iso = regina.HomGroupPresentation(iso)
		inv_iso.invert()  # from the simplified group to the unsimplified one
		gens_list = []
		for i in range(simpG.countGenerators()):
			gens_list.append(regina_term_to_Tietze(inv_iso.evaluate(i)))
		if surface_generators:
			return gens_list
		else:
			return [self.convert_to_word_embedding(gen) for gen in gens_list]

	def regina_group(self):
		"""
		Returns the fundamental group of this surface as a Regina GroupPresentation object.

		The fundamental group of boundary torus of the exterior of the figure 8 knot as a Regina GroupPresentation object
		>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
		>>> S.regina_group()
		<regina.GroupPresentation: < a b c d e | c^-1 b^-1, a^-1 e^-1 d, a b, d^-1 c e >>
		"""
		num_gens = len(self.fundamental_group_generators(return_tree=False))
		relations = self.relations()
		if num_gens < 27:
			G = regina.GroupPresentation(num_gens, [Tietze_to_string(rel) for rel in relations])
		else:
			G = regina.GroupPresentation(num_gens, [Tietze_to_long_string(rel) for rel in relations])
		return G

	def convert_to_word_embedding(self, word):
		"""
		Given a word in the fundamental group of the surface (represented as a list of numbers, the so-called Tietze list representation),
		gets the same word as an element of the fundamental group of the surrounding manifold by applying the homomorphism induced from the
		embedding of the surface into the manifold.
		For example, the element aBAC would be represented as [1, -2, -1, -3].
		"""
		embedding = self.fundamental_group_embedding()
		embedded_word = []
		for elt in word:
			if elt > 0:
				embedded_word += embedding[elt - 1]
			else:
				embedded_word += [-num for num in embedding[-elt - 1][::-1]]
		return embedded_word

	def sage_group(self, simplified = True):
		"""
		Returns a presentation of the fundamental group of this surface as a Sagemath finitely presented group object
		By default, returns a simplified presentation of the fundamental group, if simplified is False, returns the presentation
		coming from the generators and relations returned by fundamental_group_generators and surface_relations.

		The fundamental group of boundary torus of the exterior of the figure 8 knot as a Sagemath group object
		>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
		>>> S.sage_group()
		Finitely presented group < x0, x3 | x3^-1*x0*x3*x0^-1 >

		Unsimplified version of the above group
		>>> S.sage_group(simplified=False)
		Finitely presented group < x0, x1, x2, x3, x4 | x2^-1*x1^-1, x0^-1*x4^-1*x3, x0*x1, x3^-1*x2*x4 >
		"""
		generators = self.fundamental_group_generators()
		relations = self.relations(True)
		F = FreeGroup(len(generators))
		sage_relations = [F(relation) for relation in relations]
		G = F/sage_relations
		if simplified:
			return G.simplification_isomorphism().codomain()
		else:
			return G

	def surface_generator_of_edge(self, initial, end, faces):
		"""
		Given two normal discs and the face of the initial tetrahedron which they are glued across, finds the index of an edge in the list of edges
		of the normal surface that count as generators of its fundamental group.
		The normal discs can either be given as Polygon instances or by their ID numbers.
		The face should be given as its index inside the triangulation.
		"""
		if isinstance(initial, Polygon):
			initial = initial.get_id_numbers()
		if isinstance(end, Polygon):
			end = end.get_id_numbers()
		this_face, other_face = faces
		generators, tree = self.fundamental_group_generators(True)

		for index, gen in enumerate(generators):
			if initial == gen[0] and end == gen[1] and this_face == gen[2][0]:
				return index + 1
			elif initial == gen[1] and end == gen[0] and other_face == gen[2][0]:
				return -(index + 1)
		for gen in tree.edges(keys=True):
			if initial == gen[0] and end == gen[1] and this_face == gen[2][0]:
				return 0
			elif initial == gen[1] and end == gen[0] and other_face == gen[2][0]:
				return 0
		raise RuntimeError("An edge was given which was not in the surface")

	def plot_limit_set(self, name=None, num_points = 10000):
		"""
		Plots the limit set of the normal surface and saves it with the given file name.
		If a file name is not given the plot will be saved as 'limit_set.png'.
		By default, 10,000 randomly chosen points are plotted. This number can be changed by setting the num_points argument.
		In case the generators of the normal surface are all trivial we print the statement 'all generators = I' and no plot is saved.
		"""
		gens = self.fundamental_group_embedding()
		G = self.manifold.fundamental_group(simplify_presentation=False)
		gens_matrix = [Tietze_to_matrix(gen, G) for gen in gens]
		gens_matrix = gens_matrix + [mat.inverse() for mat in gens_matrix]
		gens_excludeI = []
		I = matrix.identity(CC, 2)
		for gen in gens_matrix:
			if (gen - I).norm() > 0.01:
				gens_excludeI.append(gen)

		# checks for the case where generators of the normal surface are all trivial
		if len(gens_excludeI) == 0:
			print('all generators = I')
		else:
			points_real = []
			points_complex = []
			for i in range(num_points):
				pt = vector(CC, [1, 0])
				n = random.randint(1000, 2000)
				for k in range(n):
					mat = random.choice(gens_excludeI)
					pt = mat * pt
				points_real.append((pt[0] / pt[1]).real())
				points_complex.append((pt[0] / pt[1]).imag())

			fig, ax = plt.subplots()
			ax.plot(points_real, points_complex, 'bo', markersize=0.5)
			ax.set_aspect('equal', 'box')
			ax.set_xticks(np.linspace(min(points_real), max(points_real), 5))
			ax.set_yticks(np.linspace(min(points_complex), max(points_complex), 5))
			if name is None:
				fig.savefig('limit_set')
			elif isinstance(name, str):
				fig.savefig(name)
			else:
				raise TypeError('Name must be a string')


class Polygon:
	"""
	A class used to store information about normal discs.
	"""
	def __init__(self, manifold):
		# Note that most of the indexing here follows the conventions from Regina, but this object stores a reference
		# to a Snappy manifold (which are not always completely compatible).

		# Stored in edge order based on Regina's permuation system
		self.adjacent_discs = []  # are our Polygon classes
		# The indices of the tetrahedron that is adjacent to this one across the edge in edges
		self.adjacent_tets = []
		# The actual faces that the edges lie on in the given triangulation
		self.faces = []
		# The index of the tetrahedron that this normal disk sits inside
		self.tetrahedron = None
		# The (Snappy) manifold that this normal disc lies in
		self.manifold = manifold
		# The disc type (a number from 0-6 inclusive, no octs or tubes!)
		self.disc_type = None
		# The index of the disc (there can be multiple discs which have the same type and tet number)
		self.disc_index = None

	def get_id_numbers(self):
		"""
		Regina indexes normal discs as tuples of three numbers, the first entry is the index of the tetrahedron it lies in,
		the second is the disc type, the third is the index of the disc among the same types.
		This function allows us to access this information for NormalSurface instances.
		"""
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


class Triangle(Polygon):
	"""
	A subclass of Polygon which is a triangle.
	"""
	def __init__(self, manifold):
		super().__init__(manifold)
		self.adjacent_discs = [None] * 3
		self.adjacent_tets = [None] * 3
		self.faces = [None] * 3

	def is_triangle(self):
		return True


class Quad(Polygon):
	"""
	A subclass of Polygon which is a quad.
	"""
	def __init__(self, manifold):
		super().__init__(manifold)
		self.adjacent_discs = [None] * 4
		self.adjacent_tets = [None] * 4
		self.faces = [None] * 4

	def is_quad(self):
		return True



def id_numbers_from_DiscSpec(disc_spec):
	"""
	Retrieves the ID numbers (index of tetrahedron, disc type, and index of disc among same types) of a Regina DiscSpec object.
	"""
	return disc_spec.tetIndex, disc_spec.type, disc_spec.number

def vec_to_NormalSurface(vector, M, coord=regina.NS_STANDARD):
	"""
	Returns a NormalSurface object from its vector, given as a list of integers, and the Snappy manifold it lies in.

	Construct the exterior of the figure 8 knot as a NormalSurface then display its polygons.
	>>> S = vec_to_NormalSurface([1, 1, 1, 1, 0, 0, 0]*2, snappy.Manifold('4_1'))
	>>> S.polygons_list
	[tri:0.0.0, tri:0.1.0, tri:0.2.0, tri:0.3.0, tri:1.0.0, tri:1.1.0, tri:1.2.0, tri:1.3.0]
	"""
	MR = regina.Triangulation3(M)
	SR = regina.NormalSurface(MR, coord, regina.VectorLarge(vector))
	S = from_regina_normal_surface(SR, M)
	return S

def from_regina_normal_surface(surface, manifold):
	"""
	Given a Regina normal surface and a Snappy manifold, returns a NormalSurface object.
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
			glue_triangles(DSS, disc, face_gluings, face_list, our_surface)
		elif disc.is_quad():
			glue_quads(DSS, disc, face_gluings, face_list, our_surface)
	return our_surface

def glue_quads(DSS, disc, face_gluings, face_list, our_surface):
	"""
	A helper function used to construct gluings of quads in NormalSurface objects.
	"""
	discSpec = regina.DiscSpec(disc.tetrahedron, disc.disc_type, disc.disc_index)
	for edge_index, perm in enumerate(regina.quadDiscArcs[disc.disc_type - 4]):
		adj_disc, other_perm = DSS.adjacentDisc(discSpec, perm)
		face_index = perm[3]
		for i, embeddings in enumerate(face_gluings):
			for embedding in embeddings:
				if embedding.face() == face_index and embedding.simplex().index() == disc.tetrahedron:
					disc.faces[edge_index] = face_list[i]
					break
		disc.adjacent_discs[edge_index] = our_surface.get_polygon(adj_disc.tetIndex, adj_disc.type, adj_disc.number)
		disc.adjacent_tets[edge_index] = adj_disc.tetIndex

def glue_triangles(DSS, disc, face_gluings, face_list, our_surface):
	"""
	A helper function used to construct gluings of triangles in NormalSurface objects.
	"""
	discSpec = regina.DiscSpec(disc.tetrahedron, disc.disc_type, disc.disc_index)
	for edge_index, perm in enumerate(regina.triDiscArcs[disc.disc_type]):
		adj_disc, other_perm = DSS.adjacentDisc(discSpec, perm)
		face_index = perm[3]
		for i, embeddings in enumerate(face_gluings):
			for embedding in embeddings:
				if embedding.face() == face_index and embedding.simplex().index() == disc.tetrahedron:
					disc.faces[edge_index] = face_list[i]
					break
		disc.adjacent_discs[edge_index] = our_surface.get_polygon(adj_disc.tetIndex, adj_disc.type, adj_disc.number)
		disc.adjacent_tets[edge_index] = adj_disc.tetIndex

def Tietze_to_string(word):
	"""
	From the Tietze list of a word in a group, returns the Snappy-like string associated to it.
	Note: this will not work for more than 26 generators, use Tietze_to_long_string instead.
	>>> Tietze_to_string([1,-3,4,2,-1])
	'aCdbA'
	"""
	alphabet = 'abcdefghijklmnopqrstuvwxyz'
	our_alphabet = '?' + alphabet + alphabet.upper()[::-1]
	return ''.join([our_alphabet[index] for index in word])

def Tietze_to_long_string(elt, regina_conventions=True):
	"""
	From the Tietze list of a word in a group, returns the Regina-like string associated to it by default.
	Note: unlike Tietze_to_string this works for any number of generators.
	>>> Tietze_to_long_string([1,-3,4,2,-1])
	'g0g2^-1g3g1g0^-1'

	For a Tietze list of a word in a group with more than 26 generators this function can be used to find the Snappy-like string
	by setting the argument regina_conventions to True.
	Note: in the case with less than 26 generators, use Tietze_to_string
	>>> Tietze_to_long_string([1,-27,10,35,-2], regina_conventions=False)
	'x1X27x10x35X2'
	"""
	word = ''
	if regina_conventions:
		# Using Regina conventions
		for num in elt:
			if num > 0:
				word += 'g' + str(num - 1)
			else:
				word += 'g' + str(abs(num) - 1) + '^-1'
	else:
		# Using Snappy conventions
		for num in elt:
			if num > 0:
				word += 'x' + str(num)
			else:
				word += 'X' + str(abs(num))
	return word

def Tietze_to_matrix(word, G):
	"""
	Given a word in the Snappy fundamental group G (coming from a call to fundamental_group on a manifold M), gives the
	SL2(C) matrix corresponding to it.

	Matrix corresponding to the word [1, -2] in the fundamental group of the exterior of the figure 8 knot
	>>> G = snappy.Manifold('4_1').fundamental_group()
	>>> Tietze_to_matrix([1, -2], G)
	[ 6.66000000000000e-16 + 3.46410161513775*I    -0.866025403784438 - 2.50000000000000*I]
	[    2.59807621135331 + 0.499999999999999*I -2.00000000000000 - 6.66000000000000e-16*I]
	"""
	return G.SL2C(Tietze_to_string(word))

def regina_term_to_Tietze(elt):
	"""
	Given a Regina GroupExpression, for example, g0g1^-1 returns the Tietze representation of the term as a list, in this
	case, [1, -2].
	"""
	Tietze = []
	for term in elt.terms():
		gen = term.generator
		exp = term.exponent
		if exp > 0:
			Tietze.extend([gen + 1] * exp)
		else:
			Tietze.extend([-(gen + 1)] * abs(exp))
	return Tietze

def get_exact_holonomy_matrices(M, size = 40, try_higher = True):
	"""
	Given a manifold M, this returns the exact form of the holonomy matrices.
	The size parameter defines an upper bound for the degree of the field extension that the entries of the holonomy matrices live in.
	If try_higher is True, this will increase the size parameter (by doubling the current size repeatedly) until it finds the right field.
	If try_higher is False and the field is not found, this returns None.
	"""
	# Gets the field information (the field as an extension of the rationals and what each matrix entry is in that field)
	group_field_info = None
	while group_field_info is None:
		entries = M.holonomy_matrix_entries([False, True, True, True])
		group_field_info = entries.find_field(10000, size)
		if not try_higher and group_field_info is None:
			return None
		size = size * 2

	F = group_field_info[0]
	entries = group_field_info[2]
	# Snappy returns the matrix entries as a list, so we reformat this as a list of matrices
	mat_space = MatrixSpace(F, 2)
	assert len(entries) % 4 == 0
	num_mats = len(entries) // 4
	exact_matrices = []
	for i in range(num_mats):
		mat = mat_space(group_field_info[2][4 * i:4 * (i + 1)])
		exact_matrices.append(mat)
	return exact_matrices

def get_exact_generators(M, S_vec):
	"""
	Given a Snappy manifold and vector of a normal surface (as a list of integers) returns the exact holonomy matrices of the fundamental
	group of the surface.
	Prints information about the returned matrices when run.
	"""
	S = vec_to_NormalSurface(S_vec, M)
	gens = S.fundamental_group_embedding()
	matrices = get_exact_holonomy_matrices(M)
	print('matrix generators for the manifold group')
	for mat in matrices:
		print(mat)
	field = matrices[0].base_ring()
	print('field information')
	print(field)
	surface_matrices = []
	for gen in gens:
		mat = matrix.identity(field, 2)
		for num in gen:
			if num < 0:
				mat = mat * matrices[-num - 1].inverse()
			else:
				mat = mat * matrices[num - 1]
		surface_matrices.append(mat)
	print('matrix generators for the surface')
	for mat in surface_matrices:
		print(mat)
		print('trace')
		print(mat.trace())
	return surface_matrices

def get_surface_trace_field(M, S_vec):
	"""
	Given a Snappy manifold and vector of a normal surface (as a list of integers) returns the exact trace field of the fundamental
	group of the surface.
	"""
	mats = get_exact_generators(M, S_vec)
	F = mats[0].base_ring()
	N = len(mats)
	non_invariant_gens = [mat.trace() for mat in mats]
	for i in range(N):
		for j in range(i, N):
			mat = mats[i]*mats[j]
			non_invariant_gens.append(mat.trace())
	for i in range(N):
		for j in range(i, N):
			for k in range(j, N):
				mat = mats[i]*mats[j]*mats[k]
				non_invariant_gens.append(mat.trace())
	trace_field = F.subfield_from_elements(non_invariant_gens)
	return trace_field

def is_obviously_a_surface_group(G):
	"""
	Tests whether the given presentation is transparently one of the fundamental group of an orientable surface.
	"""
	# Code provided by Nathan Dunfield

	# First, we check there is only one relation
	# and that every generator g appears exactly
	# once as g and once as g^-1.
	n, rels = G.ngens(), G.relations()
	if len(rels) > 1:
		return False
	R = rels[0].Tietze()
	if sorted(R) != list(range(-n, 0)) + list(range(1, n+1)):
		return False

	# Now we make sure that we take the boundary
	# of the relator and identify sides accordingly
	# the result is actually a surface.
	perms = []
	for g in range(1, n+1):
		a0, b0 = R.index(g)+1, R.index(-g)+1
		a1 = a0 + 1 if a0 < 2*n else 1
		b1 = b0 + 1 if b0 < 2*n else 1
		perms += [ [(a0, b1)], [(a1, b0)] ]

	return len(PermutationGroup(perms).orbits()) == 1

def find_surfaces(vertex_surfaces, e):
	"""
	Given a list of Regina normal surfaces that are vertex surfaces and a certain Euler characteristic, returns all normal surfaces of
	the given Euler characteristic (not genus).
	"""
	NS = vertex_surfaces
	NS_nscomplex = [surfaces.NormalSurface(S, i) for i, S in enumerate(NS) if S.eulerChar() < 0]
	all_faces = faces.admissible_faces(NS_nscomplex)
	all_surfaces = []
	for AF in all_faces:
		all_surfaces = all_surfaces + AF.surfaces_euler_in_interior(e)
	return all_surfaces