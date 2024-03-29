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

	TODO: 2/26 add doctests for the following (without a star) and continue with documentation (line 412)
	list of functions to have doctests
		get_vector *
		dual_graph (some features) *
		fundamental_group_generators *
		fundamental_group_embedding *
		relations *
		simplified_generators *
		regina_group *
		sage_group *

	"""
	def __init__(self, surface, manifold):
		"""
		Creates a normal surfaces from a regina normal surface and snappy manifold
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
		(this follows regina's normal disc indexing conventions)
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
		Note that the indices for the edges in the regina and Snappy triangulations do NOT match.
		If return_vertex_pairs is set to True then it returns the edge on the tetrahedron as a list of pairs of vertices (0, 1, 2, 3) that are
		its endpoints.
		Note further that the vertex pairs of the edges in Snappy and regina DO match up.
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
			# However, we have included the corresponding line for the regina edge indices as a comment below.
			# edge_indices.append(Tr.tetrahedron(tet_num).edge(*e).index()) # REGINA ONLY
			edge_indices.append(T.Tetrahedra[tet_num].Class[snappy.snap.t3mlite.bitmap(e)].Index)
		return edge_indices

	def relations(self, surface_relations = True):
		"""
		Returns the relations of the presentation of the fundamental group of this surface.
		If surface_relations is True, this returns the relations in terms of the generators of the surface corresponding
		to the edges outside the dual spanning tree.
		If surface_relations is False, this returns the relations written in terms of the generators of the fundamental
		group of the surrounding snappy manifold.

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
					next_our_disc = self.polygons[next_disc.tetIndex][next_disc.type][next_disc.number]  # our Polygon class instead of a regina one
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
		Returns the fundamental group of this surface as a regina GroupPresentation object.

		The fundamental group of boundary torus of the exterior of the figure 8 knot as a regina GroupPresentation object
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

	# Should not be included in final code
	def surface_relations_as_holonomy_matrices(self):
		relations = self.get_embedded_relations()
		G = self.manifold.fundamental_group(simplify_presentation=False)
		for embedded_rel in relations:
			mat = Tietze_to_matrix(embedded_rel, G)
			Id = matrix.identity(CC, 2)
			if not ((Id - mat).norm() < .01 or (Id + mat).norm() < .01):
				print(embedded_rel)
				print(mat)
				raise RuntimeError('One of the relations was not plus or minus the identity')

	# Should not be included in final code
	def relations_as_holonomy_matrices(self):
		relations = self.relations(surface_relations=False)
		G = self.manifold.fundamental_group(simplify_presentation=False)
		for relation in relations:
			mat = Tietze_to_matrix(relation, G)
			Id = matrix.identity(CC, 2)
			print(mat)
			if not ((Id - mat).norm() < .01 or (Id + mat).norm() < .01):
				print(relation)
				print(mat)
				raise RuntimeError('One of the relations was not plus or minus the identity')

	# Should not be included in final code
	def get_embedded_relations(self):
		surface_relations = self.relations(True)
		embedding = self.fundamental_group_embedding()
		embedded_relations = []
		for i in range(len(surface_relations)):

			surf_rel = surface_relations[i]
			# print('new relation')
			# print('surf_rel', surf_rel)
			# print('embedding', embedding)
			embedded_rel = []
			for elt in surf_rel:
				if elt > 0:
					embedded_rel += embedding[elt - 1]
					# print(elt, embedding[elt - 1])
				else:
					embedded_rel += [-num for num in embedding[-elt - 1][::-1]]
					# print(elt, [-num for num in embedding[-elt - 1][::-1]])
			# print('embedded_rel', embedded_rel)
			embedded_relations.append(embedded_rel)
		return embedded_relations

	# Should not be included in final code
	def relation_check(self):
		"""
		Writes the relations for the surface fundamental group in terms of the fundamental group of the manifold
		"""
		embedded_relations = self.get_embedded_relations()
		relations = self.relations(False)
		assert len(relations) == len(embedded_relations)
		F = FreeGroup(len(self.fundamental_group_generators()))
		sage_embedded_relations = [F(rel) for rel in embedded_relations]
		sage_relations = [F(rel) for rel in relations]
		return zip(sage_embedded_relations, sage_relations)


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

	def plot_limit_set(self, name=None, simplify_presentation=True, num_points = 10000):
		"""
		Plots the limit set of the normal surface and saves it with the given file name.
		If a file name is not given the plot will be saved as limit_set.png.
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

		# in case the generators of the normal surface are all trivial we print this statement and no plot is saved
		# (should be checked for totally geodesic surfaces at some point, once done this part can be removed)
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

	# Not doable as a small side project, will be a project on its own in the future (maybe?)
	def plot_finer_limit_set(self, name=None, max_dist = .1):
		graph = nx.Graph()
		graph.add_node(tuple())
		simpG = self.regina_group()
		print(simpG)
		iso = simpG.intelligentSimplify()  # from the unsimplified group to the simplified one
		inv_iso = regina.HomGroupPresentation(iso)
		inv_iso.invert()  # from the simplified group to the unsimplified one
		print(inv_iso)
		print(simpG)
		assert len(simpG.relations()) == 1
		relation = simpG.relations()[0]
		word = regina_term_to_Tietze(relation)
		assert len(word) % 2 == 0
		points = [[letter] for letter in word]
		for point in points:
			graph.add_edge(tuple(), tuple(point))
		current_max_dist = float('inf')
		# the function used to convert a Tietze representation into a string which is digestible by regina
		if simpG.countGenerators() > 26:
			Tietze_converter = Tietze_to_long_string
		else:
			Tietze_converter = Tietze_to_string
		snappyG = self.manifold.fundamental_group(simplify_presentation=False)
		origin = vector(CC, [0, 1])
		pairs_to_be_done = []
		for i, point in enumerate(points):
			# Here we're setting the assumption that we're taking the cyclic ordering to be the ordering given by
			# going along the initial word from left to right (i.e. in the word abABcdCD, the ordering goes from a to b
			# to A to B etc.
			pairs_to_be_done.append((point, points[(i+1) % len(points)]))
		while len(pairs_to_be_done) > 0:
			a_word, b_word = pairs_to_be_done.pop(0)
			# convert from the Tietze words in the simplified group to a matrix
			a = Tietze_to_matrix(self.convert_to_word_embedding(regina_term_to_Tietze(inv_iso.evaluate(regina.GroupExpression(Tietze_converter(a_word))))), snappyG)
			b = Tietze_to_matrix(self.convert_to_word_embedding(regina_term_to_Tietze(inv_iso.evaluate(regina.GroupExpression(Tietze_converter(b_word))))), snappyG)
			# Apply the matrices on a point
			pta = a * origin
			ptb = b * origin
			# Make the points be complex numbers
			pta = pta[0] / pta[1]
			ptb = ptb[0] / ptb[1]
			# If the points are close enough together we don't need to add points between
			if (pta - ptb).abs() < max_dist:
				continue
			# If the points are too far from each other, we add points between

			# Get the common suffix between the two words
			# Given a = [1,37,3], b = [5,37,3], the common_suffix_length should be 2
			length = min(len(a), len(b))
			common_suffix_length = 0
			for i in range(length):
				if a[-i - 1] != b[-i - 1]:
					common_suffix_length = i
					break
			common_suffix = a[len(a) - common_suffix_length:]
			a_prefix = a[0 : len(a) - common_suffix_length]
			b_prefix = b[0 : len(b) - common_suffix_length]
			if len(a_prefix) == 0 or len(b_prefix == 0):
				# Cases 3 or 4, they lie on the same "octagon" but don't branch
				# TODO: Case 3, there is some space between them, but on still on same octagon (Fill in points in between)
				# Make sure to know which way around the octagon to go
				# TODO: Case 4, no space between them, very bad (Enter a new octagon)
				pass
			else:
				# Check if they are on the same octagon
				a_letter = a_prefix[-1]
				b_letter = b_prefix[-1]
				a_index = word.index(a_letter)
				b_index = word.index(b_letter)
				a_backwards_index = word.index(-a_letter)
				b_backwards_index = word.index(-b_letter)
				if (a_backwards_index + 1) % len(word) == b_index or (b_backwards_index + 1) % len(word) == a_index:
					# The case where we're in the octagon
					a_prefix_in_octagon = []
					for letter in a_prefix:
						pass
					# they are in the octagon but one is not a subset of the other (fill in points in between)

				else:
					# There is no common octagon
					# This is the easiest case
					if a_index < b_index:
						letters = word[a_index + 1: b_index]
					else:
						letters = word[a_index + 1:] + word[:b_index]
					new_words = []
					for letter in letters:
						new_words.append(common_suffix + [letter])
					pair_words = [a] + new_words + [b]
					for i in range(len(pair_words) - 1):
						pairs_to_be_done.append((pair_words[i], pair_words[i + 1]))



			# if b_branch_letter is directly to the right of a_branch_letter in the circular ordering
			a_branch_index = word.index(a_branch_letter)
			b_branch_index = word.index(b_branch_letter)
			# make the point in the middle
			# if when the paths branch, the two branches are adjacent to each other
			if b_branch_index - a_branch_index in [1, -(len(word) - 1)]:
				# choose the shorter branch and turn it towards the other branch
				if len(a) < len(b):
					middle = [word[(a_branch_index + len(word) // 2 + 1) % len(word)], a_branch_letter] + common_suffix
				else:
					middle = [word[(b_branch_index + len(word) // 2 - 1) % len(word)], b_branch_letter] + common_suffix
			else:
				# if the paths branch with something between them, choose the branch in between
				if a_branch_index < b_branch_index:
					middle = [word[ceil((a_branch_index + b_branch_index) // 2)]] + common_suffix
				else:
					# we want the middle in the circular orderin which is on the outside of a and b in the list
					middle = [word[ceil((a_branch_index + b_branch_index + len(word)) / 2) % len(word)]] + common_suffix
			print('middle')
			print(middle)
			points[max_dist_i] = [word[(word.index(a[0]) + len(word) // 2) % len(word)]] + a
			points[max_dist_i + 1] = [word[(word.index(b[0]) + len(word) // 2) % len(word)]] + b
			points = points[:max_dist_i + 1] + [middle] + points[max_dist_i + 1:]


		fig, ax = plt.subplots()
		complex_points = []
		print(points)
		for point in points:
			# print(point)
			# print(Tietze_converter(point))
			# print(regina.GroupExpression(Tietze_converter(point)))
			# print(inv_iso.evaluate(regina.GroupExpression(Tietze_converter(point))))
			# print(regina_term_to_Tietze(inv_iso.evaluate(regina.GroupExpression(Tietze_converter(point)))))
			# print(regina_term_to_Tietze(inv_iso.evaluate(regina.GroupExpression(Tietze_converter(point)))))
			# print(snappyG)


			mat = Tietze_to_matrix(self.convert_to_word_embedding(regina_term_to_Tietze(inv_iso.evaluate(regina.GroupExpression(Tietze_converter(point))))),
			                     snappyG)
			projective_point = mat * origin
			complex_point = projective_point[0] / projective_point[1]
			complex_points.append(complex_point)
		points_real = []
		points_imaginary = []
		for pt in complex_points:
			points_real.append((pt[0] / pt[1]).real())
			points_imaginary.append((pt[0] / pt[1]).imag())
		ax.plot(points_real, points_imaginary, 'bo')

		if name is None:
			fig.savefig('limit_set')
		elif isinstance(name, str):
			fig.savefig(name)
		else:
			raise TypeError('Name must be a string')

		return points, complex_points

		# TODO: When running the code, we get to a point where we have two branches, a and b.
		# When we find the middle point, it ends up being a subset of one of a or b.
		# Then when we extend, we end up getting a or b exactly afterwards. This causes the same point to be added forever
		# In this process, the two points are getting slightly farther and farther away each time actually (not good)

		# Two fixes:
		# Add in some randomness?
		# Try and come up with a deterministic way to stop this from happening



class Polygon:
	def __init__(self, manifold):
		# Note that most of the indexing here follows the conventions from Regina, but this object stores a reference
		# to a snappy manifold (which are not always completely compatible).

		# Stored in edge order based on regina's permuation system
		self.adjacent_discs = []  # are our Polygon classes
		# The indices of the tetrahedron that is adjacent to this one across the edge in edges
		self.adjacent_tets = []
		# The actual faces that the edges lie on in the given triangulation
		self.faces = []
		# The index of the tetrahedron that this normal disk sits inside
		self.tetrahedron = None
		# The (snappy) manifold that this normal disc lies in
		self.manifold = manifold
		# The disc type (a number from 0-6 inclusive, no octs or tubes!)
		self.disc_type = None
		# The index of the disc (there can be multiple discs which have the same type and tet number)
		self.disc_index = None

	def get_id_numbers(self):
		"""
		Regina indexes normal discs as tuples of three numbers, the first entry is the index of the tetrahedron it lies in,
		the second is the disc type, the third is the index of the disc among the same types. This function allows us to access this information
		for NormalSurface instances.
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

def id_numbers_from_DiscSpec(disc_spec):
	return disc_spec.tetIndex, disc_spec.type, disc_spec.number


def vec_to_NormalSurface(vector, M, coord=regina.NS_STANDARD):
	'''
	vector is a list of integers, M is snappy manifold
	returns our normal surface from the vector
	'''
	MR = regina.Triangulation3(M)
	SR = regina.NormalSurface(MR, coord, regina.VectorLarge(vector))
	S = from_regina_normal_surface(SR, M)
	return S

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
			glue_triangles(DSS, disc, face_gluings, face_list, our_surface)
		elif disc.is_quad():
			glue_quads(DSS, disc, face_gluings, face_list, our_surface)
	return our_surface


def glue_quads(DSS, disc, face_gluings, face_list, our_surface):
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
	discSpec = regina.DiscSpec(disc.tetrahedron, disc.disc_type, disc.disc_index)
	for edge_index, perm in enumerate(regina.triDiscArcs[disc.disc_type]):
		adj_disc, other_perm = DSS.adjacentDisc(discSpec, perm)
		face_index = perm[3]
		for i, embeddings in enumerate(face_gluings):
			for embedding in embeddings:
				if embedding.face() == face_index and embedding.simplex().index() == disc.tetrahedron:
					disc.faces[edge_index] = face_list[i]
					break
		# WHY BEN?!
		disc.adjacent_discs[edge_index] = our_surface.get_polygon(adj_disc.tetIndex, adj_disc.type, adj_disc.number)
		disc.adjacent_tets[edge_index] = adj_disc.tetIndex


def Tietze_to_string(word):
	"""
	From the Tietze list of a word in a group, returns the snappy-like string associated to it. For example, for the
	word [1,-3,4,2,-1], this returns the string 'aCdbA'
	Note will not work for more than 26 generators
	"""
	alphabet = 'abcdefghijklmnopqrstuvwxyz'
	our_alphabet = '?' + alphabet + alphabet.upper()[::-1]
	return ''.join([our_alphabet[index] for index in word])

def Tietze_to_long_string(elt, regina_conventions=True):
	word = ''
	if regina_conventions:
		# Using regina conventions
		for num in elt:
			if num > 0:
				word += 'g' + str(num - 1)
			else:
				word += 'g' + str(abs(num) - 1) + '^-1'
	else:
		# Using snappy conventions
		for num in elt:
			if num > 0:
				word += 'x' + str(num)
			else:
				word += 'X' + str(abs(num))
	return word

def Tietze_to_matrix(word, G):
	"""
	Given a word in the snappy fundamental group G (coming from a call to fundamental_group on a manifold M), gives the
	SL2(C) matrix corresponding to it
	"""
	return G.SL2C(Tietze_to_string(word))

def regina_term_to_Tietze(elt):
	"""
	Given a regina GroupExpression, for example, g0g1^-1 returns the Tietze representation of the term as a list, in this
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
	Given a manifold M, this returns the holonomy matrices. The size parameter defines an upper bound for the degree of
	the field extension that the entries of the holonomy matrices live in. If try_higher is true, this will increase the
	size parameter (by doubling the current size repeatedly) until it finds the right field. If try_higher is false and
	the field is not found, this returns None

	TODO: Maybe add a check to make sure this is indeed a faithful representation of the manifold group?
	(Checking it's a representation should be easy, checking that it's faithful will not be easy)
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
	mat_space = MatrixSpace(F, 2)
	assert len(entries) % 4 == 0
	num_mats = len(entries) // 4
	# assert len(entries) // 4 == len(M.fundamental_group().generators())
	exact_matrices = []
	for i in range(num_mats):
		mat = mat_space(group_field_info[2][4 * i:4 * (i + 1)])
		exact_matrices.append(mat)
	return exact_matrices

def get_exact_generators(M, S_vec):
	S = vec_to_NormalSurface(S_vec, M)
	gens = S.fundamental_group_embedding()
	matrices = get_exact_holonomy_matrices(M)
	print('matrix generators for the manifold group')
	for mat in matrices:
		print(mat)
	print('field information')
	field = matrices[0].base_ring()
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
		print(mat.trace(), mat.trace())
	return surface_matrices

def get_surface_trace_field(M, S_vec):
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




def surface_generators(surface, manifold, SL2C=True):
	our_surface = from_regina_normal_surface(surface, manifold)
	generators = our_surface.fundamental_group_embedding()
	if SL2C:
		G = M.fundamental_group(simplify_presentation=False)
		return [Tietze_to_matrix(word, G) for word in generators]
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

def is_obviously_a_surface_group(G):
	"""
	Tests whether the given presentation is transparently
	one of the fundamental group of an orientable surface.
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

def find_surfaces(vertex_surfaces, e):
	'''
	vertex_surfaces: list of regina normal surfaces that are vertex surfaces
	modified to return all surfaces of a given euler characteristic (not genus)
	'''
	NS = vertex_surfaces
	NS_nscomplex = [surfaces.NormalSurface(S, i) for i, S in enumerate(NS) if S.eulerChar() < 0]
	all_faces = faces.admissible_faces(NS_nscomplex)
	all_surfaces = []
	for AF in all_faces:
		all_surfaces = all_surfaces + AF.surfaces_euler_in_interior(e)
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
	G = M.fundamental_group(simplify_presentation=False)
	gens_matrix = [Tietze_to_matrix(gen, G) for gen in gens]
	gens_string = [Tietze_to_string(gen) for gen in gens]
	# print(len(boundary_surface.polygons_list))

	# for gen in boundary_surface.fundamental_group_generators():
	# 	print(boundary_surface.surface_generator_of_edge(*gen))

	# print(boundary_surface.polygons_list)

	for rel in relations:
		print(Tietze_to_string(rel))
		G = M.fundamental_group(simplify_presentation=False)
		print(Tietze_to_matrix(rel, G))
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
		G = M.fundamental_group(simplify_presentation=False)
		gens_matrix = [Tietze_to_matrix(gen, G) for gen in gens]
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
	n = 9
	SR = regina.NormalSurface(MR, regina.NS_STANDARD, regina.VectorLarge(surface_list[n]))
	S = from_regina_normal_surface(SR, M)
	name = 'limit set-surface' + str(n)
	S.plot_limit_set(name, num_points=100_000)
	print('surface', n, SR.eulerChar())

		# print('surface', n)
		# print('euler char:', SR.eulerChar())
		# print('compressible:', not SR.isIncompressible())
		# print('connected:', SR.isConnected())
		# print('orientable:', SR.isOrientable())
		# print('embedded:', SR.embedded())
		# index of sfces with euler char < -2: [15, 17, 21, 22, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 85, 86, 87, 88, 89, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265]

def main7():
	"""
	Find all the bad manifolds (Manifolds with tetrahedra which are glued to themselves)
	"""
	bad_manifolds = []
	for M in snappy.HTLinkExteriors(alternating = False)[7.2:]:
		T = snappy.snap.t3mlite.Mcomplex(M)
		for tet in T.Tetrahedra:
			if tet in tet.Neighbor.values():
				bad_manifolds.append(M.name())
				break

	with open('bad_manifolds.pickle', 'wb') as file:
		pickle.dump(bad_manifolds, file)

def main8():
	"""
	Print exact stuff for a given manifold and normal surface
	"""
	M = regina.Triangulation3.tightDecoding(""":("*"/"3")"+"."2"\'","1"5"&"-"0"4"6"86.,7"96/,8"760,9"661,,2-862-2,87282924(325(223(522(4266768696""")
	M = snappy.Manifold(M)
	S_vec = (0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0)
	return get_surface_trace_field(M, S_vec)

def main9():
	"""
	Check our new limit set plot function
	"""
	vector = (0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
	0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
	2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0,
	2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0)
	T = regina.Triangulation3.tightDecoding('1%")"-"(*&"*"+")*\'"(","**/".,+,0"/,-,."0,,,."-*/",*0"+*.*/*0*')
	M = snappy.Manifold(T)
	S = vec_to_NormalSurface(vector, M)
	points, complex_points = S.plot_finer_limit_set('test_limit_set')
	print(points)
	print(complex_points)

def main10():
	"""
	Check simplified_generators
	"""
	M = snappy.HTLinkExteriors(alternating=False)[7.2:][120000]  # this was K15n8371(0,0) and had 18 tetrahedra
	vec = [1, 1, 1, 1, 0, 0, 0]*18
	S = vec_to_NormalSurface(vec, M)
	G = M.fundamental_group(simplify_presentation=False)
	print(G)
	gens = S.simplified_generators(False)
	gens_matrix = [Tietze_to_matrix(gen, G) for gen in gens]

def test_new_relations():
	# TODO: m412 breaks fundamental_group_embedding
	M1 = snappy.Manifold('m004')
	M2 = snappy.Manifold('m006')
	M3 = snappy.Manifold('m413')
	Ms = [M1]
	for M in Ms:
		T = regina.Triangulation3(M)
		print(M)
		# print(T.detail())
		vec = [1,1,1,1,0,0,0]*M.num_tetrahedra()
		S = vec_to_NormalSurface(vec, M)

		print('Surface relations')
		print(S.relations(True))
		for rel in S.surface_relations():
			print(rel)
		# print('Manifold relations')
		# print(S.relations_version_2(False))
		# for rel in S.relations():
		# 	print(rel)
		print(S.sage_group())
		print(S.surface_relations_as_holonomy_matrices())

def print_all_information():
	M = snappy.Manifold('m004')
	S = vec_to_NormalSurface([1,1,1,1,0,0,0]*M.num_tetrahedra(), M)
	for disc in S.polygons_list:
		print('disc', disc)
		print('adjacent faces', disc.adjacent_discs)
		print('regina face indices', disc.faces)
	T = regina.Triangulation3(M)
	print(T.detail())
	gens, tree = S.fundamental_group_generators(return_tree=True)
	print('gens')
	for gen in gens:
		print(gen)
	print('tree')
	for edge in tree.edges():
		print(edge)
	print('embedding info')
	print(S.fundamental_group_embedding())

def surface_group_is_fine():
	M = snappy.Manifold('K15n1234')
	S = vec_to_NormalSurface([1,1,1,1,0,0,0]*M.num_tetrahedra(), M)
	G = S.sage_group(False)
	gens, tree = S.fundamental_group_generators(return_tree=True)
	print('gens', gens)
	print('tree edges', tree.edges(keys= True))
	print(G)
	Gsimp = G.simplification_isomorphism().codomain()
	print(Gsimp)
	print(is_obviously_a_surface_group(Gsimp))
	for pair in S.relation_check():
		print(pair)
	print(S.surface_relations_as_holonomy_matrices())
	print(S.relations_as_holonomy_matrices())




if __name__ == '__main__':
	main10()

# TODO: Get rid of print statements
# TODO: Make limit set plot better
# TODO: Rerun all the walkofshame things
# TODO: Tidy up code
# TODO: Fork and publish somewhere
# TODO: Merge branches