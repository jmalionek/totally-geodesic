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

	def get_vector(self):
		"""
		Returns the normal coordinate of the normal surface as a tuple
		"""
		V = self.surface.vector()
		V_int = [int(n.stringValue()) for n in V]
		surface_vec = tuple(V_int)
		return surface_vec

	def add_disc(self, disc):
		self.polygons_list.append(disc)
		self.polygons[disc.tetrahedron][disc.disc_type].append(disc)

	def get_polygon(self, tet_number, disc_type, disc_index):
		return self.polygons[tet_number][disc_type][disc_index]

	def dual_graph(self):
		"""
		Given a NormalSurface object from this package, returns the dual graph as a networkx MultiDiGraph object
		The vertices correspond to normal discs in the normal surface and are indexed by the id numbers (see Polygon.get_id_numbers for more information)
		The edges correspond to edges of the normal surface triangulation and are indexed by the faces of the manifold triangulation on which they lie
		"""
		G = nx.MultiDiGraph()
		for polygon in self.polygons_list:
			G.add_node(polygon.get_id_numbers())
		for polygon in self.polygons_list:
			for index, adj_disc in enumerate(polygon.adjacent_discs):
				# TODO: SKETCHY since the face indices are not the same for regina and snappy
				if not G.has_edge(adj_disc.get_id_numbers(), polygon.get_id_numbers(), key=polygon.faces[index].index()):
					G.add_edge(polygon.get_id_numbers(), adj_disc.get_id_numbers(), key=polygon.faces[index].index())
		return G

	def fundamental_group_generators(self, return_tree=False):
		"""
		Retrieves the edges that lie outside a maximal tree of the dual graph of the normal surface.
		These edges correspond to a set of generators of the fundamental group of the normal surface.
		This function optionally returns the maximal tree also if return_tree is True.
		"""
		undirected_dual_graph = self.dual_graph().to_undirected()
		tree_graph = nx.minimum_spanning_tree(undirected_dual_graph)
		diff = nx.difference(undirected_dual_graph, tree_graph)
		nx.draw(tree_graph)
		if return_tree:
			return list(diff.edges(keys=True)), tree_graph
		return diff.edges

	def fundamental_group_embedding(self):
		"""
		Finds a set of generators of the fundamental group of the normal surface written in terms of the generators of the fundamental group
		of the manifold it lies in. These generators are given as numbers that come from the unreduced presentation of the fundamental group
		of the manifold computed by Snappy.
		"""
		self.manifold._choose_generators(False, False)
		gen_info = self.manifold._choose_generators_info()
		basepoint_tet = -1
		# choose normal disc that will be used as the basepoint of the fundamental group of the normal surface
		# try to choose normal disc inside basepoint tetrahedron of manifold, if there are no discs inside this tetrahedron
		# just choose first normal disc in list of discs of normal surface
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
		verbose = True
		if verbose:
			print('basepoint')
			print(self.basepoint)
		cycles = []
		node_paths = []
		generators, tree = self.fundamental_group_generators(True)
		# add paths from/to basepoint to given edge (generator corresponds to a single edge outside the tree) and creates a loop in the fundamental group
		for edge in generators:
			if verbose:
				print(edge)
			path_to_edge = nx.shortest_path(tree, basepoint_disc.get_id_numbers(), edge[0])
			path_from_edge = nx.shortest_path(tree, edge[1], basepoint_disc.get_id_numbers())
			edge_path = []
			for i in range(len(path_to_edge) - 1):
				edge_path.append((path_to_edge[i], path_to_edge[i+1], list(tree.get_edge_data(path_to_edge[i], path_to_edge[i+1]).keys())[0]))
			if verbose:
				print('edge_path_to_generator')
				print(edge_path)
			edge_path.append(edge)
			for i in range(len(path_from_edge) - 1):
				edge_path.append((path_from_edge[i], path_from_edge[i+1], list(tree.get_edge_data(path_from_edge[i], path_from_edge[i+1]).keys())[0]))
			cycles.append(edge_path)
			if verbose:
				print('edge_path')
				print(edge_path)
			node_paths.append(path_to_edge + path_from_edge)
			if verbose:
				print('node_path')
				print(path_to_edge + path_from_edge)
		gens_in_M = []
		for i, cycle in enumerate(cycles):
			gens_in_M.append(self.find_manifold_generator_in_tree(node_paths[i], cycle))
		return gens_in_M

	def find_manifold_generator(self, edge):
		"""
		Given an edge in the dual graph of the normal surface returns the corresponding generator in the fundamental group of the manifold.
		This generator is oriented following the orientation on the edge.
		"""
		T = self.surface.triangulation()
		self.manifold._choose_generators(True, True)
		info = self.manifold._choose_generators_info()

		starting_tet = self.get_polygon(*edge[0]).tetrahedron
		ending_tet = self.get_polygon(*edge[1]).tetrahedron
		starting_disc = self.get_polygon(*edge[0])
		ending_disc = self.get_polygon(*edge[1])
		# index = starting_disc.adjacent_discs.index(ending_disc)
		triangle_number = edge[2]
		regina_triangle = T.triangle(triangle_number)


		if starting_disc.is_triangle():
			disc_type = starting_disc.disc_type
			face_index = regina.triDiscArcs[disc_type][index][3]  # last digit of permutation gives the index of the face inside the tetrahedron (0 ~ 3)
		else:
			disc_type = starting_disc.disc_type - 4
			face_index = regina.quadDiscArcs[disc_type][index][3]  # last digit of permutation gives the index of the face inside the tetrahedron (0 ~ 3)
		gen = info[starting_tet]['generators'][face_index]
		return gen

	def find_manifold_generator_in_tree(self, node_path, edge_path):
		"""
		Given a path written as a list of consecutive edges and consecutive vertices returns a list of the corresponding generators
		in the fundamental group of the manifold.
		These generators are oriented consistently with the orientation of the nodes on the path.
		"""
		# Given a path in a graph retrieves the list of edges and list of nodes, then compares these two lists to find the generator
		# with the right orientation
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

	def intersecting_edges(self, normal_disc, return_vertex_pairs=False):
		"""
		Given a normal disc (Polygon object) returns the list of the indices of the edges inside the snappy triangulation that the normal disc intersects.
		If return_vertex_pairs is set to True then it returns the edge on the tetrahedron as a list of pairs of vertices (0, 1, 2, 3) that are
		its endpoints.
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
			# TODO: SOURCE OF SKETCHINESS. BELOW IS ONLY FOR REGINA. FIX TO RETURN SNAPPY STUFF?
			# edge_indices.append(Tr.tetrahedron(tet_num).edge(*e).index()) # REGINA ONLY
			# Below should work for snappy now
			edge_indices.append(T.Tetrahedra[tet_num].Class[snappy.snap.t3mlite.bitmap(e)].Index)
		return edge_indices

	def relations(self):
		"""
		Computes the relations of the fundamental group of the surface written in terms of the generators of the manifold fundamental group.
		Erring on the side of caution the result includes a lot of relations that are redundant.
		"""
		# TODO: APPARENTLY, snappy and regina do not always have the same information coming from the same manifold???
		# Look at the output of this code on the surface [1,1,1,1,0,0,0]*3 for the knot 5_2.
		# Known disparities: The index of the triangles and edges in snappy and regina do not line up (again, see 5_2)
		# TODO: Because of this, we need to examine every place where we use a correspondence between regina and the snappy Mcomplex/Triangulation

		# DONE: We wrote code which checks if the tetrahedra have the same internal ordering on the vertices and ran it on a ton of things
		# def check_the_same(M):
		# 	Mr = regina.Triangulation3(M)
		# 	T = Mcomplex(M)
		# 	for i in range(M.num_tetrahedra()):
		# 		tet_snappy = T.Tetrahedra[i]
		# 		tet_regina = Mr.tetrahedron(i)
		# 		for face_index in range(4):
		# 			snappy_perm = tet_snappy.Gluing[15 - 2 ** face_index]
		# 			regina_perm = tet_regina.adjacentGluing(face_index)
		# 			regina_perm = [regina_perm[i] for i in range(4)]
		# 			print(snappy_perm, regina_perm)


		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
		T_regina = regina.Triangulation3(self.manifold)
		# We need this to determine the basepoint
		self.fundamental_group_embedding()
		relators = []

		self.manifold._choose_generators(True, True)
		info = self.manifold._choose_generators_info()

		for normal_disc in self.polygons_list:
			for i, edge_index in enumerate(self.intersecting_edges(normal_disc)):

				# will contain indices of faces in the triangulation (where the normal discs are glued across) that correspond to a single relation
				relator = []
				# the edge on the tetrahedron that intersect our normal disc as a pair of endpoint vertices
				edge_pair = self.intersecting_edges(normal_disc, True)[i]

				# actual snappy edge
				edge = T.Edges[edge_index]

				# loops over the triangles glued to our edge until it finds the right arrow corresponding to our normal disc
				# needed because snappy has an internal ordering on these arrows
				arrow = edge.get_arrow()

				while True:
					if arrow.Edge == snappy.snap.t3mlite.simplex.bitmap(edge_pair) and arrow.Tetrahedron.Index == normal_disc.tetrahedron:
						break
					else:
						arrow = arrow.next()
				current_disc = normal_disc
				current_arrow = arrow
				print(edge.valence())
				for n in range(edge.valence()):
					# we look at which face we are gluing our normal disc across
					# we take the index of this face (0, 1, 2, 3) with respect to the tetrahedron that it lies in
					# corresponds to information stored in 'arrow'
					face_index_in_tetrahedron = snappy.snap.t3mlite.simplex.FaceIndex[arrow.Face]  # is a decimal, not binary!
					face = T_regina.tetrahedron(arrow.Tetrahedron.Index).triangle(face_index_in_tetrahedron)
					face_index = face.index()
					# find index of tetrahedron on which our normal disc lies in
					tetrahedron_index = arrow.Tetrahedron.Index

					print(current_disc, face_index)
					# find generator
					gen = info[tetrahedron_index]['generators'][face_index_in_tetrahedron]

					# stuff needed to find the next disc
					current_disc_face_index = [face.index() for face in current_disc.faces]
					next_index = current_disc_face_index.index(face_index)
					current_disc = current_disc.adjacent_discs[next_index]
					current_arrow = current_arrow.next()



					if gen != 0:
						relator.append(gen)

					# in case the normal disc we are at is the last one check that it glues back to the first disc
					if n == edge.valence() - 1:
						assert current_disc == normal_disc

				# find paths from/to basepoint and append to relators
				tree = self.fundamental_group_generators(return_tree=True)[1]
				node_path = nx.shortest_path(tree, self.basepoint.get_id_numbers(), normal_disc.get_id_numbers())  # list of nodes in tree
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
		Get the relations in a presentation of the fundamental group of the surface when used with the generators
		from the fundamental_group_generators function.
		TODO: Run the relations checkers on a bunch of things
		TODO: It seems like relations is working (mostly) fine, but surface_relations is bad
		TODO: Rewrite relations :( :( :( BUT ACTUALLY THIS ONE 8/31/2023
		TODO: There can be two discs which are glued to eachother across the same face MULTIPLE times (see m006 boundary torus) 8/29/2023
		TODO: Use m004 as an example to debug 8/31/2023
		"""
		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
		T_regina = regina.Triangulation3(self.manifold)
		# We need this to determine the basepoint
		self.fundamental_group_embedding()
		relators = []
		verbose = False

		for normal_disc in self.polygons_list:
			for i, edge_index in enumerate(self.intersecting_edges(normal_disc)):

				if verbose:
					print('normal_disc', normal_disc)

					print('starting edge', edge_index)
					corner = T.Edges[edge_index].Corners[0]
					print('corner', corner)


				# will contain indices of faces in the triangulation (where the normal discs are glued across) that correspond to a single relation
				relator = []
				# the edge on the tetrahedron that intersect our normal disc as a pair of endpoint vertices
				edge_embedding = self.intersecting_edges(normal_disc, True)[i]
				# actual snappy edge
				edge = T.Edges[edge_index]

				# loops over the triangles glued to our edge until it finds the right arrow corresponding to our normal disc
				# needed because snappy has an internal ordering on these arrows
				arrow = edge.get_arrow()
				while True:
					if arrow.Edge == snappy.snap.t3mlite.simplex.bitmap(edge_embedding) and arrow.Tetrahedron.Index == normal_disc.tetrahedron:
						break
					else:
						arrow = arrow.next()
				current_disc = normal_disc
				current_arrow = arrow
				if verbose:
					print('starting arrow', current_arrow)
				# print(edge.valence())
				for n in range(edge.valence()):
					# we look at which face we are gluing our normal disc across
					# we take the index of this face (0, 1, 2, 3) with respect to the tetrahedron that it lies in
					# corresponds to information stored in 'arrow'
					if verbose:
						print()
						print()
						print('current_disc', current_disc)
						print('current_disc neighbors', current_disc.adjacent_discs)
						print('current_disc adjacent faces', current_disc.faces)
						print(arrow)
					face_index_in_tetrahedron = snappy.snap.t3mlite.simplex.FaceIndex[arrow.Face]  # is a decimal, not binary!
					face = T_regina.tetrahedron(arrow.Tetrahedron.Index).triangle(face_index_in_tetrahedron)
					face_index = face.index()
					if verbose:
						print('face_index', face_index)
					# find index of tetrahedron on which our normal disc lies in
					tetrahedron_index = arrow.Tetrahedron.Index

					# stuff needed to find the next disc
					current_disc_face_index = [face.index() for face in current_disc.faces]
					if verbose:
						print(current_disc_face_index)
					next_index = current_disc_face_index.index(face_index)
					next_disc = current_disc.adjacent_discs[next_index]

					# find generator
					gen = self.surface_generator_of_edge(current_disc.get_id_numbers(), next_disc.get_id_numbers(), face_index)

					if gen != 0:
						relator.append(gen)

					# in case the normal disc we are at is the last one check that it glues back to the first disc
					if n == edge.valence() - 1:
						assert next_disc == normal_disc

					current_arrow = current_arrow.next()
					current_disc = next_disc
				relators.append(relator)
		return relators

	def relations_version_2(self, surface_relations = True):
		all_relations = []  # list of all relations that will be returned

		T = snappy.snap.t3mlite.Mcomplex(self.manifold)
		Tr = regina.Triangulation3(self.manifold)
		DSS = regina.DiscSetSurface(self.surface)

		self.manifold._choose_generators(True, True)
		info = self.manifold._choose_generators_info()

		generators, tree = self.fundamental_group_generators(return_tree=True)

		num_edges = len(Tr.edges())

		swap_last = regina.Perm4(0, 1, 3, 2)  # switches last two digits of a permutation when multiplied on the right
		swap_both = regina.Perm4(1, 0, 3, 2)  # switches first two digits and last two digits of a permutation when multiplied on the right,
											  # e.g. swaps (1, 3, 2, 0) to (3, 1, 0, 2)

		# Goes from an edge in the triangulation to pairs (disc, disc corner) (and by disc corner we mean the corner on
		# the disc). The disc corner is recorded as the unique edge of the tetrahedron (given as a pair of endpoints)
		# which, touches that specific corner.

		# From these endpoints (and information about the normal disc), we can get the arcs on the normal disc which
		# would be glued across as you go along getting the relation.
		edge_disc_list = [[] for i in range(num_edges)]

		for disc in self.polygons_list:
			tet = disc.tetrahedron
			disc_type = disc.disc_type
			edge_pairs = self.intersecting_edges(disc, return_vertex_pairs=True)
			edge_indices = [Tr.tetrahedron(tet).edge(*edge).index() for edge in edge_pairs]
			for i, edge_index in enumerate(edge_indices):
				edge_disc_list[edge_index].append((disc, set(edge_pairs[i])))

		for edge in range(num_edges):
			edge_valence = Tr.edge(edge).degree()
			disc_list = edge_disc_list[edge]
			while len(disc_list) > 0:
				# disc, corner = disc_list.pop(0)
				disc, corner = disc_list[0]
				# find the two arcs on normal disc that intersect the given edge, start_arc will be the one where we find the next disc
				# at the end of finding a cycle of discs we check whether the last disc glues back to the end_arc
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
						regina_face = regina_tet.triangle(next_arc[3])
						gen = self.surface_generator_of_edge(current_disc, next_our_disc, regina_face.index())
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
			# TODO: Deal with basepoints for manifold relations
		return all_relations

	def sage_group(self, simplified = True):
		"""
		Returns a presentation of the fundamental group of this surface as a Sagemath finitely presented group object
		By default, returns a simplified presentation of the fundamental group, if simplified is False, returns the presentation
		coming from the generators and relations returned by fundamental_group_generators and surface_relations.
		"""
		generators = self.fundamental_group_generators()
		relations = self.relations_version_2(True)
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
		for embedded_rel in relations:
			mat = Tietze_to_matrix(embedded_rel, self.manifold)
			Id = matrix.identity(CC, 2)
			if not ((Id - mat).norm() < .01 or (Id + mat).norm() < .01):
				print(embedded_rel)
				print(mat)
				raise RuntimeError('One of the relations was not plus or minus the identity')

	# Should not be included in final code
	def relations_as_holonomy_matrices(self):
		relations = self.relations()
		for relation in relations:
			mat = Tietze_to_matrix(relation, self.manifold)
			Id = matrix.identity(CC, 2)
			print(mat)
			if not ((Id - mat).norm() < .01 or (Id + mat).norm() < .01):
				print(relation)
				print(mat)
				raise RuntimeError('One of the relations was not plus or minus the identity')

	# Should not be included in final code
	def get_embedded_relations(self):
		surface_relations = self.relations_version_2(True)
		embedding = self.fundamental_group_embedding()
		embedded_relations = []
		for i in range(len(surface_relations)):

			surf_rel = surface_relations[i]
			print('new relation')
			print('surf_rel', surf_rel)
			print('embedding', embedding)
			embedded_rel = []
			for elt in surf_rel:
				if elt > 0:
					embedded_rel += embedding[elt - 1]
					print(elt, embedding[elt - 1])
				else:
					embedded_rel += [-num for num in embedding[-elt - 1][::-1]]
					print(elt, [-num for num in embedding[-elt - 1][::-1]])
			print('embedded_rel', embedded_rel)
			embedded_relations.append(embedded_rel)
		return embedded_relations

	# Should not be included in final code
	def relation_check(self):
		"""
		Writes the relations for the surface fundamental group in terms of the fundamental group of the manifold
		"""
		embedded_relations = self.get_embedded_relations()
		relations = self.relations()
		assert len(relations) == len(embedded_relations)
		F = FreeGroup(len(self.fundamental_group_generators()))
		sage_embedded_relations = [F(rel) for rel in embedded_relations]
		sage_relations = [F(rel) for rel in relations]
		return zip(sage_embedded_relations, sage_relations)


	def surface_generator_of_edge(self, initial, end, face):
		"""
		Given two normal discs and the faces along which they are glued across, finds the index of an edge in the list of edges
		of the normal surface that count as generators of its fundamental group.
		The normal discs can either be given as Polygon instances or by their ID numbers.
		The face should be given as its index inside the triangulation.
		"""
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

		return iso.codomain()

		#goal: find simplified_group = []
		pass

	def plot_limit_set(self, name=None, simplify_presentation=True, num_points = 10000):
		"""
		Plots the limit set of the normal surface and saves it with the given file name.
		If a file name is not given the plot will be saved as limit_set.png.
		"""
		gens = self.fundamental_group_embedding()
		gens_matrix = [Tietze_to_matrix(gen, self.manifold) for gen in gens]
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
				# This is the stuff to change?
	return our_surface


def Tietze_to_string(word):
	alphabet = 'abcdefghijklmnopqrstuvwxyz'
	our_alphabet = '?' + alphabet + alphabet.upper()[::-1]
	return ''.join([our_alphabet[index] for index in word])


def Tietze_to_matrix(word, M):
	return M.fundamental_group(simplify_presentation=False).SL2C(Tietze_to_string(word))


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

def find_surfaces(vertex_surfaces, g):
    '''
    vertex_surfaces: list of regina normal surfaces that are vertex surfaces
    '''
    NS = vertex_surfaces
    NS_nscomplex = [surfaces.NormalSurface(S, i) for i, S in enumerate(NS) if S.eulerChar() < 0]
    all_faces = faces.admissible_faces(NS_nscomplex)
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
		print(S.relations_version_2(True))
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
	M = snappy.Manifold('m004')
	S = vec_to_NormalSurface([1,1,1,1,0,0,0]*M.num_tetrahedra(), M)
	G = S.sage_group()
	Gsimp = G.simplification_isomorphism().codomain()
	print(is_obviously_a_surface_group(Gsimp))



if __name__ == '__main__':
	surface_group_is_fine()
