import snappy
import nscomplex
M = snappy.Manifold('6_1')
# The Mcomplex triangulation fo the manifold
T = snappy.snap.t3mlite.Mcomplex(M)
cs = nscomplex.ConnectedSurfaces(M, -10)

S = cs.normal[0]
F = cs.incompressible[0]

# The triangulation of the original manifold
T = F.triangulation()
# The regina surface
RF = F.surface


#Example with totally geodesic surface(?)
[t12052(0,0)(0,0),
 10^2_104(0,0)(0,0),
 L10n36(0,0)(0,0),
 ooct02_00006(0,0)(0,0)]
 #volume = 7.3277247534
 #num_tetrahedra = 8
 
 
 #another example with totally geodesic surface
 #stored as isometry signature
-cGcvLvvvLLLvvwvAvAvMvLvLwLvwMvPwvzAQvvzvQLzPLLwAzvPAQPAAzPLzLzMMwAPwvPLAzPQQQwMQPwQzAMvQPQMMQzAMPQPQQQAQMAQQQkdagapatavaoaraEaqasaCaAazaRaXaDaEa1a+aHaYabbNaSaPa3aMa3a+aTaabtbWaXarb2aBb0a6a8aEb6aIbxbmbabmbpbFbhb-agbubgbAbfbBbFbRblbjbkbobnbMbobDbLb0bybtbOb-bxbKbQbdcfcEbCbGbQbSbLbGbNbPb1bJbKb9b7bicNbkcWb-bRbVb3b3b2bocYbpcYbrcmcacjc5b-bhcackcncjcicrchcbcdc8bwcjczcbcucxcgcBceclcBcmcpcxcCctcwcvcmcscwcycAcqcCcscxcAcvcFcFcEcDcBcDcFcEcEcofaaaoaaoooooaafoaaoaaffoofofofafoafaofoaoaaafaqafooqoqoaaafoofooofoqfofaoooaafoooffoofqfoqqaoqoqofaaqafaoaafofoffooooooffoqfafqqfaooqooqqfaaooaaaoqoaqfoooooofoo(0,0)(0,0)
#volume = 22.8293234089
 #num_tetrahedra = 26
 
 
 #another another example
 [L14n30011(0,0)(0,0)]
 #volume = 17.9499597719
 #num_tetrahedra = 19
 
 
 #anotha one
 [L12n910(0,0)(0,0)]
  #volume = 15.0686956377
 #num_tetrahedra = 17
 
 
 #exhibit A
 [10^4_8(0,0)(0,0)(0,0)(0,0),
 L10n101(0,0)(0,0)(0,0)(0,0),
 otet10_00014(0,0)(0,0)(0,0)(0,0),
 ocube02_00026(0,0)(0,0)(0,0)(0,0)]
 #volume = 10.1494160641
 #num_tetrahedra = 10
 
 
 
# I switched names lol
Sr = RF
# The bizarre way that regina stores the gluings of adjacent edges
DSS = regina.DiscSetSurface(Sr)
# arguments are tet number, disc type, disc index
discSpec = regina.DiscSpec(0, 1, 0)
# permutation specifies which arc you are looking at:
# where 0 is sent specifies which vertex of the tetrahedron is "cut off" byu the arc
# Then the next two numbers in the tetrahedron specify which edges of the tetrahedron the arc goes between
perm = regina.Perm4(1,0,2,3)
# returns the adjacent disc and the arc (given as a permutation) 
DSS.adjacentDisc(discSpec, perm)

#Look at regina.triDiscArcs for looping over the arcs of a given triangle


F = snappy.snap.fundamental_polyhedron.FundamentalPolyhedronEngine.from_manifold_and_shapes(M, M.tetrahedra_shapes('rect'))

# <Gluing info on figure8 knot complement>
# tet0    [tet1, tet1, tet1, tet1]
#         [(0, 1, 3, 2), (1, 3, 0, 2), (1, 0, 2, 3), (2, 0, 3, 1)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e0 (int)     E02 : e0 (int)     E12 : e1 (int)
#                E03 : e1 (int)     E31 : e1 (int)     E23 : e0 (int)
# tet1    [tet0, tet0, tet0, tet0]
#         [(0, 1, 3, 2), (1, 3, 0, 2), (1, 0, 2, 3), (2, 0, 3, 1)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e0 (int)     E02 : e0 (int)     E12 : e1 (int)
#                E03 : e1 (int)     E31 : e1 (int)     E23 : e0 (int)
#
# Edges:
# e0 (int)         Edge of valence 6      Endpoints [v0 (int) , v0 (int) ]
#         < E01 | F2 | tet0 >  < E01 | F3 | tet1 >  < E02 | F3 | tet0 >
#         < E23 | F0 | tet1 >  < E23 | F1 | tet0 >  < E02 | F1 | tet1 >
# e1 (int)         Edge of valence 6      Endpoints [v0 (int) , v0 (int) ]
#         < E12 | F3 | tet0 >  < E03 | F2 | tet1 >  < E31 | F0 | tet0 >
#         < E12 | F3 | tet1 >  < E03 | F2 | tet0 >  < E31 | F0 | tet1 >

# [{'index': 0,
#   'generators': (0, -3, -2, -1),
#   'neighbors': (1, 1, 1, 1),
#   'gluings': ((0, 1, 3, 2), (1, 3, 0, 2), (1, 0, 2, 3), (2, 0, 3, 1)),
#   'corners': (0.000000000000000,
#    0.000000000000000,
#    0.000000000000000,
#    0.000000000000000),
#   'generator_path': 0},
#  {'index': 1,
#   'generators': (0, 1, 2, 3),
#   'neighbors': (0, 0, 0, 0),
#   'gluings': ((0, 1, 3, 2), (1, 3, 0, 2), (1, 0, 2, 3), (2, 0, 3, 1)),
#   'corners': (0.000000000000000,
#    0.000000000000000,
#    0.000000000000000,
#    0.000000000000000),
#   'generator_path': -1}]

