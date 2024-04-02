============================================================================
Totally Geodesic Surfaces in Hyperbolic 3-Manifolds: Algorithms and Examples
============================================================================

This is the code and data which accompany the paper "Totally Geodesic Surfaces in Hyperbolic 3-Manifolds: Algorithms
and Examples" by Brannon Basilio, Chaeryn Lee, Joseph Malionek.

Data
====

The data which accompanies the paper is contained in the zipped csv file all_data.csv.gz.
The uncompressed version of this file is around 4GB, so it is recommended that this file be opened using a terminal or command prompt (i.e. not with a graphical user interface).
To unzip, on linux and MacOS machines, the command

  gunzip all_data.csv

should uncompress the file into its full 4GB version.
On Windows machines, this file can be extracted with several external utilities including 7-Zip, WinRAR, WinZip, or Gzip.

The uncompressed version of this file is a "csv" with columns being delimited by the ; character.
This file has 9 columns containing data and a leading column which contains the row number.
These columns are:

1. manifold: This is the Regina tight encoding string which is a compact way to store information about a triangulation.
   When in a python terminal, the Regina Triangulation object can be recovered with the command
     regina.Triangulation3.tightDecoding(manifold)

2. manifold_name: The name of the manifold in Snappy's OrientableCuspedCensus or HTLinkExteriors.
   For a link, this is useful because it allows you to directly get the Snappy Manifold object with the code
     snappy.Manifold(manifold_name)
   For a cover of a cusped manifold, the manifold_name is a string created by Snappy when taking the cover of a manifold.
   It consists of the name of the base manifold followed by the type (either 'cyc' for cyclic or 'irr' for irregular) of the covering.
   An example is 't00143~irr', where the base manifold is the manifold t00143, followed by '~irr' denoting that the cover was irregular.

3. runtime_vertex_surfaces: The amount of time it took regina to find the vertex surfaces of the given manifold.

4. runtime_enumerate_surfaces: The amount of time it took nscomplex_tg to find the surfaces of the given manifold up to the correct genus bound, given the vertex surfaces.

5. runtime_tot_geo: The collective amount of time it took to determine whether each surface up to the genus bound was totally geodesic.

6. vertex_surfaces_vec: A list of all of the vertex surfaces for the manifold as found by regina.
   This is stored as a string representing a python list containing a sequence of tuples.
   In python, to obtain the list of vectors, run the command

     list_of_vecs = eval(vertex_surfaces_vec)

   Once this list has been obtained, if M is the variable containing the Snappy manifold object, and vec is an element of list_of_vecs, each surface can be instantiated with the command

     normal_surfaces.vec_to_NormalSurface(vec, M)

7. all_surfaces_vec: This is the list of all surfaces found up to the requisite Euler characteristic bound.
   This is formatted and can be accessed in python as above.

8. tot_geo: This is the list of all surfaces which were found to be totally geodesic.
   This is formatted and can be accessed in python as in 6.

9. potential_tot_geo: This is the list of all surfaces which were found to be Fuchsian but were orientable and so further analysis needs to be done to determine whether or not they are actually totally geodesic.
   (Note that for all of the surfaces found in this column in the database, all are orientable double covers of surfaces found in the tot_geo column.)
   This is formatted and can be accessed in python as in 6.

For some entries, the vertex_surfaces_vec column may end up being empty because we only decided to store this information partway through our experiments.

An example of how to open our file (after uncompressing):

   >>> import pandas as pd
   >>> df = pd.read_csv('all_data.csv', delimiter = ';')
   >>> df.loc[5]
   Unnamed: 0                                                                    5
   manifold                          -#"$"%"&"$2("'")"*"'"(0+")4+"'.,"*,,0*8,"+&,4
   manifold_name                                                          L12n1344
   runtime_vertex_surfaces                                                0.367275
   runtime_enumerate_surfaces                                             0.276709
   runtime_tot_geo                                                        0.143235
   vertex_surfaces_vec           [(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0,...
   all_surfaces_vec              [(1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0,...
   tot_geo                                                                      []
   potential_tot_geo                                                            []
   Name: 5, dtype: object