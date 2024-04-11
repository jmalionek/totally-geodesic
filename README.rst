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

::

  gunzip all_data.csv

should uncompress the file into its full 4GB version.
On Windows machines, this file can be extracted with several external utilities including 7-Zip, WinRAR, WinZip, or Gzip.

The uncompressed version of this file is a "csv" with columns being delimited by the ; character.
This file has 9 columns containing data and a leading column which contains the row number.
These columns are:

1. manifold: This is the Regina tight encoding string which is a compact way to store information about a triangulation.
   When in a python terminal, the Regina Triangulation object can be recovered with the command
   ::

     regina.Triangulation3.tightDecoding(manifold)

2. manifold_name: The name of the manifold in Snappy's OrientableCuspedCensus or HTLinkExteriors.
   For a link, this is useful because it allows you to directly get the Snappy Manifold object with the code
   ::

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
   ::

     list_of_vecs = eval(vertex_surfaces_vec)

   Once this list has been obtained, if M is the variable containing the Snappy manifold object, and vec is an element of list_of_vecs, each surface can be instantiated with the command
   ::

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

Code Quick Start
================
In order to use our software, some additional packages and software are necessary.
SageMath, SnapPy, and Regina are necessary to work with this software.
In order to install these, a helpful place to look is
  https://snappy.computop.org/installing.html#kitchen-sink.
Here you can find detailed instructions on how to install SageMath, SnapPy, and Regina in such a way that they can work together.
Additionally, in order to enumerate surfaces in manifolds with multiple cusps you will specifically need Nathan Dunfield's fork of Regina.
This fork has two branches: the multicusp_closed branch and the multicusp_closed_service-7.3 branch.
You will need to clone this fork, checkout one of these branches and then compile from the source code.
Instructions for how to do this can be found at
  https://github.com/NathanDunfield/regina/tree/multicusp_closed.
The versions of the packages used for the experiments accompanying the paper were Python 3.7.8, SageMath 9.2, Regina 7.2, and SnapPy 3.1.
However, other versions have also been found to work with our code.


Example
=======

Here is an example of how our code could be used::

  sage: import regina, snappy
  sage: import normal_surfaces as ns
  sage: import pandas as pd
  sage: df = pd.read_csv('all_data.csv', delimiter = ';')
  sage: row = df.loc[int(5)] #pandas does not work well with sagemath integers
  sage: row
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
  sage: M_regina = regina.Triangulation3.tightDecoding(row.loc['manifold'])
  sage: M_snappy = snappy.Manifold(M_regina)
  sage: vs_vectors = eval(row.loc['vertex_surfaces_vec'])
  sage: example_vector = vs_vectors[0]
  sage: S = ns.vec_to_NormalSurface(example_vector, M_snappy)
  sage: S.surface.eulerChar()
  -2
  sage: S.surface.isOrientable()
  True
  sage: S.sage_group()
  Finitely presented group < x1, x2, x3, x9 | x3^-1*x9^-1*x3*x9*x2*x1^-1*x2^-1*x1 >


Code Details
============

Here we have an explanation of our code files.
Details about all of the functions here can be found in the documentation for the code.
For a specific function calling

::

  help(function_name)

in python or

::

  function_name?

in sage will give you the documentation for the function.

------------------
normal_surfaces.py
------------------

This file contains the main portion of our code.
It mainly contains the class NormalSurface which stores the corresponding Regina NormalSurface object and further information constructed from the Regina NormalSurface.
This class has methods concerning calculation of the fundamental group of the surface and how it relates to the fundamental group of the manifold.
There are also some supplementary/helper functions about converting between Regina, Snappy and our code.
Some functions which may be especially useful are vec_to_normal_surface which converts the vector for a normal surface into a NormalSurface object and find_surfaces which finds all surfaces in a given 3-manifold up to a specific Euler characteristic bound.

-----------------
detect_tot_geo.py
-----------------

This file contains the functions which are directly related to algorithms 2 and 3 in the associated paper.

-------
test.py
-------

This is a file which runs all of the doctests in normal_surfaces.py and detect_tot_geo.py adapted from similar code by Nathan Dunfield.
The file can be run with the command

::

   python test.py

------------
nscomplex_tg
------------

This is a folder of code taken from https://doi.org/10.7910/DVN/FZIHMB which is data associated to the paper "Counting essential surfaces in 3-manifolds" by Nathan Dunfield, Stavros Garoufalidis, and Hyam Rubinstein.
It has been modified by us slightly so as to extend some of its functionality.

License
=======

All code and data herein is hereby released into the public domain by
its above-named authors, as per CC0::

  https://creativecommons.org/publicdomain/zero/1.0/legalcode