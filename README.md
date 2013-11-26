About:
	MPI-SCVT is a Parallel SCVT Generator written in C++.
	It makes use of boost's mpi and serialization libraries.
	It was written by Douglas Jacobsen as part of his 
	Dissertation work while he was a Ph.D. student at Florida State University.

Compiling:
	To compile, you need the boost mpi and serialization libraries
	to be in your LD_LIBRARY_PATH, and the header files in your PATH variables.
	Afterwards, compilation is as simple as typing "make", a debug version can be compiled
	using "make debug" as well. These commands create MpiScvt.x which can be run from the
	command line.

Boost:
	Boost can be downloaded from:
	http://www.boost.org/

	The only libraries that need to be setup at mpi, and serialization. Please refer to boost's documentation
	to set up correctly.

Running:
	To run, simply type:
		./MpiScvt.x
	
	To run using mpi, simply type:
		mpirun -n P ./MpiScvt.x
	Where P is replaced with the number of processors used. The number of Processors has to be either
	1 or the number of points in RegionList as described below.

Input:
	To run MPI-SCVT 4 Files can be used. Each file has a label saying if it is required or not.

	Params: (Required)
		This file contains parameters required to run MPI-SCVT. If the file does not exist
		a run of MPI-SCVT will create this file, with default values.


	RegionList: (Required)
		This file contains a list of region centers, for parallelization. Some example lists are provided in the
		original source package. Each region center needs to be given as x, y, z with one region per line.
		One region should be provided for each processor to be used. A simple example two processor setup (the 
		smallest number MPI-SCVT is capable of using) is the north and south pole, with a RegionList file containing
		the following data:
			0 0 1
			0 0 -1

	RegionTriangulation: (Required)
		This file contains the connectivity of the regions in RegionList. Each line contains one "triangle" given as 
		three indices, where each index relates to a point number from the RegionList file. Contradictory to the name, 
		the connectivity does not need to be a Triangulation, though a Triangulation is easy to come by. As an example, 
		the connectivity for the example provided for RegionList would be:
			0 1 0

		However, a triangulation (in triangles.dat) as output from MPI-SCVT can be used as well

	SaveVertices: (Optional)
		This file contains an initial point set for iterating on. These points should be given as x, y, z
		with one point per line. This is required if the method of generating points is set to 0 in Params.

	SaveBoundaries: (Optional)
		This file lists all points which lie along some boundary segment or segments. Generator points can
		be projected onto these segments during generation, but do not have to be.

	SaveLoopCounts: (Optional)
		This file defines the continous segements of boundary lines. The first
		column is a beginning index, and the second column is the number of
		points to read from SaveBoundaries to create a loop.

Output:
	After a successful run of MPI-SCVT, Three files are created.

	end_points.dat:
		This file contains the ending pointset from MPI-SCVT, given as x, y, z with one point per line

	triangles.dat:
		This file contains the final triangulation of the points in end_points.dat, given as one triangle
		per line, with each number being the index of a point in end_points.dat.

	bisected_points.dat:
		This file contains the bisected pointset from end_points.dat. For a full sphere grid, it should have
		4*n - 6 points, where n is the number of points in end_points.dat. 

[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/douglasjacobsen/mpi-scvt/trend.png)](https://bitdeli.com/free "Bitdeli Badge")
