Copyright (C) since 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington

RelaxMap is a parallel community detection algorithm to optimize flow-based information-theoretic objective function, called the map equation. 
RelaxMap is under GNU General Public License, detailed information is in LICENSE.txt.



1. How to compile

RelaxMap is implemented in C++ and uses OpenMP for shared-memory parallelism,
so your C++ compiler should contain OpenMP library.

You can compile by simple make command:

	> make


2. How to run RelaxMap

	[Usage] >./ompRelaxmap <seed> <network data> <# threads> <# attempts> <threshold> <vThresh> <maxIter> <outDir> <prior/normal>

The required arguments are following:

	1) seed: this is for random seed value for generating random sequential order of vertices during the iterative optimization process.

	2) network data: RelaxMap currently supports two different types of the input data format. 1) pajek format (.net) and 2) edge list format (.txt)
				You can find an example of each format in the data/ directory.
				You can also find more datasets from Stanford Network Analysis Project (SNAP) network data repository.
				You should note that RelaxMap assumes the input graph is directed graph.

	3) # thread: the number of threads 

	4) # attempts: the number of attempts (this is not applied yet, so it only return with 1 attempt.)

	5) threshold: the stop condition threshold  (recommended 1e-3 or 1e-4)

	6) vThresh: the threshold value for each vertex movement (recommended 0.0)

	7) maxIter: the number of maximum iteration for each super-step. (recommended 10)

	8) outDir: the directory where the output files will be located.

	9) prior/normal flag: apply the prioritized search for efficient runs (prior) or not (normal).  (recommended prior)


3. Reference

If you would like to add a reference for this application in documents, please put the following bibliography information:

	Seung-Hee Bae, Daniel Halperin, Jevin West, Martin Rosvall, and Bill Howe, 
	"Scalable Flow-Based Community Detection for Large-Scale Network Analysis,"
	In Proceedings of IEEE 13th International Conference on Data Mining Workshop (ICDMW), 2013


4. Bug Report

If you find a bug, please send a bug report to 
	Seung-Hee Bae: shbae@cs.washington.edu
	Bill Howe: billhowe@cs.washington.edu
