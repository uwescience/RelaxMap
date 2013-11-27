This version (v0.9.2) separate active Modules in two groups: small active modules and large-active-modules.
The two different groups are dealt in different way for better efficiency.


--------------------------------------

## Version 0.9 ##

Version 0.9 adds the two extensions in the paper, for improving the quality of the found communities, to the version (v0.7 - HogWild-Infomap).
	- Single-node movements
	- Sub-Module movements


## Version 0.7 ##

Version 0.7 implements infomap with OpenMP in HogWild fashion which is immediately update the decision,
so each decision will be made based on information at the point of access each module information.

The module information is possible to be stale information since module information updated in parallel.
However, if the graph (G) is sparse, it could happen rarely.  
Even if the information is stale, it could be similar value due to short time difference, it will be not the best decision but still could be valid decision,
which makes some improvement of output quality in high probability.
Possibly, some nodes' movements will make worse the quality, but in overall the progress will be made in each loop (a.k.a. iteration).

Other minor change is that this version make specifies output directory and the outputs will be located to the specified directory.
We can choose either fixed threshold value or adaptive threshold value in this version.


## Version 0.6 ##

Version 0.6 implements infomap with OpenMP in "Move-If-Improved" (MII) fashion.
In this version, we added adaptive threshold feature of the stop condition with respect to the number of supernodes at each super-step.


## Version 0.5 ##

Version 0.5 is actually optimized version of version 0.4.
It tries to optimize running performance of convertModuleToSPNode() method, which will be a huge bottle-neck with larger data (e.g. soc-LiveJournal1).



## Version 0.4 ##

We found that the v0.2 and v0.3 takes large amount of time for the converting modules to superNode procedure.
Thus, we optimized the convertModuleToSuperNode() function as following:
	- Instead of looking at outLinks and inLinks of each node, we first investigate outLinks of the member nodes in each module in parallel 
		to generate outLinks of the corresponding SuperNode.  
	- After we add all the outLinks of all superNodes in parallel, we just add inLinks of superNodes based on the corresponding outLinks of superNodes, in sequential fashion.
	- Since the number of nodes are reduced about an order of magnitude, it will reduce a lot works even run in sequential.
	- The false sharing might be the reason of the poor parallel performance of this, so it would be changed the MPI version with a large number of parallelism.
		- In other words, this feature might not work well in MPI version with 100's of processes or more.

	- So far, with test of 4 core or 8 core machine with OpenMP implementation, this added feature make some benefit of performance.


## Version 0.3 ##

v0.3 version tries to monotonic quality improvement by examining whether the suggested movement actually decreases MDL or not.
Then, it only executes the suggestion when MDL improved.  Otherwise, just ignores it.
Because we investigate the actual quality gain of each movement so that it prevent circular movement - If one move occurred the corresponding move will not be occurred.
Thus, we don't need to constrain the direction of new module, i.e. forward or backward in order to prevent circular movement of nodes. 
This results in faster convergence than the previous direction-constrained version.

In addition, we don't need to call updateCodeLength() function (AND updateMembersInModule() function) in this version, because during the examination we can directly 
use the necessary values which updated during the examination, such as exitPr, sumPr, and codeLength, etc.
The drawback of this version is that it requires recompute the flow information (outFlow, inFlow for oldMod and newMod),
which is also computationally expensive and sequentially implemented.  However, it will reduce the unnecessary movements happened in v0.2.
Therefore, finally it will make benefits in computationally and qualitatively.


## Version 0.2 ##

v0.2 moves all the decided movements of each node without regard to actual MDL improvement,
so it could make a movement which increases MDL. Also, it should call updateCodeLength() function which is O(E).

