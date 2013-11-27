#include <iostream>
#include <omp.h>
#include <cstdlib>
#if defined __GNUCC__ || defined __APPLE__
#include <ext/hash_map>
#else
#include <hash_map>
#endif
#include <ctime>
#include <sys/time.h>
#include "Module.h"
#include "Node.h"
#include "timing.h"


//typedef __gnu_cxx::hash_map<int, double> flowmap;
typedef map<int, double> flowmap;
typedef map<int, pair<double, double> > modInfo;	// <modID, <exitPr, sumPr> >


void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID);

//#define INIT_LENGTH 10000
const double InitLength = 10000.0;

using namespace std;

struct MoveSummary {
	double diffCodeLen;			// delta(L(M))
	int newModule;
	double sumPr1, sumPr2;		// updated sumPr1 and sumPr2.  1 --> oldM, 2 --> newM.
	double exitPr1, exitPr2;	// updated exitPr1 and exitPr2.
	//double sumDangling1;		// updated sumDangling1.
	//double sumDangling2;		// updated sumDangling2.
	double newSumExitPr;			// SUM(q_i) = this.sumAllExitPr + q_1_new + q_2_new - q_1_old - q_2_old.
	//double delta_sum_exit_log_exit;		// DELTA [sum {q_i log (q_i)}]
	//double delta_sum_stay_log_stay;		// DELTA [sum {(q_i+sumPr_i) log (q_i+sumPr_i)}]
};

// Constructors of Module class.
Module::Module()
:index(0),
 exitPr(0.0),
 stayPr(0.0),
 sumPr(0.0),
 sumTPWeight(0.0),
 sumDangling(0.0),
 numMembers(1)
{}


Module::Module(int idx, double exitPr, double sumPr)
:index(idx),
 exitPr(exitPr),
 stayPr(exitPr + sumPr),
 sumPr(sumPr),
 sumTPWeight(0.0),
 sumDangling(0.0),
 numMembers(1)
{}


//Module::Module(int idx, int nNode, Node nd)
Module::Module(int idx, Node * nd)
:index(idx),
 exitPr(0.0),
 stayPr(0.0),
 sumPr(nd->Size()),
 sumTPWeight(nd->TeleportWeight()),
 sumDangling(0.0),
 numMembers(1)
{
	// sum (p_a * w_ab)
	double sumExitFlow = 0.0;
	int nOutLinks = nd->outLinks.size();

//#pragma omp parallel for reduction (+:sumExitFlow) --> This Constructor will be called in parallel per node.
	for (int i = 0; i < nOutLinks; i++)
		sumExitFlow += nd->outLinks[i].second;	// after call initiate(), w_ab is updated as p_a * w_ab.
		
	// exitPr = tau * (1 - n_i/n)*p_a + (1-tau) * sumExitFlow.
	//exitPr = Network::alpha * sumPr * (nNode - 1) / nNode + Network::beta * sumExitFlow;
	
	//if (nd->IsSuper() || nd->IsModule()) {
	if (nd->IsSuper()) {	// If SuperNode, we already calculated exitPr in the corresponding module.
		exitPr = nd->ExitPr();
		sumDangling = nd->DanglingSize();
	}
	// exitPr = tau * (1 - sum(tau_a))*sum(p_a) + (1-tau) * sumExitFlow.
	else if (!nd->IsDangling()) {
		//exitPr = Network::alpha * (1.0 - nd.NodeWeight()) * sumPr + Network::beta * sumExitFlow;
		exitPr = Network::alpha * (1.0 - nd->TeleportWeight()) * sumPr + Network::beta * sumExitFlow;
		nd->setExitPr(exitPr);
	}
	else {
		//exitPr = (1.0 - nd.NodeWeight()) * sumPr;
		exitPr = (1.0 - nd->TeleportWeight()) * sumPr;
		nd->setExitPr(exitPr);
		sumDangling = nd->Size();
	}

	stayPr = exitPr + sumPr;	// exitPr + sum (P_a).
	members.push_back(nd);
}



SubModule::SubModule()
: modIdx(0), 
  numMembers(0), 
  sumPr(0.0), 
  sumTPWeight(0.0), 
  sumDangling(0.0) 
{}


/**
 *	Module mod : module of newNetwork (network based on a module in original network.)
 */
SubModule::SubModule(Module& mod, map<int, int>& origNodeID, int modIndex) {
	numMembers = mod.NumMembers();
	sumPr = mod.SumPr();
	sumTPWeight = mod.SumTPWeight();
	sumDangling = mod.SumDangling();
	modIdx = modIndex;

	//for (int i = 0; i < numMembers; i++) {
	//	members.push_back(origNodeID[mod.members->ID()]);
	for (vector<Node*>::iterator it = mod.members.begin(); it != mod.members.end(); it++) {
		members.push_back(origNodeID[(*it)->ID()]);
	}

	// numMembers should be equal to members.size() ...
	//	assert(numMembers == members.size());
	if (numMembers != members.size())
		cout << "SOMETHING WRONG!!! -- numMembers != members.size() in SubModule()..." << endl;
}

SubModule::SubModule(Module& mod) {
	modIdx = mod.Index();
	numMembers = 1;
	sumPr = mod.SumPr();
	sumTPWeight = mod.SumTPWeight();
	sumDangling = mod.SumDangling();

	members.push_back(mod.members[0]->ID());
}



// Constructors of Network class.
Network::Network()
:level(0),
 codeLength(InitLength),
 nNode(0),
 nEdge(0),
 nEmptyMod(0),
 nModule(0),
 totNodeWeights(0.0),
 nDanglings(0),
 //alpha(0.15),
 //beta(0.85),
 allNodes_log_allNodes(0.0),
 sumAllExitPr(0.0)
{}


Network::Network(int numNode, int numEdge, int level)
:level(level),
 codeLength(InitLength),
 nNode(numNode),
 nEdge(numEdge),
 nEmptyMod(0),
 nModule(numNode),
 totNodeWeights(0.0),
 nDanglings(0),
 //alpha(0.15),
 //beta(0.85),
 allNodes_log_allNodes(0.0),
 sumAllExitPr(0.0)
{
	emptyModules.reserve(numNode);
}



Network::Network(int numNode, int numEdge, int level, double codeLen)
:level(level),
 codeLength(codeLen),
 nNode(numNode),
 nEdge(numEdge),
 nEmptyMod(0),
 nModule(numNode),
 totNodeWeights(0.0),
 nDanglings(0),
 //alpha(0.15),
 //beta(0.85),
 allNodes_log_allNodes(0.0),
 sumAllExitPr(0.0)
{
	emptyModules.reserve(numNode);
}


/*
 * Other Member functions of Network class.
 */

void Network::findDanglingNodes() {

	nDanglings = 0;		// reset the number of dangling nodes.

	for (int i = 0; i < nNode; i++) {
		//if (nodes[i].outLinks.empty() && (nodes[i].selfLink <= 0.0)) {
		if (nodes[i].outLinks.empty()) {
			danglings.push_back(i);
			nDanglings++;
			nodes[i].setIsDangling(true);
		}
	}
	
}

//void Network::initiate() {
void Network::initiate(int numTh) {

	//nDanglings = 0;		// reset the number of dangling nodes.
	int nDangNodes = 0;

	struct timeval startT, endT;

	omp_set_num_threads(numTh);

	gettimeofday(&startT, NULL);

	//#pragma omp parallel for reduction (+:nDangNodes)
	#pragma omp parallel reduction (+:nDangNodes) 
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		//for (int i = 0; i < nNode; i++) {
		for (int i = start; i < end; i++) {
			//if (nodes[i].outLinks.empty() && (nodes[i].selfLink <= 0.0)) {
			if (nodes[i].outLinks.empty()) {
				#pragma omp critical (updateDang)
				{
					danglings.push_back(i);
				}
				//nDanglings++;
				nDangNodes++;
				nodes[i].setIsDangling(true);
			}
			else {	// normal nodes. --> Normalize edge weights.
				int nOutLinks = nodes[i].outLinks.size();
				//double sum = nodes[i].selfLink;	// don't support selfLink yet.
				double sum = 0.0;
				for (int j = 0; j < nOutLinks; j++)
					sum += nodes[i].outLinks[j].second;

				//nodes[i].selfLink /= sum;
				for (int j = 0; j < nOutLinks; j++)
					nodes[i].outLinks[j].second /= sum;
			}
		}
	}

	gettimeofday(&endT, NULL);

	nDanglings = nDangNodes;

	cout << "Level " << level << ": the number of dangling nodes = " << nDanglings << endl;
	cout << "Time for finding dangling nodes : " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	gettimeofday(&startT, NULL);
	//calculateSteadyState();		//eigenvector();
	calculateSteadyState(numTh);		//eigenvector();
	gettimeofday(&endT, NULL);
	cout << "Time for calculating steady state of nodes (eigenvector): " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;


	gettimeofday(&startT, NULL);

	// Update edges to represent flow based on calculated steady state (aka size).
	//#pragma omp parallel for
	#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		//for (int i = 0; i < nNode; i++) {
		for (int i = start; i < end; i++) {
			//nodes[i].selfLink = beta * nodes[i].Size() * nodes[i].selfLink;


			if (!nodes[i].IsDangling()) {
				int nOutLinks = nodes[i].outLinks.size();
				for (int j = 0; j < nOutLinks; j++) {
					//nodes[i].outLinks[j].second = beta * nodes[i].Size() * nodes[i].outLinks[j].second;
					nodes[i].outLinks[j].second = nodes[i].Size() * nodes[i].outLinks[j].second;

					// lazy update of inLinks after we find out final flow of this networks.
					//#pragma omp critical
					//nodes[nodes[i].outLinks[j].first].inLinks.push_back(make_pair(i,nodes[i].outLinks[j].second));
				}
			}
			else {
				nodes[i].setDanglingSize(nodes[i].Size());
			}
		}
	}

	gettimeofday(&endT, NULL);
	cout << "Time for updating edge weights based on flow information: " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;
	
	// Calculate SUM of p_log_p over all nodes.
	//allNodes_log_allNodes = 0.0;
	double allNds_log_allNds = 0.0;

	//#pragma omp parallel for reduction (+:allNds_log_allNds)
	#pragma omp parallel reduction (+:allNds_log_allNds)
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);
	
		//for (int i = 0; i < nNode; i++)
		for (int i = start; i < end; i++)
			//allNodes_log_allNodes += pLogP(nodes[i].Size());
			allNds_log_allNds += pLogP(nodes[i].Size());
	}

	allNodes_log_allNodes = allNds_log_allNds;

	// DEBUG PURPOSE: Print nodeSize (P_a) for each node.
	//cout << "Node Size of each node...\n";
	//for (int i = 0; i < nNode; i++)
	//	cout << "Node[" << i << "]:\t" << nodes[i].Size() << endl;

	/////////////////////////////
	// Make modules from nodes //
	/////////////////////////////

	// The exit flow from each node at initiation.
	// each node == each module at initial time.
	//for (int i = 0; i < nNode; i++)
	//	nodes[i].exit = nodes[i].Size() - (alpha * nodes[i].Size() + beta * nodes[i].danglingSize) * nodes[i].TeleportWeight() - nodes[i].selfLink;

	//#pragma omp parallel for
	#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		//for (int i = 0; i < nNode; i++) {
		for (int i = start; i < end; i++) {
			modules[i] = Module(i, &nodes[i]);	// Assign each Node to a corresponding Module object.
												// Initially, each node is in its own module.
			nodes[i].setModIdx(i);
		}
	}
	nModule = nNode;

	gettimeofday(&startT, NULL);
	// calibrate();
	calibrate(numTh);
	gettimeofday(&endT, NULL);
	cout << "Time for calculating of initial code length: " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;

}


// calculating steady state of nodes by Power Iteration Method.
// same as Greedy::eigenvector() in Infomap implementation.
// Modify for the design of this implementation.

//void Network::calculateSteadyState() {
void Network::calculateSteadyState(int numTh) {
	// initial probability distribution = 1 / N.
	vector<double> size_tmp = vector<double>(nNode, 1.0/nNode);
	//double nodeSizes[nNode]; // temporary vector for nodeSize to use locks.

	//omp_lock_t* nodeSizes_locks = new omp_lock_t[nNode];

	int iter = 0;
	double danglingSize = 0.0;
	double sqdiff = 1.0;
	//double sqdiff_old = 1.0;
	double sum = 0.0;

	// initialize lock for each node.
	//#pragma omp parallel for 
	//for (int i = 0; i < nNode; i++)
	//	omp_init_lock(&nodeSizes_locks[i]);

	// Generate addedSize array per each thread, so that we don't need to use lock for each nodeSize addition.
	double** addedSize = new double*[numTh];
	for (int i = 0; i < numTh; i++)
		addedSize[i] = new double[nNode];
			
	do {
		// calculate sum of the size of dangling nodes.
		danglingSize = 0.0;
	
		//#pragma omp parallel for reduction (+:danglingSize)
		#pragma omp parallel reduction (+:danglingSize)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nDanglings, nTh, myID);
			
			//for (int i = 0; i < nDanglings; i++)
			for (int i = start; i < end; i++)
				danglingSize += size_tmp[danglings[i]];
		}

		// flow via teleportation.
		//#pragma omp parallel for
		#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			//for (int i = 0; i < nNode; i++) 
			for (int i = start; i < end; i++) 
				// alpha * nodes[i].TeleporWeight() + beta*danglingSize*nodes[i].TeleportWeight();
				nodes[i].setSize( (alpha + beta*danglingSize) * nodes[i].TeleportWeight());
				//nodeSizes[i] = (alpha + beta*danglingSize) * nodes[i].TeleportWeight();
		}


		int realNumTh = 0;

		// flow from network steps via following edges.
		//#pragma omp parallel for
		#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			#pragma omp master
			{
				realNumTh = nTh;
			}

			double *myAddedSize = addedSize[myID];
			for (int i = 0; i < nNode; i++)
				myAddedSize[i] = 0.0;	// initialize for the temporary addedSize array.

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			//for (int i = 0; i < nNode; i++) {
			for (int i = start; i < end; i++) {
				//nodes[i].addSize(beta * nodes[i].selfLink * size_tmp[i]);
				int nOutLinks = nodes[i].outLinks.size();
				for (int j = 0; j < nOutLinks; j++) {
					//nodes[nodes[i].outLinks[j].first].addSize( beta * nodes[i].outLinks[j].second * size_tmp[i]);
					myAddedSize[nodes[i].outLinks[j].first] +=  beta * nodes[i].outLinks[j].second * size_tmp[i];

					//int target = nodes[i].outLinks[j].first;
					//omp_set_lock(&nodeSizes_locks[target]);
					//nodeSizes[target] += beta * nodes[i].outLinks[j].second * size_tmp[i];
					//omp_unset_lock(&nodeSizes_locks[target]);
				}
			}
		}
		//set nodes[i].size by added addedSize[] values.
		#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			//for (int i = 0; i < nNode; i++) {
			for (int i = start; i < end; i++) {
				for (int j = 0; j < realNumTh; j++)
					nodes[i].addSize(addedSize[j][i]);
			}
		}

		// set nodes[i].size value from nodeSizes[i]
		//#pragma omp parallel for
		//for (int i = 0; i < nNode; i++)
		//	nodes[i].setSize(nodeSizes[i]);
			


		// Normalize of node size.
		sum = 0.0;
		//#pragma omp parallel for reduction (+:sum)
		#pragma omp parallel reduction (+:sum)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			//for (int i = 0; i < nNode; i++) 
			for (int i = start; i < end; i++) 
				sum += nodes[i].Size();
		}
		//sqdiff_old = sqdiff;
		sqdiff = 0.0;

		//#pragma omp parallel for reduction (+:sqdiff)
		#pragma omp parallel reduction (+:sqdiff)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			//for (int i = 0; i < nNode; i++) {
			for (int i = start; i < end; i++) {
				nodes[i].setSize(nodes[i].Size()/sum);
				sqdiff += fabs(nodes[i].Size() - size_tmp[i]);
				size_tmp[i] = nodes[i].Size();
			}
		}

		iter++;

		// I AM NOT SURE BELOW PART SHOULD BE NECESSARY.
		/*if (sqdiff == sqdiff_old) {
			alpha += 1.0e-10;
			beta = 1.0 - alpha;
		}
		*/

	} while ((iter < 200) && (sqdiff > 1.0e-15 || iter < 50));

	// deallocate 2D array.
	for (int i = 0; i < numTh; i++)
		delete [] addedSize[i];
	delete [] addedSize;

	cout << "Calculating flow done in " << iter << " iterations!" << endl;

	// destroy locks.
	//#pragma omp parallel for
	//for (int i = 0; i < nNode; i++)
	//	omp_destroy_lock(&nodeSizes_locks[i]);

	//delete nodeSizes_locks;
}


// This function calculate current codeLength.
// This implementation is modified version of infomap implementation.
void Network::calibrate(int numTh) {
	//This is the calculation of Equation (4) in the paper.
	double sum_exit_log_exit = 0.0;
	double sum_stay_log_stay = 0.0;
	double sumExit = 0.0;

	omp_set_num_threads(numTh);

	//#pragma omp parallel for reduction (+:sum_exit_log_exit, sum_stay_log_stay, sumExit)
	#pragma omp parallel reduction (+:sum_exit_log_exit, sum_stay_log_stay, sumExit)
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nModule, nTh, myID);
			
		//for (unsigned int i = 0; i < nModule; i++) {
		for (unsigned int i = start; i < end; i++) {
			sum_exit_log_exit += pLogP(modules[i].ExitPr());
			sum_stay_log_stay += pLogP(modules[i].StayPr());
			sumExit += modules[i].ExitPr();
		}
	}

	sumAllExitPr = sumExit;
	double sumExit_log_sumExit = pLogP(sumExit);

	//cout << "sumExit_log_sumExit = " << sumExit_log_sumExit << endl;
	//cout << "sum_size_log_size = " << allNodes_log_allNodes << endl;
	//cout << "sum_exit_log_exit = " << sum_exit_log_exit << endl;
	//cout << "sum_stay_log_stay = " << sum_stay_log_stay << endl;

	codeLength = sumExit_log_sumExit - 2.0 * sum_exit_log_exit + sum_stay_log_stay - allNodes_log_allNodes;
}

/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in SeqInfomap.cpp.
 *
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 */
bool Network::move() {

	// Generate random sequential order of nodes.
	vector<int> randomOrder(nNode);
	for (int i = 0; i < nNode; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nNode; i++) {
		int target = R->randInt(nNode - 1);
		//int target = R.randInt(nNode - 1);
		//int target = rand() % nNode;

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	bool movedAny = false;

	// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nNode; i++) {
		//int nd = randomOrder[i];		// look at i_th Node of the random sequential order.
		Node& nd = nodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		//double outFlowNotToOldMod = 0.0;	// Sum_b (P_nd * w_nd_b) where b not in oldMod.
		//double inFlowFromOldMod = 0.0;		// SUM_a (P_a * w_a_nd)  where a in oldMod.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		// During this counting, aggregate outFlowNotToOldMod and inFlowFromOldMod values.
		//int nOutEdges = nd.outLinks.size();
		//for (int j = 0; j < nOutEdges; j++) {
		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
			//int newMod = nodes[nd.outLinks[j].first].ModIdx();
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				//outFlowToMod[newMod] += beta * nd.outLinks[j].second;
				outFlowToMod[newMod] += beta * linkIt->second;
			}
			else {
				//outFlowToMod[newMod] = beta * nd.outLinks[j].second;	// initialization of the outFlow of the current newMod.
				outFlowToMod[newMod] = beta * linkIt->second;	// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		//int nInEdges = nd.inLinks.size();
		//for (int j = 0; j < nInEdges; j++) {
		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
			//int newMod = nodes[nd.inLinks[j].first].ModIdx();
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				//inFlowFromMod[newMod] += beta * nd.inLinks[j].second;
				inFlowFromMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = 0.0;
				//inFlowFromMod[newMod] = beta * nd.inLinks[j].second;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size() in Network::move()." << endl;

		// Debug purpose...
		//cout << "Node[" << nd << "]: oldMod = " << oldMod << ", nOutEdges = " << nOutEdges << ", nInEdges = " << nInEdges << ", nModLinks = " << nModLinks << endl;



		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();					// p_nd.
		double ndTPWeight = nd.TeleportWeight();		// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].ExitPr();
		double oldSumPr1 = modules[oldMod].SumPr();
		double oldSumDangling1 = modules[oldMod].SumDangling();
		double oldModTPWeight = modules[oldMod].SumTPWeight();
		
		double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			}
			else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
				inFlowFromMod[newMod] += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}


		//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			int emptyTarget = emptyModules.back();
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;



		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;	// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		// We don't need to check newMod != oldMod, since all newMod values are not equal to oldMod.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].ExitPr();
				double oldSumPr2 = modules[newMod].SumPr();
				//double oldSumDangling2 = modules[newMod].SumDangling();
				
				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1 + currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1) + pLogP(currentResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(currentResult.exitPr1 + currentResult.sumPr1) + pLogP(currentResult.exitPr2 + currentResult.sumPr2) \
											- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
					bestResult.sumPr1 = currentResult.sumPr1;
					bestResult.sumPr2 = currentResult.sumPr2;
					bestResult.exitPr1 = currentResult.exitPr1;
					bestResult.exitPr2 = currentResult.exitPr2;
					bestResult.newSumExitPr = currentResult.newSumExitPr;
				}
			}
		}



		// Make best possible move for the current node nd.
		if (bestResult.diffCodeLen < 0.0) {
		//if (bestResult.diffCodeLen < -1.0e-10) {
			// update related to newMod...
			int newMod = bestResult.newModule;

			//cout << "NewModule: " << newMod << endl;
			//if (modules[newMod].members.size() == 0) {
			if (modules[newMod].NumMembers() == 0) {
				newMod = emptyModules.back();
				emptyModules.pop_back();
				nEmptyMod--;
				nModule++;
			}
			nd.setModIdx(newMod);

			//modules[newMod].members.push_back(&nd);
			modules[newMod].increaseNumMembers();
			modules[newMod].setExitPr(bestResult.exitPr2);
			modules[newMod].setSumPr(bestResult.sumPr2);
			modules[newMod].setStayPr(bestResult.exitPr2 + bestResult.sumPr2);
			modules[newMod].addSumTPWeight(ndTPWeight);

			if (nd.IsDangling()) {
				modules[newMod].addSumDangling(ndSize);
				modules[oldMod].minusSumDangling(ndSize);
			}

			// update related to the oldMod...
			modules[oldMod].decreaseNumMembers();
			modules[oldMod].setExitPr(bestResult.exitPr1);
			modules[oldMod].setSumPr(bestResult.sumPr1);
			modules[oldMod].setStayPr(bestResult.exitPr1 + bestResult.sumPr1);
			modules[oldMod].minusSumTPWeight(ndTPWeight);

			if (modules[oldMod].NumMembers() == 0) {
				nEmptyMod++;
				nModule--;
				emptyModules.push_back(oldMod);
				//cout << "Number of Empty Module: " << nEmptyMod << endl;
			}

			// remove the current node from oldMod
/*			int nodeID = nd.ID();
			int nMembers = modules[oldMod].members.size();

			if (nMembers == 1) {
				modules[oldMod].members.pop_back();
			}
			else {
				for (int n = 0; n < nMembers; n++) {
					if (modules[oldMod].members[n]->ID() == nodeID) {
						// swap the found node element and the back() element,
						// so that we can reduce the time complexity for erasing the node to the constant.
						// compare vector::erase() to vector::pop_back().
						Node * tempND = modules[oldMod].members[n];
						modules[oldMod].members[n] = modules[oldMod].members[nMembers - 1];
						modules[oldMod].members[nMembers-1] = tempND;

						modules[oldMod].members.pop_back();		// remove the current node which swapped to the end of the vector.

						break;	// We don't need to search any more.
					}
				}
			}
*/

			sumAllExitPr = bestResult.newSumExitPr;

			codeLength += bestResult.diffCodeLen;
			movedAny = true;
		}
	}

	/////////// NEED TO CONCATENATE MODULE VECTOR ///////////////
/*	int idx = 0;
	int backIdx = modules.size() - 1;

	while (movedAny && (idx <= backIdx)) {
		if (modules[idx].NumMembers() == 0) {
			// find non-empty modules from the end and swap it with this.
			while (modules[backIdx].NumMembers() == 0) {
				backIdx--;
				modules.pop_back();
			}
			
			if ( idx < backIdx) {
				//modules[idx] = modules[backIdx];	// copy a non-empty module at 'backIdx' to 'idx'
				//vector<Node *>().swap(modules[idx].members);
				//for (int k = 0; k < modules[backIdx].members.size(); k++)
				//	modules[idx].members.push_back(modules[backIdx].members[k]);

				//copyModule(*(modules[idx]), *(modules[backIdx]));	// copy a non-empty module from 'backIdx' to 'idx'.
				copyModule(idx, backIdx);	// copy a non-empty module from 'backIdx' to 'idx'.
				backIdx--;
				modules.pop_back();	// Didn't swap it, because the empty module should be removed.

				modules[idx].setIndex(idx);

				// update modIdx of each node in modules[idx] which was at modules[backIdx].
				int memSize = modules[idx].members.size();
				for (int i = 0; i < memSize; i++)
					modules[idx].members[i]->setModIdx(idx);
			}
		}
		idx++;
	}
*/

	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod." << endl;
	}

	return movedAny;
}



/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in OmpInfomap.cpp.
 *
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 *
 *	+ This is a parallel implementation of Network::move() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */
bool Network::parallelMove(int numTh, double& tSequential) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	// Generate random sequential order of nodes.
	vector<int> randomOrder(nNode);
	for (int i = 0; i < nNode; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nNode; i++) {
		int target = R->randInt(nNode - 1);
		//int target = R.randInt(nNode - 1);
		//int target = rand() % nNode;

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	bool movedAny = false;

	omp_set_num_threads(numTh);

	// generate vectors which contains node-movements results from each thread.
	typedef pair<int, int> moveinfo;	// first = nodeIndex, second = new module index.
	vector<vector<moveinfo> > movements(numTh);

	const int emptyTarget = nNode + 1;	// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);
	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
	// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nNode; i++) {
		//int nd = randomOrder[i];		// look at i_th Node of the random sequential order.
		Node& nd = nodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;	// The number of links to/from between the current node and other modules.


		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		//modInfo modSnapShot;		// <modID, <exitPr, sumPr> > for neighbor modules

		// count other modules that are connected to the current node.
		// During this counting, aggregate outFlowNotToOldMod and inFlowFromOldMod values.
		//int nOutEdges = nd.outLinks.size();
		//for (int j = 0; j < nOutEdges; j++) {
		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
			//int newMod = nodes[nd.outLinks[j].first].ModIdx();
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = beta * linkIt->second;	// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		//int nInEdges = nd.inLinks.size();
		//for (int j = 0; j < nInEdges; j++) {
		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
			//int newMod = nodes[nd.inLinks[j].first].ModIdx();
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size() in Network::move()." << endl;


		/***********************************************************
		 * Here we store snapshot of related modules' information.
		 ***********************************************************/
/*		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int nMod = it->first;
			// This is a kind of snapshot of module information at the time of access each module. NOW...
			// Since the module information can be changed any time (but rarely happen for each module due to sparsity of the Graph),
			// instead of using lock for correct information, we decide to use a snapshot information which could be stale info.
			modSnapShot[nMod] = make_pair(modules[nMod].ExitPr(), modules[nMod].SumPr());
		}
*/

		// Debug purpose...
		//cout << "Node[" << nd << "]: oldMod = " << oldMod << ", nOutEdges = " << nOutEdges << ", nInEdges = " << nInEdges << ", nModLinks = " << nModLinks << endl;



		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();					// p_nd.
		double ndTPWeight = nd.TeleportWeight();		// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		// Below values are related current module information, but it could be changed in the middle of this decision process.
		// However, we used this snapshot for finding next module of current node.
		// These will be a guideline of decision process and the correct values will be calculated at the end of this iteration.
		double oldExitPr1 = modules[oldMod].ExitPr();
		double oldSumPr1 = modules[oldMod].SumPr();
		double oldSumDangling1 = modules[oldMod].SumDangling();
		double oldModTPWeight = modules[oldMod].SumTPWeight();
		
		double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			}
			else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
				inFlowFromMod[newMod] += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}


		//////////////////// THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			//int emptyTarget = nNode + 1;
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;



		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;	// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		// We don't need to check newMod != oldMod, since all newMod values are not equal to oldMod.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].ExitPr();
				double oldSumPr2 = modules[newMod].SumPr();
				//double oldSumDangling2 = modules[newMod].SumDangling();
				
				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1 + currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1) + pLogP(currentResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(currentResult.exitPr1 + currentResult.sumPr1) + pLogP(currentResult.exitPr2 + currentResult.sumPr2) \
											- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
					//bestResult.sumPr1 = currentResult.sumPr1;
					//bestResult.sumPr2 = currentResult.sumPr2;
					//bestResult.exitPr1 = currentResult.exitPr1;
					//bestResult.exitPr2 = currentResult.exitPr2;
					//bestResult.newSumExitPr = currentResult.newSumExitPr;
				}
			}
		}

		// store the best possilbe movement information if necessary.
		// In this version, we are trying to move as soon as possible, i.e. right after decision making...
		if (bestResult.diffCodeLen < 0.0) {
			//int thID = omp_get_thread_num();
			//movements[thID].push_back(make_pair(randomOrder[i], bestResult.newModule));
			
			bool isEmptyTarget = false;
			bool validMove = true;		// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;

			#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);
			
				// if newMod == emptyTarget, it indicates moves to empty module.
				if ( (nEmptyMod > 0) && (newMod == emptyTarget) && (modules[oldMod].NumMembers() > 1) ) {
					newMod = emptyModules.back();
					//emptyModules.pop_back();
				
					//nEmptyMod--;
					//nModule++;
					isEmptyTarget = true;
				}
				else if (newMod == emptyTarget) {
					//continue;
					validMove = false;
				}
				//else if (modules[oldMod].NumMembers() == 1 && modules[newMod].NumMembers() == 0) {
				//	continue;	// This case don't need to move at all.
				//}
				else if (modules[newMod].NumMembers() == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there. 
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}


				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					//int nOutEdges = nd.outLinks.size();
					//int nInEdges = nd.inLinks.size();
	
					if (!isEmptyTarget) {
						for (link_iterator linkIt =nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();
	
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
							else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt =nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();
	
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
							else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					}
					else {
						for (link_iterator linkIt =nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt =nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}


					oldExitPr1 = modules[oldMod].ExitPr();
					oldSumPr1 = modules[oldMod].SumPr();
					oldSumDangling1 = modules[oldMod].SumDangling();
					oldModTPWeight = modules[oldMod].SumTPWeight();
		
					//double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
					//double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
					inFlowFromNewMod += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;


					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					//MoveSummary moveResult;

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;


					// copy module specific values...
					double oldExitPr2 = modules[newMod].ExitPr();
					double oldSumPr2 = modules[newMod].SumPr();
					//double oldSumDangling2 = modules[newMod].SumDangling();
				
					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1 + moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1) + pLogP(moveResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(moveResult.exitPr1 + moveResult.sumPr1) + pLogP(moveResult.exitPr2 + moveResult.sumPr2) \
												- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

//					if (moveResult.diffCodeLen >= 0.0)
//						//continue;	// ignore this movement option.
//						validMove = false;
//				}

//				if (validMove) {

					if (isEmptyTarget) {
						//newMod = emptyModules.back();		// Done before...
						emptyModules.pop_back();
				
						nEmptyMod--;
						nModule++;
					}


					nd.setModIdx(newMod);

					//modules[newMod].members.push_back(&nodes[ndIdx]);
					modules[newMod].increaseNumMembers();
					modules[newMod].setExitPr(moveResult.exitPr2);
					modules[newMod].setSumPr(moveResult.sumPr2);
					modules[newMod].setStayPr(moveResult.exitPr2 + moveResult.sumPr2);
					modules[newMod].addSumTPWeight(ndTPWeight);

					if (nd.IsDangling()) {
						modules[newMod].addSumDangling(ndDanglingSize);
						modules[oldMod].minusSumDangling(ndDanglingSize);
					}

					// update related to the oldMod...
					modules[oldMod].decreaseNumMembers();
					modules[oldMod].setExitPr(moveResult.exitPr1);
					modules[oldMod].setSumPr(moveResult.sumPr1);
					modules[oldMod].setStayPr(moveResult.exitPr1 + moveResult.sumPr1);
					modules[oldMod].minusSumTPWeight(ndTPWeight);

					if (modules[oldMod].NumMembers() == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;
			
					movedAny = true;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

		}

	}	// END parallel for


	/////// SEQUENTIAL PART STARTING HERE!!! ///////////
//	gettimeofday(&tStart, NULL);


	// Check if there any node is moved to a new module...
//	for (int i = 0; i < numTh; i++) {
//		if (movements[i].size() > 0) {
//			movedAny = true;
//			break;
//		}
//	}


	/*
	 *	This part is the actual movement of the node-movement decision from each threads.
	 *	We will actually move nodes if the decision actually improves the MDL, and ignore, otherwise.
	 *	In this way, we may recompute the exitPr of the movement, 
	 *	but we don't need to call updateCodeLength, since we can update codeLength as same as in sequential mode.
	 */
/*	for (int i = 0; i < numTh; i++) {
		int numMove = movements[i].size();
		for (int j = 0; j < numMove; j++) {
			// update related to newMod...
			int ndIdx = movements[i][j].first;
			int oldMod = nodes[ndIdx].ModIdx();
			int newMod = movements[i][j].second;

			double ndSize = nodes[ndIdx].Size();
			double ndTPWeight = nodes[ndIdx].TeleportWeight();
			double ndDanglingSize = nodes[ndIdx].DanglingSize();

			bool isEmptyTarget = false;

			// if newMod == emptyTarget, it indicates moves to empty module.
			if ( (nEmptyMod > 0) && (newMod == emptyTarget) && (modules[oldMod].NumMembers() > 1) ) {
				newMod = emptyModules.back();
				//emptyModules.pop_back();
				
				//nEmptyMod--;
				//nModule++;
				isEmptyTarget = true;
			}
			else if (newMod == emptyTarget) {
				continue;
			}
			//else if (modules[oldMod].NumMembers() == 1 && modules[newMod].NumMembers() == 0) {
			//	continue;	// This case don't need to move at all.
			//}
			else if (modules[newMod].NumMembers() == 0) {
				// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
				// but the nodes are all moved to other modules so there are no nodes in there. 
				// Thus, we don't know whether generating a new module will be better or not.
				// Therefore, don't move this option.
				continue;
			}

			////////////////////////////////////////////////////////////////////
			/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
			////////////////////////////////////////////////////////////////////

			double outFlowToOldMod = 0.0;
			double outFlowToNewMod = 0.0;
			double inFlowFromOldMod = 0.0;
			double inFlowFromNewMod = 0.0;

			int nOutEdges = nodes[ndIdx].outLinks.size();
			int nInEdges = nodes[ndIdx].inLinks.size();

			if (!isEmptyTarget) {
				for (int j = 0; j < nOutEdges; j++) {
					int toMod = nodes[nodes[ndIdx].outLinks[j].first].ModIdx();

					if (toMod == oldMod) {
						outFlowToOldMod += beta * nodes[ndIdx].outLinks[j].second;
					}
					else if (toMod == newMod) {
						outFlowToNewMod += beta * nodes[ndIdx].outLinks[j].second;
					}
				}

				for (int j = 0; j < nInEdges; j++) {
					int fromMod = nodes[nodes[ndIdx].inLinks[j].first].ModIdx();

					if (fromMod == oldMod) {
						inFlowFromOldMod += beta * nodes[ndIdx].inLinks[j].second;
					}
					else if (fromMod == newMod) {
						inFlowFromNewMod += beta * nodes[ndIdx].inLinks[j].second;
					}
				}
			}
			else {
				for (int j = 0; j < nOutEdges; j++) {
					int toMod = nodes[nodes[ndIdx].outLinks[j].first].ModIdx();
					if (toMod == oldMod) 
						outFlowToOldMod += beta * nodes[ndIdx].outLinks[j].second;
				}

				for (int j = 0; j < nInEdges; j++) {
					int fromMod = nodes[nodes[ndIdx].inLinks[j].first].ModIdx();
					if (fromMod == oldMod) 
						inFlowFromOldMod += beta * nodes[ndIdx].inLinks[j].second;
				}
			}



			double oldExitPr1 = modules[oldMod].ExitPr();
			double oldSumPr1 = modules[oldMod].SumPr();
			double oldSumDangling1 = modules[oldMod].SumDangling();
			double oldModTPWeight = modules[oldMod].SumTPWeight();
		
			//double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
			//double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

			// For teleportation and danling nodes.
			outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
			inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;
			outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
			inFlowFromNewMod += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;


			if (isEmptyTarget) {
				outFlowToNewMod = 0.0;
				inFlowFromNewMod = 0.0;
			}

			MoveSummary moveResult;

			moveResult.exitPr1 = oldExitPr1 - nodes[ndIdx].ExitPr() + outFlowToOldMod + inFlowFromOldMod;


			// copy module specific values...
			double oldExitPr2 = modules[newMod].ExitPr();
			double oldSumPr2 = modules[newMod].SumPr();
			//double oldSumDangling2 = modules[newMod].SumDangling();
				
			// Calculate status of current investigated movement of the node nd.
			moveResult.newModule = newMod;
			moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
			moveResult.sumPr2 = oldSumPr2 + ndSize;
			moveResult.exitPr2 = oldExitPr2 + nodes[ndIdx].ExitPr() - outFlowToNewMod - inFlowFromNewMod;

			moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1 + moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

			// Calculate delta_L(M) = L(M)_new - L(M)_old
			double delta_allExit_log_allExit = pLogP(moveResult.newSumExitPr) - pLogP(sumAllExitPr);
			double delta_exit_log_exit = pLogP(moveResult.exitPr1) + pLogP(moveResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
			double delta_stay_log_stay = pLogP(moveResult.exitPr1 + moveResult.sumPr1) + pLogP(moveResult.exitPr2 + moveResult.sumPr2) \
										- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

			// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
			moveResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

			////////////////////////////////////
			////// THE END OF EXAMINATION //////
			////////////////////////////////////

			if (moveResult.diffCodeLen >= 0.0)
				continue;	// ignore this movement option.


			if (isEmptyTarget) {
				//newMod = emptyModules.back();		// Done before...
				emptyModules.pop_back();
				
				nEmptyMod--;
				nModule++;
			}


			nodes[ndIdx].setModIdx(newMod);

			//modules[newMod].members.push_back(&nodes[ndIdx]);
			modules[newMod].increaseNumMembers();
			modules[newMod].setExitPr(moveResult.exitPr2);
			modules[newMod].setSumPr(moveResult.sumPr2);
			modules[newMod].setStayPr(moveResult.exitPr2 + moveResult.sumPr2);
			modules[newMod].addSumTPWeight(ndTPWeight);

			if (nodes[ndIdx].IsDangling()) {
				modules[newMod].addSumDangling(ndDanglingSize);
				modules[oldMod].minusSumDangling(ndDanglingSize);
			}

			// update related to the oldMod...
			modules[oldMod].decreaseNumMembers();
			modules[oldMod].setExitPr(moveResult.exitPr1);
			modules[oldMod].setSumPr(moveResult.sumPr1);
			modules[oldMod].setStayPr(moveResult.exitPr1 + moveResult.sumPr1);
			modules[oldMod].minusSumTPWeight(ndTPWeight);

			if (modules[oldMod].NumMembers() == 0) {
				nEmptyMod++;
				nModule--;
				emptyModules.push_back(oldMod);
			}

			sumAllExitPr = moveResult.newSumExitPr;
			codeLength += moveResult.diffCodeLen;
		}
	}
*/


	/////////// NEED TO CONCATENATE MODULE VECTOR ///////////////
/*	int idx = 0;
	int backIdx = modules.size() - 1;

	while (movedAny && (idx <= backIdx)) {
		if (modules[idx].NumMembers() == 0) {
			// find non-empty modules from the end and swap it with this.
			while (modules[backIdx].NumMembers() == 0) {
				backIdx--;
				modules.pop_back();
			}
			
			if ( idx < backIdx) {
				//modules[idx] = modules[backIdx];	// copy a non-empty module at 'backIdx' to 'idx'
				//vector<Node *>().swap(modules[idx].members);
				//for (int k = 0; k < modules[backIdx].members.size(); k++)
				//	modules[idx].members.push_back(modules[backIdx].members[k]);

				//copyModule(*(modules[idx]), *(modules[backIdx]));	// copy a non-empty module from 'backIdx' to 'idx'.
				copyModule(idx, backIdx);	// copy a non-empty module from 'backIdx' to 'idx'.
				backIdx--;
				modules.pop_back();	// Didn't swap it, because the empty module should be removed.

				modules[idx].setIndex(idx);

				// update modIdx of each node in modules[idx] which was at modules[backIdx].
				int memSize = modules[idx].members.size();
				for (int i = 0; i < memSize; i++)
					modules[idx].members[i]->setModIdx(idx);
			}
		}
		idx++;
	}
*/

	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod." << endl;
	}

	//gettimeofday(&tEnd, NULL);

	//tSequential += elapsedTimeInSec(tStart, tEnd);
	
	//gettimeofday(&tStart, NULL);
	//updateMembersInModule();
	//updateCodeLength(numTh, false);
	//gettimeofday(&tEnd, NULL);

	//tSequential += elapsedTimeInSec(tStart, tEnd);

	return movedAny;
}










void Network::removeEmptyModules() {

	/////////// NEED TO CONCATENATE MODULE VECTOR ///////////////
	int idx = 0;
	int backIdx = modules.size() - 1;

	while (idx <= backIdx) {
		if (modules[idx].NumMembers() == 0) {
			// find non-empty modules from the end and swap it with this.
			while (modules[backIdx].NumMembers() == 0) {
				backIdx--;
				modules.pop_back();
			}
			
			if ( idx < backIdx) {
				//modules[idx] = modules[backIdx];	// copy a non-empty module at 'backIdx' to 'idx'
				//vector<Node *>().swap(modules[idx].members);
				//for (int k = 0; k < modules[backIdx].members.size(); k++)
				//	modules[idx].members.push_back(modules[backIdx].members[k]);

				//copyModule(*(modules[idx]), *(modules[backIdx]));	// copy a non-empty module from 'backIdx' to 'idx'.
				copyModule(idx, backIdx);	// copy a non-empty module from 'backIdx' to 'idx'.
				backIdx--;
				modules.pop_back();	// Didn't swap it, because the empty module should be removed.

				modules[idx].setIndex(idx);

				// update modIdx of each node in modules[idx] which was at modules[backIdx].
				int memSize = modules[idx].members.size();
				for (int i = 0; i < memSize; i++)
					modules[idx].members[i]->setModIdx(idx);
			}
		}
		idx++;
	}

	// the following should be true: modules.size() == nModule.
	//cout << "modules.size() = " << modules.size() << ", nModule = " << nModule << endl;

	if (modules.size() != nModule) {
		cout << "Something wrong!! modules.size() != nModule." << endl;
	}
}





/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in SeqInfomap.cpp.
 *
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 */

bool Network::moveSuperNodes() {

	int nSuperNodes = superNodes.size();

//	cout << "Num of SuperNode: " << nSuperNodes << endl;

//	for (int k = 0; k < modules.size(); k++) {
//		cout << "Module[" << k << "]: NumMembers =  " << modules[k].NumMembers();
//		cout << ", members.size() = " << modules[k].members.size() << endl;
//	}


	// Generate random sequential order of nodes.
	vector<int> randomOrder(nSuperNodes);	
	for (int i = 0; i < nSuperNodes; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nSuperNodes; i++) {
		int target = R->randInt(nSuperNodes - 1);
		//int target = rand() % nSuperNodes;

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	bool movedAny = false;

	// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nSuperNodes; i++) {
		//int nd = randomOrder[i];		// look at i_th Node of the random sequential order.
		SuperNode& nd = superNodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();
		//int oldMod = nd.members[0]->ModIdx();

		//cout << i+1 << ": Start looking for new mod of SuperNode[" << nd.ID() << "] - oldMod = " << oldMod << endl;

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		//double outFlowNotToOldMod = 0.0;	// Sum_b (P_nd * w_nd_b) where b not in oldMod.
		//double inFlowFromOldMod = 0.0;		// SUM_a (P_a * w_a_nd)  where a in oldMod.

		//map<int, double> outFlowToMod;	// <modID, flow> for outFlow...
		//map<int, double> inFlowFromMod;		// <modID, flow> for inFlow...
		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		// During this counting, aggregate outFlowNotToOldMod and inFlowFromOldMod values.
		//int nOutEdges = nd.outLinks.size();
		//for (int j = 0; j < nOutEdges; j++) {
		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
			//int newMod = superNodes[nd.outLinks[j].first].ModIdx();
			int newMod = superNodes[linkIt->first].ModIdx();
			//int newMod = superNodes[nd.outLinks[j].first].members[0]->ModIdx();

			//outFlowNotToOldMod += nd.outLinks[j].second;
			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = beta * linkIt->second;	// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		//int nInEdges = nd.inLinks.size();
		//for (int j = 0; j < nInEdges; j++) {
		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
			//int newMod = superNodes[nd.inLinks[j].first].ModIdx();
			int newMod = superNodes[linkIt->first].ModIdx();
			//int newMod = superNodes[nd.inLinks[j].first].members[0]->ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// Debug purpose...
		//cout << "Node[" << nd << "]: oldMod = " << oldMod << ", nOutEdges = " << nOutEdges << ", nInEdges = " << nInEdges << ", nModLinks = " << nModLinks << endl;


		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();					// p_nd.
		double ndTPWeight = nd.TeleportWeight();		// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].ExitPr();
		double oldSumPr1 = modules[oldMod].SumPr();
		double oldSumDangling1 = modules[oldMod].SumDangling();
		double oldModTPWeight = modules[oldMod].SumTPWeight();
		
		double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			}
			else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
				inFlowFromMod[newMod] += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > ndSize && emptyModules.size() > 0) {
			int emptyTarget = emptyModules.back();
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}


		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;

		//double newExitPr1 = oldExitPr1 + Network::alpha * (ndTPWeight * oldSumPr1 - (1.0 - modules[oldMod].SumTPWeight() + ndTPWeight) * ndSize);
		//newExitPr1 += Network::beta * (inFlowFromOldMod - outFlowNotToOldMod + ndTPWeight * oldSumDangling1);
		//if (nd.DanglingSize() > 0.0)	// Need to adjust teleport flow for dangling node.
		//	newExitPr1 -= Network::beta * (1.0 - modules[oldMod].SumTPWeight() + ndTPWeight) * nd.DanglingSize();
		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;	// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		// We don't need to check newMod != oldMod, since all newMod values are not equal to oldMod.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].ExitPr();
				double oldSumPr2 = modules[newMod].SumPr();
				//double oldSumDangling2 = modules[newMod].SumDangling();
				
				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				//currentResult.exitPr2 = oldExitPr2 + Network::alpha * ((1.0 - modules[newMod].SumTPWeight() - ndTPWeight) * ndSize - ndTPWeight * oldSumPr2);
				//currentResult.exitPr2 += Network::beta * (it->second - inFlowFromNewMod[newMod] - ndTPWeight * oldSumDangling2);
				//if (nd.DanglingSize() > 0.0)	// Need to adjust teleport flow for dangling node.
				//	currentResult.exitPr2 += Network::beta * (1.0 - modules[newMod].SumTPWeight() - ndTPWeight) * nd.DanglingSize();
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1 + currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1) + pLogP(currentResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(currentResult.exitPr1 + currentResult.sumPr1) + pLogP(currentResult.exitPr2 + currentResult.sumPr2) \
											- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
					bestResult.sumPr1 = currentResult.sumPr1;
					bestResult.sumPr2 = currentResult.sumPr2;
					bestResult.exitPr1 = currentResult.exitPr1;
					bestResult.exitPr2 = currentResult.exitPr2;
					bestResult.newSumExitPr = currentResult.newSumExitPr;
				}
			}
		}


		// Make best possible move for the current node nd.
		if (bestResult.diffCodeLen < 0.0) {
			// update related to newMod...
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

			//if (modules[newMod].members.size() == 0) {
			if (modules[newMod].NumMembers() == 0) {
				newMod = emptyModules.back();
				emptyModules.pop_back();
				nEmptyMod--;
				nModule++;
			}

			nd.setModIdx(newMod);

			for (int j = 0; j < spMembers; j++) {
				nd.members[j]->setModIdx(newMod);
				
				//modules[newMod].members.push_back(nd.members[j]);
			}

			modules[newMod].increaseNumMembers(spMembers);
			modules[newMod].setExitPr(bestResult.exitPr2);
			modules[newMod].setSumPr(bestResult.sumPr2);
			modules[newMod].setStayPr(bestResult.exitPr2 + bestResult.sumPr2);
			modules[newMod].addSumTPWeight(ndTPWeight);

			double spDanglingSize = nd.DanglingSize();
			if (spDanglingSize > 0.0) {
				modules[newMod].addSumDangling(spDanglingSize);
				modules[oldMod].minusSumDangling(spDanglingSize);
			}

			// update related to the oldMod...
			modules[oldMod].decreaseNumMembers(spMembers);
			modules[oldMod].setExitPr(bestResult.exitPr1);
			modules[oldMod].setSumPr(bestResult.sumPr1);
			modules[oldMod].setStayPr(bestResult.exitPr1 + bestResult.sumPr1);
			modules[oldMod].minusSumTPWeight(ndTPWeight);


			if (modules[oldMod].NumMembers() == 0) {
				nEmptyMod++;
				nModule--;
				emptyModules.push_back(oldMod);
				//cout << "Number of Empty Module: " << nEmptyMod << endl;
			}

			// remove the member nodes of the current superNode from oldMod
/*			int nMembers = modules[oldMod].members.size();
			if (nMembers == spMembers) {
				//for (int j = 0; j < nMembers; j++)
				//	modules[oldMod].members.pop_back();
				modules[oldMod].members.clear();
			}
			else {
				///////// DEBUG PURPOSE ////////
				//cout << "SuperNodes[" << nd.ID() << "] - remove member nodes from OldMod." << endl;

				set<int> removeCand;
				for (int j = 0; j < spMembers; j++)
					removeCand.insert(nd.members[j]->ID());

				int removeCount = 0;
				int backIdx = nMembers - 1;
				for (int n = 0; n < nMembers; n++) {
					if (removeCand.count(modules[oldMod].members[n]->ID()) > 0) {
						// swap the found node element and the back() element,
						// so that we can reduce the time complexity for erasing the node to the constant.
						// compare vector::erase() to vector::pop_back().
						Node * tempND = modules[oldMod].members[n];

						while ((n < backIdx) && (removeCand.count(modules[oldMod].members[backIdx]->ID()) > 0)) {
							modules[oldMod].members.pop_back();
							backIdx--;
							removeCount++;
						}

						if (removeCand.count(modules[oldMod].members[backIdx]->ID()) == 0)
							modules[oldMod].members[n] = modules[oldMod].members[backIdx];
							
						modules[oldMod].members[backIdx] = tempND;

						modules[oldMod].members.pop_back();		// remove the current node which swapped to the end of the vector.
						backIdx--;

						removeCount++;	
						if (removeCount == spMembers)
							break;
					}
				}

				if (removeCount != spMembers)
					cout << "SuperNode[" << nd << "] - removeCount != spMembers..." << endl;
			}
*/
			sumAllExitPr = bestResult.newSumExitPr;

			codeLength += bestResult.diffCodeLen;
			movedAny = true;
		}
		//cout << "Finish looking for new mod of SuperNode[" << nd.ID() << "]!!\n";
	}

	/////////// NEED TO CONCATENATE MODULE VECTOR ///////////////
/*	int idx = 0;
	int emptyIdx = modules.size() - 1;


	while (movedAny && (idx <= emptyIdx)) {

		if (modules[idx].NumMembers() == 0) {

			// find non-empty modules from the end and swap it with this.
			while (modules[emptyIdx].NumMembers() == 0) {
				emptyIdx--;
				modules.pop_back();
			}

			if (idx < emptyIdx) {
				//modules[idx] = modules[emptyIdx];	// copy a non-empty module at 'emptyIdx' to 'idx'
				//copyModule(*(modules[idx]), *(modules[emptyIdx]));
				copyModule(idx, emptyIdx);
				emptyIdx--;
				modules.pop_back();	// Didn't swap it, because the empty module should be removed.

				modules[idx].setIndex(idx);
				
				// update modIdx of each node in modules[idx] which was at modules[emptyIdx].
				int memSize = modules[idx].members.size();

				for (int i = 0; i < memSize; i++)
					modules[idx].members[i]->setModIdx(idx);
			}
		}
		idx++;

		if (nModule == modules.size())
			break;
	}

	for (int i = 0; i < nSuperNodes; i++) {
		int updatedMod = superNodes[i].members[0]->ModIdx();
		superNodes[i].setModIdx(updatedMod);
	}
*/
	// the following should be true: modules.size() == nModule.
	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod." << endl;
	}

	return movedAny;
}




/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in OmpInfomap.cpp.
 *
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 */

bool Network::parallelMoveSuperNodes(int numTh, double& tSequential) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	int nSuperNodes = superNodes.size();


	// Generate random sequential order of nodes.
	vector<int> randomOrder(nSuperNodes);	
	for (int i = 0; i < nSuperNodes; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nSuperNodes; i++) {
		int target = R->randInt(nSuperNodes - 1);
		//int target = rand() % nSuperNodes;

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	bool movedAny = false;

	omp_set_num_threads(numTh);

	// generate vectors which contains node-movements results from each thread.
	typedef pair<int, int> moveinfo;	// first = nodeIndex, second = new module index.
	vector<vector<moveinfo> > movements(numTh);

	const int emptyTarget = nNode + 1;	// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);

	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
	// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nSuperNodes; i++) {
		//int nd = randomOrder[i];		// look at i_th Node of the random sequential order.
		SuperNode& nd = superNodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		//double outFlowNotToOldMod = 0.0;	// Sum_b (P_nd * w_nd_b) where b not in oldMod.
		//double inFlowFromOldMod = 0.0;		// SUM_a (P_a * w_a_nd)  where a in oldMod.

		//map<int, double> outFlowToMod;	// <modID, flow> for outFlow...
		//map<int, double> inFlowFromMod;		// <modID, flow> for inFlow...
		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		// During this counting, aggregate outFlowNotToOldMod and inFlowFromOldMod values.
		//int nOutEdges = nd.outLinks.size();
		//for (int j = 0; j < nOutEdges; j++) {
		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
			//int newMod = superNodes[nd.outLinks[j].first].ModIdx();
			int newMod = superNodes[linkIt->first].ModIdx();
			//int newMod = superNodes[nd.outLinks[j].first].members[0]->ModIdx();

			//outFlowNotToOldMod += nd.outLinks[j].second;
			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = beta * linkIt->second;	// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		//int nInEdges = nd.inLinks.size();
		//for (int j = 0; j < nInEdges; j++) {
		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
			//int newMod = superNodes[nd.inLinks[j].first].ModIdx();
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// Debug purpose...
		//cout << "Node[" << nd << "]: oldMod = " << oldMod << ", nOutEdges = " << nOutEdges << ", nInEdges = " << nInEdges << ", nModLinks = " << nModLinks << endl;


		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();					// p_nd.
		double ndTPWeight = nd.TeleportWeight();		// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].ExitPr();
		double oldSumPr1 = modules[oldMod].SumPr();
		double oldSumDangling1 = modules[oldMod].SumDangling();
		double oldModTPWeight = modules[oldMod].SumTPWeight();
		
		double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			}
			else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
				inFlowFromMod[newMod] += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].SumPr() > ndSize && emptyModules.size() > 0) {
			//int emptyTarget = emptyModules.back();
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}


		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;	// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		// We don't need to check newMod != oldMod, since all newMod values are not equal to oldMod.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].ExitPr();
				double oldSumPr2 = modules[newMod].SumPr();
				//double oldSumDangling2 = modules[newMod].SumDangling();
				
				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1 + currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1) + pLogP(currentResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(currentResult.exitPr1 + currentResult.sumPr1) + pLogP(currentResult.exitPr2 + currentResult.sumPr2) \
											- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
					//bestResult.sumPr1 = currentResult.sumPr1;
					//bestResult.sumPr2 = currentResult.sumPr2;
					//bestResult.exitPr1 = currentResult.exitPr1;
					//bestResult.exitPr2 = currentResult.exitPr2;
					//bestResult.newSumExitPr = currentResult.newSumExitPr;
				}
			}

		}

		// store the best possilbe movement information if necessary.
		if (bestResult.diffCodeLen < 0.0) {
			//int thID = omp_get_thread_num();
			//movements[thID].push_back(make_pair(randomOrder[i], bestResult.newModule));

			bool isEmptyTarget = false;
			bool validMove = true;		// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

			#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);
			
				// if newMod == emptyTarget, it indicates moves to empty module.
				if ( (nEmptyMod > 0) && (newMod == emptyTarget) && (modules[oldMod].NumMembers() > 1) ) {
					newMod = emptyModules.back();
					//emptyModules.pop_back();
				
					//nEmptyMod--;
					//nModule++;
					isEmptyTarget = true;
				}
				else if (newMod == emptyTarget) {
					//continue;
					validMove = false;
				}
				//else if (modules[oldMod].NumMembers() == 1 && modules[newMod].NumMembers() == 0) {
				//	continue;	// This case don't need to move at all.
				//}
				else if (modules[newMod].NumMembers() == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there. 
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}


				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					//int nOutEdges = nd.outLinks.size();
					//int nInEdges = nd.inLinks.size();
	
					if (!isEmptyTarget) {
						for (link_iterator linkIt =nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();
	
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
							else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt =nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();
	
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
							else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					}
					else {
						for (link_iterator linkIt =nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt =nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}


					oldExitPr1 = modules[oldMod].ExitPr();
					oldSumPr1 = modules[oldMod].SumPr();
					oldSumDangling1 = modules[oldMod].SumDangling();
					oldModTPWeight = modules[oldMod].SumTPWeight();
		
					//double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
					//double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
					inFlowFromNewMod += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;


					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					//MoveSummary moveResult;

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;


					// copy module specific values...
					double oldExitPr2 = modules[newMod].ExitPr();
					double oldSumPr2 = modules[newMod].SumPr();
					//double oldSumDangling2 = modules[newMod].SumDangling();
				
					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1 + moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1) + pLogP(moveResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(moveResult.exitPr1 + moveResult.sumPr1) + pLogP(moveResult.exitPr2 + moveResult.sumPr2) \
												- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

//					if (moveResult.diffCodeLen >= 0.0)
//						//continue;	// ignore this movement option.
//						validMove = false;
//				}

//				if (validMove) {

					if (isEmptyTarget) {
						//newMod = emptyModules.back();		// Done before...
						emptyModules.pop_back();
				
						nEmptyMod--;
						nModule++;
					}


					nd.setModIdx(newMod);
					for (int k = 0; k < spMembers; k++) {
						nd.members[k]->setModIdx(newMod);
					}

					//modules[newMod].members.push_back(&nodes[ndIdx]);
					modules[newMod].increaseNumMembers(spMembers);
					modules[newMod].setExitPr(moveResult.exitPr2);
					modules[newMod].setSumPr(moveResult.sumPr2);
					modules[newMod].setStayPr(moveResult.exitPr2 + moveResult.sumPr2);
					modules[newMod].addSumTPWeight(ndTPWeight);

					if (ndDanglingSize > 0.0) {
						modules[newMod].addSumDangling(ndDanglingSize);
						modules[oldMod].minusSumDangling(ndDanglingSize);
					}

					// update related to the oldMod...
					modules[oldMod].decreaseNumMembers(spMembers);
					modules[oldMod].setExitPr(moveResult.exitPr1);
					modules[oldMod].setSumPr(moveResult.sumPr1);
					modules[oldMod].setStayPr(moveResult.exitPr1 + moveResult.sumPr1);
					modules[oldMod].minusSumTPWeight(ndTPWeight);

					if (modules[oldMod].NumMembers() == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;
			
					movedAny = true;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical





		}

	}	// END parallel for


	/////////// SEQUENTIAL PART STARTING FROM HERE !!!! ////////////
	//gettimeofday(&tStart, NULL);

	
	// Check if there any node is moved to a new module...
/*	for (int i = 0; i < numTh; i++) {
		if (!movedAny && movements[i].size() > 0) {
			movedAny = true;
			break;
		}
	}
*/

	/*
	 *	This part is the actual movement of the node-movement decision from each threads.
	 *	We will actually move nodes if the decision actually improves the MDL, and ignore, otherwise.
	 *	In this way, we may recompute the exitPr of the movement, 
	 *	but we don't need to call updateCodeLength, since we can update codeLength as same as in sequential mode.
	 */
/*	for (int i = 0; i < numTh; i++) {
		int numMove = movements[i].size();
		for (int j = 0; j < numMove; j++) {
			// update related to newMod...
			int ndIdx = movements[i][j].first;
			int oldMod = superNodes[ndIdx].ModIdx();
			int newMod = movements[i][j].second;
			int spMembers = superNodes[ndIdx].members.size();

			double ndSize = superNodes[ndIdx].Size();
			double ndTPWeight = superNodes[ndIdx].TeleportWeight();
			double spDanglingSize = superNodes[ndIdx].DanglingSize();

			bool isEmptyTarget = false;


			// if newMod == emptyTarget, it indicates moves to empty module.
			if ( (nEmptyMod > 0) && (newMod == emptyTarget) && (modules[oldMod].SumPr() > ndSize) ) {
				newMod = emptyModules.back();
				//emptyModules.pop_back();
				
				//nEmptyMod--;
				//nModule++;
				isEmptyTarget = true;
			}
			else if (newMod == emptyTarget) {
				continue;
			}
			//else if (modules[oldMod].SumPr() == ndSize && modules[newMod].NumMembers() == 0) {
			//	continue;	// This case don't need to move at all.
			//}
			else if (modules[newMod].NumMembers() == 0) {
				// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
				// but the nodes are all moved to other modules so there are no nodes in there. 
				// Thus, we don't know whether generating a new module will be better or not.
				// Therefore, don't move this option.
				continue;
			}

			////////////////////////////////////////////////////////////////////
			/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
			////////////////////////////////////////////////////////////////////

			double outFlowToOldMod = 0.0;
			double outFlowToNewMod = 0.0;
			double inFlowFromOldMod = 0.0;
			double inFlowFromNewMod = 0.0;

			int nOutEdges = superNodes[ndIdx].outLinks.size();
			int nInEdges = superNodes[ndIdx].inLinks.size();

			if (!isEmptyTarget) {
				for (int j = 0; j < nOutEdges; j++) {
					int toMod = superNodes[superNodes[ndIdx].outLinks[j].first].ModIdx();

					if (toMod == oldMod) {
						outFlowToOldMod += beta * superNodes[ndIdx].outLinks[j].second;
					}
					else if (toMod == newMod) {
						outFlowToNewMod += beta * superNodes[ndIdx].outLinks[j].second;
					}
				}

				for (int j = 0; j < nInEdges; j++) {
					int fromMod = superNodes[superNodes[ndIdx].inLinks[j].first].ModIdx();

					if (fromMod == oldMod) {
						inFlowFromOldMod += beta * superNodes[ndIdx].inLinks[j].second;
					}
					else if (fromMod == newMod) {
						inFlowFromNewMod += beta * superNodes[ndIdx].inLinks[j].second;
					}
				}
			}
			else {
				for (int j = 0; j < nOutEdges; j++) {
					int toMod = superNodes[superNodes[ndIdx].outLinks[j].first].ModIdx();
					if (toMod == oldMod) 
						outFlowToOldMod += beta * superNodes[ndIdx].outLinks[j].second;
				}

				for (int j = 0; j < nInEdges; j++) {
					int fromMod = superNodes[superNodes[ndIdx].inLinks[j].first].ModIdx();
					if (fromMod == oldMod) 
						inFlowFromOldMod += beta * superNodes[ndIdx].inLinks[j].second;
				}
			}


			double oldExitPr1 = modules[oldMod].ExitPr();
			double oldSumPr1 = modules[oldMod].SumPr();
			double oldSumDangling1 = modules[oldMod].SumDangling();
			double oldModTPWeight = modules[oldMod].SumTPWeight();
		
			//double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
			//double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

			// For teleportation and danling nodes.
			outFlowToOldMod += (alpha * ndSize + beta * spDanglingSize) * (oldModTPWeight - ndTPWeight);
			inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - spDanglingSize)) * ndTPWeight;
			outFlowToNewMod += (alpha * ndSize + beta * spDanglingSize) * modules[newMod].SumTPWeight();
			inFlowFromNewMod += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;


			if (isEmptyTarget) {
				outFlowToNewMod = 0.0;
				inFlowFromNewMod = 0.0;
			}

			MoveSummary moveResult;

			moveResult.exitPr1 = oldExitPr1 - superNodes[ndIdx].ExitPr() + outFlowToOldMod + inFlowFromOldMod;


			// copy module specific values...
			double oldExitPr2 = modules[newMod].ExitPr();
			double oldSumPr2 = modules[newMod].SumPr();
			//double oldSumDangling2 = modules[newMod].SumDangling();
				
			// Calculate status of current investigated movement of the node nd.
			moveResult.newModule = newMod;
			moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
			moveResult.sumPr2 = oldSumPr2 + ndSize;
			moveResult.exitPr2 = oldExitPr2 + superNodes[ndIdx].ExitPr() - outFlowToNewMod - inFlowFromNewMod;

			moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1 + moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

			// Calculate delta_L(M) = L(M)_new - L(M)_old
			double delta_allExit_log_allExit = pLogP(moveResult.newSumExitPr) - pLogP(sumAllExitPr);
			double delta_exit_log_exit = pLogP(moveResult.exitPr1) + pLogP(moveResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
			double delta_stay_log_stay = pLogP(moveResult.exitPr1 + moveResult.sumPr1) + pLogP(moveResult.exitPr2 + moveResult.sumPr2) \
										- pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

			// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
			moveResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

			////////////////////////////////////
			////// THE END OF EXAMINATION //////
			////////////////////////////////////

			if (moveResult.diffCodeLen >= 0.0)
				continue;	// ignore this movement option, since it makes worse.


			if (isEmptyTarget) {
				//newMod = emptyModules.back();		// Done before...
				emptyModules.pop_back();
				
				nEmptyMod--;
				nModule++;
			}




			superNodes[ndIdx].setModIdx(newMod);

			for (int k = 0; k < spMembers; k++) {
				superNodes[ndIdx].members[k]->setModIdx(newMod);
				
				//modules[newMod].members.push_back(superNodes[ndIdx].members[k]);
			}

			modules[newMod].increaseNumMembers(spMembers);
			modules[newMod].setExitPr(moveResult.exitPr2);
			modules[newMod].setSumPr(moveResult.sumPr2);
			modules[newMod].setStayPr(moveResult.exitPr2 + moveResult.sumPr2);
			modules[newMod].addSumTPWeight(ndTPWeight);

			if (spDanglingSize > 0.0) {
				modules[newMod].addSumDangling(spDanglingSize);
				modules[oldMod].minusSumDangling(spDanglingSize);
			}

			// update related to the oldMod...
			modules[oldMod].decreaseNumMembers(spMembers);
			modules[oldMod].setExitPr(moveResult.exitPr1);
			modules[oldMod].setSumPr(moveResult.sumPr1);
			modules[oldMod].setStayPr(moveResult.exitPr1 + moveResult.sumPr1);
			modules[oldMod].minusSumTPWeight(ndTPWeight);


			if (modules[oldMod].NumMembers() == 0) {
				nEmptyMod++;
				nModule--;
				emptyModules.push_back(oldMod);
			}

			sumAllExitPr = moveResult.newSumExitPr;
			codeLength += moveResult.diffCodeLen;
		}
	}
*/


	/////////// NEED TO CONCATENATE MODULE VECTOR ///////////////
/*	int idx = 0;
	int emptyIdx = modules.size() - 1;


	while (movedAny && (idx <= emptyIdx)) {

		if (modules[idx].NumMembers() == 0) {

			// find non-empty modules from the end and swap it with this.
			while (modules[emptyIdx].NumMembers() == 0) {
				emptyIdx--;
				modules.pop_back();
			}

			if (idx < emptyIdx) {
				//modules[idx] = modules[emptyIdx];	// copy a non-empty module at 'emptyIdx' to 'idx'
				//copyModule(*(modules[idx]), *(modules[emptyIdx]));
				copyModule(idx, emptyIdx);
				emptyIdx--;
				modules.pop_back();	// Didn't swap it, because the empty module should be removed.

				modules[idx].setIndex(idx);
				
				// update modIdx of each node in modules[idx] which was at modules[emptyIdx].
				int memSize = modules[idx].members.size();

				for (int i = 0; i < memSize; i++)
					modules[idx].members[i]->setModIdx(idx);
			}
		}
		idx++;

		if (nModule == modules.size())
			break;
	}

	for (int i = 0; i < nSuperNodes; i++) {
		int updatedMod = superNodes[i].members[0]->ModIdx();
		superNodes[i].setModIdx(updatedMod);
	}
*/
	// the following should be true: modules.size() == nModule.
	//if (nodes.size() != nModule + nEmptyMod) {
	//	cout << "Something wrong!! nodes.size() != nModule + nEmptyMod." << endl;
	//}


	//gettimeofday(&tStart, NULL);
	//updateSPMembersInModule();
	//updateCodeLength(numTh, true);
//	gettimeofday(&tEnd, NULL);

//	tSequential += elapsedTimeInSec(tStart, tEnd);

	return movedAny;
}

/**
 *	This function will update members vector in modules correspondingly.
 */
void Network::updateMembersInModule() {
	int numModules = modules.size();
	//vector<int>().swap(activeModules);
	//activeModules.reserve(nModule);
	vector<int>().swap(smActiveMods);
	smActiveMods.reserve(nModule);
	vector<int>().swap(lgActiveMods);
	lgActiveMods.reserve(nModule);

	for (int i = 0; i < numModules; i++) {
		vector<Node*>().swap(modules[i].members);
		//if(modules[i].NumMembers() > 0)
		//	activeModules.push_back(i);
		if(modules[i].NumMembers() > 1000)
			lgActiveMods.push_back(i);
		else if(modules[i].NumMembers() > 0)
			smActiveMods.push_back(i);
	}

	for (int i = 0; i < nNode; i++)
		modules[nodes[i].ModIdx()].members.push_back(&nodes[i]);
}

void Network::updateSPMembersInModule() {
	int numModules = modules.size();

	for (int i = 0; i < numModules; i++)
		vector<Node*>().swap(modules[i].members);

	int numSPNodes = superNodes.size();
	for (int i = 0; i < numSPNodes; i++)
		modules[superNodes[i].ModIdx()].members.push_back(&superNodes[i]);
}



// calculate exit-probability based on node information and current module assignment.
void Network::updateCodeLength(int numTh, bool isSPNode) {

	// calculate exit-probability for each module.
	//vector<double> calculatedExitProb(nModule, 0.0);
	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;

	//updateMembersInModule();
	omp_set_num_threads(numTh);

#pragma omp parallel for reduction (+:tempSumAllExit, exit_log_exit, stay_log_stay)
	//for (unsigned int i = 0; i < nModule; i++) {
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].NumMembers() > 0) {
			// exitPr w.r.t. teleportation.
			double exitPr = Network::alpha * (1.0 - modules[i].SumTPWeight()) * modules[i].SumPr();
			int curMod = modules[i].Index();

			double sumOutFlow = 0.0;
			int nMembers = modules[i].members.size();
			for (int j = 0; j < nMembers; j++) {
				//int nOutLinks = modules[i].members[j]->outLinks.size();
				Node * nd = modules[i].members[j];
				int nOutLinks = nd->outLinks.size();

				for (int k = 0; k < nOutLinks; k++) {
					if (!isSPNode && nodes[nd->outLinks[k].first].ModIdx() != curMod)
						sumOutFlow += nd->outLinks[k].second;
					else if (isSPNode && superNodes[nd->outLinks[k].first].ModIdx() != curMod)
						sumOutFlow += nd->outLinks[k].second;
					/*if (!isSPNode) {
						if (nodes[nd->outLinks[k].first].ModIdx() != curMod)
							sumOutFlow += nd->outLinks[k].second;
					}
					else {
						if (superNodes[nd->outLinks[k].first].ModIdx() != curMod)
							sumOutFlow += nd->outLinks[k].second;
					}*/
				}
			}

			exitPr += Network::beta * (sumOutFlow + (1.0 - modules[i].SumTPWeight()) * modules[i].SumDangling());

			modules[i].setExitPr(exitPr);
			modules[i].setStayPr(exitPr + modules[i].SumPr());

			tempSumAllExit += exitPr;
			exit_log_exit += pLogP(exitPr);
			stay_log_stay += pLogP(modules[i].StayPr());

			//calculatedExitProb[i] = exitPr;
		}
	}

	sumAllExitPr = tempSumAllExit;

	//codeLength = pLogP(tempSumAllExit) - 2.0 * exit_log_exit + stay_log_stay - allNodes_log_allNodes;
	codeLength = pLogP(sumAllExitPr) - 2.0 * exit_log_exit + stay_log_stay - allNodes_log_allNodes;
}


// calculate exit-probability based on node information and current module assignment.
double Network::calculateCodeLength() {

	// calculate exit-probability for each module.
	//vector<double> calculatedExitProb(nModule, 0.0);
	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;


	//for (unsigned int i = 0; i < nModule; i++) {
#pragma omp parallel for reduction (+:tempSumAllExit, exit_log_exit, stay_log_stay)
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].NumMembers() > 0) {
			// exitPr w.r.t. teleportation.
			double exitPr = Network::alpha * (1.0 - modules[i].SumTPWeight()) * modules[i].SumPr();
			int curMod = modules[i].Index();

			double sumOutFlow = 0.0;
			int nMembers = modules[i].members.size();
			for (int j = 0; j < nMembers; j++) {
				int nOutLinks = modules[i].members[j]->outLinks.size();
				for (int k = 0; k < nOutLinks; k++) {
					if (nodes[modules[i].members[j]->outLinks[k].first].ModIdx() != curMod)
						sumOutFlow += modules[i].members[j]->outLinks[k].second;
				}
			}

			exitPr += Network::beta * (sumOutFlow + (1.0 - modules[i].SumTPWeight()) * modules[i].SumDangling());

			tempSumAllExit += exitPr;
			exit_log_exit += pLogP(exitPr);
			stay_log_stay += pLogP(exitPr + modules[i].SumPr());

			//calculatedExitProb[i] = exitPr;
		}
	}

	double computedCodeLength = pLogP(tempSumAllExit) - 2.0 * exit_log_exit + stay_log_stay - allNodes_log_allNodes;

	return computedCodeLength;
}


// similar function of Greedy::level() in the original infomap_dir implementation.
// make a group of nodes as a Node (here, SuperNode)
// Granularity of the network is same but a group of nodes move together instead of moving individually.
void Network::convertModulesToSuperNodes(int numTh) {
	//initialize superNodes vector for updating w.r.t. the current module status.
	vector<SuperNode>().swap(superNodes);
	superNodes.reserve(nModule);

	vector<unsigned int> modToSPNode(nNode);	// indicate corresponding superNode ID from each module.

	int idx = 0;
	//for (unsigned int i = 0; i < nModule; i++) {
	for (unsigned int i = 0; i < nNode; i++) {
		//superNodes.push_back(SuperNode(modules[i], i));
		if (modules[i].NumMembers() > 0) {
			superNodes.push_back(SuperNode(modules[i], idx));
			modToSPNode[i] = idx;
			idx++;
		}
	}

	/*
	 *	Calculate outLinks and inLinks between superNodes...
	 */
	int numSPNode = superNodes.size();

	omp_set_num_threads(numTh);

//#pragma omp parallel for
//	for (int i = 0; i < numSPNode; i++) {
	#pragma omp parallel
	{
		//double* newOutLinks = new double[numSPNode];
		//double* newInLinks = new double[numSPNode];


		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, numSPNode, nTh, myID);

		for (int i = start; i < end; i++) {

			int numNodesInSPNode = superNodes[i].members.size();

			//for (int j = 0; j < numSPNode; j++) {
			//	newOutLinks[j] = 0.0;
			//	//newInLinks[j] = 0.0;
			//}

			//set<int> toList;	// set of nodes reached to by outLinks. 
			////set<int> fromList;	// set of nodes connected by inLinks.
			typedef map<int, double> EdgeMap;
			EdgeMap newOutLinks;

			for (int j = 0; j < numNodesInSPNode; j++) {
				// Calculate newOutLinks from a superNode to other superNodes.
				Node* nd = superNodes[i].members[j];
				int nOutEdge = nd->outLinks.size();

				for (int k = 0; k < nOutEdge; k++) {
					int toSPNode = modToSPNode[nodes[nd->outLinks[k].first].ModIdx()];
			
					if (toSPNode != i) {	// outLinks to another superNode...
						//newOutLinks[toSPNode] += nd->outLinks[k].second;
						//toList.insert(toSPNode);
						pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(make_pair(toSPNode, nd->outLinks[k].second));
						if (!ret.second)
							ret.first->second += nd->outLinks[k].second;
					}
				}

				// Calculate newInLinks to a superNode from other superNodes.
/*				int nInEdge = nd->inLinks.size();

				for (int k = 0; k < nInEdge; k++) {
					int fromSPNode = modToSPNode[nodes[nd->inLinks[k].first].ModIdx()];

					if (fromSPNode != i) {	// inLinks from another superNode...
						newInLinks[fromSPNode] += nd->inLinks[k].second;
						fromList.insert(fromSPNode);
					}
				}
*/
			}
		
			superNodes[i].outLinks.reserve(newOutLinks.size());
			//superNodes[i].outLinks.reserve(toList.size());
			//superNodes[i].inLinks.reserve(fromList.size());

			//for (set<int>::iterator it = toList.begin(); it != toList.end(); it++)
			//	superNodes[i].outLinks.push_back(make_pair(*it, newOutLinks[*it]));
			for (EdgeMap::iterator it = newOutLinks.begin(); it != newOutLinks.end(); it++) {
				superNodes[i].outLinks.push_back(make_pair(it->first, it->second));
				//#pragma omp critical
				//superNodes[it->first].inLinks.push_back(make_pair(i, it->second));
			}

			/*#pragma omp critical
			{
				for (EdgeMap::iterator it = newOutLinks.begin(); it != newOutLinks.end(); it++) {
					superNodes[it->first].inLinks.push_back(make_pair(i, it->second));
				}
			}*/

			//for (set<int>::iterator it = fromList.begin(); it != fromList.end(); it++)
			//	superNodes[i].inLinks.push_back(make_pair(*it, newInLinks[*it]));

			/*for (int j = 0; j < numSPNode; j++) {
				if (newOutLinks[j] > 0.0)
					superNodes[i].outLinks.push_back(make_pair(j, newOutLinks[j]));

				if (newInLinks[j] > 0.0)
					superNodes[i].inLinks.push_back(make_pair(j, newInLinks[j]));
			}*/
		}

		//delete [] newOutLinks;
		//delete [] newInLinks;
	}

	// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) {
		int nOutLinks = superNodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
			superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(make_pair(i,superNodes[i].outLinks[j].second));
	}

}




//void Network::generateSuperNodesFromSubModules() {
void Network::generateSuperNodesFromSubModules(int numTh) {
	//initialize superNodes vector for updating w.r.t. the current module status.
	vector<SuperNode>().swap(superNodes);

	int nSubMods = subModules.size();
	superNodes.reserve(nSubMods);

	//vector<unsigned int> modToSPNode(nNode);	// indicate corresponding superNode ID from each module.

	for (int i = 0; i < nSubMods; i++) {
		superNodes.push_back(SuperNode(subModules[i], i, *this));
		//modToSPNode[i] = i;

		/// DEBUG PURPOSE...
/*		cout << "SubModule[" << i << "]: members - ";
		for (vector<int>::iterator it = subModules[i].members.begin(); it != subModules[i].members.end(); it++)
			cout << (*it) << ", ";
		cout << endl;
*/
	}

	/*
	 *	Calculate outLinks and inLinks between superNodes...
	 */
	int numSPNode = superNodes.size();

	//cout << "nSubMods == " << nSubMods << ", numSPNode == " << numSPNode << endl;

	typedef map<int, double> EdgeMap;

	omp_set_num_threads(numTh);

	//#pragma omp parallel
	//{
	//	int myID = omp_get_thread_num();
	//	int nTh = omp_get_num_threads();

	//	int start, end;
	//	findAssignedPart(&start, &end, numSPNode, nTh, myID);

	//	for (int i = start; i < end; i++) {
	#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < numSPNode; i++) {
			//cout << "Converting sub-Module[" << i << "] to superNode starting... " << flush;

			SuperNode* spNode = &(superNodes[i]);

			//int numNodesInSPNode = superNodes[i].members.size();
			int numNodesInSPNode = spNode->members.size();
			//cout << "# nodes in SPNode = " << numNodesInSPNode << " ... " << flush;

			EdgeMap newOutLinks;

			//cout << "before inner for-loop ..." << endl;
			for (int j = 0; j < numNodesInSPNode; j++) {
				// Calculate newOutLinks from a superNode to other superNodes.
				//Node* nd = superNodes[i].members[j];
				Node* nd = spNode->members[j];
				int nOutEdge = nd->outLinks.size();

				//cout << "node[" << nd->ID() <<"] - # outEdges = " << nOutEdge << ": " << flush;

				for (int k = 0; k < nOutEdge; k++) {
					//int toSPNode = modToSPNode[nodes[nd->outLinks[k].first].ModIdx()];
					int toSPNode = ndToSubMod[nd->outLinks[k].first];
					//cout << toSPNode << " " << flush;
			
					if (toSPNode != i) {	// outLinks to another superNode...
						pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(make_pair(toSPNode, nd->outLinks[k].second));
						if (!ret.second)
							ret.first->second += nd->outLinks[k].second;
					}
				}
				//cout << endl;
			}
			//cout << "after inner for-loop ..." << flush;
		
			//superNodes[i].outLinks.reserve(newOutLinks.size());
			spNode->outLinks.reserve(newOutLinks.size());

			double sumExitFlow = 0.0;

			for (EdgeMap::iterator it = newOutLinks.begin(); it != newOutLinks.end(); it++) {
				//superNodes[i].outLinks.push_back(make_pair(it->first, it->second));
				spNode->outLinks.push_back(make_pair(it->first, it->second));
				sumExitFlow += it->second;
			}

			double teleportExitFlow = 0.0;
			for (int j = 0; j < numNodesInSPNode; j++) {
				Node* nd = spNode->members[j];
				teleportExitFlow += (alpha * nd->Size() + beta * nd->DanglingSize()) * (1.0 - nd->TeleportWeight());
			}


			if (sumExitFlow == 0.0) {
				//double exitProb = (1.0 - spNode->TeleportWeight()) * spNode->Size();
				//spNode->setExitPr(exitProb);
				spNode->setExitPr(teleportExitFlow);
			} else {
				//double exitProb = alpha * (1.0 - spNode->TeleportWeight()) * spNode->Size() + beta * sumExitFlow;
				double exitProb = teleportExitFlow + beta * sumExitFlow;
				spNode->setExitPr(exitProb);
			}

			//cout << "Converting sub-Module[" << i << "] to superNode done!" << endl;
		}
	//}

	// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) {
		int nOutLinks = superNodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
			superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(make_pair(i,superNodes[i].outLinks[j].second));
	}
}









void Network::copyModule(Module * newM, Module * oldM) {
	newM->setIndex(oldM->Index());
	newM->setExitPr(oldM->ExitPr());
	newM->setStayPr(oldM->StayPr());
	newM->setSumPr(oldM->SumPr());
	newM->setSumTPWeight(oldM->SumTPWeight());
	newM->setSumDangling(oldM->SumDangling());
	newM->setNumMembers(oldM->NumMembers());

	//vector<Node *>().swap(newM->members);
	newM->members.clear();
	int nMembers = oldM->members.size();
	for (int i = 0; i < nMembers; i++)
		newM->members.push_back(oldM->members[i]);
}

void Network::copyModule(int newM, int oldM) {
	modules[newM].setIndex(modules[oldM].Index());
	modules[newM].setExitPr(modules[oldM].ExitPr());
	modules[newM].setStayPr(modules[oldM].StayPr());
	modules[newM].setSumPr(modules[oldM].SumPr());
	modules[newM].setSumTPWeight(modules[oldM].SumTPWeight());
	modules[newM].setSumDangling(modules[oldM].SumDangling());
	modules[newM].setNumMembers(modules[oldM].NumMembers());

	//vector<Node *>().swap(modules[newM].members);
	modules[newM].members.clear();
	int nMembers = modules[oldM].members.size();
	for (int i = 0; i < nMembers; i++)
		modules[newM].members.push_back(modules[oldM].members[i]);
}


// Miscellaneous function..
void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID) {
	if (numNodes % numTh == 0) {
		*start = (numNodes / numTh) * myID;
		*end = (numNodes/ numTh) * (myID + 1);
	}
	else {
		int block = numNodes / numTh;
		int modular = numNodes % numTh;

		if (myID < modular) {
			*start = (block + 1) * myID;
			*end = (block + 1) * (myID + 1);
		}
		else {
			*start = block * myID + modular;
			*end = block * (myID + 1) + modular;
		}
	}
}
