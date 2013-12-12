/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Dec. 2013
 *	Copyright (C) 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

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


typedef map<int, double> flowmap;
typedef map<int, pair<double, double> > modInfo;	// <modID, <exitPr, sumPr> >


void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID);

const double InitLength = 10000.0;

using namespace std;

struct MoveSummary {
	double diffCodeLen;			// delta(L(M))
	int newModule;
	double sumPr1, sumPr2;		// updated sumPr1 and sumPr2.  1 --> oldM, 2 --> newM.
	double exitPr1, exitPr2;	// updated exitPr1 and exitPr2.
	double newSumExitPr;			// SUM(q_i) = this.sumAllExitPr + q_1_new + q_2_new - q_1_old - q_2_old.
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
	double sumExitFlow = 0.0;
	int nOutLinks = nd->outLinks.size();

	for (int i = 0; i < nOutLinks; i++)
		sumExitFlow += nd->outLinks[i].second;	// after call initiate(), w_ab is updated as p_a * w_ab.
		
	if (nd->IsSuper()) {	// If SuperNode, we already calculated exitPr in the corresponding module.
		exitPr = nd->ExitPr();
		sumDangling = nd->DanglingSize();
	}
	// exitPr = tau * (1 - sum(tau_a))*sum(p_a) + (1-tau) * sumExitFlow.
	else if (!nd->IsDangling()) {
		exitPr = Network::alpha * (1.0 - nd->TeleportWeight()) * sumPr + Network::beta * sumExitFlow;
		nd->setExitPr(exitPr);
	}
	else {
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

	for (vector<Node*>::iterator it = mod.members.begin(); it != mod.members.end(); it++) {
		members.push_back(origNodeID[(*it)->ID()]);
	}

	// numMembers should be equal to members.size() ...
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
		if (nodes[i].outLinks.empty()) {
			danglings.push_back(i);
			nDanglings++;
			nodes[i].setIsDangling(true);
		}
	}
	
}

//void Network::initiate() {
void Network::initiate(int numTh) {

	int nDangNodes = 0;

	struct timeval startT, endT;

	omp_set_num_threads(numTh);

	gettimeofday(&startT, NULL);

	#pragma omp parallel reduction (+:nDangNodes) 
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		for (int i = start; i < end; i++) {
			if (nodes[i].outLinks.empty()) {
				#pragma omp critical (updateDang)
				{
					danglings.push_back(i);
				}
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
	calculateSteadyState(numTh);
	gettimeofday(&endT, NULL);
	cout << "Time for calculating steady state of nodes (eigenvector): " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;


	gettimeofday(&startT, NULL);

	// Update edges to represent flow based on calculated steady state (aka size).
	#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		for (int i = start; i < end; i++) {

			if (!nodes[i].IsDangling()) {
				int nOutLinks = nodes[i].outLinks.size();
				for (int j = 0; j < nOutLinks; j++) {
					nodes[i].outLinks[j].second = nodes[i].Size() * nodes[i].outLinks[j].second;
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
	double allNds_log_allNds = 0.0;

	#pragma omp parallel reduction (+:allNds_log_allNds)
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);
	
		for (int i = start; i < end; i++)
			allNds_log_allNds += pLogP(nodes[i].Size());
	}

	allNodes_log_allNodes = allNds_log_allNds;

	/////////////////////////////
	// Make modules from nodes //
	/////////////////////////////
	#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		for (int i = start; i < end; i++) {
			modules[i] = Module(i, &nodes[i]);	// Assign each Node to a corresponding Module object.
												// Initially, each node is in its own module.
			nodes[i].setModIdx(i);
		}
	}
	nModule = nNode;

	gettimeofday(&startT, NULL);
	calibrate(numTh);
	gettimeofday(&endT, NULL);
	cout << "Time for calculating of initial code length: " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;

}


// calculating steady state of nodes by Power Iteration Method.
// same as Greedy::eigenvector() in Infomap implementation.
// Modify for the design of this implementation.

void Network::calculateSteadyState(int numTh) {
	// initial probability distribution = 1 / N.
	vector<double> size_tmp = vector<double>(nNode, 1.0/nNode);

	int iter = 0;
	double danglingSize = 0.0;
	double sqdiff = 1.0;
	double sum = 0.0;

	// Generate addedSize array per each thread, so that we don't need to use lock for each nodeSize addition.
	double** addedSize = new double*[numTh];
	for (int i = 0; i < numTh; i++)
		addedSize[i] = new double[nNode];
			
	do {
		// calculate sum of the size of dangling nodes.
		danglingSize = 0.0;
	
		#pragma omp parallel reduction (+:danglingSize)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nDanglings, nTh, myID);
			
			for (int i = start; i < end; i++)
				danglingSize += size_tmp[danglings[i]];
		}

		// flow via teleportation.
		#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			for (int i = start; i < end; i++) 
				nodes[i].setSize( (alpha + beta*danglingSize) * nodes[i].TeleportWeight());
		}


		int realNumTh = 0;

		// flow from network steps via following edges.
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
			
			for (int i = start; i < end; i++) {
				int nOutLinks = nodes[i].outLinks.size();
				for (int j = 0; j < nOutLinks; j++) {
					myAddedSize[nodes[i].outLinks[j].first] +=  beta * nodes[i].outLinks[j].second * size_tmp[i];

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
			
			for (int i = start; i < end; i++) {
				for (int j = 0; j < realNumTh; j++)
					nodes[i].addSize(addedSize[j][i]);
			}
		}


		// Normalize of node size.
		sum = 0.0;
		#pragma omp parallel reduction (+:sum)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			for (int i = start; i < end; i++) 
				sum += nodes[i].Size();
		}
		sqdiff = 0.0;

		#pragma omp parallel reduction (+:sqdiff)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);
			
			for (int i = start; i < end; i++) {
				nodes[i].setSize(nodes[i].Size()/sum);
				sqdiff += fabs(nodes[i].Size() - size_tmp[i]);
				size_tmp[i] = nodes[i].Size();
			}
		}

		iter++;

	} while ((iter < 200) && (sqdiff > 1.0e-15 || iter < 50));

	// deallocate 2D array.
	for (int i = 0; i < numTh; i++)
		delete [] addedSize[i];
	delete [] addedSize;

	cout << "Calculating flow done in " << iter << " iterations!" << endl;

}


// This function calculate current codeLength.
// This implementation is modified version of infomap implementation.
void Network::calibrate(int numTh) {
	//This is the calculation of Equation (4) in the Map Equation paper.
	double sum_exit_log_exit = 0.0;
	double sum_stay_log_stay = 0.0;
	double sumExit = 0.0;

	omp_set_num_threads(numTh);

	#pragma omp parallel reduction (+:sum_exit_log_exit, sum_stay_log_stay, sumExit)
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nModule, nTh, myID);
			
		for (unsigned int i = start; i < end; i++) {
			sum_exit_log_exit += pLogP(modules[i].ExitPr());
			sum_stay_log_stay += pLogP(modules[i].StayPr());
			sumExit += modules[i].ExitPr();
		}
	}

	sumAllExitPr = sumExit;
	double sumExit_log_sumExit = pLogP(sumExit);

	codeLength = sumExit_log_sumExit - 2.0 * sum_exit_log_exit + sum_stay_log_stay - allNodes_log_allNodes;
}



/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in OmpRelaxmap.cpp.
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

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	bool movedAny = false;

	// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nNode; i++) {
		Node& nd = nodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;	// The number of links to/from between the current node and other modules.


		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
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

		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
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


		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			int emptyTarget = emptyModules.back();
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

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

			if (modules[newMod].NumMembers() == 0) {
				newMod = emptyModules.back();
				emptyModules.pop_back();
				nEmptyMod--;
				nModule++;
			}
			nd.setModIdx(newMod);

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
			}

			sumAllExitPr = bestResult.newSumExitPr;

			codeLength += bestResult.diffCodeLen;
			movedAny = true;
		}
	}


	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod." << endl;
	}

	return movedAny;
}



/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in OmpRelaxmap.cpp.
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
		Node& nd = nodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;	// The number of links to/from between the current node and other modules.


		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
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

		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
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
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;



		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;	// This is the default value. If we can't find diffCodeLen < 0, then don't move the node.

		// We don't need to check newMod != oldMod, since all newMod values are not equal to oldMod.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].ExitPr();
				double oldSumPr2 = modules[newMod].SumPr();
				
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

				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}
		}

		// store the best possilbe movement information if necessary.
		// In this version, we are trying to move as soon as possible, i.e. right after decision making...
		if (bestResult.diffCodeLen < 0.0) {
			
			bool isEmptyTarget = false;
			bool validMove = true;		// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;

			#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);
			
				if ( (nEmptyMod > 0) && (newMod == emptyTarget) && (modules[oldMod].NumMembers() > 1) ) {
					newMod = emptyModules.back();
				
					isEmptyTarget = true;
				}
				else if (newMod == emptyTarget) {
					validMove = false;
				}
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
		
					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
					inFlowFromNewMod += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;


					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;


					// copy module specific values...
					double oldExitPr2 = modules[newMod].ExitPr();
					double oldSumPr2 = modules[newMod].SumPr();
				
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

					moveResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();
				
						nEmptyMod--;
						nModule++;
					}


					nd.setModIdx(newMod);

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


	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod." << endl;
	}

	return movedAny;
}







bool Network::moveSuperNodes() {

	int nSuperNodes = superNodes.size();

	// Generate random sequential order of nodes.
	vector<int> randomOrder(nSuperNodes);	
	for (int i = 0; i < nSuperNodes; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nSuperNodes; i++) {
		int target = R->randInt(nSuperNodes - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	bool movedAny = false;

	// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nSuperNodes; i++) {
		SuperNode& nd = superNodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();


		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = beta * linkIt->second;	// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
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

		if (modules[oldMod].members.size() > ndSize && emptyModules.size() > 0) {
			int emptyTarget = emptyModules.back();
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}


		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;	

		// We don't need to check newMod != oldMod, since all newMod values are not equal to oldMod.
		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].ExitPr();
				double oldSumPr2 = modules[newMod].SumPr();
				
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

			if (modules[newMod].NumMembers() == 0) {
				newMod = emptyModules.back();
				emptyModules.pop_back();
				nEmptyMod--;
				nModule++;
			}

			nd.setModIdx(newMod);

			for (int j = 0; j < spMembers; j++) {
				nd.members[j]->setModIdx(newMod);
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
			}

			sumAllExitPr = bestResult.newSumExitPr;

			codeLength += bestResult.diffCodeLen;
			movedAny = true;
		}
	}

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
		SuperNode& nd = superNodes[randomOrder[i]];		// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			}
			else {
				outFlowToMod[newMod] = beta * linkIt->second;	// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) {
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

		if (modules[oldMod].SumPr() > ndSize && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}


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

				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}

		}

		// store the best possilbe movement information if necessary.
		if (bestResult.diffCodeLen < 0.0) {

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
					isEmptyTarget = true;
				}
				else if (newMod == emptyTarget) {
					validMove = false;
				}
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
		
					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].SumTPWeight();
					inFlowFromNewMod += (alpha * modules[newMod].SumPr() + beta * modules[newMod].SumDangling()) * ndTPWeight;


					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;


					// copy module specific values...
					double oldExitPr2 = modules[newMod].ExitPr();
					double oldSumPr2 = modules[newMod].SumPr();
				
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

					moveResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();
				
						nEmptyMod--;
						nModule++;
					}


					nd.setModIdx(newMod);
					for (int k = 0; k < spMembers; k++) {
						nd.members[k]->setModIdx(newMod);
					}

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

	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;

	omp_set_num_threads(numTh);

#pragma omp parallel for reduction (+:tempSumAllExit, exit_log_exit, stay_log_stay)
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].NumMembers() > 0) {
			// exitPr w.r.t. teleportation.
			double exitPr = Network::alpha * (1.0 - modules[i].SumTPWeight()) * modules[i].SumPr();
			int curMod = modules[i].Index();

			double sumOutFlow = 0.0;
			int nMembers = modules[i].members.size();
			for (int j = 0; j < nMembers; j++) {
				Node * nd = modules[i].members[j];
				int nOutLinks = nd->outLinks.size();

				for (int k = 0; k < nOutLinks; k++) {
					if (!isSPNode && nodes[nd->outLinks[k].first].ModIdx() != curMod)
						sumOutFlow += nd->outLinks[k].second;
					else if (isSPNode && superNodes[nd->outLinks[k].first].ModIdx() != curMod)
						sumOutFlow += nd->outLinks[k].second;
				}
			}

			exitPr += Network::beta * (sumOutFlow + (1.0 - modules[i].SumTPWeight()) * modules[i].SumDangling());

			modules[i].setExitPr(exitPr);
			modules[i].setStayPr(exitPr + modules[i].SumPr());

			tempSumAllExit += exitPr;
			exit_log_exit += pLogP(exitPr);
			stay_log_stay += pLogP(modules[i].StayPr());
		}
	}

	sumAllExitPr = tempSumAllExit;

	codeLength = pLogP(sumAllExitPr) - 2.0 * exit_log_exit + stay_log_stay - allNodes_log_allNodes;
}


// calculate exit-probability based on node information and current module assignment.
double Network::calculateCodeLength() {

	// calculate exit-probability for each module.
	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;


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
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].NumMembers() > 0) {
			superNodes.push_back(SuperNode(modules[i], idx));
			modToSPNode[i] = idx;
			idx++;
		}
	}

	int numSPNode = superNodes.size();

	omp_set_num_threads(numTh);

	#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, numSPNode, nTh, myID);

		for (int i = start; i < end; i++) {

			int numNodesInSPNode = superNodes[i].members.size();

			typedef map<int, double> EdgeMap;
			EdgeMap newOutLinks;

			for (int j = 0; j < numNodesInSPNode; j++) {
				// Calculate newOutLinks from a superNode to other superNodes.
				Node* nd = superNodes[i].members[j];
				int nOutEdge = nd->outLinks.size();

				for (int k = 0; k < nOutEdge; k++) {
					int toSPNode = modToSPNode[nodes[nd->outLinks[k].first].ModIdx()];
			
					if (toSPNode != i) {	// outLinks to another superNode...
						pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(make_pair(toSPNode, nd->outLinks[k].second));
						if (!ret.second)
							ret.first->second += nd->outLinks[k].second;
					}
				}
			}
		
			superNodes[i].outLinks.reserve(newOutLinks.size());

			for (EdgeMap::iterator it = newOutLinks.begin(); it != newOutLinks.end(); it++) {
				superNodes[i].outLinks.push_back(make_pair(it->first, it->second));
			}
		}
	}

	// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) {
		int nOutLinks = superNodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
			superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(make_pair(i,superNodes[i].outLinks[j].second));
	}

}




void Network::generateSuperNodesFromSubModules(int numTh) {
	//initialize superNodes vector for updating w.r.t. the current module status.
	vector<SuperNode>().swap(superNodes);

	int nSubMods = subModules.size();
	superNodes.reserve(nSubMods);

	for (int i = 0; i < nSubMods; i++) {
		superNodes.push_back(SuperNode(subModules[i], i, *this));
	}

	int numSPNode = superNodes.size();

	typedef map<int, double> EdgeMap;

	omp_set_num_threads(numTh);

	#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < numSPNode; i++) {
			SuperNode* spNode = &(superNodes[i]);

			int numNodesInSPNode = spNode->members.size();

			EdgeMap newOutLinks;

			for (int j = 0; j < numNodesInSPNode; j++) {
				Node* nd = spNode->members[j];
				int nOutEdge = nd->outLinks.size();

				for (int k = 0; k < nOutEdge; k++) {
					int toSPNode = ndToSubMod[nd->outLinks[k].first];
			
					if (toSPNode != i) {	// outLinks to another superNode...
						pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(make_pair(toSPNode, nd->outLinks[k].second));
						if (!ret.second)
							ret.first->second += nd->outLinks[k].second;
					}
				}
			}
		
			spNode->outLinks.reserve(newOutLinks.size());

			double sumExitFlow = 0.0;

			for (EdgeMap::iterator it = newOutLinks.begin(); it != newOutLinks.end(); it++) {
				spNode->outLinks.push_back(make_pair(it->first, it->second));
				sumExitFlow += it->second;
			}

			double teleportExitFlow = 0.0;
			for (int j = 0; j < numNodesInSPNode; j++) {
				Node* nd = spNode->members[j];
				teleportExitFlow += (alpha * nd->Size() + beta * nd->DanglingSize()) * (1.0 - nd->TeleportWeight());
			}


			if (sumExitFlow == 0.0) {
				spNode->setExitPr(teleportExitFlow);
			} else {
				double exitProb = teleportExitFlow + beta * sumExitFlow;
				spNode->setExitPr(exitProb);
			}
		}

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
