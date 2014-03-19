/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Mar. 2014
 *	Copyright (C) since 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <ctime>
#include <sys/time.h>
#include <omp.h>
#include "MersenneTwister.h"
#include "Node.h"
#include "Module.h"
#include "FileIO.h"
#include "timing.h"

using namespace std;

unsigned stou(char *s) {
	return strtoul(s, (char **)NULL, 10);
}


void stochastic_greedy_partition(Network &network, int numTh, double threshold, double vThresh, int maxIter, bool prior, bool fineTune, bool fast);
void partition_module_network(Network &network, int numTh, double threshold, int maxIter, bool fast);
void generate_sub_modules(Network &network, int numTh, double threshold, int maxIter);
void generate_network_from_module(Network &newNetwork, Module* mod, map<int, int>& origNodeID, int numTh);
void generate_network_from_module(Network &newNetwork, Module* mod, map<int, int> &origNodeID);
void print_twoLevel_Cluster(Network network, string networkName, string outDir);

void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID);


int main(int argc, char *argv[]) {

	if( argc < 10){ 
		cout << "Call: ./ompRelaxmap <seed> <network.net> <# threads> <# attempts> <threshold> <vThresh> <maxIter> <outDir> <prior/normal>  [selflinks]" << endl;
		exit(-1);
	}
	
	string outDir = string(argv[8]);
	int maxIter = atoi(argv[7]);	// Set the maximum number of iteration..
	int Ntrials = atoi(argv[4]);  // Set number of partition attempts
	int numThreads = atoi(argv[3]);	// Set number of threads...
	string line;
	string buf;
  
	MTRand *R = new MTRand(stou(argv[1]));

	string infile = string(argv[2]);
	string networkFile = string(argv[2]);
	string networkName(networkFile.begin() + networkFile.find_last_of("/"), networkFile.begin() + networkFile.find_last_of("."));
	string networkType(infile.begin() + infile.find_last_of("."), infile.end());

	double threshold = atof(argv[5]);
	double vThresh = atof(argv[6]);		// vertex-threshold: threshold for each vertex-movement.

	cout << "Threshold = " << threshold << ", Vertex-Threshold = " << vThresh << endl;
	vThresh *= -1;	// change the threshold value for negative.

	string priorFlag = string(argv[9]);

	bool prior = false;
	if (priorFlag == "prior")
		prior = true;

	bool includeSelfLinks = false;
	if(argc == 11) {
		string selfLinks(argv[10]);
		if(selfLinks == "selflinks")
			includeSelfLinks = true;
	}

	Network origNetwork;	// default constructor is called.

	origNetwork.R = R;
	
	// time values for measuring elapsed times for each step...
	struct timeval allStart, allEnd;
	struct timeval noIOstart, noIOend;
	struct timeval start, end;

	gettimeofday(&allStart, NULL);
	gettimeofday(&start, NULL);

	if(networkType == ".net"){
		load_pajek_format_network(networkFile, origNetwork);    
	}
	else{
		load_linkList_format_network(networkFile, origNetwork); 
	}

	gettimeofday(&end, NULL);

	cout << "Time for reading input data : " << elapsedTimeInSec(start, end) << " (sec)" << endl;

	gettimeofday(&noIOstart, NULL);

	int nNode = origNetwork.NNode();
	double totNodeWeights = origNetwork.TotNodeWeights();
  
	cout << "total Node weights = " << totNodeWeights << endl;

	gettimeofday(&start, NULL);

	for (int i = 0; i < nNode; i++) {
		origNetwork.nodes[i].setTeleportWeight(origNetwork.nodes[i].NodeWeight()/totNodeWeights);
	}

	int NselfLinks = 0;
	for(map<pair<int,int>,double>::iterator it = origNetwork.Edges.begin(); it != origNetwork.Edges.end(); it++){
        
		int from = it->first.first;
		int to = it->first.second;
		double weight = it->second;
		if(weight > 0.0){
			if(from == to){
				NselfLinks++;
                //if(includeSelfLinks)
                //	origNetwork.nodes[from]->selfLink += weight;
			}   
			else{
				origNetwork.nodes[from].outLinks.push_back(make_pair(to,weight));
				// we will going to update inLinks, after we got final flow of the network.
				//origNetwork.nodes[to].inLinks.push_back(make_pair(from,weight));
			}   
		}   
	}
  

    if(includeSelfLinks)
		//cout << "including " <<  NselfLinks << " self link(s)." << endl;  
		cout << "current version always excludes self links.\nignoring " << NselfLinks << " self link(s)." << endl;
    else
        cout << "ignoring " <<  NselfLinks << " self link(s)." << endl;

	//Swap vector to free memory
	map<pair<int,int>,double>().swap(origNetwork.Edges);
        
	cout << "DONE: Parsing the given network  ..." << endl;

	gettimeofday(&end, NULL);
	cout << "Time for parsing the given network : " << elapsedTimeInSec(start, end) << " (sec)" << endl;


	gettimeofday(&start, NULL);
	// Initialization.
	origNetwork.initiate(numThreads);

	// Now update inLinks..
	for (int i = 0; i < nNode; i++) {
		int nOutLinks = origNetwork.nodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
			origNetwork.nodes[origNetwork.nodes[i].outLinks[j].first].inLinks.push_back(make_pair(i,origNetwork.nodes[i].outLinks[j].second));
	}

	gettimeofday(&end, NULL);

	cout << "DONE: Initiate() ... in " << elapsedTimeInSec(start, end) << " (sec)" << endl;
	cout << "Initial Code Length: " << origNetwork.CodeLength()/log(2.0) << " in " << origNetwork.NModule() << " modules.\n";

	
	// copy size of each node for print in order.
	vector<double> nodeSize(nNode);
	for (int i = 0; i < nNode; i++)
		nodeSize[i] = origNetwork.nodes[i].Size();

	
	cout << "Now partition the network starts...\n";
	gettimeofday(&start, NULL);


	
	bool fineTune = true;
	bool fast = false;		// This will be true only for sub-module partitioning...

	int step = 1;

	// Initial SuperStep running ...
	double oldCodeLength = origNetwork.CodeLength();

	stochastic_greedy_partition(origNetwork, numThreads, threshold, vThresh, maxIter, prior, fineTune, fast);
	cout << "SuperStep [" << step << "] - codeLength = " << origNetwork.CodeLength()/log(2.0) << " in " << origNetwork.NModule() << " modules." << endl;

	bool nextIter = true;

	if ((oldCodeLength - origNetwork.CodeLength())/log(2.0) < threshold)
		nextIter = false;

	struct timeval subStart, subEnd;

	while (nextIter) {	
		oldCodeLength = origNetwork.CodeLength();

		stochastic_greedy_partition(origNetwork, numThreads, threshold, vThresh, maxIter, prior, fineTune, fast);

		step++;
		cout << "SuperStep [" << step << "] - codeLength = " << origNetwork.CodeLength()/log(2.0) << " in " << origNetwork.NModule() << " modules." << endl;

		if ((oldCodeLength - origNetwork.CodeLength())/log(2.0) < threshold)
			nextIter = false;

		fineTune = !fineTune;	// fine-tune and coarse-tune will be done alternatively.

		if (nextIter && !fineTune) {
			// Next iteration will be Coarse Tune.
			gettimeofday(&subStart, NULL);
			generate_sub_modules(origNetwork, numThreads, threshold, maxIter);
			gettimeofday(&subEnd, NULL);
			cout << "Time for finding sub-modules: " << elapsedTimeInSec(subStart, subEnd) << " (sec)" << endl;
		}
	}
	
	gettimeofday(&end, NULL);
	cout << "Time for partitioning : " << elapsedTimeInSec(start, end) << " (sec)" << endl;

	cout << "DONE: Code Length = " << origNetwork.CodeLength()/log(2.0) << " in ";
	cout << origNetwork.NModule() << " modules, with " << nNode << " nodes.\n" << endl;

	gettimeofday(&noIOend, NULL);
	gettimeofday(&allEnd, NULL);
	cout << "Overall Elapsed Time for Module Detection (w/o file IO): " << elapsedTimeInSec(noIOstart, noIOend) << " (sec)" << endl;
	cout << "Overall Elapsed Time for Module Detection (w/ file Reading): " << elapsedTimeInSec(allStart, allEnd) << " (sec)" << endl;

	
	cout << "\nComputed Code Length = " << origNetwork.calculateCodeLength()/log(2.0) << endl;

	//Print two-level clustering result in .tree file
	print_twoLevel_Cluster(origNetwork, networkName, outDir);


	// Print partition in Pajek's .clu format
	ofstream outFile;
	ostringstream oss;

	oss.str("");
	oss << outDir << "/" << networkName << ".clu";
	outFile.open(oss.str().c_str());
	outFile << "*Vertices " << nNode << "\x0D\x0A";
	for(int i=0;i<nNode;i++)
		outFile << origNetwork.nodes[i].ModIdx() + 1 << "\x0D\x0A";
	outFile.close();


}


/*
 * Procedure will be following:
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 *	2) repeated 1) procedure, each time in a new random sequential order, until no move generates a gain of the map EQ.
 *
 *	The 1) procedure is implemented in Network::move() function.
 */
void stochastic_greedy_partition(Network &network, int numTh, double threshold, double vThresh, int maxIter, bool prior, bool fineTune, bool fast) {

	double oldCodeLength = network.CodeLength();
	int iter = 0;
	bool stop = false;

	struct timeval outer_T1, outer_T2;
	struct timeval inner_T1, inner_T2;
	struct timeval seq_T1, seq_T2;
	struct timeval convert_T1, convert_T2;

	double tSequential = 0.0;

	gettimeofday(&outer_T1, NULL);

	int nActiveUnits = (fineTune) ? network.NNode() : network.superNodes.size();
	cout << nActiveUnits << ", ";


	// set initial active nodes list ...
	vector<char>(nActiveUnits).swap(network.isActives);
	vector<int>(nActiveUnits).swap(network.activeNodes);
	for (int i = 0; i < nActiveUnits; i++) {
		network.activeNodes[i] = i;
		network.isActives[i] = 0;	// initially set inactive nodes.
	}


	int numMoved = 0;

	while (!stop && iter < maxIter) {
		gettimeofday(&inner_T1, NULL);

		oldCodeLength = network.CodeLength();

		if (fineTune) {
			if (numTh == 1) {
				if (prior)
					numMoved = network.prioritize_move(vThresh);
				else
					numMoved = network.move();
			}
			else { 
				if (prior)
					numMoved = network.prioritize_parallelMove(numTh, tSequential, vThresh);
				else
					numMoved = network.parallelMove(numTh, tSequential);
			}
		} 
		else {
			if (numTh == 1) {
				if (prior)
					numMoved = network.prioritize_moveSPnodes(vThresh);
				else 
					numMoved = network.moveSuperNodes();		// If at least one node is moved, return true. Otherwise, return false.
			}
			else {
				if (prior) 
					numMoved = network.prioritize_parallelMoveSPnodes(numTh, tSequential, vThresh);
				else
					numMoved = network.parallelMoveSuperNodes(numTh, tSequential);
			}
		}
		iter++;

		if (oldCodeLength - network.CodeLength() >= 0 && (oldCodeLength - network.CodeLength())/log(2.0) < threshold)
			stop = true;	//moved = false;

		gettimeofday(&inner_T2, NULL);

		// Print code length per iteration for DEBUG purpose.
		cout << "Iteration " << iter << ": code length =\t" << network.CodeLength()/log(2.0) << "\t, ";
		cout << "elapsed time:\t" << elapsedTimeInSec(inner_T1, inner_T2) << "\t(sec), ";
		cout << "accumulated time:\t" << elapsedTimeInSec(outer_T1, inner_T2) << "\t(sec)\t";
		cout << "sumExitPr = " << network.SumAllExitPr() << "\t";
		cout << "numMoved:\t" << numMoved << endl;
	}

	gettimeofday(&seq_T1, NULL);

	int outerLoop = 1;

	network.updateMembersInModule();
	
	gettimeofday(&seq_T2, NULL);
	tSequential += elapsedTimeInSec(seq_T1, seq_T2);

	if (fast)
		return;

	double tConvert = 0.0;

	do {
		oldCodeLength = network.CodeLength();

		stop = false;

		gettimeofday(&convert_T1, NULL);
		network.convertModulesToSuperNodes(numTh);
		gettimeofday(&convert_T2, NULL);

		tConvert += elapsedTimeInSec(convert_T1, convert_T2);

		nActiveUnits = network.superNodes.size();
		cout << nActiveUnits << ", ";


		// set initial active nodes list ...
		vector<char>(nActiveUnits).swap(network.isActives);
		vector<int>(nActiveUnits).swap(network.activeNodes);
		for (int i = 0; i < nActiveUnits; i++) {
			network.activeNodes[i] = i;	// initially all active units are active.
			network.isActives[i] = 0;	// initially set inactive nodes.
		}


		int spIter = 0;
		while (!stop && spIter < maxIter) {
			gettimeofday(&inner_T1, NULL);

			double innerOldCodeLength = network.CodeLength();

			if (numTh == 1) {
				if (prior)
					numMoved = network.prioritize_moveSPnodes(vThresh);
				else 
					numMoved = network.moveSuperNodes();
			}
			else { 
				if (prior)
					numMoved = network.prioritize_parallelMoveSPnodes(numTh, tSequential, vThresh);
				else
					numMoved = network.parallelMoveSuperNodes(numTh, tSequential);
			}

			spIter++;

			if (innerOldCodeLength - network.CodeLength() >= 0.0 && (innerOldCodeLength - network.CodeLength())/log(2.0) < threshold)
				stop = true;	//moved = false;

			gettimeofday(&inner_T2, NULL);

			// Print code length per spIter for DEBUG purpose.
			cout << "SuperIteration " << outerLoop << "-" << spIter << ": code length =\t" << network.CodeLength()/log(2.0) << "\t, ";
			cout << "elapsed time:\t" << elapsedTimeInSec(inner_T1, inner_T2) << "\t(sec), ";
			cout << "accumulated time:\t" << elapsedTimeInSec(outer_T1, inner_T2) << "\t(sec)\t";
			cout << "sumExitPr = " << network.SumAllExitPr() << "\t";
			cout << "numMoved:\t" << numMoved << endl;
		}

		gettimeofday(&seq_T1, NULL);

		network.updateMembersInModule();

		outerLoop++;
		
		gettimeofday(&seq_T2, NULL);
		tSequential += elapsedTimeInSec(seq_T1, seq_T2);

	} while ((oldCodeLength - network.CodeLength())/log(2.0) > threshold);
	
	gettimeofday(&outer_T2, NULL);
	
	cout << "Sequential running time for partition: " << tSequential << " (sec)" << endl;
	cout << "Time for converting Module to SuperNode: " << tConvert << " (sec)" << endl;
	cout << "Overall time for partition: " << elapsedTimeInSec(outer_T1, outer_T2) << "\t(sec)" << endl;
}


/**
 *	This function will be called for partitioning sub-Module of each module of the original graph.
 *	Thus, we would like to reduce printing from this function for providing high-level log.
 */
void partition_module_network(Network &network, int numTh, double threshold, int maxIter, bool fast) {

	double oldCodeLength = network.CodeLength();
	int iter = 0;
	bool stop = false;

	double tSequential;

	int numMoved = 0;

	while (!stop && iter < maxIter) {
		oldCodeLength = network.CodeLength();

		if (numTh == 1) {
			numMoved = network.move();
		}
		else { 
			numMoved = network.parallelMove(numTh, tSequential);
		}

		iter++;

		if ((oldCodeLength - network.CodeLength())/log(2.0) < threshold)
			stop = true;
	}

	int outerLoop = 1;

	network.updateMembersInModule();
	
	if (fast)
		return;

	do {
		oldCodeLength = network.CodeLength();

		stop = false;
		network.convertModulesToSuperNodes(numTh);

		int spIter = 0;
		while (!stop && spIter < maxIter) {
			double innerOldCodeLength = network.CodeLength();

			if (numTh == 1) {
				numMoved = network.moveSuperNodes();
			}
			else {
				numMoved = network.parallelMoveSuperNodes(numTh, tSequential);
			}

			spIter++;

			if ((innerOldCodeLength - network.CodeLength())/log(2.0) < threshold)
				stop = true;
		}

		network.updateMembersInModule();

		outerLoop++;

	} while ((oldCodeLength - network.CodeLength())/log(2.0) > threshold);

}



void generate_sub_modules(Network &network, int numTh, double threshold, int maxIter) {
	int numNodes = network.NNode();

	struct timeval t1, t2;

	gettimeofday(&t1, NULL);

	vector<SubModule*>().swap(network.subModules);	// swap network.subModules vector with the empty vector.
	network.subModules.reserve(numNodes);

	vector<int>(numNodes).swap(network.ndToSubMod);

	omp_set_num_threads(numTh);

	vector<vector<SubModule> > tmpSubModList(numTh);	// generate 'numTh' of vector<SubModule>.
	for(int i = 0; i < numTh; i++) {
		tmpSubModList[i].reserve(numNodes);
	}

	gettimeofday(&t2, NULL);
	cout << "initialization time for generate_sub_modules(): " << elapsedTimeInSec(t1, t2) << endl;


	struct timeval tPar1, tPar2;

	gettimeofday(&tPar1, NULL);

	MTRand *Rand = new MTRand();

	int numSmallMods = network.smActiveMods.size();
	cout << "number of small modules: " << numSmallMods << endl;

#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < numSmallMods; i++) {
		int myID = omp_get_thread_num();	// get my thread ID.
		
		Module* mod = &(network.modules[network.smActiveMods[i]]);

		// check whether the current module has more than one node or not.
		if (mod->NumMembers() > 1) {
			int modIdx = mod->Index();

			map<int, int> origNodeID;	//map from newNodeID to original Node ID. a.k.a. <newNodeID, origNodeID>

			Network newNetwork;
			generate_network_from_module(newNetwork, mod, origNodeID);
			newNetwork.R = Rand;

			partition_module_network(newNetwork, 1, threshold, maxIter, true);		// fast = true..

			int nActiveMods = newNetwork.smActiveMods.size();
			// Adding sub-modules from a new network of the corresponding module to the list of subModules...
			for (int j = 0; j < nActiveMods; j++) {
				SubModule subMod(newNetwork.modules[newNetwork.smActiveMods[j]], origNodeID, modIdx);
				tmpSubModList[myID].push_back(subMod);
			}
		}
		else {
			// This is the special case that the module has ONLY ONE member.
			SubModule subMod(*mod);
			tmpSubModList[myID].push_back(subMod);
		}
	} // End of for
//}	// End of parallel.
	gettimeofday(&tPar2, NULL);

	cout << "Time for parallel for loop for SMALL-MODULES:\t" << elapsedTimeInSec(tPar1, tPar2) << " (sec)" << endl;


	///////////////////////////////
	// Larger-Modules

	gettimeofday(&tPar1, NULL);


	int numLargeMods = network.lgActiveMods.size();
	cout << "number of large modules: " << numLargeMods << endl;

	for (int i = 0; i < numLargeMods; i++) {
		Module* mod = &(network.modules[network.lgActiveMods[i]]);

		// NO-NEED to check the size of the current module.
		int modIdx = mod->Index();

		map<int, int> origNodeID;	//map from newNodeID to original Node ID. a.k.a. <newNodeID, origNodeID>
			
		Network newNetwork;
		generate_network_from_module(newNetwork, mod, origNodeID, numTh);
		newNetwork.R = Rand;

		partition_module_network(newNetwork, numTh, threshold, maxIter, true);		// fast = true..

		// Adding sub-modules from a new network of the corresponding module to the list of subModules...
		int nActiveMods = newNetwork.smActiveMods.size();
		#pragma omp parallel for //schedule(dynamic)
		for (int j = 0; j < nActiveMods; j++) {
			SubModule subMod(newNetwork.modules[newNetwork.smActiveMods[j]], origNodeID, modIdx);
			tmpSubModList[omp_get_thread_num()].push_back(subMod);
		}

		nActiveMods = newNetwork.lgActiveMods.size();
		#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < nActiveMods; j++) {
			SubModule subMod(newNetwork.modules[newNetwork.lgActiveMods[j]], origNodeID, modIdx);
			tmpSubModList[omp_get_thread_num()].push_back(subMod);
		}
	}	// End of parallel for.

	gettimeofday(&tPar2, NULL);

	cout << "Time for parallel for loop for LARGE-MODULES:\t" << elapsedTimeInSec(tPar1, tPar2) << " (sec)" << endl;


	gettimeofday(&t1, NULL);

	int numSubMods = 0;
	for (int i = 0; i < numTh; i++) {
		for (vector<SubModule>::iterator it = tmpSubModList[i].begin(); it != tmpSubModList[i].end(); it++) {
			network.subModules.push_back(&(*it));
			for (vector<int>::iterator ndIt = it->members.begin(); ndIt != it->members.end(); ndIt++)
				network.ndToSubMod[*ndIt] = numSubMods;

			numSubMods++;
		}
	}

	gettimeofday(&t2, NULL);
	cout << "sequential subModules push_back() time:\t" << elapsedTimeInSec(t1, t2) << " (sec)" << endl;


	gettimeofday(&t1, NULL);
	network.generateSuperNodesFromSubModules(numTh);
	gettimeofday(&t2, NULL);

	cout << "generateSuperNodesFromSubModules() time:\t" << elapsedTimeInSec(t1, t2) << " (sec)" << endl;

}




void generate_network_from_module(Network &newNetwork, Module* mod, map<int, int> &origNodeID) {
	int numMembers = mod->NumMembers();
	newNetwork.modules = vector<Module>(numMembers);

	map<int, int> newNodeID;	// key = origNodeID --> value =  newNodeID.

	int newIdx = 0;
	for (vector<Node *>::iterator it = mod->members.begin(); it != mod->members.end(); it++) {
		newNodeID[(*it)->ID()] = newIdx;
		origNodeID[newIdx] = (*it)->ID();

		Node nd(newIdx, (*it)->Size());
		nd.setNodeWeight((*it)->NodeWeight());
		nd.setTeleportWeight((*it)->TeleportWeight());
		nd.setDanglingSize((*it)->DanglingSize());
		nd.setIsDangling((*it)->IsDangling());

		newNetwork.nodes.push_back(nd);
		newIdx++;	// newIdx is equal to the number of nodes in this module (or network.)
	}

	// add outLinks within the newNetwork.
	for (int i = 0; i < numMembers; i++) {
		Node* it = mod->members[i];
		int nid = newNodeID[it->ID()];
		Node* nd_ptr = &(newNetwork.nodes[nid]);

		for (link_iterator link_it = it->outLinks.begin(); link_it != it->outLinks.end(); link_it++) {
			// check whether the edge within the module or not.
			map<int, int>::iterator m_it = newNodeID.find(link_it->first);
			if (m_it != newNodeID.end()) {
				nd_ptr->outLinks.push_back(make_pair(m_it->second, link_it->second));
			}
		}
	}

	// add inLinks within the newNetwork based on the generated outLinks above.
	for (vector<Node>::iterator it = newNetwork.nodes.begin(); it != newNetwork.nodes.end(); it++) {
		for (link_iterator l_it = it->outLinks.begin(); l_it != it->outLinks.end(); l_it++) {
			newNetwork.nodes[l_it->first].inLinks.push_back(make_pair(it->ID(), l_it->second));
		}
	}

	double sum_size_log_size = 0.0;

	for (vector<Node>::iterator it = newNetwork.nodes.begin(); it != newNetwork.nodes.end(); it++)
		sum_size_log_size += pLogP(it->Size());

	newNetwork.setAllLogAll(sum_size_log_size);

	for (int i = 0; i < newIdx; i++) {
		newNetwork.modules[i] = Module(i, &newNetwork.nodes[i]);
		newNetwork.nodes[i].setModIdx(i);
	}


	newNetwork.setNModule(newIdx);
	newNetwork.setNNode(newIdx);

	newNetwork.calibrate(1);	// This function is run in sequential.

}








void generate_network_from_module(Network &newNetwork, Module* mod, map<int, int> &origNodeID, int numTh) {
	int numMembers = mod->NumMembers();
	newNetwork.modules = vector<Module>(numMembers);

	map<int, int> newNodeID;	// key = origNodeID --> value =  newNodeID.

	int newIdx = 0;
	for (vector<Node *>::iterator it = mod->members.begin(); it != mod->members.end(); it++) {
		newNodeID[(*it)->ID()] = newIdx;
		origNodeID[newIdx] = (*it)->ID();

		Node nd(newIdx, (*it)->Size());
		nd.setNodeWeight((*it)->NodeWeight());
		nd.setTeleportWeight((*it)->TeleportWeight());
		nd.setDanglingSize((*it)->DanglingSize());
		nd.setIsDangling((*it)->IsDangling());

		newNetwork.nodes.push_back(nd);
		newIdx++;	// newIdx is equal to the number of nodes in this module (or network.)
	}

	omp_set_num_threads(numTh);

	double sum_size_log_size = 0.0;

#pragma omp parallel
{
	// add outLinks within the newNetwork.
	#pragma omp for //nowait
	for (int i = 0; i < numMembers; i++) {
		Node* it = mod->members[i];
		int nid = newNodeID[it->ID()];
		Node* nd_ptr = &(newNetwork.nodes[nid]);

		for (link_iterator link_it = it->outLinks.begin(); link_it != it->outLinks.end(); link_it++) {
			// check whether the edge within the module or not.
			map<int, int>::iterator m_it = newNodeID.find(link_it->first);
			if (m_it != newNodeID.end()) {
				nd_ptr->outLinks.push_back(make_pair(m_it->second, link_it->second));
			}
		}
	}

	#pragma omp master
	// add inLinks within the newNetwork based on the generated outLinks above.
	for (vector<Node>::iterator it = newNetwork.nodes.begin(); it != newNetwork.nodes.end(); it++) {
		for (link_iterator l_it = it->outLinks.begin(); l_it != it->outLinks.end(); l_it++) {
			newNetwork.nodes[l_it->first].inLinks.push_back(make_pair(it->ID(), l_it->second));
		}
	}

	#pragma omp for reduction(+:sum_size_log_size)
	for (int i = 0; i < numMembers; i++) {
		sum_size_log_size += pLogP(newNetwork.nodes[i].Size());
	}

	#pragma omp master
	newNetwork.setAllLogAll(sum_size_log_size);

	#pragma omp for
	for (int i = 0; i < newIdx; i++) {
		newNetwork.modules[i] = Module(i, &newNetwork.nodes[i]);
		newNetwork.nodes[i].setModIdx(i);
	}
}

	newNetwork.setNModule(newIdx);
	newNetwork.setNNode(newIdx);

	newNetwork.calibrate(numTh);
}






void print_twoLevel_Cluster(Network network, string networkName, string outDir) {
	ofstream outFile;
	ostringstream oss;
	oss << outDir << "/" << networkName << ".tree";
	outFile.open(oss.str().c_str());

	outFile << "# Code length " << network.CodeLength()/log(2.0) << " in " << network.NModule() << " modules." << endl;

	int nModules = network.modules.size();
	int modIdx = 0;
	for (int i = 0; i < nModules; i++) {

		int nMembers = network.modules[i].NumMembers();
		if (nMembers > 0)
			modIdx++;

		for (int j = 0; j < nMembers; j++) {
			outFile << modIdx << ":" << j+1 << " " << network.modules[i].members[j]->Size() << " \"" << network.modules[i].members[j]->Name() << "\"" << endl;
		}
	}

	outFile.close();
}



