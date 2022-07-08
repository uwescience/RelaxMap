/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Mar. 2014
 *	Copyright (C) since 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

#ifndef MODULE_H
#define MODULE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include "Node.h"
#include "MersenneTwister.h"


//class Module;
//class Network;

// define inline function of p_log_p.
inline double pLogP(double pr) {
	return (pr > 0.0 ? pr * log(pr) : 0.0);
}


using namespace std;


class Module {
	int index;			// Index number of this module.
	double exitPr;		// exit probability of this module.
	double stayPr;		// stay probability of this module, which is sum of p_alpha + exit probability.
	double sumPr;		// sum of p_alpha.  alpha is in this module i.
	double sumTPWeight;	// sum of teleport weight.   SUM (tau_a).
	double sumDangling;	// sum of dangling nodes weight.
	int numMembers;		// number of members of this Module.

  public:
	vector<Node *> members;	// member nodes in this module.

	// Constructors and member functions
	Module();
	Module(int idx, double exitPr, double sumPr);
	Module(int idx, Node * nd);	// used for Node to Module transformation.

	// Getter -- Setter
	int Index() { return index; }
	void setIndex(int idx) { index = idx; }

	double ExitPr() { return exitPr; }
	void setExitPr(double exitProb) { exitPr = exitProb; }

	double StayPr() { return stayPr; }
	void setStayPr(double stayProb) { stayPr = stayProb; }

	double SumPr() { return sumPr; }
	void setSumPr(double sumProb) { sumPr = sumProb; }
	void addSumPr(double ndSize) { sumPr += ndSize; }
	void minusSumPr(double ndSize) { sumPr -= ndSize; }

	double SumTPWeight() { return sumTPWeight; }
	void setSumTPWeight(double sumTPW) { sumTPWeight = sumTPW; }
	void addSumTPWeight(double addedTPW) { sumTPWeight += addedTPW; }
	void minusSumTPWeight(double minusTPW) { sumTPWeight -= minusTPW; }

	double SumDangling() { return sumDangling; }
	void setSumDangling(double sumTPW) { sumDangling = sumTPW; }
	void addSumDangling(double addedWeight) { sumDangling += addedWeight; }
	void minusSumDangling(double minusWeight) { sumDangling -= minusWeight; }

	int NumMembers() { return numMembers; }
	void setNumMembers(int nMember) { numMembers = nMember; }
	void increaseNumMembers() { numMembers++; }
	void increaseNumMembers(int n) { numMembers += n; }
	void decreaseNumMembers() { numMembers--; }
	void decreaseNumMembers(int n) { numMembers -= n; }

};


struct SubModule {
	int modIdx;				// the original module index.
	int numMembers;			// the number of members in this sub-module.
	double sumPr;			// sum of p_alpha.
	double sumTPWeight;		// sum of the teleport weight. SUM(tau_a).
	double sumDangling;		// sum of the weight of dangling nodes. a.k.a. dangling-weight.

	vector<int> members;	// the original node-IDs of the members in this sub-module. 

	SubModule();
	SubModule(Module& mod, map<int, int>& origNodeID, int modIndex);
	SubModule(Module& mod);
};



class Network {
	int level;				// current level of hierarchical mapping.
	double codeLength;		// current best codeLength value of this network.
	int nNode;				// number of Nodes in this Network.
	int nEdge;				// number of Edges in this Network.
	int nEmptyMod;			// number of empty modules.
	unsigned int nModule;			// number of non-empty modules.
	double totNodeWeights;	// total node weights.

	int nDanglings;			// number of dangling nodes.

	double allNodes_log_allNodes;	// SUM (p * log(p))
	double sumAllExitPr;			// SUM (exitPr_i)  i = 1 .. nModule

  public:
	static constexpr double alpha = 0.15;
	//static constexpr double beta = 1.0 - alpha;
	static constexpr double beta = 1.0 - 0.15;

	MTRand *R;

	vector<Module> modules;
	vector<Node> nodes;
	vector<SuperNode> superNodes;
	vector<int> emptyModules;		// keep the list of the indices of empty modules.
	vector<int> smActiveMods;		// keep the list of the indices of small active modules.
	vector<int> lgActiveMods;		// keep the list of the indices of large active modules.
	map<pair<int, int>, double> Edges;	// pair<int, int> will be src and desc of an edge and double is the corresponding weight.
	
	vector<int> danglings;	// node index of dangling nodes in nodes vector.

	vector<int> ndToSubMod;		// This will store subModule ID for each node.
	vector<SubModule*> subModules;	// This will store SubModules of this network
									// for coarse-tune optimization..

	// two-vectors for maintaining activeNodes for prioritizing.
	vector<char> isActives;		// 0 - inactive, 1 - active nodes. Working as a boolean array.
	vector<int> activeNodes;	// Actual node IDs for active nodes.

	// Constructors and member functions
	Network();
	Network(int numNode, int numEdge, int level);
	Network(int numNode, int numEdge, int level, double codeLen);

	//Getter -- Setter
	int Level() { return level; }
	void setLevel(int l) { level = l; }

	double CodeLength() { return codeLength; }
	void setCodeLength(double cdLength) { codeLength = cdLength; }

	int NNode() { return nNode; }
	void setNNode(int nNd) { nNode = nNd; }

	int NEdge() { return nEdge; }
	void setNEdge(int nEd) { nEdge = nEd; }

	int NEmptyMod() { return nEmptyMod; }
	void setNEmptyMod(int nEmpty) { nEmptyMod = nEmpty; }
	void increaseNEmptyMod() { nEmptyMod++; }

	int NModule() { return nModule; }
	void setNModule(int nMod) { nModule = nMod; }

	double TotNodeWeights() { return totNodeWeights; }
	void setTotNodeWeights(double tWeight) { totNodeWeights = tWeight; }

	int NDanglings() { return nDanglings; }
	void setNDanglings(int nDangle) { nDanglings = nDangle; }

	double AllLogAll() { return allNodes_log_allNodes; }
	void setAllLogAll(double all) { allNodes_log_allNodes = all; }

	double SumAllExitPr() { return sumAllExitPr; }
	void setSumAllExitPr(double sumExitPr) { sumAllExitPr = sumExitPr; }

	// other Member functions ...
	void findDanglingNodes();
	//void initiate();
	void initiate(int numTh);
	//void calculateSteadyState();	// eigenvector();
	void calculateSteadyState(int numTh);	// eigenvector();
	void calibrate(int numTh);
	int move();
	int prioritize_move(double vThresh);
	int parallelMove(int numTh, double & tUpdate);
	int prioritize_parallelMove(int numTh, double & tUpdate, double vThresh);
	int moveSuperNodes();
	int prioritize_moveSPnodes(double vThresh);
	int parallelMoveSuperNodes(int numTh, double & tUpdate);
	int prioritize_parallelMoveSPnodes(int numTh, double & tUpdate, double vThresh);
	void convertModulesToSuperNodes(int numTh);
	//void generateSuperNodesFromSubModules();	// generate SuperNodes from SubModule..
	void generateSuperNodesFromSubModules(int numTh);	// parallel version of the above function.
	
	double calculateCodeLength();
	void updateMembersInModule();
	void updateSPMembersInModule();
	void updateCodeLength(int numTh, bool isSPNode);
	void copyModule(Module * newM, Module * oldM);	// copy from oldM to newM.
	void copyModule(int newM, int oldM);	// copy modules[oldM] to modules[newM].
};

#endif
