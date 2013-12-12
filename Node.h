/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Dec. 2013
 *	Copyright (C) 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>
#include <map>
class Node;
class SuperNode;
#include "Module.h"

class Module;
class Network;
struct SubModule;

using namespace std;

typedef vector<pair<int, double> >::iterator	link_iterator;	// outLinks, inLinks iterator definition.

class Node {
  protected:
	int id;					// node ID number
	string name;			// node name.	// initialized with empty string via string().
	double size;			// p_alpha.  
	double exitPr;			// exitPr for SuperNode or Node from original module.
	double nodeWeight;		// node weight for non-uniform teleport weight.
	int modIdx;				// module index number.
	bool isModule;			// true: if it is a module in the previous level
							// false: otherwise
	bool isSuper;			// true: if it is superNode (it is a group of nodes in the current level)
							// false: otherwise
	Module *orig_module;	// pointer to the original module of this node, NULL if isModule == false.

	double teleportWeight;
	
	double danglingSize;	// danglingSize will be equivalent to sumDangliing of corresponding module.
	bool isDangling;		// true: if it is a dangling node, false: otherwise.

  public:
	vector<pair<int, double> > outLinks;	// outlinks: <nodeIDto, weight>
	vector<pair<int, double> > inLinks;		// inlinks: <nodeIDfrom, weight>	


	// Constructors and functions
	Node();
	Node(int id);
	Node(int id, double size);
	Node(int id, double size, int module, bool isModule, Module *original_module);


	// Getter -- Setter
	int ID() { return id; }
	void setID(int id) { this->id = id; }

	string Name() { return name; }
	void setName(string str) { name = str; }

	double Size() { return size; }
	void setSize(double sz) { size = sz; }
	void addSize(double val) { size += val; }

	double ExitPr() { return exitPr; }
	void setExitPr(double exit) { exitPr = exit; }

	double NodeWeight() { return nodeWeight; }
	void setNodeWeight(double weight) { nodeWeight = weight; }

	int ModIdx() { return modIdx; }
	void setModIdx(int mod) { modIdx = mod; }

	bool IsModule() { return isModule; }
	void setIsModule(bool module) { isModule = module; }

	bool IsSuper() { return isSuper; }
	void setIsSuper(bool super) { isSuper = super; }

	Module* OrigModule() { return orig_module; }
	void setOrigModule(Module* origMod) { orig_module = origMod; }

	double TeleportWeight() { return teleportWeight; }
	void setTeleportWeight(double tpweight) { teleportWeight = tpweight; }

	double DanglingSize() { return danglingSize; }
	void setDanglingSize(double dSize) { danglingSize = dSize; }

	bool IsDangling() { return isDangling; }
	void setIsDangling(bool dangle) { isDangling = dangle; }

};


class SuperNode : public Node {
  public: 
	vector<Node *> members;
	
	// Constructor
	SuperNode();
	SuperNode(Module module, int newID);	// will call default constructor of Node class.
	SuperNode(SubModule* subMod, int newID, Network& network);
};


#endif
