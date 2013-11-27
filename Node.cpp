#include <iostream>
#include "Node.h"

using namespace std;

Node::Node()
:id(0),
 size(0.0),
 exitPr(0.0),
 nodeWeight(1.0),
 modIdx(0),
 isModule(false),
 isSuper(false),
 orig_module(NULL),
 teleportWeight(0.0),
 danglingSize(0.0),
 isDangling(false)
{}


Node::Node(int id)
:id(id),
 size(0.0),
 exitPr(0.0),
 nodeWeight(1.0),
 modIdx(0),
 isModule(false),
 isSuper(false),
 orig_module(NULL),
 teleportWeight(0.0),
 danglingSize(0.0),
 isDangling(false)
{}


Node::Node(int id, double size)
:id(id),
 size(size),
 exitPr(0.0),
 nodeWeight(1.0),
 modIdx(0),
 isModule(false),
 isSuper(false),
 orig_module(NULL),
 teleportWeight(0.0),
 danglingSize(0.0),
 isDangling(false)
{}


Node::Node(int id, double size, int module, bool isModule, Module *original_module)
:id(id),
 size(size),
 exitPr(0.0),
 nodeWeight(1.0),
 modIdx(module),
 isModule(isModule),
 isSuper(false),
 orig_module(original_module),
 teleportWeight(0.0),
 danglingSize(0.0),
 isDangling(false)
{}


// Constructor of SuperNode class which is a derived class of Node.
SuperNode::SuperNode() {
}

SuperNode::SuperNode(Module module, int newID)
{
	// initialize member variables w.r.t. the given Module object 'module.'
	id = newID;
	name = "SuperNode";
	size = module.SumPr();
	exitPr = module.ExitPr();
	nodeWeight = 0.0;	// later adding nodeWeight of each member node.
	modIdx = module.Index();
	isModule = false;
	isSuper = true;
	orig_module = NULL;
	teleportWeight = module.SumTPWeight();
	danglingSize = module.SumDangling();
	isDangling = false;		// NEED TO BE CONFIRMED !!!

	int numMembers = module.members.size();
	if (numMembers > 0) {
		for (int i = 0; i < numMembers; i++) {
			//module.members[i]->setModIdx(newID);
			members.push_back(module.members[i]);
			nodeWeight += module.members[i]->NodeWeight();
		}
	}

	if (numMembers != module.NumMembers())
		cout << "SuperNode[" << id << "] members.size() != numMembers" << endl;

	if (members.size() != module.members.size())
		cout << "SuperNode[" << id << "] members.size() != module.members.size()" << endl;
}

SuperNode::SuperNode(SubModule* subMod, int newID, Network& network)
{
	id = newID;
	name = "SuperNode";
	size = subMod->sumPr;
	exitPr = 0.0;
	nodeWeight = 0.0;
	modIdx = subMod->modIdx;
	isModule = false;
	isSuper = true;
	orig_module = NULL;
	teleportWeight = subMod->sumTPWeight;
	danglingSize = subMod->sumDangling;
	isDangling = false;

	//size = 0.0;
	//teleportWeight = 0.0;
	//danglingSize = 0.0;

	int numMembers = subMod->members.size();	// numMembers > 0 always.
	for (int i = 0; i < numMembers; i++) {
		members.push_back( &(network.nodes[subMod->members[i]]) );
		nodeWeight += members[i]->NodeWeight();
		//size += members[i]->Size();
		//teleportWeight += members[i]->TeleportWeight();
		//danglingSize += members[i]->DanglingSize();
	}

	// outLinks and inLinks vectors will be added later...
}
