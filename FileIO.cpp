/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Mar. 2014
 *	Copyright (C) since 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include "Node.h"
#include "Module.h"
#include "FileIO.h"

using namespace std;

template <class T>
inline string to_string(const T& t) {
	stringstream ss;
	ss << t;
	return ss.str();
}

/*
 * Modified from the Infomap implementation.
 * Reading the given network data in Pajek format.
 */
void load_pajek_format_network(string fName, Network &network) {

  string line;
  string buf;
  string nameBuf;
  
  /* Read network in Pajek format with nodes ordered 1, 2, 3, ..., N,            */
  /* each directed link occurring only once, and link weights > 0.               */
  /* For more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/.   */
  /* Node weights are optional and sets the relative proportion to which         */
  /* each node receives teleporting random walkers. Default value is 1.          */
  /* Example network with three nodes and four directed and weighted links:      */
  /* *Vertices 3                                                                 */
  /* 1 "Name of first node" 1.0                                                  */
  /* 2 "Name of second node" 2.0                                                 */
  /* 3 "Name of third node" 1.0                                                  */
  /* *Arcs 4                                                                     */
  /* 1 2 1.0                                                                     */
  /* 1 3 1.7                                                                     */
  /* 2 3 2.0                                                                     */
  /* 3 2 1.2                                                                     */
  
  cout << "Reading network " << fName << " file\n" << flush;
  ifstream net(fName.c_str());
  network.setNNode(0);

  istringstream ss;
  while(network.NNode() == 0){ 
    if(!getline(net,line)){
      cout << "the network file is not in Pajek format...exiting" << endl;
      exit(-1);
    }
    else{
      ss.clear();
      ss.str(line);
      ss >> buf;
      if(buf == "*Vertices" || buf == "*vertices" || buf == "*VERTICES"){
        ss >> buf;
        network.setNNode(atoi(buf.c_str()));
      }
      else{
        cout << "the network file is not in Pajek format...exiting" << endl;
        exit(-1);
      }
    }
  }
  
  int nNode = network.NNode();

  network.modules = vector<Module>(nNode);
  network.nodes = vector<Node>(nNode);
  double totNodeWeights = 0.0;
  
  // Read node names, assuming order 1, 2, 3, ...
  for(int i = 0; i < nNode; i++){
    getline(net,line);

	// Assume Vertices format in "id name weight"
	ss.clear();
	ss.str(line);
	ss >> buf;		// read first item, which is ID.
	//network.nodes[i].setID(atoi(buf.c_str()));
	network.nodes[i].setID(i);

    int nameStart = line.find_first_of("\"");
    int nameEnd = line.find_last_of("\"");
    if(nameStart < nameEnd){
      network.nodes[i].setName(string(line.begin() + nameStart + 1,line.begin() + nameEnd));
      line = string(line.begin() + nameEnd + 1, line.end());
      ss.clear();
      ss.str(line);
    }
    else{
      ss.clear();
      ss.str(line);		// After reading ID...
      //ss >> buf; 
      //ss >> network.nodeNames[i];
	  ss >> buf;
	  ss >> buf;
	  network.nodes[i].setName(buf);
    }
    
    buf = "1";
    ss >> buf;
    double nodeWeight = atof(buf.c_str());
    if(nodeWeight <= 0.0)
      nodeWeight = 1.0;
    network.nodes[i].setNodeWeight(nodeWeight);
    totNodeWeights += nodeWeight; 
  }

  network.setTotNodeWeights(totNodeWeights);
  
  // Read the number of links in the network
  getline(net,line);
  ss.clear();
  ss.str(line);
  ss >> buf;
  
  if(buf != "*Edges" && buf != "*edges" && buf != "*Arcs" && buf != "*arcs"){
    cout << endl << "Number of nodes not matching, exiting" << endl;
    exit(-1);
  }

  double newLinkWeight;
  int NdoubleLinks = 0;
  
  // Read links in format "from to weight", for example "1 3 0.7"
  while(getline(net,line)){
    ss.clear();
    ss.str(line);
    ss >> buf;
    int linkEnd1 = atoi(buf.c_str());
    ss >> buf;
    int linkEnd2 = atoi(buf.c_str());
    buf.clear();
    ss >> buf;
    double linkWeight;
    if(buf.empty()) // If no information 
      linkWeight = 1.0;
    else
      linkWeight = atof(buf.c_str());
    linkEnd1--; // Nodes start at 1, but C++ arrays at 0.
    linkEnd2--;	// As of now, these values are indices of corresponding nodes in network.nodes vector.
     
    newLinkWeight = network.Edges[make_pair(linkEnd1,linkEnd2)] += linkWeight;
    if(newLinkWeight > linkWeight)
      NdoubleLinks++;
  }

  net.close();

  network.setNEdge(network.Edges.size());
  
  cout << "done! (found " << network.NNode() << " nodes and " << network.NEdge() << " edges.)";
  if(NdoubleLinks > 0)
    cout << ", aggregated " << NdoubleLinks << " link(s) defined more than once" << endl;
  else
	cout << endl;
}

/*
 *	Read network in the format "FromNodeID	ToNodeID	[Weight]"
 *	not assuming a complete list of nodes, 1..maxnode.
 */
void load_linkList_format_network(string fName, Network &network) {

	string line;
	string buf;
	istringstream ss;

	cout << "Reading network " << fName << " file\n" << flush;
	ifstream net(fName.c_str());
	network.setNNode(0);

	int nDoubleLinks = 0;
	double newLinkWeight = 0.0;

	set<long> Nodes;


	// generate temporary Edge map with long value for large id graphs..
	map<pair<long,long>, double> tempEdges;

	// Read links in format "from to [weight]"
	while (getline(net,line)) {
		if (line[0] != '#') {
			ss.clear();
			ss.str(line);

			ss >> buf;
			long linkEnd1 = atol(buf.c_str());
			
			ss >> buf;
			long linkEnd2 = atol(buf.c_str());
			
			buf.clear();
			ss >> buf;
			double linkWeight;
			if (buf.empty()) // If no linkWeight information..
				linkWeight = 1.0;
			else
				linkWeight = atof(buf.c_str());

			// To keep track of all nodes for renaming 0..N-1
			Nodes.insert(linkEnd1);
			Nodes.insert(linkEnd2);

			newLinkWeight = tempEdges[make_pair(linkEnd1, linkEnd2)] += linkWeight;
			if (newLinkWeight > linkWeight)
				nDoubleLinks++;
		}
	}

	net.close();

	network.setNNode(Nodes.size());
	network.setTotNodeWeights(Nodes.size());	// Since no node weight info given, assume uniform weight.
	network.setNEdge(network.Edges.size());
	
	network.modules = vector<Module>(Nodes.size());
	network.nodes = vector<Node>(Nodes.size());

	// Renaming all nodes if necessary...
	vector<string> nodeNames = vector<string>(Nodes.size());

	long nodeCounter = 0;
	map<long, long> renumber;
	bool renum = false;

	for (set<long>::iterator it = Nodes.begin(); it != Nodes.end(); it++) {
		nodeNames[nodeCounter] = to_string((*it));
		renumber.insert(make_pair((*it), nodeCounter));

		network.nodes[nodeCounter].setID( (int) nodeCounter);
		network.nodes[nodeCounter].setName(nodeNames[nodeCounter]);

		if (nodeCounter != (*it))
			renum = true;
		nodeCounter++;
	}

	for (map<pair<long, long>, double>::iterator it = tempEdges.begin(); it != tempEdges.end(); it++)
		network.Edges.insert(make_pair(make_pair((int)renumber.find(it->first.first)->second, (int)renumber.find(it->first.second)->second), it->second));


	cout << "done! (found " << network.NNode() << " nodes and " << network.NEdge() << " edges.)";
	if(nDoubleLinks > 0)
		cout << ", aggregated " << nDoubleLinks << " link(s) defined more than once" << endl;
	else
		cout << endl;

	cout << "done! (found " << network.nodes.size() << " nodes and " << network.Edges.size() << " edges.)" << endl;
}

void write_cluster_result(string outFile, Network &finalNetwork) {

}
