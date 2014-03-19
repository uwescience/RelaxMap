/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Mar. 2014
 *	Copyright (C) since 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

#ifndef FILEIO_H
#define FILEIO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include "Node.h"
#include "Module.h"

using namespace std;

// declaration of functions related to file IO.
void load_pajek_format_network(string fName, Network &network);
void load_linkList_format_network(string fName, Network &network);
void write_cluster_result(string outFile, Network &finalNetwork);

#endif
