/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Dec. 2013
 *	Copyright (C) 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>
#include <ctime>

float gettime();
double elapsedTimeInSec(struct timeval start, struct timeval end);

#endif
