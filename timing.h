#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>
#include <ctime>

float gettime();
double elapsedTimeInSec(struct timeval start, struct timeval end);

#endif
