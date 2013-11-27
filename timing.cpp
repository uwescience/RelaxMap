#include <iostream>
#include "timing.h"
#include <sys/times.h>
#include <math.h>
#include <unistd.h>
#define DEC 10000

using namespace std;

float gettime(){
  float clockticks;
  struct tms Time;

  times(&Time);

  clockticks = (double) sysconf(_SC_CLK_TCK);

  return floor(DEC*Time.tms_utime/clockticks)/DEC;
}

double elapsedTimeInSec(struct timeval start, struct timeval end) {
	return ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.0e6;
}

