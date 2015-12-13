#ifndef CPU_H

#include <unistd.h>

// TODO ifdef linux...
int cpu_count()
{
	// http://stackoverflow.com/questions/4586405/get-number-of-cpus-in-linux-using-c
	// XXX ^^^ ought to use /proc/cpuinfo?
	return sysconf(_SC_NPROCESSORS_ONLN);
}

#define CPU_H
#endif
