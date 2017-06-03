#ifndef CPU_H

#ifdef __linux__
#include <unistd.h>
int cpu_count()
{
	// http://stackoverflow.com/questions/4586405/get-number-of-cpus-in-linux-using-c
	// XXX ^^^ ought to use /proc/cpuinfo?
	return sysconf(_SC_NPROCESSORS_ONLN);
}
#elif __APPLE__
#define u_long unsigned long
#define u_int unsigned int
#define u_short unsigned short
#define u_char unsigned char
#include <sys/types.h>
#include <sys/sysctl.h>
//#include <sys/param.h>
int cpu_count()
{
	int mib[2];
	int ncpu = 0;
	size_t len = sizeof ncpu;
	mib[0] = CTL_HW;
	mib[1] = HW_AVAILCPU;
	sysctl(mib, 2, &ncpu, &len, NULL, 0);
	if (ncpu < 1) {
		mib[1] = HW_NCPU;
		sysctl(mib, 2, &ncpu, &len, NULL, 0);
		if (ncpu < 1) ncpu = 1;
	}
	return ncpu;
}
#else
#error "no platform for cpu_count()"
#endif

#define CPU_H
#endif
