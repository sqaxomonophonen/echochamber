#ifndef SYS_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <errno.h>

static inline void* calloc_or_die(size_t nmemb, size_t size)
{
	void* ptr = calloc(nmemb, size);
	if (ptr == NULL) {
		fprintf(stderr, "calloc: %s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	return ptr;
}

static inline void wrong(char* msg)
{
	fprintf(stderr, msg);
	exit(EXIT_FAILURE);
}

static inline ssize_t writen(int fd, const void* vptr, size_t n)
{
	size_t nleft = n;
	const char* ptr = vptr;

	while (nleft > 0) {
		ssize_t nwritten = write(fd, ptr, nleft);
		if (nwritten < 0) {
			if (errno == EINTR) {
				nwritten = 0;
			} else {
				return -1;
			}
		}
		nleft -= nwritten;
		ptr += nwritten;
	}

	return n;
}

static inline ssize_t writen_or_die(int fd, const void* vptr, size_t n)
{
	ssize_t e = writen(fd, vptr, n);
	if (e < 0) {
		fprintf(stderr, "failed to write %zd bytes: %s\n", e, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return e;
}

#define SYS_H
#endif
