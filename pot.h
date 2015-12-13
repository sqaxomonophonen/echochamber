#ifndef POT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct pot {
	int index;
	int argc;
	char** argv;

	int is_swtch;
	int is_opt;
	int is_arg;

	char swtch;
	int swtch_idx;
	char* opt;
	char* arg;

	char* value;
};

static inline void pot_init(struct pot* pot, int argc, char** argv)
{
	memset(pot, 0, sizeof(*pot));
	pot->argc = argc;
	pot->argv = argv;
}

static inline char* _pot_following_value(struct pot* pot)
{
	if ((pot->index + 1) < pot->argc && (strlen(pot->argv[pot->index+1]) == 1 || pot->argv[pot->index+1][0] != '-')) {
		return pot->argv[pot->index + 1];
	} else {
		return NULL;
	}
}

static inline int pot_next(struct pot* pot)
{
	if (pot->index >= (pot->argc - 1)) return 0;

	pot->is_swtch = 0;
	pot->is_opt = 0;
	pot->is_arg = 0;
	pot->swtch = 0;
	pot->opt = NULL;
	pot->arg = NULL;
	pot->value = NULL;

	if (pot->swtch_idx > 0) {
		char* rest = pot->argv[pot->index] + pot->swtch_idx + 1;
		size_t rl = strlen(rest);
		if (rl > 0) {
			pot->is_swtch = 1;
			pot->swtch = rest[0];
			pot->value = (rl == 1) ? _pot_following_value(pot) : rest + 1;
			return 1;
		}
		pot->swtch_idx = 0;
	}

	char* arg = pot->argv[++pot->index];
	size_t l = strlen(arg);

	if (l >= 2 && arg[0] == '-' && arg[1] != '-') {
		pot->is_swtch = 1;
		pot->swtch = arg[1];
		pot->value = (l == 2) ? _pot_following_value(pot) : arg + 2;
	} else if (l >= 3 && arg[0] == '-' && arg[1] == '-') {
		pot->is_opt = 1;
		pot->opt = arg;
		pot->value = _pot_following_value(pot);
	} else {
		pot->is_arg = 1;
		pot->arg = arg;
	}

	return 1;
}

static inline void pot_reject(struct pot* pot)
{
	if (pot->is_arg) {
		fprintf(stderr, "Unexpected argument: %s\n", pot->arg);
	} else if (pot->is_opt) {
		fprintf(stderr, "Unknown option: %s\n", pot->opt);
	} else if (pot->is_swtch) {
		fprintf(stderr, "Unknown switch: -%c\n", pot->swtch);
	} else {
		fprintf(stderr, "pot_reject: invalid state");
	}
	exit(EXIT_FAILURE);
}

static inline void pot_swtch_no_arg(struct pot* pot)
{
	if (!pot->is_swtch) {
		fprintf(stderr, "pot_swtch_no_arg: invalid state; not a switch");
		exit(EXIT_FAILURE);
	}
	pot->swtch_idx++;
}

static inline char* pot_value(struct pot* pot)
{
	if (pot->value == NULL) {
		if (pot->is_swtch) {
			fprintf(stderr, "expected argument for switch -%c\n", pot->swtch);
		} else if (pot->is_opt) {
			fprintf(stderr, "expected argument for option %s\n", pot->opt);
		} else {
			fprintf(stderr, "pot_value: no value + invalid state\n");
		}
	} else {
		if (pot->is_swtch || pot->is_opt) {
			return pot->value;
		} else {
			fprintf(stderr, "pot_value: value + invalid state\n");
		}
	}
	exit(EXIT_FAILURE);
}

static inline int pot_int(struct pot* pot)
{
	char* value = pot_value(pot);
	return atoi(value);
	// XXX TODO fuck atoi; do your own
}

/*
static inline float pot_float(struct pot* pot)
{
	char* value = pot_value(pot);
	size_t l = strlen(value);
}
*/


#endif
