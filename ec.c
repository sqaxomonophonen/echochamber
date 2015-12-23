#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <getopt.h>
#include <string.h>

#include "sys.h"
#include "rng.h"
#include "pot.h"
#include "ecs.h"
#include "ec_run.h"


static void usage(char* prg, int full, int exit_status)
{
	FILE* out = (exit_status == EXIT_SUCCESS) ? stdout : stderr;

	fprintf(out, "usage: %s <command> [<args>]\n\n", prg);
	if (full) {
		fprintf(out, "commands:\n");
		fprintf(out, "   init <input.ecs>     Initialize echochamber scene\n");
		fprintf(out, "   run <input.ecs>      Run path tracer\n");
		fprintf(out, "   mixdown <input.ecs>  Create a mixdown\n");
		fprintf(out, "   help [<command>]     Get help!\n");
	}

	exit(exit_status);
}

static char* scene_file(struct pot* pot)
{
	while (pot_next(pot)) {
		if (!pot->is_arg) pot_reject(pot);
		return pot->arg;
	}
	fprintf(stderr, "No scene file specified!\n");
	exit(EXIT_FAILURE);
}

static void cmd_init(struct pot* pot)
{
	char* path = scene_file(pot);

	int force = 0;

	int sample_rate = 48000;
	float speed_of_sound_units_per_second = 343.2f;
	int impulse_length_samples = 5 * sample_rate;
	int fir_length = 64;
	enum ecs_fir_window_function fir_window_function = ECS_FIR_WINDOW_FN_KAISER_BESSEL;
	int indirect_only = 0;

	while (pot_next(pot)) {
		if ((pot->is_swtch && pot->swtch == 'f') || (pot->is_opt && strcmp(pot->opt, "--force") == 0)) {
			force = 1;
		}
		else if ((pot->is_swtch && pot->swtch == 'r') || (pot->is_opt && strcmp(pot->opt, "--sample-rate") == 0)) {
			sample_rate = pot_int(pot);
		}
		else if ((pot->is_swtch && pot->swtch == 's') || (pot->is_opt && strcmp(pot->opt, "--sound-speed") == 0)) {
			speed_of_sound_units_per_second = pot_int(pot); // XXX pot_float
		}
		else if ((pot->is_swtch && pot->swtch == 'l') || (pot->is_opt && strcmp(pot->opt, "--impulse-length") == 0)) {
			impulse_length_samples = pot_int(pot) * sample_rate;
		}
		else if ((pot->is_swtch && pot->swtch == 'I') || (pot->is_opt && strcmp(pot->opt, "--indirect-only") == 0)) {
			indirect_only = 1;
		}
		else {
			pot_reject(pot);
		}
	}

	ecs_init(
		path,
		force,
		sample_rate,
		speed_of_sound_units_per_second,
		impulse_length_samples,
		fir_length,
		fir_window_function,
		indirect_only
	);

	ecs_info(path);
}

static void cmd_info(struct pot* pot)
{
	char* ecs = scene_file(pot);
	while (pot_next(pot)) pot_reject(pot);
	ecs_info(ecs);
}

static void cmd_run(struct pot* pot)
{
	char* ecs = scene_file(pot);
	while (pot_next(pot)) pot_reject(pot);

	ec_run(ecs, -1);
}

static void cmd_mixdown(struct pot* pot)
{
	char* ecs = scene_file(pot);
	ecs_mixdown(ecs);
}

static void cmd_help(struct pot* pot)
{
}

static void cmd(struct pot* pot)
{
	while (pot_next(pot)) {
		if (!pot->is_arg) pot_reject(pot);

		char* cmd = pot->arg;

		if (strcmp(cmd, "init") == 0) {
			cmd_init(pot);
		} else if(strcmp(cmd, "info") == 0) {
			cmd_info(pot);
		} else if(strcmp(cmd, "run") == 0) {
			cmd_run(pot);
		} else if(strcmp(cmd, "mixdown") == 0) {
			cmd_mixdown(pot);
		} else if(strcmp(cmd, "help") == 0) {
			cmd_help(pot);
		} else {
			fprintf(stderr, "'%s' is not a valid command\n", cmd);
			exit(EXIT_FAILURE);
		}

		exit(EXIT_SUCCESS);
	}
}

int main(int argc, char** argv)
{
	char* prg = argv[0];

	struct pot pot;
	pot_init(&pot, argc, argv);

	cmd(&pot);

	usage(prg, 1, EXIT_FAILURE);
	return 1;
}
