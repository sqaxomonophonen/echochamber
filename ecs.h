#ifndef ECS_H

#include <stdlib.h>
#include <stdint.h>

enum ecs_fir_window_function {
	ECS_FIR_WINDOW_FN_KAISER_BESSEL = 1
};

static inline const char* ecs_fir_window_function_string(enum ecs_fir_window_function fn)
{
	switch (fn) {
		case ECS_FIR_WINDOW_FN_KAISER_BESSEL: return "Kaiser-Bessel";
	}
	return "?";
}

struct ecs_worker_stats {
	uint64_t n_rays;
	uint64_t n_bounces;
	uint64_t n_hits;
	uint64_t n_rays_escaped;
	uint64_t n_commits;
	double units_traveled;
};

static inline void ecs_worker_stats_accumulate(struct ecs_worker_stats* target, struct ecs_worker_stats* source)
{
	target->n_rays += source->n_rays;
	target->n_hits += source->n_hits;
	target->n_bounces += source->n_bounces;
	target->n_rays_escaped += source->n_rays_escaped;
	target->n_commits += source->n_commits;
	target->units_traveled += source->units_traveled;
}


struct ecs_global_stats {
	uint64_t render_time_ms;
};

struct ecs_initialization {
	int sample_rate;
	float speed_of_sound_units_per_second;
	int impulse_length_samples;
	int fir_length;
	enum ecs_fir_window_function fir_window_function;
	int indirect_only;

	struct ecs_worker_stats worker_stats;
	struct ecs_global_stats global_stats;
};

struct ecs {
	int initialized;
	size_t sz;
	void* ptr;

	/* mandatory blocks */
	void* _flt0_block;
	size_t _flt0_sz;

	void* _emg0_block;
	size_t _emg0_sz;

	void* _mic0_block;
	size_t _mic0_sz;

	void* _mat0_block;
	size_t _mat0_sz;

	void* _ply0_block;
	size_t _ply0_sz;

	struct ecs_initialization* initialization;

	size_t* accumulator_sz;
	float** accumulators;
};

struct ecs_emission_group {
	char name[64];
};

struct ecs_microphone {
	char name[64];
	float position[3];
};

struct ecs_material {
	int32_t emission_group_index;
	int32_t diffuse_filter_index;
	int32_t specular_filter_index;
	int32_t hardness;
};


void ecs_open(struct ecs* ecs, char* path);
void ecs_close(struct ecs* ecs);
void ecs_info(char* path);
void ecs_init(
	char* path,
	int force,
	int sample_rate,
	float speed_of_sound_units_per_second,
	int impulse_length_samples,
	int fir_length,
	enum ecs_fir_window_function fir_window_function,
	int indirect_only
);

void ecs_mixdown(char* path);

char* ecs_get_filter_strings(struct ecs* ecs, size_t* sz);

int ecs_get_emission_group_count(struct ecs* ecs);
struct ecs_emission_group* ecs_get_emission_group(struct ecs* ecs, int i);

int ecs_get_microphone_count(struct ecs* ecs);
struct ecs_microphone* ecs_get_microphone(struct ecs* ecs, int i);

int ecs_get_material_count(struct ecs* ecs);
struct ecs_material* ecs_get_material(struct ecs* ecs, int i);


struct ecs_poly {
	uint32_t material_index;
	uint32_t n_vertices;
	float* vertex_data;
};

typedef size_t ecs_poly_iterator;
int ecs_poly_iterator_next(struct ecs* ecs, ecs_poly_iterator* it, struct ecs_poly* poly);


#define ECS_H

#endif
