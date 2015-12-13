#ifndef ECS_H

#include <stdlib.h>
#include <stdint.h>

struct ecs_initialization {
	int sample_rate;
	float speed_of_sound_units_per_second;
	int impulse_length_samples;
	float attenuation_product_threshold;
	// TODO stats!
};

struct ecs {
	int initialized;
	size_t sz;
	void* ptr;

	/* mandatory blocks */
	void* _bqd0_block;
	size_t _bqd0_sz;

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

struct ecs_biquad {
	float a0, a1, a2, b1, b2;
};

struct ecs_emission_group {
	char name[64];
};

struct ecs_microphone {
	char name[64];
	float position[3];
};

struct ecs_material {
	uint32_t emission_group_index;
	uint32_t biquad_index;
	float hardness;
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
	float attenuation_product_threshold
);

void ecs_mixdown(char* path);

int ecs_get_biquad_count(struct ecs* ecs);
struct ecs_biquad* ecs_get_biquad(struct ecs* ecs, int i);

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
