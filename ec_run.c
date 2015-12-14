#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <math.h>
#include <signal.h>

#include "ec_run.h"
#include "ecs.h"
#include "rng.h"
#include "sys.h"
#include "cpu.h"

#define EC_MAX_BOUNCES (150)

#ifndef M_PI
#define M_PI (3.141592653589793)
#endif

union vec3 {
	float s[3];
	struct { float x,y,z; };
};

static inline void vec3_dump(union vec3 v)
{
	printf(" [%.3f %.3f %.3f]\n", v.x, v.y, v.z);
}

static inline union vec3 vec3_add(union vec3 a, union vec3 b)
{
	union vec3 r;
	for (int i = 0; i < 3; i++) {
		r.s[i] = a.s[i] + b.s[i];
	}
	return r;
}


static inline union vec3 vec3_sub(union vec3 a, union vec3 b)
{
	union vec3 r;
	for (int i = 0; i < 3; i++) {
		r.s[i] = a.s[i] - b.s[i];
	}
	return r;
}

static inline union vec3 vec3_scale(union vec3 v, float scalar)
{
	for (int i = 0; i < 3; i++) v.s[i] *= scalar;
	return v;
}

static inline float vec3_dot(union vec3 a, union vec3 b)
{
	float sum = 0;
	for (int i = 0; i < 3; i++) sum += a.s[i] * b.s[i];
	return sum;
}



static inline union vec3 vec3_cross(union vec3 a, union vec3 b)
{
	union vec3 r = {
		.x = a.y*b.z - a.z*b.y,
		.y = a.z*b.x - a.x*b.z,
		.z = a.x*b.y - a.y*b.x
	};
	return r;
}

static inline union vec3 vec3_unit(union vec3 v)
{
	return vec3_scale(v, 1.0 / sqrtf(vec3_dot(v, v)));
}

struct biquad_setup {
	float a0, a1, a2, b1, b2;
	//float attenuation; // XXX must be derived
};

struct biquad_state {
	float a0, a1, a2, b1, b2, z0, z1;
};

struct material {
	int biquad_setup_index;
	/* attenuation is included in the biquad as a0 I guess? */
	float hardness;
	int emission_group_index; // -1 -> not an emitter
};

#if 0
struct emission_group {
	char* name;
};
#endif

struct microphone {
	//char* name;
	union vec3 position;
};

struct session {
	struct biquad_setup* biquad_setups;
	struct material* materials;

	int n_emission_groups;
	//struct emission_group* emission_groups;

	int n_microphones;
	struct microphone* microphones;

	int n_polys;
	uint8_t* poly_vertex_counts;
	uint8_t* poly_material_indices;
	union vec3* poly_vertex_cords;

	pthread_mutex_t commit_mutex;
	float** accumulators;
	int* accumulator_fds;

	/* setup */
	float sample_rate;
	float speed_of_sound_units_per_second;
	int impulse_length_samples;
	float attenuation_product_threshold;
	int n_rays_per_commit;

	/* derived from setup */
	float distance_units_to_samples;
	float ray_max_distance;
	float impulse_length_seconds;
};


struct worker_stats {
	uint64_t n_rays;
	uint64_t n_bounces;
	uint64_t n_hits;
	uint64_t n_rays_escaped;
	uint64_t n_commits;
};

static void worker_stats_accumulate(struct worker_stats* target, struct worker_stats* source)
{
	target->n_rays += source->n_rays;
	target->n_hits += source->n_hits;
	target->n_bounces += source->n_bounces;
	target->n_rays_escaped += source->n_rays_escaped;
	target->n_commits += source->n_commits;
}

static void worker_stats_print(struct worker_stats* st)
{
	printf("%ld rays\t%ld hits\t%ld bounces\t%ld escaped\t%ld commits\n",
		st->n_rays,
		st->n_hits,
		st->n_bounces,
		st->n_rays_escaped,
		st->n_commits);
}

struct worker {
	uint32_t index;

	struct rng rng;

	float** accumulators;

	struct worker_stats stats;

	int quit;
	pthread_mutex_t mutex;
};

struct worker_thread_args {
	pthread_t tid;
	struct session* session;
	struct worker* worker;
};


static inline void biquad_chain_impulse_run(
	int n_biquads,
	struct biquad_setup* setups,
	int* setup_indices,
	struct biquad_state* states,
	int n_samples,
	float* out)
{
	float y0 = 1.0f;
	for (int bi = 0; bi < n_biquads; bi++) {
		struct biquad_setup* setup = &setups[setup_indices[bi]];
		struct biquad_state* state = &states[bi];
		float x = y0;
		y0 *= setup->a0;
		state->a0 = setup->a0;
		state->a1 = setup->a1;
		state->a2 = setup->a2;
		state->b1 = setup->b1;
		state->b2 = setup->b2;
		state->z0 = x * setup->a1 + y0 * setup->b1;
		state->z1 = x * setup->a2 + y0 * setup->b2;
	}
	out[0] = y0;

	for (int si = 1; si < n_samples; si++) {
		float y = 0.0f;
		for (int bi = 0; bi < n_biquads; bi++) {
			struct biquad_state* state = &states[bi];
			float x = y;
			y *= state->a0;
			y += state->z0;
			state->z0 = x * state->a1 + y * state->b1 + state->z1;
			state->z1 = x * state->a2 + y * state->b2;
		}
		out[si] = y;
	}
}

static void worker_init(struct worker* w, struct session* s, uint32_t index)
{
	memset(w, 0, sizeof(*w));

	w->index = index;

	struct timespec tspec;
	clock_gettime(CLOCK_REALTIME, &tspec);
	uint32_t seed = tspec.tv_nsec + index;

	rng_seed(&w->rng, seed);
	pthread_mutex_init(&w->mutex, NULL);

	int n_accumulators = s->n_emission_groups * s->n_microphones;
	w->accumulators = calloc_or_die(n_accumulators, sizeof(*w->accumulators));
	for (int i = 0; i < n_accumulators; i++) {
		w->accumulators[i] = calloc_or_die(s->impulse_length_samples, sizeof(**w->accumulators));
	}
}

static inline float ray_poly_intersection(union vec3 p0, union vec3 dir, int nv, union vec3* vs, union vec3* x)
{
	union vec3 normal = vec3_unit(vec3_cross(vec3_sub(vs[1], vs[0]), vec3_sub(vs[2], vs[0])));
	float dot = vec3_dot(dir, normal);
	if (dot >= 0.0f) return -1;
	float t = vec3_dot(vec3_sub(vs[0], p0), normal) / dot;

	*x = vec3_add(p0, vec3_scale(dir, t));

	int previ = nv-1;
	for (int i = 0; i < nv; i++) {
		union vec3 a = vs[previ];
		union vec3 b = vs[i];
		if (vec3_dot(vec3_sub(*x, a), vec3_cross(normal, vec3_sub(b, a))) < 0) {
			return -2 - i;
		}
		previ = i;
	}

	return t;
}


static inline union vec3 vec3_random_unit(struct rng* rng)
{
	/*
	(33) Generate random point on sphere
	http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf
	*/

	float r1 = rng_float(rng);
	float r2 = rng_float(rng);

	float phi = 2 * M_PI * r1;

	float s = sqrtf(r2 * (1 - r2));

	return (union vec3) {
		.x = 2.0 * cosf(phi) * s,
		.y = 2.0 * sinf(phi) * s,
		.z = 1.0 - 2.0 * r2
	};
}

static inline union vec3 vec3_random_hemisphere_cosine(struct rng* rng)
{
	/*
	(19a) Polar map
	http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf
	*/

	float r1 = rng_float(rng);
	float r2 = rng_float(rng);

	float r = sqrtf(r1);
	float phi = 2 * M_PI * r2;

	/*
	(35) Generate random direction on unit hemisphere proportional to
	cosine-weight solid angle
	http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf
	*/

	return (union vec3) {
		.x = r * cosf(phi),
		.y = r * sinf(phi),
		.z = sqrtf(1 - r*r)
	};
}

static void ray_pew_pew(
	struct session* restrict session,
	struct worker* restrict worker,
	struct worker_stats* restrict stats,
	int micidx)
{
	stats->n_rays++;

	//float attenuation_product = 1.0f;
	float distance = 0.0f;
	float ray_max_distance = session->ray_max_distance;

	struct microphone* mic = &session->microphones[micidx];

	// initial ray
	union vec3 ray_origin = mic->position;
	union vec3 ray_direction = vec3_random_unit(&worker->rng);

	int biquad_stack[EC_MAX_BOUNCES];
	int biquad_stack_top = 0;

	do {
		stats->n_bounces++;

		int nearest_poly_index = -1;
		float nearest_t = 0;
		union vec3 nearest_position;
		int nearest_voff;

		int voff = 0;
		for (int i = 0; i < session->n_polys; i++) {
			int nv = session->poly_vertex_counts[i];
			union vec3 x;
			float t = ray_poly_intersection(
				ray_origin,
				ray_direction,
				nv,
				&session->poly_vertex_cords[voff],
				&x);
			voff += nv;

			if (t < 0) continue;

			if (nearest_poly_index == -1 || t < nearest_t) {
				nearest_poly_index = i;
				nearest_t = t;
				nearest_position = x;
				nearest_voff = voff-nv;
			}
		}

		if (nearest_poly_index < 0) {
			stats->n_rays_escaped++;
			return;
		}

		distance += nearest_t;
		// stats->total_distance += nearest_t?

		struct material* m = &session->materials[session->poly_material_indices[nearest_poly_index]];

		if (m->emission_group_index < 0) {
			/* hit a reflective surface; bounce */

			int biquad_setup_index = m->biquad_setup_index;

			/* FIXME doesn't make sense to calculate and act on the
			 * attenuation product before the biquad attenuation
			 * has a meaningful value */
			#if 0
			attenuation_product *= session->biquad_setups[biquad_setup_index].attenuation;
			if (attenuation_product < session->attenuation_product_threshold) {
				/* signal got too weak; abandoning ray */
				return;
			}
			#endif

			biquad_stack[biquad_stack_top++] = biquad_setup_index;

			// calculate new ray, and continue
			ray_origin = nearest_position;
			{
				union vec3* vs = &session->poly_vertex_cords[nearest_voff];
				union vec3 bx = vec3_unit(vec3_sub(vs[1], vs[0]));
				union vec3 bz = vec3_unit(vec3_cross(bx, vec3_sub(vs[2], vs[0])));
				union vec3 by = vec3_cross(bz, bx);
				union vec3 d = vec3_random_hemisphere_cosine(&worker->rng);
				ray_direction = vec3_add(
					vec3_scale(bx, d.x),
					vec3_add(
						vec3_scale(by, d.y),
						vec3_scale(bz, d.z)
					)
				);
			}
		} else {
			/* hit an emitter; calculate contribution */

			int acci = micidx * session->n_microphones + m->emission_group_index;
			float* acc = worker->accumulators[acci];

			int soff = distance * session->distance_units_to_samples + rng_float(&worker->rng);
			int n = 64;
			if (soff + n > session->impulse_length_samples) {
				n = session->impulse_length_samples - soff;
				if (n < 2) return;
			}

			float* accoff = acc + soff;

			struct biquad_state states[EC_MAX_BOUNCES];
			biquad_chain_impulse_run(
				biquad_stack_top,
				session->biquad_setups,
				biquad_stack,
				states,
				n,
				accoff);

			stats->n_hits++;
			return;
		}
	} while(distance < ray_max_distance && biquad_stack_top < EC_MAX_BOUNCES);
}

static void commit(
	struct session* restrict session,
	struct worker* restrict worker)
{
	int n_accumulators = session->n_emission_groups * session->n_microphones;
	int impulse_length_samples = session->impulse_length_samples ;

	/* offload local accumulators to global one */
	pthread_mutex_lock(&session->commit_mutex);
	for (int i = 0; i < n_accumulators; i++) {
		float* src = worker->accumulators[i];
		float* dst = session->accumulators[i];
		for (int j = 0; j < impulse_length_samples; j++) {
			dst[j] += src[j];
		}
	}
	pthread_mutex_unlock(&session->commit_mutex);

	/* clear local accumulators */
	for (int i = 0; i < n_accumulators; i++) {
		memset(worker->accumulators[i], 0, sizeof(**worker->accumulators) * impulse_length_samples);
	}

	worker->stats.n_commits++;
}

static void run(
	struct session* restrict session,
	struct worker* restrict worker)
{
	const int n_rays_per_microphone = 1024;
	int n_microphones = session->n_microphones;

	int n_rays_since_commit = 0;
	for (;;) {
		for (int micidx = 0; micidx < n_microphones; micidx++) {
			struct worker_stats stats;
			memset(&stats, 0, sizeof(stats));
			for (int j = 0; j < n_rays_per_microphone; j++) {
				ray_pew_pew(session, worker, &stats, micidx);
			}
			n_rays_since_commit += n_rays_per_microphone;

			int quitting = 0;
			pthread_mutex_lock(&worker->mutex);
			worker_stats_accumulate(&worker->stats, &stats);
			quitting = worker->quit;
			pthread_mutex_unlock(&worker->mutex);

			int must_commit = n_rays_since_commit > session->n_rays_per_commit;

			if (quitting || must_commit) commit(session, worker);
			if (must_commit) n_rays_since_commit = 0;
			if (quitting) return;
		}
	}
}

static void* worker_thread_start(void* arg)
{
	struct worker_thread_args* args = arg;
	run(args->session, args->worker);
	return NULL;
}


static volatile int master_quit = 0;
void sigint_handler(int sig)
{
	printf("quitting...\n");
	master_quit = 1;
}

void ec_run(char* path, int n_workers)
{

	ecs_info(path);

	struct ecs ecs;
	ecs_open(&ecs, path);

	struct session* s = calloc_or_die(1, sizeof(*s));

	{
		struct ecs_initialization* ei = ecs.initialization;

		if (ei == NULL) {
			fprintf(stderr, "ecs file is not initialized, use ec init\n");
			exit(EXIT_FAILURE);
		}

		s->sample_rate = ei->sample_rate;
		s->speed_of_sound_units_per_second = ei->speed_of_sound_units_per_second;
		s->impulse_length_samples = ei->impulse_length_samples;
		s->attenuation_product_threshold = ei->attenuation_product_threshold;

		/* derived units */
		s->distance_units_to_samples = s->sample_rate / s->speed_of_sound_units_per_second;
		s->impulse_length_seconds = s->impulse_length_samples / s->sample_rate;
		s->ray_max_distance = s->speed_of_sound_units_per_second * s->impulse_length_seconds;
		s->impulse_length_samples = s->impulse_length_seconds * s->sample_rate;
	}

	{
		int n = ecs_get_biquad_count(&ecs);
		s->biquad_setups = calloc_or_die(n, sizeof(struct biquad_setup));
		for (int i = 0; i < n; i++) {
			struct biquad_setup* bqs = &s->biquad_setups[i];
			struct ecs_biquad* bq = ecs_get_biquad(&ecs, i);
			bqs->a0 = bq->a0;
			bqs->a1 = bq->a1;
			bqs->a2 = bq->a2;
			bqs->b1 = bq->b1;
			bqs->b2 = bq->b2;
			/* I need to calculate the frequency response of the
			 * biquad (using z-transform?) and figure out which
			 * frequency is attenuated the least. i.e. find the
			 * peak of the frequency response curve; this is the
			 * attenuation value */
			//bqs->attenuation = 0.99f; // XXX XXX
		}
	}

	{
		int n = ecs_get_emission_group_count(&ecs);
		s->n_emission_groups = n;
	}

	{
		int n = ecs_get_microphone_count(&ecs);
		s->n_microphones = n;
		s->microphones = calloc_or_die(n, sizeof(struct microphone));
		for (int i = 0; i < n; i++) {
			struct microphone* m = &s->microphones[i];
			struct ecs_microphone* em = ecs_get_microphone(&ecs, i);
			for (int i = 0; i < 3; i++) m->position.s[i] = em->position[i];
		}
	}

	{
		int n = ecs_get_material_count(&ecs);
		s->materials = calloc_or_die(n, sizeof(struct material));
		for (int i = 0; i < n; i++) {
			struct material* m = &s->materials[i];
			struct ecs_material* em = ecs_get_material(&ecs, i);
			m->biquad_setup_index = em->biquad_index;
			m->emission_group_index = em->emission_group_index;
			m->hardness = em->hardness;
		}
	}

	{
		struct ecs_poly poly;

		ecs_poly_iterator it = 0;

		int n = 0;
		int nv = 0;
		while (ecs_poly_iterator_next(&ecs, &it, &poly)) {
			n++;
			nv += poly.n_vertices;
		}

		s->n_polys = n;
		s->poly_material_indices = calloc_or_die(n, sizeof(*s->poly_material_indices));
		s->poly_vertex_counts = calloc_or_die(n, sizeof(*s->poly_vertex_counts));
		s->poly_vertex_cords = calloc_or_die(nv, sizeof(*s->poly_vertex_cords));

		it = 0;
		int index = 0;
		int o = 0;
		while (ecs_poly_iterator_next(&ecs, &it, &poly)) {
			s->poly_material_indices[index] = poly.material_index;
			s->poly_vertex_counts[index] = poly.n_vertices;
			for (int i = 0; i < poly.n_vertices; i++) {
				for (int j = 0; j < 3; j++) {
					s->poly_vertex_cords[o].s[j] = poly.vertex_data[i*3 + j];
				}
				o++;
			}
			index++;
		}
	}

	s->accumulators = ecs.accumulators;

	int use_cpu_count = n_workers <= 0;
	n_workers = use_cpu_count ? cpu_count() : n_workers;

	s->n_rays_per_commit = 1<<20;

	printf("starting path tracing with these parameters:\n");
	printf("  worker count: %d%s\n", n_workers, use_cpu_count ? " (ncpus)" : "");
	printf("  rays per commit: %d\n", s->n_rays_per_commit);

	struct worker* workers = calloc_or_die(n_workers, sizeof(*workers));
	struct worker_thread_args* thread_args = calloc_or_die(n_workers, sizeof(*thread_args));

	pthread_attr_t attr;
	pthread_attr_init(&attr);

	for (int i = 0; i < n_workers; i++) {
		struct worker_thread_args* args = &thread_args[i];
		worker_init(&workers[i], s, i);
		args->session = s;
		args->worker = &workers[i];
		pthread_create(&args->tid, &attr, worker_thread_start, args);
	}

	master_quit = 0;
	signal(SIGINT, sigint_handler);

	struct worker_stats stats;

	while (!master_quit) {
		sleep(1);
		memset(&stats, 0, sizeof(stats));
		for (int i = 0; i < n_workers; i++) {
			struct worker* worker = &workers[i];
			pthread_mutex_lock(&worker->mutex);
			worker->quit = master_quit;
			worker_stats_accumulate(&stats, &worker->stats);
			pthread_mutex_unlock(&worker->mutex);

		}
		worker_stats_print(&stats);
	}

	printf("waiting for workers to quit...\n");
	memset(&stats, 0, sizeof(stats));
	for (int i = 0; i < n_workers; i++) {
		struct worker* worker = &workers[i];
		pthread_join(thread_args[i].tid, NULL);
		worker_stats_accumulate(&stats, &worker->stats);
	}
	worker_stats_print(&stats);

	printf("done\n");

	ecs_close(&ecs);
}
