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
#include "flt0.h"
#include "flt0.yy.h"

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

struct material {
	int filter_index;
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
	float* filter_coefficients;
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
	int fir_length;
	int n_rays_per_commit;
	int indirect_only;

	/* derived from setup */
	float distance_units_to_samples;
	float ray_max_distance;
	float impulse_length_seconds;
	float* fir_window;
};

static void stats_print(struct ecs_initialization* init)
{
	struct ecs_global_stats* gst = &init->global_stats;
	struct ecs_worker_stats* wst = &init->worker_stats;
	printf("%lds\t%ld rays\t%ld hits\t%ld bounces\t%ld escaped\t%ld commits\t%.2e units traveled\n",
		gst->render_time_ms / 1000,
		wst->n_rays,
		wst->n_hits,
		wst->n_bounces,
		wst->n_rays_escaped,
		wst->n_commits,
		wst->units_traveled);
}

struct worker {
	uint32_t index;

	struct rng rng;
	float* filter_coefficients;

	float** accumulators;

	float* coefficient_product;

	struct ecs_worker_stats stats;

	int quit;
	pthread_mutex_t mutex;
};

struct worker_thread_args {
	pthread_t tid;
	struct session* session;
	struct worker* worker;
};


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

	w->coefficient_product = calloc_or_die(s->fir_length, sizeof(*w->coefficient_product));
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
	struct ecs_worker_stats* restrict stats,
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

	int fir_length = session->fir_length;
	float* coefficient_product = worker->coefficient_product;

	/* reset all coefficients to 1+0j */
	for (int i = 0; i < fir_length; i++) {
		coefficient_product[i] = i&1 ? 0.0f : 1.0f;
	}

	int bounces = 0;

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
		stats->units_traveled += (double)nearest_t;

		struct material* m = &session->materials[session->poly_material_indices[nearest_poly_index]];

		if (m->emission_group_index < 0) {
			/* hit a reflective surface; bounce */

			bounces++;

			float* filter_coefficients = session->filter_coefficients + fir_length * m->filter_index;
			for (int i = 0; i < fir_length/2; i++) {
				/* (a+bi)*(c+di) = (a*c - b*d) + (b*c + a*d)i */
				int idx = i<<1;
				coefficient_product[idx] =
					coefficient_product[idx] * filter_coefficients[idx]
					- coefficient_product[idx+1] * filter_coefficients[idx+1];
				coefficient_product[idx+1] =
					coefficient_product[idx+1] * filter_coefficients[idx]
					+ coefficient_product[idx] * filter_coefficients[idx+1];
			}


			// calculate new ray, and continue
			ray_origin = nearest_position;
			#if 1
			{
				/* diffuse reflection */
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
			#else
			{
				/* perfect reflection */
				union vec3* vs = &session->poly_vertex_cords[nearest_voff];
				union vec3 normal = vec3_unit(vec3_cross(vec3_sub(vs[1], vs[0]), vec3_sub(vs[2], vs[0])));
				union vec3 incident = vec3_unit(ray_direction);
				ray_direction = vec3_unit(vec3_sub(incident, vec3_scale(normal, 2.0 * vec3_dot(normal, incident))));
			}
			#endif
		} else {
			/* hit an emitter; calculate contribution */

			if (session->indirect_only && bounces == 0) {
				return;
			}

			float fsoff = distance * session->distance_units_to_samples;
			int soff = fsoff;

			int N = session->fir_length;
			int s0 = soff - N/2;
			float s0f = distance * session->distance_units_to_samples - N/2;
			int s1 = soff + N/2;

			int impulse_length_samples = session->impulse_length_samples;

			int s0s = s0 < 0 ? 0 : s0;
			int s1s = s1 > impulse_length_samples-1 ? impulse_length_samples-1 : s1;

			if (s1s-s0s <= 0) {
				/* outside range; cull */
				return;
			}

			int acci = micidx * session->n_microphones + m->emission_group_index;
			float* acc = worker->accumulators[acci];

			int Nhalf = N>>1;
			float C = (2 * M_PI) / (float)N;
			float* fir_window = session->fir_window;
			for (int i = s0s; i < s1s; i++) {
				float k = (float)i - s0f;
				float Ck = C * k;

				float y = 1;

				for (int n = 0; n < Nhalf; n++) {
					float phi0 = Ck * (float)(n+1);
					float phi1 = Ck * (float)(N-n-1);
					/*
					s x*cos(phi) + x*sin(phi)*i
					(a+bi)*(c+di) = (a*c - b*d) + (b*c + a*d)i
					*/
					float impfft = (n&1) ? 1 : -1;
					y +=
						impfft * coefficient_product[n<<1] * (cosf(phi0) + cosf(phi1)) -
						impfft * coefficient_product[(n<<1)+1] * (sinf(phi0) - sinf(phi1));
				}
				acc[i] = y * fir_window[i - s0];
			}

			stats->n_hits++;
			return;
		}
	} while(distance < ray_max_distance && bounces < EC_MAX_BOUNCES);
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
			struct ecs_worker_stats stats;
			memset(&stats, 0, sizeof(stats));
			for (int j = 0; j < n_rays_per_microphone; j++) {
				ray_pew_pew(session, worker, &stats, micidx);
			}
			n_rays_since_commit += n_rays_per_microphone;

			int quitting = 0;
			pthread_mutex_lock(&worker->mutex);
			ecs_worker_stats_accumulate(&worker->stats, &stats);
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

static uint64_t timespec_diff_ms(struct timespec* t0, struct timespec* t1)
{
	uint64_t t0ms = t0->tv_sec * 1000L + t0->tv_nsec / 1000000L;
	uint64_t t1ms = t1->tv_sec * 1000L + t1->tv_nsec / 1000000L;
	return t1ms - t0ms;
}

static int bz0_fp_frq_sort(const void* va, const void* vb)
{
	const struct flt0_bz0_point* a = va;
	const struct flt0_bz0_point* b = vb;
	return a->frq - b->frq;
}

static float bz0_crv(float t, float y0, float y1)
{
	float t1 = 1-t;
	return
		y0 * t1*t1*t1 +
		y0 * 3*t * t1*t1 +
		y1 * 3*t*t * t1 +
		y1 * t*t*t;
}


static float bessel_I0(float x)
{
	float d = 0;
	float ds = 1;
	float s = 1;
	do {
		d += 2;
		ds *= (x*x) / (d*d);
		s += ds;
	} while (ds > s*1e-6);
	return s;
}

static void calc_kaiser_bessel_window(struct session* s, float att)
{
	float alpha = 0;
	if (att > 50) {
		alpha = 0.1102 * (att - 8.7);
	} else if (att >= 21) {
		alpha = 0.5842 * powf(att - 21, 0.4) + 0.07886 * (att - 21);
	}

	float I0alpha = bessel_I0(alpha);
	int n = s->fir_length;
	for (int i = 0; i < n; i++) {
		float x = ((float)i - ((float)n / 2.0f)) * 2.0f;
		s->fir_window[i] = bessel_I0(alpha * sqrtf(1 - ((x*x) / (float)(n*n)))) / I0alpha;
	}
}
void ec_run(char* path, int n_workers)
{

	ecs_info(path);

	struct ecs ecs;
	ecs_open(&ecs, path);

	struct session* session = calloc_or_die(1, sizeof(*session));

	{
		struct ecs_initialization* ei = ecs.initialization;

		if (ei == NULL) {
			fprintf(stderr, "ecs file is not initialized, use ec init\n");
			exit(EXIT_FAILURE);
		}

		session->sample_rate = ei->sample_rate;
		session->speed_of_sound_units_per_second = ei->speed_of_sound_units_per_second;
		session->impulse_length_samples = ei->impulse_length_samples;
		session->fir_length = ei->fir_length;
		session->indirect_only = ei->indirect_only;

		/* derived units */
		session->distance_units_to_samples = session->sample_rate / session->speed_of_sound_units_per_second;
		session->impulse_length_seconds = session->impulse_length_samples / session->sample_rate;
		session->ray_max_distance = session->speed_of_sound_units_per_second * session->impulse_length_seconds;
		session->impulse_length_samples = session->impulse_length_seconds * session->sample_rate;

		session->fir_window = calloc_or_die(session->fir_length, sizeof(*session->fir_window));

		switch (ei->fir_window_function) {
			case ECS_FIR_WINDOW_FN_KAISER_BESSEL:
				calc_kaiser_bessel_window(session, 40.0f);
				break;
			default:
				fprintf(stderr, "invalid window function %d\n", ei->fir_window_function);
		}
	}

	{
		size_t sz;
		char* ss = ecs_get_filter_strings(&ecs, &sz);

		/* count strings */
		int n = 0;
		for (size_t i = 0; i < sz; i++) if (ss[i] == 0) n++;

		session->filter_coefficients = calloc_or_die(session->fir_length * n, sizeof(*session->filter_coefficients));

		int filter_index = 0;
		size_t of = 0;
		while (of < sz) {
			float* coeffs = session->filter_coefficients + filter_index * session->fir_length;

			char* s = ss + of;
			flt0_scan_string(s);

			int tok = flt0lex();

			if (tok == 0) {
				/* "NULL" material */
				memset(coeffs, 0, sizeof(*coeffs) * session->fir_length);
			} else if (tok == FLT0_BZ0) {
				tok = flt0lex();
				if (tok != FLT0_BEGIN) flt0_unexpected_token(tok);

				#define MAX_FLT0_POINTS (256)
				struct flt0_bz0_point fps[MAX_FLT0_POINTS];
				int fpi = 0;

				for (;;) {
					if (fpi >= MAX_FLT0_POINTS) {
						fprintf(stderr, "more than %d filter points!\n", MAX_FLT0_POINTS);
						exit(EXIT_FAILURE);
					}
					struct flt0_bz0_point fp;

					tok = flt0lex();
					if (tok != FLT0_FLOAT) flt0_unexpected_token(tok);
					fp.frq = flt0_get_floatval();

					tok = flt0lex();
					if (tok != FLT0_SUBSEP) flt0_unexpected_token(tok);

					tok = flt0lex();
					if (tok != FLT0_FLOAT) flt0_unexpected_token(tok);
					fp.att = flt0_get_floatval();

					tok = flt0lex();
					if (tok != FLT0_SUBSEP) flt0_unexpected_token(tok);

					tok = flt0lex();
					if (tok != FLT0_FLOAT) flt0_unexpected_token(tok);
					fp.pha = flt0_get_floatval();

					fps[fpi] = fp;

					fpi++;

					tok = flt0lex();
					if (tok == FLT0_SEP) {
						continue;
					} else if (tok == FLT0_END) {
						break;
					} else {
						flt0_unexpected_token(tok);
					}
				}

				tok = flt0lex();
				if (tok != 0) flt0_unexpected_token(tok);

				/* sort filter points by frequency */
				qsort(fps, fpi, sizeof(*fps), bz0_fp_frq_sort);

				int N = session->fir_length >> 1;
				for (int i = 0; i < N; i++) {
					float* coeff = coeffs + i*2;

					float fft_freq = ((float)(i+1) / (float)N) * 0.5f * (float)session->sample_rate;

					float att = 0;
					float pha = 0;

					for (int j = 0; j < fpi; j++) {
						struct flt0_bz0_point* fp = &fps[j];
						float fp_freq = fp->frq;

						if ((j == 0 && fft_freq < fp_freq) || (j == fpi-1 && fft_freq > fp_freq)) {
							att = fp->att;
							pha = fp->pha;
							break;
						} else if (j > 0 && fft_freq >= fps[j-1].frq  && fft_freq < fp_freq) {
							float x0 = logf(fps[j-1].frq);
							float x1 = logf(fp_freq);
							float x = logf(fft_freq);
							float t = (x-x0)/(x1-x0);
							att = bz0_crv(t, fps[j-1].att, fps[j].att);
							pha = bz0_crv(t, fps[j-1].pha, fps[j].pha);
							break;
						}
					}

					coeff[0] = att;
					/* at the nyquist frequency the value
					 * must be its own conjugate, i.e. a
					 * purely real value, so cancel out the
					 * phase shift when i is N-1 */
					coeff[1] = i==(N-1) ? 0 : pha;
				}
			} else {
				flt0_unexpected_token(tok);
			}

			of += strlen(s) + 1;

			filter_index++;
		}
	}

	{
		int n = ecs_get_emission_group_count(&ecs);
		session->n_emission_groups = n;
	}

	{
		int n = ecs_get_microphone_count(&ecs);
		session->n_microphones = n;
		session->microphones = calloc_or_die(n, sizeof(struct microphone));
		for (int i = 0; i < n; i++) {
			struct microphone* m = &session->microphones[i];
			struct ecs_microphone* em = ecs_get_microphone(&ecs, i);
			for (int i = 0; i < 3; i++) m->position.s[i] = em->position[i];
		}
	}

	{
		int n = ecs_get_material_count(&ecs);
		session->materials = calloc_or_die(n, sizeof(struct material));
		for (int i = 0; i < n; i++) {
			struct material* m = &session->materials[i];
			struct ecs_material* em = ecs_get_material(&ecs, i);
			m->filter_index = em->filter_index;
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

		session->n_polys = n;
		session->poly_material_indices = calloc_or_die(n, sizeof(*session->poly_material_indices));
		session->poly_vertex_counts = calloc_or_die(n, sizeof(*session->poly_vertex_counts));
		session->poly_vertex_cords = calloc_or_die(nv, sizeof(*session->poly_vertex_cords));

		it = 0;
		int index = 0;
		int o = 0;
		while (ecs_poly_iterator_next(&ecs, &it, &poly)) {
			session->poly_material_indices[index] = poly.material_index;
			session->poly_vertex_counts[index] = poly.n_vertices;
			for (int i = 0; i < poly.n_vertices; i++) {
				for (int j = 0; j < 3; j++) {
					session->poly_vertex_cords[o].s[j] = poly.vertex_data[i*3 + j];
				}
				o++;
			}
			index++;
		}
	}

	session->accumulators = ecs.accumulators;

	int use_cpu_count = n_workers <= 0;
	n_workers = use_cpu_count ? cpu_count() : n_workers;

	session->n_rays_per_commit = 1<<20;

	printf("starting path tracing with these parameters:\n");
	printf("  worker count: %d%s\n", n_workers, use_cpu_count ? " (ncpus)" : "");
	printf("  rays per commit: %d\n", session->n_rays_per_commit);

	struct worker* workers = calloc_or_die(n_workers, sizeof(*workers));
	struct worker_thread_args* thread_args = calloc_or_die(n_workers, sizeof(*thread_args));

	pthread_attr_t attr;
	pthread_attr_init(&attr);

	for (int i = 0; i < n_workers; i++) {
		struct worker_thread_args* args = &thread_args[i];
		worker_init(&workers[i], session, i);
		args->session = session;
		args->worker = &workers[i];
		pthread_create(&args->tid, &attr, worker_thread_start, args);
	}

	master_quit = 0;
	signal(SIGINT, sigint_handler);

	uint64_t t0 = ecs.initialization->global_stats.render_time_ms;
	struct timespec ts0;
	clock_gettime(CLOCK_MONOTONIC, &ts0);

	while (!master_quit) {
		sleep(1);
		for (int i = 0; i < n_workers; i++) {
			struct worker* worker = &workers[i];
			pthread_mutex_lock(&worker->mutex);
			worker->quit = master_quit;
			ecs_worker_stats_accumulate(&ecs.initialization->worker_stats, &worker->stats);
			memset(&worker->stats, 0, sizeof(worker->stats));
			pthread_mutex_unlock(&worker->mutex);

		}
		struct timespec ts;
		clock_gettime(CLOCK_MONOTONIC, &ts);
		uint64_t local_ms_elapsed = timespec_diff_ms(&ts0, &ts);
		ecs.initialization->global_stats.render_time_ms = t0 + local_ms_elapsed;

		stats_print(ecs.initialization);
	}

	printf("waiting for workers to quit...\n");
	for (int i = 0; i < n_workers; i++) {
		struct worker* worker = &workers[i];
		pthread_join(thread_args[i].tid, NULL);
		ecs_worker_stats_accumulate(&ecs.initialization->worker_stats, &worker->stats);
	}
	stats_print(ecs.initialization);

	printf("done\n");

	ecs_close(&ecs);
}
