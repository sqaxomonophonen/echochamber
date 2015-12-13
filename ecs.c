#include "ecs.h"

#include "sys.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>

#define ECS_BLOCK_TERMINATOR (0xdeadbeef)

static inline void _ecs_check_rw(struct ecs* ecs, size_t off, size_t sz)
{
	if (ecs->sz < (off+sz)) {
		fprintf(stderr, "invalid read of size %zd from offset %zd; ecs is only %zd bytes long\n", sz, off, ecs->sz);
		abort();
	}
}

static inline char* _ecs_cptr(struct ecs* ecs, size_t off)
{
	return ((char*)ecs->ptr) + off;
}

static inline uint32_t _ecs_read_u32(struct ecs* ecs, size_t off)
{
	_ecs_check_rw(ecs, off, sizeof(uint32_t));
	return *((uint32_t*)_ecs_cptr(ecs, off));
}

static inline float _ecs_read_f32(struct ecs* ecs, size_t off)
{
	_ecs_check_rw(ecs, off, sizeof(float));
	return *((float*)_ecs_cptr(ecs, off));
}

static inline void _ecs_write_u32(struct ecs* ecs, size_t off, uint32_t value)
{
	_ecs_check_rw(ecs, off, sizeof(uint32_t));
	*((uint32_t*)_ecs_cptr(ecs, off)) = value;
}

static inline void _ecs_write_f32(struct ecs* ecs, size_t off, float value)
{
	_ecs_check_rw(ecs, off, sizeof(float));
	*((float*)_ecs_cptr(ecs, off)) = value;
}

static inline void _assert_valid_magic(char* magic)
{
	if (strlen(magic) != 4) {
		fprintf(stderr, "magic '%s' is not 4 bytes long\n", magic);
		exit(EXIT_FAILURE);
	}
}

static inline int _ecs_match_magic(struct ecs* ecs, size_t off, char* magic)
{
	_assert_valid_magic(magic);
	char* ptr = (char*)ecs->ptr + off;
	return memcmp(ptr, magic, 4) == 0;
}


void* _ecs_find(struct ecs* ecs, char* magic, size_t* block_size)
{
	size_t off = 12;
	size_t sz = ecs->sz;

	void* result = NULL;

	while (off < sz) {
		if (sz < (off+12)) {
			fprintf(stderr, "invalid block size found at end of file (off=%zd vs sz=%zd)\n", off, sz);
			exit(EXIT_FAILURE);
		}

		int magic_match = _ecs_match_magic(ecs, off, magic);
		size_t bsz = _ecs_read_u32(ecs, off + 4);
		if (sz < (off+12+bsz)) {
			fprintf(stderr, "invalid block size found at end of file (sz=%zd)\n", bsz);
			exit(EXIT_FAILURE);
		}
		uint32_t terminator = _ecs_read_u32(ecs, off + 8 + bsz);
		if (terminator != ECS_BLOCK_TERMINATOR) {
			fprintf(stderr, "invalid block terminator %u found at offset %zd\n", terminator, off);
			exit(EXIT_FAILURE);
		}

		if (magic_match) {
			if (result != NULL) {
				fprintf(stderr, "found more than one block with magic '%s'\n", magic);
				exit(EXIT_FAILURE);
			}
			if (block_size != NULL) *block_size = bsz;
			result = (void*)((char*)ecs->ptr + off + 8);
		}

		off += bsz + 12;
	}

	return result;
}

void* _ecs_find_or_die(struct ecs* ecs, char* magic, size_t* block_size)
{
	void* ptr = _ecs_find(ecs, magic, block_size);
	if (ptr == NULL) {
		fprintf(stderr, "no block found with magic '%s'", magic);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

void _append_block(char* path, char* magic, void* block, size_t block_size)
{
	_assert_valid_magic(magic);

	int fd = open(path, O_RDWR | O_APPEND);
	if (fd == -1) {
		fprintf(stderr, "%s: %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	writen_or_die(fd, magic, 4);

	uint32_t sz = block_size;
	writen_or_die(fd, &sz, sizeof(sz));

	writen_or_die(fd, block, block_size);

	uint32_t terminator = ECS_BLOCK_TERMINATOR;
	writen_or_die(fd, &terminator, sizeof(terminator));

	close(fd);
}

#define MAX_ACC_PATH_LENGTH (1<<13)

static inline void _get_acc_path(char* acc_path, char* base, struct ecs* ecs, int mic_index, int emg_index)
{
	char* mic_name = ecs_get_microphone(ecs, mic_index)->name;
	char* emg_name = ecs_get_emission_group(ecs, emg_index)->name;
	snprintf(acc_path, MAX_ACC_PATH_LENGTH, "%s-acc.%s.%s.f32", base, mic_name, emg_name);
}


void ecs_open(struct ecs* ecs, char* path)
{
	memset(ecs, 0, sizeof(*ecs));

	int fd = open(path, O_RDWR);
	if (fd == -1) {
		fprintf(stderr, "%s: %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	struct stat st;
	if (fstat(fd, &st) == -1) {
		fprintf(stderr, "%s: %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}
	ecs->sz = st.st_size;

	if (ecs->sz < 12) {
		fprintf(stderr, "%s: not an ecs file (sz<12)\n", path);
		exit(EXIT_FAILURE);
	}

	ecs->ptr = mmap(NULL, ecs->sz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (ecs->ptr == MAP_FAILED) {
		fprintf(stderr, "%s: %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	close(fd);

	char* cs = ecs->ptr;
	if (cs[0] != 'E' || cs[1] != 'C' || cs[2] != 'S' || cs[3] != 'n') {
		fprintf(stderr, "%s: not an ecs file (invalid magic)\n", path);
		exit(EXIT_FAILURE);
	}

	uint32_t endian_magic = _ecs_read_u32(ecs, 4);
	if (endian_magic != 0xec511235) {
		fprintf(stderr, "%s: unhandled endianess (%u)\n", path, endian_magic);
		exit(EXIT_FAILURE);
	}

	uint32_t version = _ecs_read_u32(ecs, 8);
	if (version != 1) {
		fprintf(stderr, "%s: unhandled version (%u)\n", path, version);
		exit(EXIT_FAILURE);
	}

	/* look for mandatory blocks */
	ecs->_bqd0_block = _ecs_find_or_die(ecs, "BQD0", &ecs->_bqd0_sz);
	ecs->_emg0_block = _ecs_find_or_die(ecs, "EMG0", &ecs->_emg0_sz);
	ecs->_mic0_block = _ecs_find_or_die(ecs, "MIC0", &ecs->_mic0_sz);
	ecs->_mat0_block = _ecs_find_or_die(ecs, "MAT0", &ecs->_mat0_sz);
	ecs->_ply0_block = _ecs_find_or_die(ecs, "PLY0", &ecs->_ply0_sz);

	/* init block */
	size_t ini9_sz;
	ecs->initialization = _ecs_find(ecs, "INI9", &ini9_sz);
	if (ecs->initialization != NULL && ini9_sz != sizeof(*ecs->initialization)) {
		fprintf(stderr, "wrong block size of ini9, expected %zd, got %zd\n", sizeof(*ecs->initialization), ini9_sz);
		exit(EXIT_FAILURE);
	}

	/* if initialization is present, so must accumulator files; mmap them */
	if (ecs->initialization) {
		size_t sz = ecs->initialization->impulse_length_samples * sizeof(**ecs->accumulators);
		int n_mic = ecs_get_microphone_count(ecs);
		int n_emg = ecs_get_emission_group_count(ecs);

		int n = n_mic * n_emg;

		ecs->accumulator_sz = calloc_or_die(n, sizeof(*ecs->accumulator_sz));
		ecs->accumulators = calloc_or_die(n, sizeof(*ecs->accumulators));

		int i = 0;

		for (int mic_index = 0; mic_index < n_mic; mic_index++) {
			for (int emg_index = 0; emg_index < n_emg; emg_index++) {
				char acc_path[MAX_ACC_PATH_LENGTH];
				_get_acc_path(acc_path, path, ecs, mic_index, emg_index);

				int fd = open(acc_path, O_RDWR);
				if (fd == -1) {
					fprintf(stderr, "%s: %s\n", acc_path, strerror(errno));
					exit(EXIT_FAILURE);
				}

				struct stat st;
				if (fstat(fd, &st) == -1) {
					fprintf(stderr, "%s: %s\n", acc_path, strerror(errno));
					exit(EXIT_FAILURE);
				}

				if (sz != st.st_size) {
					fprintf(stderr, "%s: unexpected size %zd (expected %zd)\n", acc_path, st.st_size, sz);
					exit(EXIT_FAILURE);
				}

				void* ptr = ecs->accumulators[i] = mmap(NULL, sz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
				if (ptr == MAP_FAILED) {
					fprintf(stderr, "%s: %s\n", acc_path, strerror(errno));
					exit(EXIT_FAILURE);
				}

				close(fd);
			}
		}
	}
}

void ecs_close(struct ecs* ecs)
{
	munmap(ecs->ptr, ecs->sz);
}

static inline uint32_t _std_block_count(void* ptr)
{
	return *(uint32_t*)((char*)ptr);
}

void ecs_info(char* path)
{
	struct ecs ecs;
	ecs_open(&ecs, path);

	printf("file: %s:\n", path);

	if (ecs.initialization) {
		struct ecs_initialization* ini9 = ecs.initialization;
		printf("setup:\n");
		printf("  sample rate                 = %d hz\n", ini9->sample_rate);
		printf("  speed of sound              = %.3f u/s\n", ini9->speed_of_sound_units_per_second);
		printf("  impulse length              = %d samples\n", ini9->impulse_length_samples);
		printf("  attenuation threshold       = %.3e\n", ini9->attenuation_product_threshold);
	}

	printf("scene:\n");
	printf("  number of biquads           = %d\n", _std_block_count(ecs._bqd0_block));
	printf("  number of emission groups   = %d\n", _std_block_count(ecs._emg0_block));
	printf("  number of microphones       = %d\n", _std_block_count(ecs._mic0_block));
	printf("  number of materials         = %d\n", _std_block_count(ecs._mat0_block));
	// XXX count is implicit now...
	//printf("  number of polys             = %d\n", _std_block_count(ecs._ply0_block));

	ecs_close(&ecs);
}

void ecs_init(
	char* path,
	int force,
	int sample_rate,
	float speed_of_sound_units_per_second,
	int impulse_length_samples,
	float attenuation_product_threshold
) {
	struct ecs ecs;
	ecs_open(&ecs, path);
	int is_initialized = ecs.initialization != NULL;
	int n_emg = ecs_get_emission_group_count(&ecs);
	int n_mic = ecs_get_microphone_count(&ecs);

	if (!force && is_initialized) {
		fprintf(stderr, "%s is already initialized, use -f to force re-initialization\n", path);
		exit(EXIT_FAILURE);
	}

	if (!is_initialized) {
		/* append INI9 block */
		{
			struct ecs_initialization ini9;
			memset(&ini9, 0, sizeof(ini9));
			ini9.sample_rate = sample_rate;
			ini9.speed_of_sound_units_per_second = speed_of_sound_units_per_second;
			ini9.impulse_length_samples = impulse_length_samples;
			ini9.attenuation_product_threshold = attenuation_product_threshold;
			_append_block(path, "INI9", &ini9, sizeof(ini9));
		}
	}

	/* create empty accumulation buffers */
	{
		size_t sz = impulse_length_samples * sizeof(float);

		for (int mic_index = 0; mic_index < n_mic; mic_index++) {
			for (int emg_index = 0; emg_index < n_emg; emg_index++) {
				char acc_path[MAX_ACC_PATH_LENGTH];
				_get_acc_path(acc_path, path, &ecs, mic_index, emg_index);
				int fd = creat(acc_path, 0666);
				if (fd < 0) {
					fprintf(stderr, "%s: %s\n", acc_path, strerror(errno));
					exit(EXIT_FAILURE);
				}
				if (ftruncate(fd, sz) < 0) {
					fprintf(stderr, "%s: %s\n", acc_path, strerror(errno));
					exit(EXIT_FAILURE);
				}
				close(fd);
			}
		}
	}

	ecs_close(&ecs);
}

void ecs_mixdown(char* path)
{
	struct ecs ecs;
	ecs_open(&ecs, path);

	int n_mic = ecs_get_microphone_count(&ecs);
	int n_emg = ecs_get_emission_group_count(&ecs);
	int n_acc = n_mic * n_emg;

	if (!ecs.initialization) {
		fprintf(stderr, "%s: not initialized\n", path);
		exit(EXIT_FAILURE);
	}

	int n_samples = ecs.initialization->impulse_length_samples;

	float* mix = calloc_or_die(n_samples, sizeof(*mix));

	for (int ai = 0; ai < n_acc; ai++) {
		float* src = ecs.accumulators[ai];
		for (int s = 0; s < n_samples; s++) {
			mix[s] += src[s];
		}
	}

	float max = 0;
	int n_nan = 0;
	int n_inf = 0;
	for (int s = 0; s < n_samples; s++) {
		float v = mix[s];
		if (isnan(v)) {
			n_nan++;
			continue;
		}
		if (isinf(v)) {
			n_inf++;
			continue;
		}
		float a = fabsf(v);
		if (a > max) max = a;
	}

	printf("max value: %.5f\n", max);
	printf("nans: %d\n", n_nan);
	printf("infs: %d\n", n_inf);
	char out[1<<13];
	snprintf(out, 1<<13, "%s.wav", path);
	printf("writing %s...\n", out);

	int fd = creat(out, 0666);
	if (fd == -1) {
		fprintf(stderr, "%s: %s\n", out, strerror(errno));
		exit(EXIT_FAILURE);
	}

	/* write .wav */
	{
		#define W16(_x) do { \
			uint16_t x = _x; \
			writen_or_die(fd, &x, 2); \
		} while(0);

		#define W32(_x) do { \
			uint32_t x = _x; \
			writen_or_die(fd, &x, 4); \
		} while(0);

		writen_or_die(fd, "RIFF", 4);
		W32(36 + 2 * n_samples);

		writen_or_die(fd, "WAVEfmt ", 8);
		W32(16); // SubChunk1Size
		W16(1); // PCM format
		W16(1); // number of channels
		W32(ecs.initialization->sample_rate); // sample rate
		W32(ecs.initialization->sample_rate * 2); // byte rate
		W16(2); // block align (n_channels * bytes per sample)
		W16(16); // bits per sample

		writen_or_die(fd, "data ", 4);
		W32(2 * n_samples);
		for (int s = 0; s < n_samples; s++) {
			W16((mix[s] / max) * 32767.0f);
		}

		#undef W16
		#undef W32
	}

	close(fd);
}

static inline void* _ecs_get_std_obj(struct ecs* ecs, void* block, size_t sz, int i, const char* what)
{
	uint32_t n = _std_block_count(block);
	if (n < i) {
		fprintf(stderr, "invalid %s index %d/%d\n", what, i, n);
	}
	return (char*)block + 4 + sz*i;
}

int ecs_get_biquad_count(struct ecs* ecs)
{
	return _std_block_count(ecs->_bqd0_block);
}

struct ecs_biquad* ecs_get_biquad(struct ecs* ecs, int i)
{
	return _ecs_get_std_obj(ecs, ecs->_bqd0_block, sizeof(struct ecs_biquad), i, "biquad");
}

int ecs_get_emission_group_count(struct ecs* ecs)
{
	return _std_block_count(ecs->_emg0_block);
}

struct ecs_emission_group* ecs_get_emission_group(struct ecs* ecs, int i)
{
	return _ecs_get_std_obj(ecs, ecs->_emg0_block, sizeof(struct ecs_emission_group), i, "emission group");
}

int ecs_get_microphone_count(struct ecs* ecs)
{
	return _std_block_count(ecs->_mic0_block);
}

struct ecs_microphone* ecs_get_microphone(struct ecs* ecs, int i)
{
	return _ecs_get_std_obj(ecs, ecs->_mic0_block, sizeof(struct ecs_microphone), i, "microphone");
}

int ecs_get_material_count(struct ecs* ecs)
{
	return _std_block_count(ecs->_mat0_block);
}

struct ecs_material* ecs_get_material(struct ecs* ecs, int i)
{
	return _ecs_get_std_obj(ecs, ecs->_mat0_block, sizeof(struct ecs_material), i, "material");
}

int ecs_poly_iterator_next(struct ecs* ecs, ecs_poly_iterator* it, struct ecs_poly* poly)
{
	size_t sz = ecs->_ply0_sz;
	if (*it == sz) return 0;
	if (*it > sz) {
		fprintf(stderr, "invalid iterator\n");
		exit(EXIT_FAILURE);
	}

	void* base = ecs->_ply0_block;

	poly->material_index = *(uint32_t*)((char*)base + *it);
	poly->n_vertices = *(uint32_t*)((char*)base + *it + 4);
	poly->vertex_data = (float*)((char*)base + *it + 8);
	*it += (8 + poly->n_vertices * 3 * sizeof(*poly->vertex_data));

	return 1;
}
