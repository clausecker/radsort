#include <assert.h>
#include <stdalign.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __amd64__
#include <immintrin.h>
#endif
#ifdef __aarch64__
#include <arm_neon.h>
#endif
#ifdef __ALTIVEC__
#include <altivec.h>
#endif

#include "radixsort.h"

#define BLOCKSIZE 64
#define BLOCKALIGN 128

typedef uint32_t histogram32[4][256];
#if SIZE_MAX > UINT32_MAX
typedef uint64_t histogram64[4][256];
#endif

static void
histogram32_uint32pair(histogram32 h, const uint32pair a[restrict], size_t n)
{
	size_t i, j;

	memset(h, 0, sizeof (histogram32));

	for (i = 0; i < n; i++)
		for (j = 0; j < 4; j++)
			h[j][a[i][0] >> 8 * j & 0xff]++;
}

#if SIZE_MAX > UINT32_MAX
static void
histogram64_uint32pair(histogram64 h, const uint32pair a[restrict], size_t n)
{
	size_t i, j;

	memset(h, 0, sizeof (histogram64));

	for (i = 0; i < n; i++)
		for (j = 0; j < 4; j++)
			h[j][a[i][0] >> 8 * j & 0xff]++;
}
#endif

static void
flush_block(uint32pair out[restrict BLOCKSIZE], const uint32pair in[BLOCKSIZE])
{
	size_t k;

#ifdef __AVX__
	for (k = 0; k < BLOCKSIZE; k += sizeof(__m256)/sizeof(uint32pair)) {
		__m256 x;

		x = _mm256_load_ps((float *)(in + k));
		_mm256_stream_ps((float *)(out + k), x);
	}
#elif defined(__SSE__)
	for (k = 0; k < BLOCKSIZE; k += sizeof(__m128)/sizeof(uint32pair)) {
		__m128 x;

		x = _mm_load_ps((float *)(in + k));
		_mm_stream_ps((float *)(out + k), x);
	}
#elif defined(__aarch64__)
	for (k = 0; k < BLOCKSIZE; k += 32/sizeof(uint32pair))
		asm("ldp q0, q1, [%0]\n\tstnp q0, q1, [%1]" :: "r"(in + k), "r"(out + k) : "q0", "q1", "memory");
#elif defined(__ALTIVEC__)
	for (k = 0; k < BLOCKSIZE; k += sizeof(vector unsigned int)/sizeof(uint32pair)) {
		vector unsigned int x;

		x = vec_ld(0, (unsigned int *)(in + k));
		vec_stl(x, 0, (unsigned int *)(out + k));
	}
#else
	for (k = 0; k < BLOCKSIZE; k++) {
		out[k][0] = in[k][0];
		out[k][1] = in[k][1];
	}
#endif
}

static void
radixsort32_pairs_round(uint32pair *restrict dest, const uint32pair *src, histogram32 h, size_t n, int round)
{
	size_t i, j, k;
	uint32_t *buckets, start, next, starts[256];
	alignas(BLOCKALIGN) uint32pair temp[256][BLOCKSIZE];
	unsigned fill[256], b;

	assert((uintptr_t)dest % BLOCKALIGN == 0);

	/* compute bucket start offsets */
	buckets = h[round];
	start = 0;
	for (i = 0; i < 256; i++) {
		next = start + buckets[i];
		starts[i] = start;
		start = next;
	}

	/* compute initial bucket fills */
	for (i = 0; i < 256; i++) {
		fill[i] = starts[i] % (BLOCKALIGN / sizeof(uint32pair));
		assert(starts[i] >= fill[i]);
		starts[i] -= fill[i];
		assert((uintptr_t)(dest + starts[i]) % BLOCKALIGN == 0);
		assert(fill[i] < BLOCKSIZE);
	}

	/* sort into buckets */
	for (i = 0; i < n; i++) {
		b = src[i][0] >> 8 * round & 0xff;
		j = fill[b]++;
		temp[b][j][0] = src[i][0];
		temp[b][j][1] = src[i][1];
		if (j == BLOCKSIZE - 1) {
			flush_block(dest + starts[b], temp[b]);
			starts[b] += BLOCKSIZE;
			fill[b] = 0;
		}
	}

	/* flush final items */
	for (b = 0; b < 256; b++) {
		size_t j0;

		j0 = buckets[b] >= fill[b] ? 0 : fill[b] - buckets[b];

		for (j = j0; j < fill[b]; j++) {
			dest[starts[b] + j][0] = temp[b][j][0];
			dest[starts[b] + j][1] = temp[b][j][1];
		}
	}
}

#if SIZE_MAX > UINT32_MAX
static void
radixsort64_pairs_round(uint32pair *restrict dest, const uint32pair *src, histogram64 h, size_t n, int round)
{
	size_t i, j, k;
	uint64_t *buckets, start, next, starts[256];
	alignas(BLOCKALIGN) uint32pair temp[256][BLOCKSIZE];
	unsigned fill[256], b;

	assert((uintptr_t)dest % BLOCKALIGN == 0);

	/* compute bucket start offsets */
	buckets = h[round];
	start = 0;
	for (i = 0; i < 256; i++) {
		next = start + buckets[i];
		starts[i] = start;
		start = next;
	}

	/* compute initial bucket fills */
	for (i = 0; i < 256; i++) {
		fill[i] = starts[i] % (BLOCKALIGN / sizeof(uint32pair));
		assert(starts[i] >= fill[i]);
		starts[i] -= fill[i];
		assert((uintptr_t)(dest + starts[i]) % BLOCKALIGN == 0);
		assert(fill[i] < BLOCKSIZE);
	}

	/* sort into buckets */
	for (i = 0; i < n; i++) {
		b = src[i][0] >> 8 * round & 0xff;
		j = fill[b]++;
		temp[b][j][0] = src[i][0];
		temp[b][j][1] = src[i][1];
		if (j == BLOCKSIZE - 1) {
			flush_block(dest + starts[b], temp[b]);
			starts[b] += BLOCKSIZE;
			fill[b] = 0;
		}
	}

	/* flush final items */
	for (b = 0; b < 256; b++) {
		size_t j0;

		j0 = buckets[b] >= fill[b] ? 0 : fill[b] - buckets[b];

		for (j = j0; j < fill[b]; j++) {
			dest[starts[b] + j][0] = temp[b][j][0];
			dest[starts[b] + j][1] = temp[b][j][1];
		}
	}
}
#endif

static void
radixsort_pairs(uint32pair *src, size_t n)
{
	size_t bytes;
	uint32pair *tmp;

	bytes = n * sizeof *tmp;
	tmp = aligned_alloc(BLOCKALIGN, bytes + BLOCKALIGN - 1 & ~(size_t)(BLOCKALIGN - 1));
	if (tmp == NULL) {
		perror("radixsort_pairs");
		exit(EXIT_FAILURE);
	}

#if SIZE_MAX > UINT32_MAX
	if (n > UINT32_MAX) {
		histogram64 h;

		histogram64_uint32pair(h, src, n);
		radixsort64_pairs_round(tmp, src, h, n, 0);
		radixsort64_pairs_round(src, tmp, h, n, 1);
		radixsort64_pairs_round(tmp, src, h, n, 2);
		radixsort64_pairs_round(src, tmp, h, n, 3);
	} else
#endif
	{
		histogram32 h;

		histogram32_uint32pair(h, src, n);
		radixsort32_pairs_round(tmp, src, h, n, 0);
		radixsort32_pairs_round(src, tmp, h, n, 1);
		radixsort32_pairs_round(tmp, src, h, n, 2);
		radixsort32_pairs_round(src, tmp, h, n, 3);
	}

	free(tmp);
}

const struct radixsort_impl swc_impl = {
	"swc",
	radixsort_pairs,
	always_supported,
};
