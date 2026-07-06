#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "radixsort.h"

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
radixsort32_pairs_round(uint32pair *restrict dest, const uint32pair *src, histogram32 h, size_t n, int round)
{
	size_t i, j;
	uint32_t *buckets, start, next;

	/* compute bucket start offsets */
	buckets = h[round];
	start = 0;
	for (i = 0; i < 256; i++) {
		next = start + buckets[i];
		buckets[i] = start;
		start = next;
	}

	/* sort into buckets */
	for (i = 0; i < n; i++) {
		j = buckets[src[i][0] >> 8 * round & 0xff]++;
		dest[j][0] = src[i][0];
		dest[j][1] = src[i][1];

#if __amd64__ || __i386__
		__builtin_prefetch(dest[j + 2], 1, 2);
#else
		__builtin_prefetch(dest[j + 1], 1, 3);
#endif
	}
}

#if SIZE_MAX > UINT32_MAX
static void
radixsort64_pairs_round(uint32pair *restrict dest, const uint32pair *src, histogram64 h, size_t n, int round)
{
	size_t i, j;
	uint64_t *buckets, start, next;

	/* compute bucket start offsets */
	buckets = h[round];
	start = 0;
	for (i = 0; i < 256; i++) {
		next = start + buckets[i];
		buckets[i] = start;
		start = next;
	}

	/* sort into buckets */
	for (i = 0; i < n; i++) {
		j = buckets[src[i][0] >> 8 * round & 0xff]++;
		dest[j][0] = src[i][0];
		dest[j][1] = src[i][1];

#if __amd64__ || __i386__
		__builtin_prefetch(dest[j + 2], 1, 2);
#else
		__builtin_prefetch(dest[j + 1], 1, 3);
#endif
	}
}
#endif

static void
radixsort_pairs(uint32pair *src, size_t n)
{
	uint32pair *tmp;

	tmp = malloc(n * sizeof *tmp);
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

int
always_supported(void)
{
	return (1);
}

const struct radixsort_impl generic_impl = {
	"generic",
	radixsort_pairs,
	always_supported,
};
