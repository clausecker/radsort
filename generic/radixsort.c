#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "radixsort.h"

/*
 * Modular inverse using Hurchalla's algorithm.
 * Assumes x is odd, computes a number such that
 * modinv(a) * a == 1.
 */
static uintptr_t
modinv(uintptr_t a) {
	uintptr_t x, y;

	assert(a & 1);

	x = 3 * a ^ 2;
	y = 1 - a * x;
	x *= (1 + y);
	y *= y;
	x *= (1 + y);
	y *= y;
	x *= (1 + y);
	if (sizeof(a) > 32) {
		y *= y;
		x *= (1 + y);
	}

	assert(x * a == 1);

	return (x);
}

/*
 * Compute the total memory consumed by all data structures of a sorter.
 */
static size_t
mem_use(size_t n, size_t esize, unsigned n_thread, unsigned bshift)
{
	struct rs_sorter srt;
	size_t bytes = 0, n_scratch, n_blocks;

	n_scratch = 2 * RS_RADIX * n_thread + 1;;
	n_blocks = (n >> bshift) + n_scratch;
	bytes += sizeof(srt);
	bytes += n_scratch * (1 << bshift) * esize; /* scratch */
	bytes += 2 * n_blocks * sizeof(*srt.perm); /* perm, perm2 */
	bytes += n_blocks * sizeof(*srt.assignments);
	bytes += n_thread * RS_RADIX * sizeof(*srt.partials);
	bytes += n_thread * sizeof(struct rs_thread);

	return (bytes);
}

/*
 * Add a partial block to the end of the partials array
 */
static void
add_partial(struct rs_sorter *srt, unsigned index, unsigned length) {
	assert(srt->n_partial < RS_RADIX * srt->n_thread);
	assert(srt->n_partial == 0 || srt->partials[srt->n_partial - 1].index < index);
	srt->partials[srt->n_partial].index = index;
	srt->partials[srt->n_partial++].length = length;
}

/*
 * Find the block size that minimises the amount of memory allocated.
 */
static unsigned
find_bshift(size_t n, size_t esize, unsigned n_thread)
{
	unsigned bshift;
	size_t cursize, minsize = SIZE_MAX;

	bshift = 6;
	for (bshift = 6; bshift < 32; bshift++) {
		cursize = mem_use(n, esize, n_thread, bshift);
		if (cursize > minsize)
			return (bshift - 1);

		minsize = cursize;
	}

	return (31);
}

static void
rs_step_generic(struct rs_sorter *srt, struct rs_thread *thr, size_t keyb)
{
	unsigned i_in, i_out, i_partial;
	unsigned char **buckets;

	buckets = (unsigned char **)thr->buckets;
	i_partial = thr->i_partial;
	for (i_in = i_out = thr->i_start + RS_RADIX; i_in < thr->fill; i_in++) {
		unsigned j, len;
		unsigned char *in_blk;

		in_blk = rs_blockat(srt, srt->perm2[i_in]);
		len = srt->bsize;
		if (srt->partials[i_partial].index == i_in)
			len = srt->partials[i_partial++].length;

		for (j = 0; j < len; j++) {
			unsigned bucket;
			unsigned char *bucketp;

			bucket = in_blk[j * srt->esize + keyb];
			bucketp = buckets[bucket];
			memcpy(bucketp, in_blk + j * srt->esize, srt->esize);
			bucketp += srt->esize;

			/* bucket full? */
			buckets[bucket] = bucketp;
			if (((uintptr_t)bucketp & srt->blockmask) == 0) {
				unsigned char *out_blk;

				out_blk = rs_blockat(srt, srt->perm2[i_out]);
				buckets[bucket] = out_blk;
				srt->assignments[i_out] = bucket;
				thr->counts[bucket]++;
				i_out++;
			}
		}
	}

	thr->fill = i_out;

	
}

void
rs_free(struct rs_sorter *srt)
{
	free(srt->scratch);
	free(srt->perm);
	free(srt->perm2);
	free(srt->assignments);
	free(srt->partials);
	free(srt->threads);
	free(srt);
}

struct rs_sorter *
rs_new(rs_key *array, size_t nmemb, size_t size, rs_step_func *step)
{
	struct rs_sorter *srt;
	size_t n_head, n_tail, n_blocks_main;
	size_t i_perm, i_scratch, i_main;
	unsigned eshift;

	if (step == NULL)
		step = rs_step_generic;

	srt = calloc(sizeof *srt, 1);
	if (srt == NULL)
		return (NULL);

	srt->step = step;

#ifdef _OPENMP
	srt->n_thread = omp_get_max_threads();
	assert (0 < srt->n_thread);
#else
	srt->n_thread = 1;
#endif
	srt->n = nmemb;
	srt->esize = size;
	srt->bshift = find_bshift(srt->n, size, srt->n_thread);
	srt->bsize = 1U << srt->bshift;
	eshift = __builtin_ctz(srt->esize);

	/* iff (ptr & blockmask) == 0, ptr points to the beginning of a block */
	srt->blockmask = (uintptr_t)(srt->bsize - 1) << __builtin_ctz(srt->esize);

	/* number of elements by which main is offset from array */
	n_head = (-(uintptr_t)array * modinv(srt->esize >> eshift) & srt->blockmask) >> eshift;
	assert(n_head <= srt->n);
	n_blocks_main = srt->n - n_head / srt->bsize;
	n_tail = srt->n - n_head % srt->bsize;

	srt->array = array;
	srt->main = (unsigned char *)array + srt->esize * n_head;
	assert(((uintptr_t)srt->main & srt->blockmask) == 0);

	/* allocate storage for the various arrays */
	srt->n_scratch = 2 * RS_RADIX * srt->n_thread + 1;
	srt->scratch = calloc(srt->esize, srt->n_scratch * srt->bsize);
	srt->n_block = n_blocks_main + srt->n_scratch;
	srt->perm = calloc(srt->n_block, sizeof *srt->perm);
	srt->perm2 = calloc(srt->n_block, sizeof *srt->perm2);
	srt->assignments = calloc(srt->n_block, sizeof *srt->assignments);
	srt->partials = calloc(RS_RADIX * srt->n_thread, sizeof *srt->partials);
	srt->threads = calloc(srt->n_thread, sizeof *srt->threads);

	if (srt->scratch == NULL || srt->perm == NULL || srt->perm2 == NULL ||
	    srt->assignments == NULL || srt->partials == NULL || srt->threads == NULL) {
		rs_free(srt);
		return (NULL);
	}

	i_perm = 0;
	i_scratch = 0;

	/* set up head block */
	memcpy(rs_blockat(srt, i_scratch), srt->array, n_head * srt->esize);
	srt->perm[i_perm] = i_scratch++;
	add_partial(srt, i_perm++, n_head);

	/* set up main permutation */
	for (i_main = 0; i_main < n_blocks_main; i_perm++, i_main++)
		srt->perm[i_perm] = srt->n_scratch + i_main;

	/* set up tail block */
	memcpy(rs_blockat(srt, i_scratch), rs_blockat(srt, i_main), n_tail * srt->esize);
	srt->perm[i_perm] = i_scratch;
	add_partial(srt, i_perm++, n_tail);

	srt->fill = i_perm;

	/* initialise remaining permutation */
	for (; i_scratch < srt->n_scratch; i_perm++, i_scratch++)
		srt->perm[i_perm] = i_scratch;

	assert(i_perm == srt->n_block);

	return (srt);
}

/*
 * Divide an array of n elements fairly among n_thread threads
 * such that each thread gets roughly the same number of elements.
 * Return the number of elements allocated to thread 0 to t
 * inclusive (i.e. the boundary at which the elements to thread t
 * end).
 */
static size_t
thread_quota(size_t n, size_t n_thread, size_t t)
{
	size_t quot, rem;

	quot = n / n_thread;
	rem = n % n_thread;
	t++;

	return (t * quot + (rem <= t ? rem : t));
}

/*
 * Take the blocks in srt->perm and group them into
 * srt->n_thread groups, each with a headstart of RS_RADIX unused
 * blocks, such that each group has roughly the same amount of
 * elements to work on.  Write the results to srt->perm2 and the
 * per-thread structures in srt->thread.  Take the number of
 * allocated blocks from block_count.
 */
static void
distribute_work(struct rs_sorter *srt)
{
	size_t i_perm, i_perm2, i_partial, items;
	unsigned t;

	assert(RS_RADIX * srt->n_thread <= srt->n_block - srt->fill);

	items = 0;
	i_perm = 0;
	i_perm2 = 0;
	i_partial = 0;

	for (t = 0; t < srt->n_thread; t++) {
		struct rs_thread *thr;
		size_t quota;

		thr = srt->threads + t;
		thr->i_start = i_perm2;
		thr->i_partial = i_partial;
		quota = thread_quota(srt->n, srt->n_thread, t);

		/* assign head start */
		memcpy(srt->perm2 + i_perm2, srt->perm + srt->fill + RS_RADIX * t, RS_RADIX * sizeof *srt->perm2);
		i_perm2 += RS_RADIX;

		if (t == srt->n_thread - 1) {
			/* assign all remaining blocks to final thread */
			while (i_perm < srt->fill) {
				assert(i_partial < srt->n_partial);
				if (srt->partials[i_partial].index == i_perm) {
					srt->partials[i_partial].index = i_perm2;
					items += srt->partials[i_partial++].length;
				} else
					items += srt->bsize;

				srt->perm2[i_perm2++] = srt->perm[i_perm++];
			}
		} else {
			/* assign blocks until quota is met */
			while (items < quota) {
				assert(i_partial < srt->n_partial);
				if (srt->partials[i_partial].index == i_perm) {
					srt->partials[i_partial].index = i_perm2;
					items += srt->partials[i_partial++].length;
				} else
					items += srt->bsize;

				srt->perm2[i_perm2++] = srt->perm[i_perm++];
			}
		}

		thr->fill = i_perm2;

		assert(i_perm <= srt->fill);
	}

	/* assign remaining blocks to srt->perm2 */
	assert(i_partial == srt->n_partial);
	assert(i_perm == srt->fill);
	assert(i_perm2 == srt->fill + RS_RADIX * t);
	memcpy(srt->perm2 + i_perm2, srt->perm + i_perm2, (srt->n_block - i_perm2) * sizeof *srt->perm2);
}

/*
 * Following the sorting procedure, sort the resulting blocks
 * according to bucket.  The blocks are placed in one chunk
 * at the beginning of srt->perm, the unused blocks thereafter.
 * Take note of the partial blocks encountered in srt->partials
 * and record their number in srt->n_partial.  Write the total
 * number of allocated blocks to srt->fill.
 */
static void
sort_blocks(struct rs_sorter *srt)
{
	size_t i, j, t;
	size_t sum, block_count, free_start;
	unsigned starts[RS_RADIX];

	memset(srt->partials, -1, RS_RADIX * srt->n_thread * sizeof *srt->partials);

	/* sum up the bucket sizes of the thread-arrays */
	block_count = 0;
	for (i = 0; i < RS_RADIX; i++) {
		sum = 0;
		starts[i] = block_count;
		for (j = 0; j < srt->n_thread; j++)
			sum += srt->threads[j].counts[i];

		block_count += sum;
	}

	free_start = block_count;

	/* iterate over srt->perm2 and bring the blocks into order. */
	for (t = 0; t < srt->n_thread; t++) {
		struct rs_thread *thr;
		size_t k, bucket, free_count;

		thr = srt->threads + t;

		/* allocated blocks in the chunk of thread i */
		for (j = thr->i_start; j < thr->fill; j++) {
			const char *block;

			bucket = srt->assignments[j];
			k = starts[bucket]++;
			srt->perm[k] = srt->perm2[j];

			/* partial block? */
			block = rs_blockat(srt, srt->perm2[j]);
			if (block <= (const char *)thr->buckets[bucket]
			    && (const char *)thr->buckets[bucket] < block + srt->bsize * srt->esize) {
				srt->partials[bucket * srt->n_thread + t].index = k;
				srt->partials[bucket * srt->n_thread + t].length =
				    ((const char *)thr->buckets[bucket] - block) / srt->esize;
			}
		}

		/* unallocated blocks at the end of the chunk */
		free_count = (t < srt->n_thread - 1 ? srt->threads[t + 1].i_start : srt->n_block) - thr->fill;
		assert(free_start + free_count <= srt->n_block);
		memcpy(srt->perm + free_start, srt->perm2 + thr->fill, free_count * sizeof *srt->perm);
		free_start += free_count;
	}

	/* gather the partial blocks at the beginning of the partials array */
	srt->n_partial = 0;
	for (i = 0; i < RS_RADIX * srt->n_thread; i++)
		if (srt->partials[i].index != (unsigned)-1)
			srt->partials[srt->n_partial++] = srt->partials[i];

	assert(free_start == srt->n_block);
	srt->fill = block_count;
}

void
rs_step(struct rs_sorter *srt, size_t keyb)
{
	unsigned t;

	distribute_work(srt);

#	pragma omp parallel for schedule(static,1)
	for (t = 0; t < srt->n_thread; t++)
		srt->step(srt, srt->threads + t, keyb);

	sort_blocks(srt);
}
