#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "radixsort.h"

#define RSHIFT CHAR_BIT		/* number of bits per sort step */
#define RADIX (1<<RSHIFT)	/* number of buckets */
#define NSCRATCH (2*RADIX)	/* number of scratch blocks */
#define BLOCKSIZE 512		/* elements per block */
#define NPERM(srt) ((srt)->n_scratch + (srt)->n / BLOCKSIZE) /* number of blocks in total */

struct block {
	uint32pair *start, *end;
};

struct partial {
	unsigned	index;		/* block index */
	unsigned	length;		/* number of items in the partial block */
};

/*
 * Per-thread structure.  Unused and uninitialised outside of sorter_step().
 */
struct thread {
	size_t		 i_start;	/* first block to be processed by this thread */
	size_t		 i_partial;	/* first partial bucket in the input */
	size_t		 fill;		/* end of allocated blocks */
	struct block	 buckets[RADIX]; /* buckets (output) */
	unsigned	 counts[RADIX];	/* number of blocks per output bucket */
};

/*
 * The sorter structure represents the state of the sorting process.
 * Between sorting iterations, the following invariants hold:
 *
 *  - n holds the number of pairs in pairs.
 *  - perm[0] to perm[NSCRATCH - 1] and perm[fill] to perm[NPERM(N)]
 *    refer to unallocated blocks
 *  - perm[NSCRATCH] to perm[fill - 1] refer to allocated blocks
 *    in the of the current sorting progress.
 *  - counts[b] holds the number of allocated blocks of bucket[b] for each b
 *  - For each bucket b, there may be one partial allocated block.  If such
 *    a block exists, counts[b].start points to the first free element in the
 *    block and counts[b].end to the end of the block.  Otherwise both are
 *    null pointers.
 *  - perm and perm2 as well as buckets and buckets2 point to distinct arrays
 *    of the same length that can be swapped by the sorting step.
 */
struct sorter {
	uint32pair	 *pairs;	/* input array */
	size_t		  n;		/* number of input elements */
	int		  n_thread;	/* number of threads in use */
	size_t		  n_scratch;	/* number of scratch blocks */
	size_t		  n_partial;	/* number of partial blocks */
	size_t		  fill;		/* number of allocated blocks in perm */
	unsigned	 *perm;		/* block permutation */
	unsigned	 *perm2;	/* block permutation (auxillary) */
	unsigned char	 *assignments;	/* assignments of logical blocks to buckets */
	struct partial	 *partials;	/* partial blocks (input) */
	struct thread	 *threads;	/* per-thread structures */
	uint32pair	(*scratch)[BLOCKSIZE];
};

/*
 * Retrieve the start address at index i from srt.
 */
static inline uint32pair *
blockat(struct sorter *srt, size_t i)
{
	return (i < srt->n_scratch ? srt->scratch[i] : srt->pairs + (i - srt->n_scratch) * BLOCKSIZE);
}

/*
 * Get the bucket a given uint32pair belongs into according to shift.
 */
static inline unsigned
bucketof(uint32pair p, int shift)
{
	return ((unsigned)(p[0] >> shift & RADIX - 1));
}

static void
sorter_free(struct sorter *srt)
{
	free(srt->perm);
	free(srt->perm2);
	free(srt->assignments);
	free(srt->partials);
	free(srt->threads);
	free(srt->scratch);
	free(srt);
}

/*
 * Create a new sorter from the given pairs and length n.  Return NULL on failure.
 * The newly created sorter represents an identity permutation of the input, with
 * all blocks being in bucket 0.
 */
static struct sorter *
sorter_initialize(uint32pair *pairs, size_t n)
{
	struct sorter *srt;
	size_t i_scratch, i_perm, n_tail;

	srt = calloc(sizeof *srt, 1);
	if (srt == NULL)
		return (NULL);

	srt->pairs = pairs;
	srt->n = n;
#ifdef _OPENMP
	srt->n_thread = omp_get_max_threads();
	assert(0 < srt->n_thread);
#else
	srt->n_thread = 1;
#endif
	srt->n_scratch = NSCRATCH * srt->n_thread;
	srt->perm = calloc(sizeof *srt->perm, NPERM(srt));
	if (srt->perm == NULL)
		goto fail;
	srt->perm2 = calloc(sizeof *srt->perm2, NPERM(srt));
	if (srt->perm2 == NULL)
		goto fail;
	srt->assignments = calloc(sizeof *srt->assignments, NPERM(srt));
	if (srt->assignments == NULL)
		goto fail;
	srt->partials = calloc(sizeof *srt->partials, srt->n_thread * RADIX);
	if (srt->partials == NULL)
		goto fail;
	srt->threads = calloc(sizeof *srt->threads, srt->n_thread);
	if (srt->threads == NULL)
		goto fail;
	srt->scratch = calloc(sizeof *srt->scratch, srt->n_scratch);
	if (srt->scratch == NULL)
		goto fail;

	/* set up identity permutation */
	for (i_perm = 0; i_perm * BLOCKSIZE < n; i_perm++)
		srt->perm[i_perm] = srt->n_scratch + i_perm;

	srt->fill = i_perm;

	/* copy partial input tail to the scratchpad */
	n_tail = n % BLOCKSIZE;
	if (n_tail != 0) {
		i_perm--;
		memcpy(srt->scratch[0], pairs + i_perm * BLOCKSIZE, n_tail * sizeof *pairs);
		srt->partials[0].index = i_perm;
		srt->partials[0].length = n_tail;
		srt->n_partial++;
	} else {
		/* ensure there is at least one partial block */
		srt->partials[0].index = i_perm - 1;
		srt->partials[0].length = BLOCKSIZE;
		srt->n_partial++;
	}

	/* store scratchpad blocks in permutation */
	for (i_scratch = 0; i_scratch < srt->n_scratch; i_scratch++)
		srt->perm[i_perm++] = i_scratch;

	assert(i_perm == NPERM(srt));

	return (srt);

fail:	sorter_free(srt);
	return (NULL);
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
 * srt->n_thread groups, each with a headstart of RADIX unused
 * blocks, such that each group has roughly the same amount of
 * elements to work on.  Write the results to srt->perm2 and the
 * per-thread structures in srt->thread.  Take the number of
 * allocated blocks from block_count.
 */
static void
distribute_work(struct sorter *srt)
{
	size_t i_perm, i_perm2, i_partial, items;
	int t;

	assert(RADIX * srt->n_thread <= NPERM(srt) - srt->fill);

	items = 0;
	i_perm = 0;
	i_perm2 = 0;
	i_partial = 0;

	for (t = 0; t < srt->n_thread; t++) {
		struct thread *thr;
		size_t quota;

		thr = srt->threads + t;
		thr->i_start = i_perm2;
		thr->i_partial = i_partial;
		quota = thread_quota(srt->n, srt->n_thread, t);

		/* assign head start */
		memcpy(srt->perm2 + i_perm2, srt->perm + srt->fill + RADIX * t, RADIX * sizeof *srt->perm2);
		i_perm2 += RADIX;

		if (t  == srt->n_thread - 1) {
			/* assign all remaining blocks to final thread */
			while (i_perm < srt->fill) {
				assert(i_partial < srt->n_partial);
				if (srt->partials[i_partial].index == i_perm) {
					srt->partials[i_partial].index = i_perm2;
					items += srt->partials[i_partial++].length;
				} else
					items += BLOCKSIZE;

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
					items += BLOCKSIZE;

				srt->perm2[i_perm2++] = srt->perm[i_perm++];
			}
		}

		thr->fill = i_perm2;

		assert(i_perm <= srt->fill);
	}

	/* assign remaining blocks to srt->perm2 */
	assert(i_partial == srt->n_partial);
	assert(i_perm == srt->fill);
	assert(i_perm2 == srt->fill + RADIX * t);
	memcpy(srt->perm2 + i_perm2, srt->perm + i_perm2, (NPERM(srt) - i_perm2) * sizeof *srt->perm2);
}

/*
 * Perform the work of thread t in sorting at the given shift.
 */
static void
sort_step_thread(struct sorter *srt, int shift, int t)
{
	struct thread *thr;
	size_t i_in, i_out, i_partial;

	thr = srt->threads + t;

	/* allocate initial buckets */
	for (i_out = thr->i_start; i_out < thr->i_start + RADIX; i_out++) {
		unsigned bucket;
		uint32pair *out_blk;

		bucket = i_out - thr->i_start;
		out_blk = blockat(srt, srt->perm2[i_out]);
		thr->buckets[bucket].start = out_blk;
		thr->buckets[bucket].end = out_blk + BLOCKSIZE;
		srt->assignments[i_out] = bucket;
		thr->counts[bucket] = 1;
	}

	/* invariant: i_out <= i_in */
	i_partial = thr->i_partial;

	for (i_in = i_out; i_in < thr->fill; i_in++) {
		size_t j, len;
		uint32pair *in_blk;

		/* grab the next source block */
		assert(srt->perm[i_in] < NPERM(srt));
		in_blk = blockat(srt, srt->perm2[i_in]);
		len = BLOCKSIZE;
		if (srt->partials[i_partial].index == i_in)
			len = srt->partials[i_partial++].length;

		for (j = 0; j < len; j++) {
			unsigned bucket;
			uint32pair *bucketp;

			bucket = bucketof(in_blk[j], shift);
			bucketp = thr->buckets[bucket].start; /* avoid reload of pointer after memcpy */
			memcpy(bucketp++, in_blk + j, sizeof in_blk[j]);

			/* bucket full? */
			thr->buckets[bucket].start = bucketp;
			if (bucketp == thr->buckets[bucket].end) {
				uint32pair *out_blk;

				assert(i_out <= i_in);
				assert(srt->perm2[i_out] < NPERM(srt));
				out_blk = blockat(srt, srt->perm2[i_out]);
				thr->buckets[bucket].start = out_blk;
				thr->buckets[bucket].end = out_blk + BLOCKSIZE;
				srt->assignments[i_out] = bucket;
				thr->counts[bucket]++;
				i_out++;
			}
		}
	}

	thr->fill = i_out;
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
sort_blocks(struct sorter *srt)
{
	size_t i, j, t;
	size_t sum, block_count, free_start;
	unsigned starts[RADIX];

	memset(srt->partials, -1, RADIX * srt->n_thread * sizeof *srt->partials);

	/* sum up the bucket sizes of the thread-arrays */
	block_count = 0;
	for (i = 0; i < RADIX; i++) {
		sum = 0;
		starts[i] = block_count;
		for (j = 0; j < srt->n_thread; j++)
			sum += srt->threads[j].counts[i];

		block_count += sum;
	}

	free_start = block_count;

	/* iterate over srt->perm2 and bring the blocks into order. */
	for (t = 0; t < srt->n_thread; t++) {
		struct thread *thr;
		size_t k, bucket, free_count;

		thr = srt->threads + t;

		/* allocated blocks in the chunk of thread i */
		for (j = thr->i_start; j < thr->fill; j++) {
			bucket = srt->assignments[j];
			k = starts[bucket]++;
			srt->perm[k] = srt->perm2[j];

			/* partial block? */
			if (blockat(srt, srt->perm2[j]) + BLOCKSIZE
			    == thr->buckets[bucket].end) {
				srt->partials[bucket * srt->n_thread + t].index = k;
				srt->partials[bucket * srt->n_thread + t].length = BLOCKSIZE -
				    (thr->buckets[bucket].end - thr->buckets[bucket].start);
			}
		}

		/* unallocated blocks at the end of the chunk */
		free_count = (t < srt->n_thread - 1 ? srt->threads[t + 1].i_start : NPERM(srt)) - thr->fill;
		assert(free_start + free_count <= NPERM(srt));
		memcpy(srt->perm + free_start, srt->perm2 + thr->fill, free_count * sizeof *srt->perm);
		free_start += free_count;
	}

	/* gather the partial blocks at the beginning of the partials array */
	srt->n_partial = 0;
	for (i = 0; i < RADIX * srt->n_thread; i++)
		if (srt->partials[i].index != -1)
			srt->partials[srt->n_partial++] = srt->partials[i];

	assert(free_start == NPERM(srt));
	srt->fill = block_count;
}

/*
 * Perform one iteration of radix sort at the given shift.
 */
static void
sort_step(struct sorter *srt, int shift)
{
	int t;

	distribute_work(srt);

#	pragma omp parallel for schedule(static,1)
	for (t = 0; t < srt->n_thread; t++)
		sort_step_thread(srt, shift, t);

	sort_blocks(srt);
}

/*
 * Take a sorter for which the data structure invariant holds and move its
 * into place in srt->pairs.  This is done using a push-pull approach:
 * space is made in srt->pairs by pushing the current block away.  The
 * appropriate block is then pulled into the now free spot, filling any
 * gaps from partial blocks.
 */
static void
compact(struct sorter *srt)
{
	unsigned *perm, *pinv;

	size_t start, blocklen, i_partial;

	/* convention: i for logical blocks, j for physical blocks */
	size_t i, i_in, i_out, i_free;
	size_t j, j_in, j_out, j_free;

	/* compute the inverse permutation of srt->perm */
	perm = srt->perm;
	pinv = srt->perm2;

	for (i = 0; i < NPERM(srt); i++) {
		j = perm[i];
		pinv[j] = i;
	}

	start = 0;			/* start of next block in srt->pairs */
	i_partial = 0;			/* next partial block in srt->partials */
	i_free = srt->fill;		/* some free block */
	j_free = perm[i_free];

	/* permute the blocks overlaying srt->pairs */
	for (j_out = srt->n_scratch; j_out < NPERM(srt); j_out++) {
		i_in = j_out - srt->n_scratch;
		i_out = pinv[j_out];
		j_in = perm[i_in];

		/* j_out not the right block, but allocated? */
		if (i_out > i_in && i_out < srt->fill) {
			/* move block at j_out to j_free */
			assert(i_free >= srt->fill);
			memcpy(blockat(srt, j_free), blockat(srt, j_out), sizeof *srt->scratch);
			perm[i_free] = j_out;
			perm[i_out] = j_free;
			pinv[j_free] = i_out;
			pinv[j_out] = i_free;

			/* swap i_free and i_out */
			i_free = i_out;
			i_out = pinv[j_out];
		}

		/* move block at j_in to start, pretend it got moved to j_out */
		blocklen = BLOCKSIZE;
		if (srt->partials[i_partial].index == i_in)
			blocklen = srt->partials[i_partial++].length;

		assert(j_in < NPERM(srt));
		memmove(srt->pairs + start, blockat(srt, j_in), blocklen * sizeof *srt->pairs);
		perm[i_in] = j_out;
		perm[i_out] = j_in;
		pinv[j_in] = i_out;
		pinv[j_out] = i_in;
		start += blocklen;

		/* if we moved the block from somewhere else, there's now a new free spot */
		if (i_in != i_out) {
			j_free = j_in;
			i_free = i_out;
			assert(i_free >= srt->fill);
		}
	}

	/* permute any remaining blocks (all scratch blocks) */
	for (i_in = j_out - srt->n_scratch; i_in < srt->fill; i_in++) {
		j_in = perm[i_in];
		assert(j_in < srt->n_scratch);

		blocklen = BLOCKSIZE;
		if (srt->partials[i_partial].index == i_in)
			blocklen = srt->partials[i_partial++].length;

		assert(j_in < NPERM(srt));
		memcpy(srt->pairs + start, blockat(srt, j_in), blocklen * sizeof *srt->pairs);
		start += blocklen;
	}

	assert(start == srt->n);
	assert(i_partial == srt->n_partial);
}

extern void
parallel_sort(uint32pair *pairs, size_t n)
{
	struct sorter *srt;

	if (n == 0)
		return;

	srt = sorter_initialize(pairs, n);
	if (srt == NULL) {
		perror("sorter_initialize");
		abort();
	}

	sort_step(srt, 0 * RSHIFT);
	sort_step(srt, 1 * RSHIFT);
	sort_step(srt, 2 * RSHIFT);
	sort_step(srt, 3 * RSHIFT);
	compact(srt);

	sorter_free(srt);
}

const struct radixsort_impl parallel_impl = {
	.name = "parallel",
	.full_sort = parallel_sort,
	.supported = always_supported,
};
