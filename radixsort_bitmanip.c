#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "radixsort.h"

#define RSHIFT CHAR_BIT		/* number of bits per sort step */
#define RADIX (1<<RSHIFT)	/* number of buckets */
#define NSCRATCH (2*RADIX+1)	/* number of scratch blocks */
#define BSHIFT 9
#define BLOCKSIZE (1<<BSHIFT)	/* elements per block */

struct partial {
	unsigned index;		/* logical block index */
	unsigned length;	/* number of items in the partial block */
};

/*
 * The sorter structure represents the state of the sorting process.
 * Between sorting iterations, the following invariants hold:
 *
 *  - n holds the number of pairs in pairs.
 *  - perm[0] to perm[NSCRATCH - 1] and perm[fill] to perm[n_perm]
 *    refer to unallocated blocks
 *  - perm[NSCRATCH] to perm[fill - 1] refer to allocated blocks
 *    in the of the current sorting progress.
 *  - For each allocated block i, perm[i].bucket holds the bucket the block
 *    belongs to; allocated blocks are sorted in increasing order of buckets.
 *  - perm and perm2 point to distinct arrays of the same length that can be
 *    swapped by the sorting step.
 *  - partials holds the partial blocks in ascending order of indices,
 */
struct sorter {
	uint32pair	 *main;		/* blocks overlaying the main array */
	uint32pair	 *pairs;	/* input array */
	size_t		  n;		/* number of input elements */
	size_t		  n_perm;	/* number if blocks in total */
	size_t		  fill;		/* end of allocated blocks */
	unsigned	 *perm;		/* block permutation (input) */
	unsigned	 *perm2;	/* block permutation (output) */
	unsigned char	 *assignments;	/* assignments of logical blocks to buckets */
	uint32pair	(*scratch)[BLOCKSIZE];
	struct partial	  partials[RADIX]; /* partial input blocks */
};

static void	fix_perm(struct sorter *, unsigned[RADIX], uint32pair*[RADIX]);

/*
 * Retrieve the start address at index i from srt.
 */
static inline uint32pair *
blockat(struct sorter *srt, size_t i)
{
	return (i < NSCRATCH ? srt->scratch[i] : srt->main + (i - NSCRATCH) * BLOCKSIZE);
}

/*
 * Get the bucket a given uint32pair belongs into according to shift.
 */
static inline unsigned
bucketof(uint32pair p, int shift)
{
	return ((unsigned)(p[0] >> shift & RADIX - 1));
}

/*
 * Create a new sorter from the given pairs and length n.  Return NULL on failure.
 * The newly created sorter represents an identity permutation of the input, with
 * all blocks being in bucket 0.
 */
static struct sorter *
sorter_initialize(uint32pair *pairs, size_t n)
{
	struct sorter *srt = NULL;
	uint32pair (*scratch)[BLOCKSIZE];
	unsigned *perm = NULL, *perm2 = NULL;
	unsigned char *assignments;
	size_t i, i_scratch, n_tail, start, n_perm;

	/* first aligned address in srt->pairs */
	start = (BLOCKSIZE - (uintptr_t)pairs / sizeof **scratch) % BLOCKSIZE;
	if (n >= 2 * BLOCKSIZE)
		n_perm = NSCRATCH + (n - start) / BLOCKSIZE;
	else
		n_perm = NSCRATCH;

	/* perm is assumed to be zero-initialised by calloc() */
	srt = malloc(sizeof *srt);
	perm = calloc(sizeof *perm, n_perm);
	perm2 = calloc(sizeof *perm2, n_perm);
	assignments = calloc(sizeof *assignments, n_perm);
	scratch = aligned_alloc(sizeof *scratch, NSCRATCH * sizeof *scratch);
	if (srt == NULL || perm == NULL || perm2 == NULL || assignments == NULL || scratch == NULL) {
		free(srt);
		free(perm);
		free(perm2);
		free(assignments);
		free(scratch);

		return (NULL);
	}

	srt->pairs = pairs;
	srt->main = srt->pairs + start;
	srt->n = n;
	srt->n_perm = n_perm;
	srt->perm = perm;
	srt->perm2 = perm2;
	srt->assignments = assignments;
	srt->scratch = scratch;

	/* runahead and head: set up temp blocks */
	for (i = i_scratch = 0; i < RADIX + 1; i++)
		perm[i] = i_scratch++;

	srt->fill = RADIX + 1 + (n - start) / BLOCKSIZE + 1;

	if (n >= 2 * BLOCKSIZE) {
		/* set up identity permutation for the rest of the array */
		for (; i < srt->fill - 1; i++)
			perm[i] = i - (RADIX + 1) + NSCRATCH;

		/* if pairs is not aligned, transfer  unaligned prefix to a scratch block */
		assert((uintptr_t)srt->main / sizeof **scratch % BLOCKSIZE == 0);
		assert((uintptr_t)srt->scratch / sizeof **scratch % BLOCKSIZE == 0);
		memcpy(blockat(srt, RADIX), srt->pairs, start * sizeof *srt->pairs);
		srt->partials[0].index = RADIX;
		srt->partials[0].length = start;

		/* copy partial input tail to the scratchpad */
		n_tail = (srt->n - start) % BLOCKSIZE;
		memcpy(blockat(srt, RADIX + 1), srt->main + BLOCKSIZE * ((n - start) / BLOCKSIZE), n_tail * sizeof *pairs);
		srt->partials[1].index = srt->fill - 1;
		srt->partials[1].length = n_tail;
	} else {
		/* very short array: transfer it entirely to scratch blocks */
		memcpy(blockat(srt, RADIX), srt->pairs, n * sizeof *srt->pairs);

		srt->fill = RADIX + 1 + (n >= BLOCKSIZE);
		srt->partials[0].index = srt->fill - 1;
		srt->partials[0].length = n % BLOCKSIZE;
	}

	/* store remaining scratchpad blocks in permutation */
	for (; i < n_perm; i++)
		perm[i] = i_scratch++;

	assert(i == n_perm);
	assert(i_scratch == NSCRATCH);

	return (srt);
}

static void
sorter_free(struct sorter *srt)
{
	free(srt->perm);
	free(srt->perm2);
	free(srt->assignments);
	free(srt->scratch);
	free(srt);
}

/*
 * Perform one iteration of radix sort at the given shift.
 */
static void
sort_step(struct sorter *srt, int shift)
{
	uint32pair *buckets[RADIX];
	size_t i_in, i_out, i_partial;
	unsigned counts[RADIX];

	/* allocate initial buckets */
	for (i_out = 0; i_out < RADIX; i_out++) {
		uint32pair *out_blk;

		out_blk = blockat(srt, srt->perm[i_out]);
		buckets[i_out] = out_blk;
		srt->assignments[i_out] = i_out;
		counts[i_out] = 1;
	}

	/* invariant: i_out <= i_in */
	i_partial = 0;
	for (i_in = RADIX; i_in < srt->fill; i_in++) {
		size_t j, len;
		uint32pair *in_blk;

		/* grab the next source block */
		assert(srt->perm[i_in] < srt->n_perm);
		in_blk = blockat(srt, srt->perm[i_in]);
		len = BLOCKSIZE;
		if (srt->partials[i_partial].index == i_in)
			len = srt->partials[i_partial++].length;

		for (j = 0; j < len; j++) {
			unsigned bucket;
			uint32pair *bucketp;

			bucket = bucketof(in_blk[j], shift);
			bucketp = buckets[bucket]; /* avoid reload of pointer after memcpy */
			memcpy(bucketp++, in_blk + j, sizeof in_blk[j]);
			__builtin_prefetch(bucketp, 1, 3);

			/* bucket full? */
			buckets[bucket] = bucketp;
			if ((uintptr_t)bucketp / sizeof *buckets[bucket] % BLOCKSIZE == 0) {
				uint32pair *out_blk;

				assert(i_out <= i_in);
				assert(srt->perm[i_out] < srt->n_perm);
				out_blk = blockat(srt, srt->perm[i_out]);
				buckets[bucket] = out_blk;
				srt->assignments[i_out] = bucket;
				counts[bucket]++;
				i_out++;
			}
		}
	}

	srt->fill = i_out;
	fix_perm(srt, counts, buckets);
}

/*
 * After the main part of sort_step() has executed, srt holds blocks sorted according
 * to the current shift.  srt->perm[0] to srt->perm[srt->fill - 1] hold the blocks
 * of the newly sorted data, but with the buckets interleaved arbitrarily.  This
 * routine computes a new permutation that re-establishes the data structure invariant.
 */
static void
fix_perm(struct sorter *srt, unsigned counts[RADIX], uint32pair *buckets[RADIX])
{
	unsigned *newperm, starts[RADIX];
	size_t i;

	/* compute start offsets of the newly allocated buckets */
	starts[0] = RADIX;
	for (i = 1; i < RADIX; i++)
		starts[i] = starts[i - 1] + counts[i - 1];

	/* compute new inverse permutation and fill srt->partials */
	newperm = srt->perm2;
	for (i = 0; i < srt->fill; i++) {
		size_t j;
		int bucket;
		uint32pair *block;

		bucket = srt->assignments[i];
		j = starts[bucket]++;
		newperm[j] = srt->perm[i];

		/* partial block? */
		block = blockat(srt, srt->perm[i]);
		if (block <= buckets[bucket] && buckets[bucket] < block + BLOCKSIZE) {
			srt->partials[bucket].index = j;
			srt->partials[bucket].length = buckets[bucket] - block;
		}
	}

	/* fill in the unallocated blocks */
	assert(i + RADIX <= srt->n_perm);
	memcpy(newperm, srt->perm + i, RADIX * sizeof *srt->perm);
	i += RADIX;
	memcpy(newperm + i, srt->perm + i, (srt->n_perm - i) * sizeof *srt->perm);

	srt->fill += RADIX;
	srt->perm2 = srt->perm;
	srt->perm = newperm;
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

	for (i = 0; i < srt->n_perm; i++) {
		j = perm[i];
		pinv[j] = i;
	}

	start = 0;			/* start of next block in srt->pairs */
	i_partial = 0;			/* next partial block */
	j_free = perm[0];		/* some free block */
	i_free = 0;

	/* permute the blocks overlaying srt->pairs */
	for (j_out = NSCRATCH; j_out < srt->n_perm; j_out++) {
		i_in = j_out - (NSCRATCH - RADIX);
		i_out = pinv[j_out];
		j_in = perm[i_in];

		/* j_out not the right block, but allocated? */
		if (i_out > i_in && i_out < srt->fill) {
			/* move block at j_out to j_free */
			assert(i_free < RADIX || i_free >= srt->fill);
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

		assert(j_in < srt->n_perm);
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
			assert(i_free < RADIX || i_free >= srt->fill);
		}
	}

	/* permute any remaining blocks (all scratch blocks) */
	for (i_in = j_out - (NSCRATCH - RADIX); i_in < srt->fill; i_in++) {
		j_in = perm[i_in];
		assert(j_in < NSCRATCH);

		blocklen = BLOCKSIZE;
		if (srt->partials[i_partial].index == i_in)
			blocklen = srt->partials[i_partial++].length;

		assert(j_in < srt->n_perm);
		memcpy(srt->pairs + start, blockat(srt, j_in), blocklen * sizeof *srt->pairs);
		start += blocklen;
	}

	assert(start == srt->n);
}

extern void
bitmanip_sort(uint32pair *pairs, size_t n)
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

const struct radixsort_impl bitmanip_impl = {
	.name = "bitmanip",
	.full_sort = bitmanip_sort,
	.supported = always_supported,
};
