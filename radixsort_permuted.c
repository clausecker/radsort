#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "radixsort.h"

#define RSHIFT CHAR_BIT		/* number of bits per sort step */
#define RADIX (1<<RSHIFT)	/* number of buckets */
#define NSCRATCH (2*RADIX)	/* number of scratch blocks */
#define BSHIFT 9
#define BLOCKSIZE (1<<BSHIFT)	/* elements per block */
#define NPERM(n) (NSCRATCH + (n) / BLOCKSIZE) /* number of blocks in total */

struct bucket {
	uint32pair *start, *end;
};

struct partial {
	unsigned index;		/* logical block index */
	unsigned length;	/* number of items in the partial block */
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
 *  - For each allocated block i, perm[i].bucket holds the bucket the block
 *    belongs to; allocated blocks are sorted in increasing order of buckets.
 *  - perm and perm2 point to distinct arrays of the same length that can be
 *    swapped by the sorting step.
 *  - partials holds the partial blocks in ascending order of indices,
 */
struct sorter {
	uint32pair	*pairs;		/* input array */
	size_t		 n;		/* number of input elements */
	size_t		 fill;		/* end of allocated blocks */
	unsigned	*perm;		/* block permutation (input) */
	unsigned	*perm2;		/* block permutation (output) */
	unsigned char	*assignments;	/* assignments of logical blocks to buckets */
	struct partial	 partials[RADIX]; /* partial input blocks */
	uint32pair	 scratch[NSCRATCH][BLOCKSIZE];
};

static void	fix_perm(struct sorter *, unsigned[RADIX], struct bucket[RADIX]);

/*
 * Retrieve the start address at index i from srt.
 */
static inline uint32pair *
blockat(struct sorter *srt, size_t i)
{
	return (i < NSCRATCH ? srt->scratch[i] : srt->pairs + (i - NSCRATCH) * BLOCKSIZE);
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
 * Get the length of block i when iterating through the block set,
 * expecting the next partial block to be at srt->partials[*i_partial].
 * Increment *i_partial if a partial block was found.
 */
static inline size_t
blocksize(struct sorter *srt, size_t i, size_t *i_partial)
{
	size_t len = BLOCKSIZE;

	assert(*i_partial < RADIX);
	if (srt->partials[*i_partial].index == i) {
		len = srt->partials[*i_partial].length;
		++*i_partial;
	}

	return (len);
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
	unsigned *perm = NULL, *perm2 = NULL;
	unsigned char *assignments;
	size_t i, n_tail;

	/* perm is assumed to be zero-initialised by calloc() */
	srt = malloc(sizeof *srt);
	perm = calloc(sizeof *perm, NPERM(n));
	perm2 = calloc(sizeof *perm2, NPERM(n));
	assignments = calloc(sizeof *assignments, NPERM(n));
	if (srt == NULL || perm == NULL || perm2 == NULL || assignments == NULL) {
		free(srt);
		free(perm);
		free(perm2);
		free(assignments);

		return (NULL);
	}

	srt->pairs = pairs;
	srt->n = n;
	srt->perm = perm;
	srt->perm2 = perm2;
	srt->assignments = assignments;
	srt->fill = RADIX + n / BLOCKSIZE + 1;

	/* runahead: set up first i temp blocks */
	for (i = 0; i < RADIX; i++)
		perm[i] = i;

	/* set up identity permutation for the rest of the array */
	for (; i < srt->fill - 1; i++)
		perm[i] = i + (NSCRATCH - RADIX);

	/* store remaining scratchpad blocks in permutation */
	for (; i < NPERM(n); i++)
		perm[i] = i - (srt->fill - 1) + RADIX;

	/* copy partial input tail to the scratchpad */
	n_tail = n % BLOCKSIZE;
	memcpy(blockat(srt, RADIX), pairs + BLOCKSIZE * (n / BLOCKSIZE), n_tail * sizeof *pairs);
	srt->partials[0].index = srt->fill - 1;
	srt->partials[0].length = n_tail;

	assert(i == NPERM(n));

	return (srt);
}

static void
sorter_free(struct sorter *srt)
{
	free(srt->perm);
	free(srt->perm2);
	free(srt->assignments);
	free(srt);
}

/*
 * Perform one iteration of radix sort at the given shift.
 */
static void
sort_step(struct sorter *srt, int shift)
{
	struct bucket buckets[RADIX];
	size_t i_in, i_out, i_partial;
	unsigned counts[RADIX];

	/* allocate initial buckets */
	for (i_out = 0; i_out < RADIX; i_out++) {
		uint32pair *out_blk;

		out_blk = blockat(srt, srt->perm[i_out]);
		buckets[i_out].start = out_blk;
		buckets[i_out].end = out_blk + BLOCKSIZE;
		srt->assignments[i_out] = i_out;
		counts[i_out] = 1;
	}

	/* invariant: i_out <= i_in */
	i_partial = 0;
	for (i_in = RADIX; i_in < srt->fill; i_in++) {
		size_t j, len;
		uint32pair *in_blk;

		/* grab the next source block */
		assert(srt->perm[i_in] < NPERM(srt->n));
		in_blk = blockat(srt, srt->perm[i_in]);
		len = blocksize(srt, i_in, &i_partial);;

		for (j = 0; j < len; j++) {
			unsigned bucket;
			uint32pair *bucketp;

			bucket = bucketof(in_blk[j], shift);
			bucketp = buckets[bucket].start; /* avoid reload of pointer after memcpy */
			memcpy(bucketp++, in_blk + j, sizeof in_blk[j]);

			/* bucket full? */
			buckets[bucket].start = bucketp;
			if (bucketp == buckets[bucket].end) {
				uint32pair *out_blk;

				assert(i_out <= i_in);
				assert(srt->perm[i_out] < NPERM(srt->n));
				out_blk = blockat(srt, srt->perm[i_out]);
				buckets[bucket].start = out_blk;
				buckets[bucket].end = out_blk + BLOCKSIZE;
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
fix_perm(struct sorter *srt, unsigned counts[RADIX], struct bucket buckets[RADIX])
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

		bucket = srt->assignments[i];
		j = starts[bucket]++;
		newperm[j] = srt->perm[i];

		/* partial block? */
		if (blockat(srt, srt->perm[i]) + BLOCKSIZE == buckets[bucket].end) {
			srt->partials[bucket].index = j;
			srt->partials[bucket].length = BLOCKSIZE -
			    (buckets[bucket].end - buckets[bucket].start);
		}
	}

	/* fill in the unallocated blocks */
	assert(i + RADIX <= NPERM(srt->n));
	memcpy(newperm, srt->perm + i, RADIX * sizeof *srt->perm);
	i += RADIX;
	memcpy(newperm + i, srt->perm + i, (NPERM(srt->n) - i) * sizeof *srt->perm);

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

	for (i = 0; i < NPERM(srt->n); i++) {
		j = perm[i];
		pinv[j] = i;
	}

	start = 0;			/* start of next block in srt->pairs */
	i_partial = 0;			/* next partial block */
	j_free = perm[0];		/* some free block */
	i_free = 0;

	/* permute the blocks overlaying srt->pairs */
	for (j_out = NSCRATCH; j_out < NPERM(srt->n); j_out++) {
		i_in = j_out - RADIX;
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
		assert(j_in < NPERM(srt->n));
		blocklen = blocksize(srt, i_in, &i_partial);
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
	for (i_in = j_out - RADIX; i_in < srt->fill; i_in++) {
		j_in = perm[i_in];
		assert(j_in < NSCRATCH);

		assert(j_in < NPERM(srt->n));
		blocklen = blocksize(srt, i_in, &i_partial);
		memcpy(srt->pairs + start, blockat(srt, j_in), blocklen * sizeof *srt->pairs);
		start += blocklen;
	}

	assert(start == srt->n);
}

extern void
permuted_sort(uint32pair *pairs, size_t n)
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

const struct radixsort_impl permuted_impl = {
	.name = "permuted",
	.full_sort = permuted_sort,
	.supported = always_supported,
};
