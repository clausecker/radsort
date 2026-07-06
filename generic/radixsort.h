#ifndef RADIXSORT_H
#define RADIXSORT_H

#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>

struct rs_sorter;
struct rs_thread;

typedef void rs_key;
typedef void rs_step_func(struct rs_sorter *, struct rs_thread *, size_t);

enum { RS_RADIX = UCHAR_MAX + 1 };

struct rs_partial {
	unsigned		 index;		/* logical block index */
	unsigned		 length;	/* number of items in the partial block */
};

struct rs_thread {
	unsigned		 i_start;	/* first block in thread chunk */
	unsigned		 i_partial;	/* first partial block in thread chunk */
	unsigned		 fill;		/* end of allocated blocks in thread chunk */
	unsigned		 counts[RS_RADIX];	/* number of blocks per output bucket */
	rs_key			*buckets[RS_RADIX];	/* output buckets */
};

struct rs_sorter {
	rs_key			*main;		/* blocks overlaying the main array */
	rs_key			*scratch;	/* extra blocks */
	size_t			 n;		/* number of input array elements */
	size_t			 esize;		/* element size */
	uintptr_t		 blockmask;	/* mask identifying pointers to bucket starts */
	unsigned		*perm;		/* block permutation */
	unsigned		*perm2;		/* aux block permutation */
	unsigned char		*assignments;	/* assignments of logical blocks to buckets */
	struct rs_partial	*partials;	/* partial blocls */
	struct rs_thread	*threads;	/* per-thread structures */
	rs_key			*array;		/* main array */
	unsigned		 bsize;		/* block size */
	unsigned		 bshift;	/* log2 of block size */
	unsigned		 n_scratch;	/* number of scratch blocks */
	unsigned		 n_block;	/* number of blocks in total */
	unsigned		 n_thread;	/* max number of threads */
	unsigned		 n_partial;	/* number of partial blocks */
	unsigned		 fill;		/* number of allocated blocks in perm */
	rs_step_func		*step;
};

struct rs_sorter	*rs_new(rs_key *, size_t, size_t, rs_step_func *);
void			 rs_free(struct rs_sorter *);
rs_key			*rs_finalize(struct rs_sorter *);
void			 rs_step(struct rs_sorter *, size_t);

static inline rs_key *
rs_blockat(struct rs_sorter *srt, unsigned i)
{
	if (i < srt->n_scratch)
		return ((unsigned char *)srt->scratch + (i * srt->esize << srt->bshift));
	else
		return ((unsigned char *)srt->main + ((i - srt->n_scratch) * srt->esize << srt->bshift));
}

#define RS_MAKE_SORTER(name, type) \
void rs_step_##name(struct rs_sorter *srt, struct rs_thread *thr, size_t keyb) \
{ \
	unsigned i_in, i_out, i_partial; \
	type **buckets; \
\
	assert(sizeof(type) == srt->esize); \
	if (sizeof(type) != srt->esize) \
		__builtin_unreachable(); \
	buckets = (type **)thr->buckets; \
	i_partial = thr->i_partial; \
	for (i_in = i_out = thr->i_start + RADIX; i_in < thr->fill; i_in++) { \
		unsigned j, len; \
		type *in_blk; \
\
		in_blk = rs_blockat(srt, srt->perm2[i_in]); \
		len = srt->bsize; \
		if (srt->partials[i_partial].index == i_in) \
			len = srt->partials[i_partial++].length; \
\
		for (j = 0; j < len; j++) { \
			unsigned bucket; \
			type *bucketp; \
\
			bucket = ((unsigned char *)(in_blk + j))[keyb]; \
			bucketp = buckets[bucket]; \
			memcpy(bucketp++, in_blk + j, sizeof in_blk[j]); \
\
			/* bucket full? */ \
			buckets[bucket] = bucketp; \
			if (((uintptr_t)bucketp & srt->blockmask) == 0) { \
				type *out_blk; \
\
				out_blk = rs_blockat(srt, srt->perm2[i_out]); \
				buckets[bucket] = out_blk; \
				srt->assignments[i_out] = bucket; \
				thr->counts[bucket]++; \
				i_out++; \
			} \
		} \
	} \
\
	thr->fill = i_out; \
}

#define RS_NEW(array, nmemb, type) rs_new(array, nmemb, sizeof *(array), rs_step_##type)

#endif /* RADIXSORT_H */
