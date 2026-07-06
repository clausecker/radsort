#ifndef RADIXSORT_H
#define RADIXSORT_H

/* Radix sorting 32 bit integers */

#include <stdint.h>

typedef uint32_t uint32pair[2];

extern const struct radixsort_impl {
	const char *name;

	/*
	 * sort the given array of pairs by first entry
	 */
	void (*full_sort)(uint32pair *, size_t);

	/*
	 * Returns nonzero if this implementation is supported by the hardware.
	 * If the function pointer is NULL, it is not supported.
	 */
	int (*supported)(void);
} generic_impl, swc_impl, permuted_impl, bitmanip_impl, parallel_impl;

extern int always_supported(void);

#endif /* RADIXSORT_H */
