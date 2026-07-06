#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "benchmark.h"
#include "radixsort.h"

static	size_t nparam = 1000000;
static	uint32pair *data;
static	uint64_t cksum;
static	uint64_t seed = 0xdeadbeefc0febabe;

static const struct radixsort_impl *impls[] = {
	&generic_impl,
	NULL,
};

static void
usage(const char *argv0)
{
	fprintf(stderr, "usage: %s [-n length] [-s seed]\n", argv0);
	exit(EXIT_FAILURE);
}

/* random data generator using a xorshift rng */
static void
fill_array(uint32pair *buf, size_t n)
{
	size_t i;
	uint64_t a = seed, sum = 0;

	for (i = 0; i < n; i ++) {
		buf[i][0] = a >> 32;
		buf[i][1] = (uint32_t)a;
		sum += a;

		a ^= a << 13;
		a ^= a >>  7;
		a ^= a << 17;
	}

	cksum = sum;
}

static void
check_sort(uint32pair *buf, size_t n)
{
	size_t i;
	uint64_t sum;

	if (n == 0)
		return;

	sum = (uint64_t)buf[0][0] << 32 | buf[0][1];
	for (i = 1; i < n; i++) {
		sum += (uint64_t)buf[i][0] << 32 | buf[i][1];

		if (buf[i - 1][0] > buf[i][0]) {
			fprintf(stderr, "warning: array not sorted!\n");
			return;
		}
	}

	if (sum != cksum)
		fprintf(stderr, "warning: checksum mismatch!\n");
}

static void
run_full(struct B *b, void *payload)
{
	struct radixsort_impl *impl = payload;
	long i;

	b->bytes = nparam * sizeof *data;

	resettimer(b);

	for (i = 0; i < b->n; i++) {
		stoptimer(b);
		fill_array(data, nparam);
		starttimer(b);
		impl->full_sort(data, nparam);
	}

	stoptimer(b);

	check_sort(data, nparam);
}

int
main(int argc, char *argv[])
{
	size_t bytes;
	int ch, res;
	char dummy;

	while (ch = getopt(argc, argv, "n:s:"), ch != EOF)
		switch (ch) {
		case 'n':
			res = sscanf(optarg, "%zu %c", &nparam, &dummy);
			if (res != 1) {
				fprintf(stderr, "can't parse -n argument: %s\n", optarg);
				usage(argv[0]);
			}

			break;

		case 's':
			seed = (uint64_t)strtoull(optarg, NULL, 0);
			break;

		default:
			usage(argv[0]);
		}

	if (optind != argc)
		usage(argv[0]);

	/* align data for the swc implementation */
	bytes = nparam * sizeof *data;
	data = aligned_alloc(128, bytes + 127 & ~(size_t)127);
	if (data == NULL) {
		perror("aligned_alloc");
		return (EXIT_FAILURE);
	}

	printf("items: %zu\n", nparam);
#ifdef _OPENMP
	printf("threads: %d\n", omp_get_max_threads());
#endif

	preamble();

	runbenchmark("genericFull", run_full, (void *)&generic_impl);
	runbenchmark("swcFull", run_full, (void *)&swc_impl);
	runbenchmark("permutedFull", run_full, (void *)&permuted_impl);
	runbenchmark("bitmanipFull", run_full, (void *)&bitmanip_impl);
	runbenchmark("parallelFull", run_full, (void *)&parallel_impl);

	free(data);

	return (EXIT_SUCCESS);
}
