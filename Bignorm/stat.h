#ifndef STAT_H
#define STAT_H

#include <stdint.h>

extern uint64_t count_reads_unpaired;
extern uint64_t count_reads_paired;
extern uint64_t count_kmers;
extern uint64_t count_newkmers;
extern uint64_t count_reads_unpaired_out;
extern uint64_t count_reads_paired_out;
extern uint64_t count_reads_too_short;
extern uint64_t count_reads_too_many_N;
extern uint64_t count_kmers_with_N;
extern uint64_t in_length;

extern uint64_t waits_reader;
extern uint64_t waits_hasher;
extern uint64_t waits_decider;
extern uint64_t waits_writer;

extern uint64_t histo[];

void print_stats();

#ifndef max
#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b);\
	_a > _b ? _a : _b; })
#endif

#ifndef min
#define min(a,b) \
	({ __typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b);\
	_a < _b ? _a : _b; })
#endif

#endif
