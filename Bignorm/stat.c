#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stat.h"
#include "options.h"
#include "ringbuffer.h"

uint64_t count_reads_unpaired = 0;
uint64_t count_reads_paired = 0;
uint64_t count_kmers = 0;
uint64_t count_newkmers = 0;
uint64_t count_reads_unpaired_out = 0;
uint64_t count_reads_paired_out = 0;
uint64_t count_reads_too_short = 0;
uint64_t count_kmers_with_N = 0;
uint64_t count_reads_too_many_N = 0;
uint64_t in_length = 0;

uint64_t waits_reader = 0;
uint64_t waits_hasher = 0;
uint64_t waits_decider = 0;
uint64_t waits_writer = 0;

uint64_t histo[256] = { 0UL};

void print_stats(){
	int i;
	uint64_t n_occupied = 0;
	uint64_t sum = 0;
	float temp1, temp2;

	printf("\n\tStatistics\n\t----------\n\n");
	printf("Reads processed:\n\tunpaired: %u\n\tpaired  : %u\n\t----------\n\ttotal   : %u\n\n", count_reads_unpaired, count_reads_paired, count_reads_paired + count_reads_unpaired);
	printf("Reads shorter than k (%i): %u\n", k, count_reads_too_short);
	if(opt_N >= 0) printf("Reads containing more than %i Ns: %u\n", opt_N, count_reads_too_many_N);
	if(opt_n > 0){
		if(count_kmers_with_N > 100000000){
			printf("K-mers containing N: %u%u \n", count_kmers_with_N / 100000000, count_kmers_with_N % 100000000);
		} else {
			printf("K-mers containing N: %u\n", count_kmers_with_N);
		}
	}
	if(count_kmers > 100000000){
		printf("K-mers processed: %u%u", count_kmers/100000000, count_kmers%100000000);
		if(count_newkmers > 100000000){
 			printf(" (at least %u%u different)\n", count_newkmers/100000000, count_newkmers%100000000);
		} else {
			printf(" (at least %u different)\n", count_newkmers);
		}
	} else {
		printf("K-mers processed: %u (at least %u different)\n", count_kmers, count_newkmers);
	}
	printf("Reads accepted:\n\tunpaired: %u \n\tpaired  : %u\n\t----------\n\ttotal   : %i\n\n", count_reads_unpaired_out, count_reads_paired_out, count_reads_unpaired_out + count_reads_paired_out);
	if(verbosity > 3){
		printf("wait-cycles by:\n\treader : %u\n\thasher : %u\n\tdecider: %u\n\twriter : %u\n\t-------\n\ttotal  : %u\n", waits_reader, waits_hasher, waits_decider, waits_writer, waits_reader + waits_hasher + waits_decider + waits_writer);
	}
	if(do_hist == 1){
		printf("Histogram:\n");
		for(i=0; i<255; i++){
			if(histo[i] > 100000000){
				printf("\t%i:\t%u%u\n", i, histo[i]/100000000, histo[i]%100000000);
			} else {
				printf("\t%i:\t%u\n", i, histo[i]);
			}
		}
		if(histo[COUNTMAX] > 100000000){
			printf("\nCOUNTMAX reached %u%u times\n", histo[COUNTMAX]/100000000,histo[COUNTMAX]%100000000);
		} else {
			printf("\nCOUNTMAX reached %u times\n", histo[COUNTMAX]);
		}
		if(histo[COUNTMAX] == 0){
			for(i = COUNTMAX; i > 0; i--){
				if(histo[i] > 0){
					printf("max. count reached: %i\n", i);
					break;
				}
			}
		} 
	}
	// calc fp value
	for(i=1; i<255; i++){
		n_occupied += histo[i];
		sum += i * histo[i];
	}
	if(n_occupied > 100000000){
		if(histo[0] > 100000000){
			printf("n_occupied:\t%u%u\nhist 0:\t%u%u\ndiff:\t%u\n", n_occupied / 100000000, n_occupied % 100000000, histo[0] / 100000000, histo[0] % 100000000, memneeded - histo[0] - n_occupied);
		} else {
			printf("n_occupied:\t%u%u\nhist 0:\t%u\ndiff:\t%u\n", n_occupied / 100000000, n_occupied % 100000000, histo[0], memneeded - histo[0] - n_occupied);
		}
	} else {
		if(histo[0] > 100000000){
			printf("n_occupied:\t%u\nhist 0:\t%u%u\ndiff:\t%u\n", n_occupied, histo[0] / 100000000, histo[0] % 100000000, memneeded - histo[0] - n_occupied);
		} else {
			printf("n_occupied:\t%u\nhist 0:\t%u\ndiff:\t%u\n", n_occupied, histo[0], memneeded - histo[0] - n_occupied);
		}
	}
	printf("Sum over all counters: %lu (%.2f per hashline, %.2f per distinct k-mer)\n", sum, (float)sum / t, (float)sum / (t * count_newkmers));
	temp1 = ((float)n_occupied)/((float)m * t * 1024 * 1024);
	temp2 = temp1;
	for(i=1; i<t; i++){
		temp2 *= temp1;
	}
	printf("fp:\t%.3f (base: %.3f)\n", temp2, temp1);
	printf("reads kept: %.2f \%\n", (float) (count_reads_unpaired_out + count_reads_paired_out) * 100 / (count_reads_paired + count_reads_unpaired));		
}
