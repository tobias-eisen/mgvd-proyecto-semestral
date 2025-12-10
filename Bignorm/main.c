#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h> 

#include "options.h"
#include "ringbuffer.h"
#include "stat.h"

#define THREADS_NEEDED 4 

struct ringbuffer_t ringbuffer[RINGBUFFER_SIZE];
omp_lock_t ringlocks[RINGBUFFER_SIZE];

int main(int argc, char **argv){
	int i;
	uint64_t i2;
	struct timespec ts_start, ts_end;

	clock_gettime(0, &ts_start);
	memset(&ringbuffer,0,RINGBUFFER_SIZE * sizeof(struct ringbuffer_t));
	for(i=0; i<RINGBUFFER_SIZE; i++)
		omp_init_lock(&ringlocks[i]);

	if(options(argc,argv) == 0){
		if(verbosity > 1){
			printf("Bignorm starting\n----------------\n\n");
		}
		printf("command line:\n");
		for(i=0; i<argc; i++) printf("%s ", argv[i]);
		printf("\n\n");
		if(count_unpaired + count_paired > 0){
			if(verbosity > 2){
				printf("parameters\n\tverbosity: %i\n\tk:        %i\n", verbosity,k);
				printf("\tstrategy: %i, parameter: %i\n", d, C);
				printf("\tmemory strategy: %i, mem_frac: %f\n", m, mem_frac);
				printf("\tinput files: %i unpaired files, %i paired files:\n", count_unpaired, count_paired);
				input_akt = input_start;
				while(input_akt != NULL){
					if(input_akt->flag == 1){
						printf("\t\t%s\n",input_akt->infile1);
					} else {
						printf("\t\t%s <-> %s\n",input_akt->infile1, input_akt->infile2);
					}
					input_akt = input_akt->next;
				}
			}
			omp_set_num_threads(THREADS_NEEDED);
#pragma omp parallel sections 

{
#pragma omp section
{
// reader
	int tid = omp_get_thread_num();
	int ret = omp_get_num_threads();
	if(ret != THREADS_NEEDED){
		printf("OpenMP error: couldn't get %i threads (got %i instead)\n", THREADS_NEEDED, ret);
		//exit(-1);
	}
	printf("reader  started, id: %i\n", tid);
	reader();
	printf("reader ended\n");
}
#pragma omp section
{
        int tid = omp_get_thread_num();
        printf("hasher started, id: %i\n", tid);
        hasher();
        printf("hasher ended\n");
}
#pragma omp section
{
        int tid = omp_get_thread_num();
        printf("decider started, id: %i\n", tid);
        decider();
        printf("decider ended\n");
}
#pragma omp section
{
// writer
	int tid = omp_get_thread_num();
	printf("writer started, id: %i\n", tid);
	writer();
	printf("writer ended\n");		
}

}
			// processing ended
                        for(i2=0; i2 < memneeded; i2++){
                                histo[countmin[i2]]++;
                        }

			free(countmin);
			if(verbosity > 2)
				print_stats();
			// give statistics
			if(verbosity > 1){
				clock_gettime(0, &ts_end);
				printf("time needed: %u s\n", ts_end.tv_sec - ts_start.tv_sec);
			} 
		} else {
			printf("no input files!\nexiting now\n");
			exit(-1);
		}
	} else {
		printf("Error parsing options\n");
	}
	return 0;
}

