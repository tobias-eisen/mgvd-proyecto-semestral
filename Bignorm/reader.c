#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "options.h"
#include "ringbuffer.h"
#include "stat.h"

#include "thirdparty/kseq.h"

KSEQ_INIT(gzFile, gzread)

#define TICD 100000UL

int reader(){
	gzFile f1,f2;
	kseq_t *ks1, *ks2;
	int ret, ret2, i, j, ncount, fnr,qtr;
	uint64_t nexttic, readcount = 0UL;

	input_akt = input_start;
	i = 0;
	fnr = 0;
	qtr = 0;
	nexttic = TICD;

	while(input_akt != NULL){
		if(input_akt->flag == 1){
			// unpaired file
			if(verbosity > 3) printf("Reading unpaired file %s\n", input_akt->infile1);
			f1 = gzopen(input_akt->infile1, "r");
			if(f1 == NULL){
				perror("-----------------\nerror opening input file!\nSkipping it\n");
				printf("-----------------\n");
			} else {
				ks1 = kseq_init(f1);
				while(ret = kseq_read(ks1) != -1){
					if((ret > 0) ||(ret == -2)){
						if(verbosity > 5)
							printf("Name: %s, comment: %s\nseq: %s\n\n", ks1->name.s,ks1->comment.s, ks1->seq.s);
						if(ret == -2){
							qtr++;
							if((verbosity > 3) && (qtr < 10))
								printf("Q-String truncated for %s\n", ks1->name.s);
						}
						while(ringbuffer[i].status != STAT_READER){
							// wait for entry to become available
							usleep(WAITTIME);
							waits_reader++;
						}
						// fill read data to ringbuffer
						omp_set_lock(&ringlocks[i]);
						ringbuffer[i].flag = 1;
						ringbuffer[i].filenr = fnr;
						strncpy(ringbuffer[i].read1, ks1->seq.s, min(ks1->seq.l,R_MAX));
						strncpy(ringbuffer[i].name1, ks1->name.s, min(ks1->name.l,Name_MAX-1));
						ringbuffer[i].name1[min(ks1->name.l,Name_MAX-1)] = '\0';
						strncpy(ringbuffer[i].comment1, ks1->comment.s, min(ks1->comment.l, Comment_MAX-1));
						ringbuffer[i].comment1[min(ks1->comment.l, Comment_MAX-1)] = '\0';
						strncpy(ringbuffer[i].qual1, ks1->qual.s, min(ks1->qual.l, R_MAX));
						if(ks1->qual.l > 0) ringbuffer[i].flag |= IS_FASTQ;
						ringbuffer[i].rlen1 = min(ks1->seq.l,R_MAX);
						// len < k? -> reject
						if(ringbuffer[i].rlen1 < k){
							ringbuffer[i].flag |= READ_REJ;
							count_reads_too_short++;
						} else {
							if(opt_N >= 0){
								ncount = 0;
								for(j=0; j < ringbuffer[i].rlen1; j++){
									if((ringbuffer[i].read1[j] =='N') || (ringbuffer[i].read1[j] == 'n')){
										ncount++;
										if(ncount > opt_N){
											ringbuffer[i].flag |= READ_REJ;
											count_reads_too_many_N++;
											break;
										}
									}
								}
							}
							
						}
						ringbuffer[i].status = STAT_TRIM;
						omp_unset_lock(&ringlocks[i]);
						count_reads_unpaired++;		
					} else {
						printf("kseq_read returned %i\n", ret);
					}
#ifdef HAVE_RINGBUFFER_MASK
                			i = (i+1) & RINGBUFFER_MASK;
#else
			                i = (i+1)%RINGBUFFER_SIZE;
#endif
					if(verbosity > 2){
                                                readcount++;
                                                if(readcount >= nexttic){
                                                        nexttic += TICD;
                                                        printf("reader: %u reads processed\n", readcount);
                                                }
                                        }

				}
				kseq_destroy(ks1);
				gzclose(f1);
			}
		} else {
			// paired files
			if(verbosity > 3) printf("Reading paired files %s <-> %s\n", input_akt->infile1, input_akt->infile2);
			f1 = gzopen(input_akt->infile1, "r");
			f2 = gzopen(input_akt->infile2, "r");
			if((f1 == NULL) || (f2 == NULL)){
				perror("-----------------\nerror opening input file!\nSkipping it\n");
				printf("-----------------\n");
			} else {
				ks1 = kseq_init(f1);
				ks2 = kseq_init(f2);
				while((ret = kseq_read(ks1) != -1) && (ret2 = kseq_read(ks2) != -1)){
					if(((ret > 0) ||(ret == -2)) && ((ret2 > 0) || (ret2 == -2))){
						if(verbosity > 5){
							printf("read 1: Name: %s, comment: %s\nseq: %s\n\n", ks1->name.s,ks1->comment.s, ks1->seq.s);
							printf("read 2: Name: %s, comment: %s\nseq: %s\n\n", ks2->name.s,ks2->comment.s, ks2->seq.s);
						}
						if(ret == -2){
							qtr++;
							if((verbosity > 3) && (qtr < 10))
								printf("Q-String truncated for %s\n", ks1->name.s);
						}
						if(ret2 == -2){
							qtr++;
							if((verbosity > 3) && (qtr < 10))
								printf("Q-String truncated for %s\n", ks2->name.s);
						}
						while(ringbuffer[i].status != STAT_READER){
							// wait for entry to become available
							usleep(WAITTIME);
							waits_reader++;
						}
						// fill read data to ringbuffer
						omp_set_lock(&ringlocks[i]);
						ringbuffer[i].flag = 2;
						ringbuffer[i].filenr = fnr;
						strncpy(ringbuffer[i].read1, ks1->seq.s, min(ks1->seq.l,R_MAX));
						strncpy(ringbuffer[i].name1, ks1->name.s, min(ks1->name.l,Name_MAX));
						ringbuffer[i].name1[min(ks1->name.l,Name_MAX-1)] = '\0';
						strncpy(ringbuffer[i].comment1, ks1->comment.s, min(ks1->comment.l, Comment_MAX));
						ringbuffer[i].comment1[min(ks1->comment.l, Comment_MAX-1)] = '\0';
						strncpy(ringbuffer[i].qual1, ks1->qual.s, min(ks1->qual.l, R_MAX));
						if(ks1->qual.l > 0) ringbuffer[i].flag |= IS_FASTQ;
						strncpy(ringbuffer[i].read2, ks2->seq.s, min(ks2->seq.l,R_MAX));
						strncpy(ringbuffer[i].name2, ks2->name.s, min(ks2->name.l,Name_MAX));
						ringbuffer[i].name2[min(ks2->name.l,Name_MAX-1)] = '\0';
						strncpy(ringbuffer[i].comment2, ks2->comment.s, min(ks2->comment.l, Comment_MAX));
						ringbuffer[i].comment2[min(ks2->comment.l, Comment_MAX-1)] = '\0';
						strncpy(ringbuffer[i].qual2, ks2->qual.s, min(ks2->qual.l, R_MAX));
						ringbuffer[i].rlen1 = min(ks1->seq.l,R_MAX);
						ringbuffer[i].rlen2 = min(ks2->seq.l,R_MAX);
						// len < k? -> reject
                                                if((ringbuffer[i].rlen1 < k) || (ringbuffer[i].rlen2 < k)){
                                                        ringbuffer[i].flag |= READ_REJ;
							count_reads_too_short++;
                                                }
                                                if(opt_N >= 0){
							ncount = 0;
							for(j=0; j < ringbuffer[i].rlen1; j++){
								if((ringbuffer[i].read1[j] =='N') || (ringbuffer[i].read1[j] == 'n')){
									ncount++;
									if(ncount > opt_N){
										ringbuffer[i].flag |= READ_REJ;
										count_reads_too_many_N++;
										break;
									}
								}
							}
							if(ncount <= opt_N){
								for(j=0; j < ringbuffer[i].rlen2; j++){
									if((ringbuffer[i].read2[j] =='N') || (ringbuffer[i].read2[j] == 'n')){
										ncount++;
										if(ncount > opt_N){
											ringbuffer[i].flag |= READ_REJ;
											count_reads_too_many_N++;
											break;
										}	
									}
								}
							}
						}
						ringbuffer[i].status = STAT_TRIM;
						omp_unset_lock(&ringlocks[i]);
						count_reads_paired++;
					} else {
						printf("kseq_read returned %i / %i\n", ret, ret2);
					}
#ifdef HAVE_RINGBUFFER_MASK
			                i = (i+1) & RINGBUFFER_MASK;
#else
 			 		i = (i+1)%RINGBUFFER_SIZE;
#endif
					if(verbosity > 2){
						readcount++;
						if(readcount >= nexttic){
							nexttic += TICD;
							printf("reader: %u reads processed\n", readcount);
						}
					}
				}
				kseq_destroy(ks1);
				kseq_destroy(ks2);
				gzclose(f1);
				gzclose(f2);
			}	
		}
		input_akt = input_akt->next;
		fnr++;
	}
	while(ringbuffer[i].status != STAT_READER){
		usleep(WAITTIME);
		waits_reader++;
	}
	omp_set_lock(&ringlocks[i]);
	ringbuffer[i].flag = -1; // input ended
	ringbuffer[i].status = STAT_TRIM;
	omp_unset_lock(&ringlocks[i]);
	if(verbosity > 2){
		printf("reader: %u reads processed, finished\n", readcount);
	}
	return 0;
}

int dummy(){
	int i,s;

	i = 0;
	while(1){
		if(verbosity > 7) printf("dummy: buffer %i status %i\n", i, ringbuffer[i].status);
		while(ringbuffer[i].status != STAT_DECIDE){
			usleep(WAITTIME);
// #pragma omp flush
		}
		omp_set_lock(&ringlocks[i]);
		ringbuffer[i].status = STAT_WRITE;
		s = ringbuffer[i].flag;
		ringbuffer[i].flag |= READ_KEEP;
		if(verbosity > 7) printf("dummy: buffer %i flag %i name %s\n", i, s, ringbuffer[i].name1); 
		omp_unset_lock(&ringlocks[i]);
		if(s == -1) return 0;
#ifdef HAVE_RINGBUFFER_MASK
		i = (i+1) & RINGBUFFER_MASK;
#else
		i = (i+1)%RINGBUFFER_SIZE;
#endif
	}
}
