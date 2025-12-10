#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "options.h"
#include "ringbuffer.h"
#include "stat.h"

int writer(){
	gzFile f1,f2, r1, r2;
	struct infiles *fs;
	int i, j, len, fnr, paired;
	char filename[150];

	fs = input_start;
	i = 0;
	fnr = -1;
	f1 = NULL; f2 = NULL;
	r1 = NULL; r2 = NULL;

	while(1){
		while(ringbuffer[i].status != STAT_WRITE){
			// wait for entry to become available
			usleep(WAITTIME);
			waits_writer++;
		}
		omp_set_lock(&ringlocks[i]);
		ringbuffer[i].status = STAT_READER;
		if(verbosity > 7) printf("writer: buffer %i status %i\n", i, ringbuffer[i].status);
		if(ringbuffer[i].flag == -1){
			// end of data
			omp_unset_lock(&ringlocks[i]);
			if(f1 != NULL) gzclose(f1);
			if(f2 != NULL) gzclose(f2);
			if(r1 != NULL) gzclose(r1);
			if(r2 != NULL) gzclose(r2);
			if(verbosity > 7) printf("writer: ending\n");
			return 0;
		}
		if(fnr != ringbuffer[i].filenr){
			// we have to open new file(s)
			if(f1 != NULL) gzclose(f1);
			if(f2 != NULL) gzclose(f2);
			if(r1 != NULL) gzclose(r1);
			if(r2 != NULL) gzclose(r2);
			if(ringbuffer[i].flag & 1){
				// unpaired
				paired = 0;
				if(verbosity > 5) printf("writer: opening output file %s_keep\n", fs->infile1);
				strncpy(filename, fs->infile1, 144);
				if(compressOut == 1){
					strcat(filename,"_keep.gz");
					f1 = gzopen(filename, "w");
				} else {
					strcat(filename,"_keep");
					f1 = gzopen(filename, "wT");
				}
				if(f1 == NULL){
					perror("error opening file for reads to keep");
					exit(-1);
				}		
				if(writerej == 1){
					strncpy(filename, fs->infile1, 144);
					if(compressOut == 1){
						strcat(filename,"_rej.gz");
						r1 = gzopen(filename, "w");
					} else {
						strcat(filename,"_rej");
						r1 = gzopen(filename, "wT");
					}
					if(r1 == NULL){
						perror("error opening file for rejected reads");
						exit(-1);
					}
				}
			} else { // paired
				paired = 1;
				if(verbosity > 5) printf("writer: opening output files %s_keep and %s_keep\n", fs->infile1, fs->infile2);
				if(compressOut == 1){
					strncpy(filename, fs->infile1, 144);
					strcat(filename,"_keep.gz");
					f1 = gzopen(filename, "w");
					strncpy(filename, fs->infile2, 144);
					strcat(filename,"_keep.gz");
					f2 = gzopen(filename, "w");
				} else {
					strncpy(filename, fs->infile1, 144);
					strcat(filename,"_keep");
					f1 = gzopen(filename, "wT");
					strncpy(filename, fs->infile2, 144);
					strcat(filename,"_keep");
                                        f2 = gzopen(filename, "wT");
				}
				if((f1 == NULL) || (f2 == NULL)){
					perror("error opening file for reads to keep");
					exit(-1);
				}
				if(writerej == 1){
					if(compressOut == 1){
						strncpy(filename, fs->infile1, 144);
						strcat(filename,"_rej.gz");
						r1 = gzopen(filename, "w");
						strncpy(filename, fs->infile2, 144);
						strcat(filename,"_rej.gz");
						r2 = gzopen(filename, "w");
					} else {
						strncpy(filename, fs->infile1, 144);
						strcat(filename,"_rej");
						r1 = gzopen(filename, "wT");
						strncpy(filename, fs->infile2, 144);
						strcat(filename,"_rej");
                                                r2 = gzopen(filename, "wT");
					}
					if((r1 == NULL) || (r2 == NULL)){
						perror("error opening file for rejected reads");
						exit(-1);
					}
				}
			} // else paired	
			fnr = ringbuffer[i].filenr;
			fs = fs->next;
		} // if(fnr != ringbuffer[i].filenr)
		if(ringbuffer[i].flag & READ_KEEP){
			if(paired == 0){
				if(verbosity > 7) printf("keeping %s\n", ringbuffer[i].name1);
				if((ringbuffer[i].flag & IS_FASTQ) && (opt_f == 0)){
					gzprintf(f1,"@%s", ringbuffer[i].name1);

				} else {
					gzprintf(f1,">%s", ringbuffer[i].name1);
				} 
				if(outputdetail > 0){
					// write details to file
					if(ringbuffer[i].rlen1 > k){
						len = ringbuffer[i].rlen1 - k + 1;
					} else {
						len = 1;
					}
					gzprintf(f1," cnt:%u",ringbuffer[i].cmin1[0]);
					for(j=1; j<len; j++){
						gzprintf(f1,";%u", ringbuffer[i].cmin1[j]);
					}
				}
				if(strnlen(ringbuffer[i].comment1,Comment_MAX) > 0)
					gzprintf(f1," %s",ringbuffer[i].comment1);
				if((ringbuffer[i].flag & IS_FASTQ) && (opt_f == 0)){
					gzprintf(f1,"\n%s\n+\n%s\n", ringbuffer[i].read1, ringbuffer[i].qual1);
				} else {
					gzprintf(f1,"\n%s\n", ringbuffer[i].read1);
				}
				count_reads_unpaired_out++;	
			} else {
				if(verbosity > 7) printf("keeping %s and %s\n", ringbuffer[i].name1, ringbuffer[i].name2);
				if((ringbuffer[i].flag & IS_FASTQ) && (opt_f == 0)){
					gzprintf(f1,"@%s", ringbuffer[i].name1);
					gzprintf(f2,"@%s", ringbuffer[i].name2);
				} else {
					gzprintf(f1,">%s", ringbuffer[i].name1); 
					gzprintf(f2,">%s", ringbuffer[i].name2);
				}
				if(outputdetail > 0){
					// write details to file
					if(ringbuffer[i].rlen1 > k){
						len = ringbuffer[i].rlen1 - k + 1;
					} else {
						len = 1;
					}
					gzprintf(f1," cnt:%u",ringbuffer[i].cmin1[0]);
					for(j=1; j<len; j++){
						gzprintf(f1,";%u", ringbuffer[i].cmin1[j]);
					}
					if(ringbuffer[i].rlen2 > k){
						len = ringbuffer[i].rlen2 - k + 1;
					} else {
						len = 1;
					}
					gzprintf(f2," cnt:%u",ringbuffer[i].cmin2[0]);
					for(j=1; j<len; j++){
						gzprintf(f2,";%u", ringbuffer[i].cmin2[j]);
					}
				}
				if(strnlen(ringbuffer[i].comment1,Comment_MAX) > 0)
					gzprintf(f1," %s",ringbuffer[i].comment1);
				if(strnlen(ringbuffer[i].comment2,Comment_MAX) > 0)
					gzprintf(f2," %s",ringbuffer[i].comment2);
				if((ringbuffer[i].flag & IS_FASTQ) && (opt_f == 0)){
					gzprintf(f1,"\n%s\n+\n%s\n", ringbuffer[i].read1, ringbuffer[i].qual1);
					gzprintf(f2,"\n%s\n+\n%s\n", ringbuffer[i].read2, ringbuffer[i].qual2);
				} else {				
					gzprintf(f1,"\n%s\n", ringbuffer[i].read1);
					gzprintf(f2,"\n%s\n", ringbuffer[i].read2);
				}
				count_reads_paired_out++;
			}
		} else if(writerej == 1){
			if(paired == 0){
				if(verbosity > 7) printf("rejecting %s\n", ringbuffer[i].name1);
				gzprintf(r1,">%s", ringbuffer[i].name1);
				if(outputdetail > 0){
					// write details to file
					if(ringbuffer[i].rlen1 > k){
						len = ringbuffer[i].rlen1 - k + 1;
					} else {
						len = 1;
					}
					gzprintf(r1," cnt:%u",ringbuffer[i].cmin1[0]);
					for(j=1; j<len; j++){
						gzprintf(r1,";%u", ringbuffer[i].cmin1[j]);
					}
				}
				if(strnlen(ringbuffer[i].comment1,Comment_MAX) > 0)
					gzprintf(r1," %s",ringbuffer[i].comment1);
				gzprintf(r1,"\n%s\n", ringbuffer[i].read1);
			} else {
				if(verbosity > 7) printf("rejecting %s and %s\n", ringbuffer[i].name1, ringbuffer[i].name2);
				gzprintf(r1,">%s", ringbuffer[i].name1);
				gzprintf(r2,">%s", ringbuffer[i].name2);
				if(outputdetail > 0){
					// write details to file
					if(ringbuffer[i].rlen1 > k){
						len = ringbuffer[i].rlen1 - k + 1;
					} else {
						len = 1;
					}
					gzprintf(r1," cnt:%u",ringbuffer[i].cmin1[0]);
					for(j=1; j<len; j++){
						gzprintf(r1,";%u", ringbuffer[i].cmin1[j]);
					}
					if(ringbuffer[i].rlen2 > k){
						len = ringbuffer[i].rlen2 - k + 1;
					} else {
						len = 1;
					}
					gzprintf(r2," cnt:%u",ringbuffer[i].cmin2[0]);
					for(j=1; j<len; j++){
						gzprintf(r2,";%u", ringbuffer[i].cmin2[j]);
					}
				}
				if(strnlen(ringbuffer[i].comment1,Comment_MAX) > 0)
					gzprintf(r1," %s",ringbuffer[i].comment1);
				if(strnlen(ringbuffer[i].comment2,Comment_MAX) > 0)
					gzprintf(r2," %s",ringbuffer[i].comment2);
				gzprintf(r1,"\n%s\n", ringbuffer[i].read1);
				gzprintf(r2,"\n%s\n", ringbuffer[i].read2);
			}
		}
		// memset(&ringbuffer[i], 0, sizeof(struct ringbuffer_t)); // implicit: status = STAT_READER
		// ringbuffer[i].status = STAT_READER;
		ringbuffer[i].flag = 0;
		omp_unset_lock(&ringlocks[i]);
#ifdef HAVE_RINGBUFFER_MASK
                i = (i+1) & RINGBUFFER_MASK;
#else
                i = (i+1)%RINGBUFFER_SIZE;
#endif

	} // while
} // writer

