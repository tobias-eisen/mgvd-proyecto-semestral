#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "options.h"
#include "ringbuffer.h"
#include "stat.h"

// thirdparty code start
uint8_t qselect(uint8_t *v, int len, int k)
{
#	define SWAP(a, b) { tmp = v[a]; v[a] = v[b]; v[b] = tmp; }
	int i, st;
	uint8_t tmp;
 
	for (st = i = 0; i < len - 1; i++) {
		if (v[i] > v[len-1]) continue;
		SWAP(i, st);
		st++;
	}
 
	SWAP(len-1, st);
 
	return k == st	?v[st]
			:st > k	? qselect(v, st, k)
				: qselect(v + st, len - st, k - st);
}
// thirdparty code end

uint8_t buffer[2 * R_MAX] = { 0 };

int decidefkt1(struct ringbuffer_t *rb){ // like diginorm: compare median of read counts to constant
	int j;

	if((rb->flag & 2) == 2){ // paired read
		j = rb->hlen1 + rb->hlen2;
		memcpy(buffer, rb->cmin1, rb->hlen1);
		memcpy((buffer + rb->hlen1), rb->cmin2, rb->hlen2);

		if(qselect(buffer, j, j/2) > C){
			return READ_REJ;
		} else {
			return READ_KEEP;
		}	
	} else { // single read
		j = rb->hlen1;
		if(qselect(rb->cmin1, j, j/2) > C){
			return READ_REJ;
		} else {
			return READ_KEEP;			 
		}
	}
}

int decidefkt2(struct ringbuffer_t *rb){ // compare first quartile of read counts to constant
	int j;

        if((rb->flag & 2) == 2){ // paired read
                j = rb->hlen1 + rb->hlen2;
                memcpy(buffer, rb->cmin1, rb->hlen1);
                memcpy((buffer + rb->hlen1), rb->cmin2, rb->hlen2);

                if(qselect(buffer, j, j/4) > C){
			 return READ_REJ;
                } else {
                        return READ_KEEP;
                }
        } else { // single read
                j = rb->hlen1;
                if(qselect(rb->cmin1, j, j/4) > C){
                        return READ_REJ;
                } else {
                        return READ_KEEP;
                }
        }
}

int decidefkt5(struct ringbuffer_t *rb){
	int j, count1, count2;

	count1 = 0;
	count2 = 0;
	if((rb->flag & 2) == 2){ // paired read
		for(j=0; j<rb->hlen2; j++){
			if(rb->cmin2[j] < opt_A){
				count1++;
			} else {
				if(rb->cmin2[j] < C){
					count2 ++;
				}
			}
		}	
	}	
	for(j=0; j<rb->hlen1; j++){
		if(rb->cmin1[j] < opt_A){
			count1++;
		} else {
			if(rb->cmin1[j] < C){
				count2 ++;
			}
		}
	}
	if((count1 > k) || (count2 > opt_B)){
		return READ_KEEP;
	} else {
		return READ_REJ;
	}
}

int decidefkt6(struct ringbuffer_t *rb){
	int j, j2, count1, count2, qmin;

	count1 = 0;
	count2 = 0;
	if((rb->flag & 2) == 2){ // paired read
		for(j=0; j<rb->hlen2; j++){
			if(rb->cmin2[j] < C){
				qmin = 0;
				for(j2 = j; j2 < j+k; j2++){ // check for phred score	
					if(rb->qual2[j2] < qbase + opt_q){ // score too low
						qmin = 1;
						break;
					}
				}
				if(qmin == 0){ // phred score good
					if(rb->cmin2[j] < opt_A){
						count1++;
					} else {
						count2 ++;
					}
				}
			}
		}
		if(count1> k) return READ_KEEP;
		count1 = 0;
	}
	for(j=0; j<rb->hlen1; j++){
		if(rb->cmin1[j] < C){
			qmin = 0;
			for(j2 = j; j2 < j+k; j2++){ // check for phred score
				if(rb->qual1[j2] < qbase + opt_q){ // score too low
					qmin = 1;
					break;
				}
			}
			if(qmin == 0){ // phred score good
				if(rb->cmin1[j] < opt_A){
					count1++;
				} else {
					count2 ++;
				}
			}
		}
	}
	if((count1 > k) || (count2 > opt_B)){
		return READ_KEEP;
	} else {
		return READ_REJ;
	}
}
	
int decidefkt3(struct ringbuffer_t *rb){
	int j, count1, count2;

	count1 = 0;
	count2 = 0;

	for(j=0; j<rb->hlen1; j++){
		if(rb->cmin1[j] < opt_A){
			count1++;
		} else {
			if(rb->cmin1[j] < C){
				count2 ++;
			}
		}
	}
	if((count1 > k) || (count2 > opt_B)){
		if((rb->flag & 2) == 2){ // paired read
			count1 = 0;
			count2 = 0;
			for(j=0; j<rb->hlen2; j++){
				if(rb->cmin2[j] < opt_A){
					count1++;
				} else {
					if(rb->cmin2[j] < C){
						count2 ++;
					}
				}
			}
			if((count1 > k) || (count2 > opt_B)){
				return READ_KEEP;
			} else {
				return READ_REJ;
			}
		} else {
			return READ_KEEP;
		}
	} else {
		return READ_REJ;
	}
}

int decidefkt4(struct ringbuffer_t *rb){
	int j, count1, count2;

	count1 = 0;
	count2 = 0;

	for(j=0; j<rb->hlen1; j++){
                if(rb->cmin1[j] < opt_A){
                        count1++;
                } else {
                        if(rb->cmin1[j] < C){
                                count2 ++;
                        }
                }
        }
        if((count1 > k) || (count2 > opt_B)){
		return READ_KEEP;
	} else {
		if((rb->flag & 2) == 2){ // paired read
			count1 = 0;
			count2 = 0;
			for(j=0; j<rb->hlen2; j++){
                                if(rb->cmin2[j] < opt_A){
                                        count1++;
                                } else {
                                        if(rb->cmin2[j] < C){
                                                count2 ++;
                                        }
                                }
                        }
                        if((count1 > k) || (count2 > opt_B)){
				return READ_KEEP;
                        } else {
                                return READ_REJ;
                        }
		} else {
			return READ_REJ;
		}
	}
}


typedef int (*ptDecide)(struct ringbuffer_t *rb);

#define N_DECIDERS 6

ptDecide the_deciders[] = { &decidefkt1, &decidefkt2, &decidefkt3, &decidefkt4, &decidefkt5, &decidefkt6};

int decider(){
	int i, s, j, akt, l, ri;
	uint8_t  min;
	uint64_t m2, kmer_mask, maskR1, maskR2;

	ptDecide do_decide = NULL;	
	i = 0;

	if(opt_n > 0){
		kmer_mask = 0;
		for(j=0; j<k; j++){
			kmer_mask <<= 2;
			kmer_mask |= 3;
		}
	}
	// init: set decider function to use
	if(((d -1) > N_DECIDERS) || (d < 1)){
		printf("decider: unknown strategy for deciding given (-d %i), must be between 1 and N_DECIDERS\n", d);
		exit(-1);
	}
	do_decide = the_deciders[d - 1];
	while(m == 0)
		usleep(WAITTIME);
	m2 = m * 1024 * 1024;

	// work 
	while(1){
		while(ringbuffer[i].status != STAT_DECIDE){
			usleep(WAITTIME);
			waits_decider++;
		}
		omp_set_lock(&ringlocks[i]);
		ringbuffer[i].status = STAT_WRITE;
		s = ringbuffer[i].flag;
		if(s == -1){
			// done
			omp_unset_lock(&ringlocks[i]);
			return 0;
		} else {
			if((s & READ_REJ) == READ_REJ){
				// nothing to do...
				omp_unset_lock(&ringlocks[i]);
			} else {
				// get CM counts
				for(j=0; j<ringbuffer[i].hlen1; j++){
					ringbuffer[i].count1[j][0] = countmin[ringbuffer[i].hash1[j][0]];
					min = ringbuffer[i].count1[j][0];
					for(l=1; l<t; l++){
						ringbuffer[i].count1[j][l] = countmin[ringbuffer[i].hash1[j][l] + l * m2];
						if(min > ringbuffer[i].count1[j][l])
							min = ringbuffer[i].count1[j][l];
					}
					ringbuffer[i].cmin1[j] = min;
				}
				if((s & 2) == 2){ // paired read
					for(j=0; j<ringbuffer[i].hlen2; j++){
						ringbuffer[i].count2[j][0] = countmin[ringbuffer[i].hash2[j][0]];
						min = ringbuffer[i].count2[j][0];
						for(l=1; l<t; l++){
							ringbuffer[i].count2[j][l] = countmin[ringbuffer[i].hash2[j][l] + l * m2];
							if(min > ringbuffer[i].count2[j][l])
								min = ringbuffer[i].count2[j][l];
						}
						ringbuffer[i].cmin2[j] = min;
					}
				}
				akt = do_decide(&ringbuffer[i]);
				ringbuffer[i].flag |= akt;
				if(verbosity > 8){
					if(akt == READ_KEEP){
						printf("decider: accepting read %s\n", ringbuffer[i].read1);
					} else {
						printf("decider: rejecting read %s\n", ringbuffer[i].read1);
					}
				} 
				if(akt == READ_KEEP){ // update CM
					if(opt_n == 0){ // how to treat reads with N: opt_n == 0: no special treatment
						for(j=0; j<ringbuffer[i].hlen1; j++){
							if(ringbuffer[i].cmin1[j] == 0) count_newkmers++;
							for(l=0; l<t; l++){
								if(ringbuffer[i].count1[j][l] < COUNTMAX){
									countmin[ringbuffer[i].hash1[j][l] + l * m2] = ringbuffer[i].count1[j][l] +1;
								}
							}
						}	
						if((s & 2) == 2){ // paired read
							for(j=0; j<ringbuffer[i].hlen2; j++){
								if(ringbuffer[i].cmin2[j] == 0) count_newkmers++;
								for(l=0; l<t; l++){
									if(ringbuffer[i].count2[j][l] < COUNTMAX){
										countmin[ringbuffer[i].hash2[j][l] + l * m2] = ringbuffer[i].count2[j][l] +1;
									}
								}
							}
						}
					} else { //special treatment: for now, don't increse counter for k-mer including N
						maskR1 = 0;
						for(j=0; j< k; j++){
							maskR1 <<= 2;					
							if(ringbuffer[i].read1[j] == 'N'){
								maskR1 |= 3;
							}
						}
						maskR1 &= kmer_mask;
						ri = k;
						for(j=0; j<ringbuffer[i].hlen1; j++){
							if(maskR1 != 0){ // k-mer includes at least one N, don't increase counter
								count_kmers_with_N++;
							} else {
								if(ringbuffer[i].cmin1[j] == 0) count_newkmers++;
								for(l=0; l<t; l++){		
									if(ringbuffer[i].count1[j][l] < COUNTMAX){
										countmin[ringbuffer[i].hash1[j][l] + l * m2] = ringbuffer[i].count1[j][l] +1;
									}
								}
							}
							maskR1 <<= 2;
							if(ringbuffer[i].read1[ri] == 'N'){
								maskR1 |= 3;
							}
							maskR1 &= kmer_mask;
							ri++;
						}
						if((s & 2) == 2){ // paired read
							maskR1 = 0;
							for(j=0; j< k; j++){
								maskR1 <<= 2;
								if(ringbuffer[i].read2[j] == 'N'){
									maskR1 |= 3;
								}
							}
							maskR1 &= kmer_mask;
							ri = k;
							for(j=0; j<ringbuffer[i].hlen2; j++){
								if(maskR1 != 0){ // k-mer includes at least one N, don't increase counter
									count_kmers_with_N++;
								} else {
									if(ringbuffer[i].cmin2[j] == 0) count_newkmers++;
									for(l=0; l<t; l++){
										if(ringbuffer[i].count2[j][l] < COUNTMAX){
											countmin[ringbuffer[i].hash2[j][l] + l * m2] = ringbuffer[i].count2[j][l] + 1;
										}
									}
								}
								maskR1 <<= 2;
								if(ringbuffer[i].read1[ri] == 'N'){
									maskR1 |= 3;
								}
								maskR1 &= kmer_mask;
								ri++;
							}
						}
					}
				}
				omp_unset_lock(&ringlocks[i]);
			}
		} // else s==-1
#ifdef HAVE_RINGBUFFER_MASK
                i = (i+1) & RINGBUFFER_MASK;
#else
                i = (i+1)%RINGBUFFER_SIZE;
#endif
	} // while(1)
} //decider	
