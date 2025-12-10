#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <zlib.h>

#include "options.h"
#include "ringbuffer.h"
#include "stat.h"

#include "thirdparty/mt64.h"
#include "thirdparty/sysinfo.h"

// macOS compatibility: _SC_AVPHYS_PAGES is not available on macOS
#ifdef __APPLE__
#ifndef _SC_AVPHYS_PAGES
#define _SC_AVPHYS_PAGES _SC_PHYS_PAGES
#endif
#endif


uint64_t *A, *B;

uint64_t A_233094834[] = { (__UINT64_C(1517733481)),
        (__UINT64_C(5428660835479589707)),(__UINT64_C(4846634552415376859)),(__UINT64_C(6404847846195106823)),
        (__UINT64_C(7082504985027665327)),(__UINT64_C(7081795071714822979)),(__UINT64_C(9866633560969882985)),
        (__UINT64_C(14710122096630180867)),(__UINT64_C(4083536418104953313)),(__UINT64_C(10085560887037081519)),
        (__UINT64_C(9657684599635430077)),(__UINT64_C(5375765776848166897)),(__UINT64_C(2360808972908299819)),
        (__UINT64_C(7250681630747551001)),(__UINT64_C(11003473074184913221)),(__UINT64_C(8076466398770329967)),
        (__UINT64_C(6139719934185894223)),(__UINT64_C(11451049086794637917)),(__UINT64_C(8619122089797105967)),
        (__UINT64_C(13390601433041719853)),(__UINT64_C(8839741094390848361)),(__UINT64_C(10577078423499430681)),
        (__UINT64_C(6183755860004429897)),(__UINT64_C(3025424482010352607)),(__UINT64_C(17004177928310018557))};

uint64_t B_233094834[] = {(__UINT64_C(4856492158556586560)),
        (__UINT64_C(18150631066297409692)),(__UINT64_C(17071671312143097890)),(__UINT64_C(11246878191574502858)),
        (__UINT64_C(17338235642782036712)),(__UINT64_C(1727458701145891333)),(__UINT64_C(11294910387947023269)),
        (__UINT64_C(13949297806502574110)),(__UINT64_C(17272752438476092944)),(__UINT64_C(3248107074795813634)),
        (__UINT64_C(5922537472645189238)),(__UINT64_C(15892286188617372710)),(__UINT64_C(8318120239316119536)),
        (__UINT64_C(16343716287248999037)),(__UINT64_C(14333795491445113305)),(__UINT64_C(13398482153519806489)),
        (__UINT64_C(1409391843040956922)),(__UINT64_C(66510355693317130)),(__UINT64_C(3013922286997948648)),
        (__UINT64_C(15330831926743536898)),(__UINT64_C(6826576085057777082)),(__UINT64_C(17985603312190364493)),
        (__UINT64_C(11279862377769237622)),(__UINT64_C(13852582901409633778)),(__UINT64_C(18184946659315698415))};

uint64_t A_343559959[] = { (__UINT64_C(187121373)),
        (__UINT64_C(2511431304343562723)),(__UINT64_C(7147665428436345271)),(__UINT64_C(3721494656211659017)),
        (__UINT64_C(9028531757532806933)),(__UINT64_C(12401019029499314951)),(__UINT64_C(16080865627514826833)),
        (__UINT64_C(402776020804031653)),(__UINT64_C(8610805722305856557)),(__UINT64_C(12952807613372483653)),
        (__UINT64_C(15581649228071459537)),(__UINT64_C(14520345310724849129)),(__UINT64_C(12111067089012895411)),
        (__UINT64_C(18066810858103403623)),(__UINT64_C(2872946235864730241)),(__UINT64_C(15638625287209411079)),
        (__UINT64_C(10925177132098298387)),(__UINT64_C(7581402507076665979)),(__UINT64_C(4557512906746496033)),
        (__UINT64_C(16934785145779532849)),(__UINT64_C(2116719500795863651)),(__UINT64_C(9822085442549415677)),
        (__UINT64_C(6870107060058843499)),(__UINT64_C(10743958474768773521)),(__UINT64_C(1298506336925404333))};

uint64_t B_343559959[] = {(__UINT64_C(8971188853532042230)),
        (__UINT64_C(6325182366327784716)),(__UINT64_C(6275591155303964276)),(__UINT64_C(14715390152486958894)),
        (__UINT64_C(12702761475232146710)),(__UINT64_C(15266937932100600371)),(__UINT64_C(2408571226227901591)),
        (__UINT64_C(17325463332754345070)),(__UINT64_C(11184575537372294475)),(__UINT64_C(6060057693046191915)),
        (__UINT64_C(2351326338454614533)),(__UINT64_C(13790566351358881485)),(__UINT64_C(13542537808996347883)),
        (__UINT64_C(1139420903880027975)),(__UINT64_C(13684739206904650925)),(__UINT64_C(2647282317481559886)),
        (__UINT64_C(6641614377480153390)),(__UINT64_C(17277963026684216469)),(__UINT64_C(5711401825189432063)),
        (__UINT64_C(6722284676151193032)),(__UINT64_C(1607587368891218311)),(__UINT64_C(5054120436284156439)),
        (__UINT64_C(3251466200192997033)),(__UINT64_C(12985719191298216096)),(__UINT64_C(469490828991440354))};

uint64_t A_837778944[] = { (__UINT64_C(3440564993)),
        (__UINT64_C(6208183799242126737)),(__UINT64_C(13757271343481908477)),(__UINT64_C(8945801451070341587)),
        (__UINT64_C(16007792826405394403)),(__UINT64_C(1027394225017363123)),(__UINT64_C(3525937044538061569)),
        (__UINT64_C(6831725399721252985)),(__UINT64_C(6100030754001587537)),(__UINT64_C(15565666264822266667)),
        (__UINT64_C(13403580988819414013)),(__UINT64_C(16388442194190681527)),(__UINT64_C(15439900043822301121)),
        (__UINT64_C(14730872969127036433)),(__UINT64_C(17507893877636372741)),(__UINT64_C(17628539803969626061)),
        (__UINT64_C(16491300322101441493)),(__UINT64_C(11030294364628018399)),(__UINT64_C(10281779409445806959)),
        (__UINT64_C(12920252808930493277)),(__UINT64_C(16082838486655860947)),(__UINT64_C(16377602800963860707)),
        (__UINT64_C(10227175036805121253)),(__UINT64_C(2829561657770814737)),(__UINT64_C(4320018602004207667))};

uint64_t B_837778944[] = {(__UINT64_C(12589259896261028943)),
        (__UINT64_C(5118218068041137432)),(__UINT64_C(15881037847873053002)),(__UINT64_C(15057295301479589376)),
        (__UINT64_C(15847630078235206640)),(__UINT64_C(10318647442273996474)),(__UINT64_C(14965218606940070942)),
        (__UINT64_C(11687556776720255862)),(__UINT64_C(12692573309166735450)),(__UINT64_C(16337151505370829493)),
        (__UINT64_C(10965754989558595602)),(__UINT64_C(16687752933887670036)),(__UINT64_C(8562896084281220660)),
        (__UINT64_C(18430352550562250709)),(__UINT64_C(3716204860066294188)),(__UINT64_C(4498041028110365570)),
        (__UINT64_C(4835626070681416229)),(__UINT64_C(10521089941020280573)),(__UINT64_C(13893041809776693739)),
        (__UINT64_C(4029776852375354148)),(__UINT64_C(209526316864103497)),(__UINT64_C(14465039513016333756)),
        (__UINT64_C(14663614034840330192)),(__UINT64_C(14639539517193443665)),(__UINT64_C(1534251652891572490))};

uint64_t gcd(uint64_t a, uint64_t b){
	for(;;){
		if(a == 0)
			return b;
		b %= a;
		if(b == 0)
			return a;
		a %= b;
	}
}

int relprime(uint64_t t, uint64_t *a, int n){
	int i;

	for(i=0; i<n; i++){
		if(gcd(t,a[i]) != 1)
			return -1;
	}
	return 0;
}


int randomset(int n, uint64_t *A, uint64_t *B){
	int i;
	uint64_t t;
	struct timespec ts;

	clock_gettime(0, &ts);
	init_genrand64(ts.tv_nsec); // use nsec part of time as seed

	A[0] = genrand64_int64();   // first value for A
	A[0] |= 1;		    // ... must be odd

	for(i=1; i < n; i++){       // fill array a with n relative prime values
		t = genrand64_int64();
		t |= 1;
		while(relprime(t,A,i) == -1){
			t = genrand64_int64();
			t |= 1;
		}
		A[i] = t;
	}
	for(i=0; i < n; i++){
		B[i] = genrand64_int64(); // no constraints on B
	}
	return 0;
}

int mTO2(int m_in){
	// the m to use has to be m = 2^n for sone n
	uint64_t i, j, k;

	k = 32;
	i = 1;
	j = (uint64_t)m_in;
	while((j > 0) && (k > 0)){
		i <<= 1;
		j >>= 1;
		k--;
	}
	i >>= 1;
	return (int)i;
}


uint64_t maxmem; // max memory to use in KB
int8_t shiftcount;

void calc_mt(){
	long temp1, temp2;
	int m_temp, t_temp;

	if((t > 0) && (m > 0)){
		// user supplied t and m
		m_temp = m;
                m = mTO2(m_temp);
                if(m != m_temp){
         	       printf("changing m from %i to %i\n", m_temp, m);
                }

		maxmem = t * m * 1024;
		if(verbosity > 3){
			printf("using user supplied t (%i) and m (%i), %i MB of RAM\n", t, m, t*m);
		}
	} else {
		if(s == 0){
			// use mem_frac GB of RAM
			maxmem = (uint64_t) (mem_frac * 1024 * 1024);
			if(verbosity > 3){
				printf("using user supplied memory limif of %.1f GB\n", mem_frac);
			}
		} else {
			// automatic sizing - get RAM size / usage of system
			// first: try to use /proc
			kb_main_total = 0;
			kb_main_free = 0;
			kb_main_cached = 0;
			kb_main_buffers = 0;
			meminfo(); // use procps-ng-3.3.10 code to get free mem

			if(kb_main_total == 0){
				printf("error reading memory information using /proc, trying sysconf\n");
				temp1 = sysconf( _SC_PHYS_PAGES );
				temp2 = sysconf( _SC_PAGESIZE );
				if((temp1 <= 0) || (temp2 <= 0)){
					printf("error reading memory information using sysconf\nPlease provide memory limit os t and m values on commandline\n");
					exit(-1);
				} else {
					maxmem = (uint64_t)((float)temp1 * temp2 * mem_frac);
					if((s == 1) || (s == 2)){
						temp1 = sysconf(_SC_AVPHYS_PAGES);
						if(temp1 <= 0){
							printf("error reading available memory, using total memory instead\n");
						} else {
							maxmem = (uint64_t)((float)temp1 * temp2 * mem_frac);
						}
					}
					if(verbosity > 3){
						printf("using memory limif of %.1f GB based on sysconf\n", (((float)maxmem) /(1024 * 1024 * 1024)));
					}
				}
			} else {
				if(verbosity > 4){
					printf("memory information using /proc:\n\tmain_total (MB): %u\n", kb_main_total / 1024);
					printf("\tfree: %u (%.1f \%)\n\tcached: %u (%.1f \%)\n\tbuffers: %u (%.1f \%)\n", 
						kb_main_free / 1024, (float)kb_main_free / (float)kb_main_total * 100, 
						kb_main_cached / 1024, (float)kb_main_cached / (float)kb_main_total * 100,
						kb_main_buffers / 1024, (float)kb_main_buffers / (float)kb_main_total * 100);
				}
				switch(s){
					case 1: // default: use (kb_main_free + kb_main_cached + kb_main_buffers) * mem_frac
						maxmem = (uint64_t)((kb_main_free + kb_main_cached + kb_main_buffers) * mem_frac);
						break;
					case 2: // use kb_main_free * mem_frac
						maxmem = (uint64_t)(kb_main_free * mem_frac);
						break;
					case 3: // use kb_main_total * mem_frac
						maxmem = (uint64_t)(kb_main_total * mem_frac);
						break;
					default:
						printf("error: unknown strategy for memory allocation given (-s %i)\n", s);
						exit(-1);
				}
				if(verbosity > 3){
					printf("using a memory limit of %.1f GB based on /proc\n", (((float)maxmem) /(1024 * 1024)));
				}
			}
		}
		// now we got maxmem and need to find m and t
		if(m > 0){ // user supplied
			m_temp = m;
			m = mTO2(m);
			if(m != m_temp){
				printf("changing m from %i to %i\n", m_temp, m);
			}
			t = maxmem / (m * 1024);
			if(t < 4){
				printf("error: given m to big for available memory (t = %i)\n", t);
				exit(-1);
			}
			if(t < 7){
				printf("warning: t is very small (%i)\n", t);
			}
			if(t > 10) // t > 10 normally doesn't make sense
				t = 10;
		} else {
			if(t > 0){
				m = mTO2((int)maxmem / (t * 1024));
				if(m < 50){
					printf("error: given t to big for available memory (m = %i)\n", m);
					exit(-1);
				}
			} else { // nothing supplied
				t_temp = 10;
				m_temp = 256;
				if(((uint64_t)t_temp * m_temp) * 1024 > maxmem){
					// find a smaller solution
					while(((uint64_t)t_temp * m_temp) * 512 >= maxmem)
						m_temp = m_temp / 2;
					t_temp = maxmem / (m_temp * 1024);
				} else {
					// find a bigger solution
					while((((uint64_t)t_temp * m_temp) * 2048 <= maxmem) && (m_temp <= 1024 * 8))
						m_temp = m_temp * 2;
					t_temp = maxmem / (m_temp * 1024);
					if(t_temp > 10)
						t_temp = 10;
				}
				t = t_temp;
				if(m_temp < 4096){ // using more than 4GB per CMS hash normally doesn't make sense
					m = m_temp;
				} else {
					m = 4096;
				}
			}
		}
		if(verbosity > 3){
			printf("using m = %i MB, t = %i, resulting in %.1f GB memory usage\n", m, t, (((float)m * t))/1024);
		}
	}
}

void print_k32(uint64_t kmer){
	int i;
	char seq[33];

	seq[32]= '\0';

	for(i=31; i>=0; i--){
		switch(kmer & 3){
			case 0:
				seq[i] = 'A';
				break;
			case 1:
				seq[i] = 'C';
				break;
			case 2:
				seq[i] = 'G';
				break;
			case 3:
				seq[i] = 'T';
		}
		kmer >>= 2;
	}
	printf("%s ", seq);
}

uint64_t memneeded;
uint8_t *countmin;

inline void do_hash(uint64_t X, int index, int read, struct ringbuffer_t *rb){
	uint8_t i;
	uint64_t hash;

	count_kmers++;
	if(read == 1){
		hash = A[0] * X + B[0];
		rb->hash1[index][0] = hash >> shiftcount;
		for(i=1; i<t; i++){
			hash = A[i] * X + B[i];
			rb->hash1[index][i] = hash >> shiftcount;
		}
	} else {
		hash = A[0] * X + B[0];
		rb->hash2[index][0] = hash >> shiftcount;
		for(i=1; i<t; i++){
			hash = A[i] * X + B[i];
			rb->hash2[index][i] = hash >> shiftcount;
		}

	 }
}


int hasher(){
	int i, j, k2, s, rnum, shiftk, b_mcount;
	struct timespec ts1, ts2;
	uint64_t kmer, kmerrev, kmer_mask, kmertemp, b_mask;
	gzFile df;

	kmer = 0;
	kmerrev = 0;
	maxmem = 0;
	// calc m and t for CM-sketch
	calc_mt();
	// tell other tasks
#pragma omp flush(m,t)
	//shiftcount = 64 - (m + 20); // m in MB, 1 MB = 2^20
	i = 19;
	kmertemp = (uint64_t)m;
	while(kmertemp != 0){
		kmertemp >>= 1;
		i++;
	}
	shiftcount = 64 - (uint8_t)i;
	shiftk = 2 * k - 2;
	b_mcount = (2 * k) - (uint8_t)i;
printf("i: %i  b_mcount: %i\n", i, b_mcount);
	// choose A and B vectors for hashing
	switch(h) {
		case 0:
			// choose random hash set
			A = malloc(t * sizeof(uint64_t));
			if(A == NULL){
				perror("error allocating memory for A");
				exit(-1);
			}
			B = malloc(t * sizeof(uint64_t));
			if(B == NULL){
				perror("error allocating memory for B");
				exit(-1);
			}
			randomset(t,A,B);
			break;
		case 1:
			A = A_233094834;
			B = B_233094834;
			break;
		case 2:
			A = A_343559959;
			B = B_343559959;
			break;
		case 3:
			A = A_837778944;
			B = B_837778944;
                        break;
		} // switch

	// calc RAM needed and alloc it
	memneeded = (uint64_t)t;
	memneeded *= (uint64_t)m;
	memneeded *= 1024 * 1024;

	if(memneeded > SIZE_MAX){
		printf("hasher: memneeded > SIZE_MAX\n");
		exit(-1);
	} else {
		countmin = malloc(memneeded);
		if(countmin == NULL){
			perror("hasher: malloc countmin");
			exit(-1);
		}
		if(verbosity > 3){
			printf("hasher: got %u bytes (%.1f GB)\n", memneeded, ((float)memneeded) /(1024 * 1024 * 1024));
		}
		clock_gettime(0, &ts1);
		memset(countmin, 0, memneeded);
		clock_gettime(0, &ts2);
		if(verbosity > 5){
			printf("hasher: memset needed %i seconds\n", ts2.tv_sec - ts1.tv_sec);
		}
		i = 0;
		// generate kmer_mask
		kmer_mask = 0;
		for(j=0; j<k; j++){
			kmer_mask <<= 2;
			kmer_mask |= 3;
		}
		b_mask = 0;
		for(j=0; j< b_mcount; j++){
			b_mask <<= 1;
			b_mask |= 1;
		}
		for(j=0; j< 25; j++){
			B[j] &= b_mask;
		}
		if((verbosity > 3) ||((verbosity > 2) && (h == 0))){
				printf("b_mask: %lu\n", b_mask);
                                printf("hasher: using these values for hashing:\nA\t\tB\n");
                                for(j=0; j<t; j++){
                                        printf("%lu\t%lu\n", A[j], B[j]);
                                }
                }

		// work ...
		while(1){
			while(ringbuffer[i].status != STAT_HASH){
				usleep(WAITTIME);
				waits_hasher++;
			}
			omp_set_lock(&ringlocks[i]);
			ringbuffer[i].status = STAT_DECIDE;
			s = ringbuffer[i].flag;
			if((s > 0) && ((s & READ_REJ) == 0)){
				// new read: init kmer and kmerrev
				kmer = 0;
				kmerrev = 0;
				rnum = 0;
				k2 = min(k, ringbuffer[i].rlen1);
				for(j=0; j<k2; j++){
					switch(ringbuffer[i].read1[j]){
						case 'a':
						case 'A':
						case 'n':
						case 'N':
							kmer <<= 2;
							break;
						case 'c':
						case 'C':
							kmer <<= 2;
							kmer |= 1;
							break;
						case 'g':
						case 'G':
							kmer <<= 2;
							kmer |= 2;
							break;
						case 't':
						case 'T':
							kmer <<= 2;
							kmer |= 3;
							break;
						default:
							printf("hasher: unexpected nucleotide %c\n", ringbuffer[i].read1[j]);
							kmer <<= 2;
					} // switch
				} // for
				kmertemp = ~kmer;
				for(j=0; j<k2; j++){
					kmerrev <<= 2;
					kmerrev |= (kmertemp & 3);
					kmertemp >>= 2;
				}
				if(verbosity > 8){
					printf("hasher: kmer = %x\t kmerrev = %x\n", kmer, kmerrev);
					print_k32(kmer);print_k32(kmerrev);
					printf("\t%s\n",ringbuffer[i].read1);
				}
				if(kmer < kmerrev){ // min or max?
					do_hash(kmer, rnum, 1, &ringbuffer[i]);
				} else {
					do_hash(kmerrev, rnum, 1, &ringbuffer[i]);
				}
				rnum++;
				for(j=k; j<ringbuffer[i].rlen1; j++){
					switch(ringbuffer[i].read1[j]){
                                        	case 'a':
                                        	case 'A':
                                        	case 'n':
                                        	case 'N':
                                                	kmer <<= 2;
							kmerrev >>= 2;
							kmerrev |= (3UL << shiftk);
                                               		break;
                                        	case 'c':
                                        	case 'C':
                                                	kmer <<= 2;
                                                	kmer |= 1;
							kmerrev >>= 2;
							kmerrev |= (2UL << shiftk);
                                                	break;
                                        	case 'g':
                                        	case 'G':
                                                	kmer <<= 2;
                                                	kmer |= 2;
							kmerrev >>= 2;
							kmerrev |= (1UL << shiftk);
                                                	break;
                                        	case 't':
                                        	case 'T':
                                                	kmer <<= 2;
                                                	kmer |= 3;
							kmerrev >>= 2;
                                                	break;
                                        	default:
                                                	printf("hasher: unexpected nucleotide %c\n", ringbuffer[i].read1[j]);
                                                	kmer <<= 2;
							kmerrev >>= 2;				
					} // switch
					kmer &= kmer_mask;
					kmerrev &= kmer_mask;
					if(verbosity > 8){
	                                        printf("hasher: kmer = %x\t kmerrev = %x\n", kmer, kmerrev);
        	                                print_k32(kmer);print_k32(kmerrev);
                        	        }

					if(kmer > kmerrev){
						do_hash(kmer, rnum, 1, &ringbuffer[i]);
					} else {
						do_hash(kmerrev, rnum, 1, &ringbuffer[i]);
					}
					rnum++;
				} // for
				ringbuffer[i].hlen1 = rnum;
				// paired?
				if((ringbuffer[i].flag & 2) ==  2){
					// the same as above for second part
					kmer = 0;
					kmerrev = 0;
					rnum=0;
					k2 = min(k, ringbuffer[i].rlen2);
					for(j=0; j<k2; j++){
						switch(ringbuffer[i].read2[j]){
							case 'a':
							case 'A':
							case 'n':
							case 'N':
								kmer <<= 2;
								break;
							case 'c':
							case 'C':
								kmer <<= 2;
								kmer |= 1;
								break;
							case 'g':
							case 'G':
								kmer <<= 2;
								kmer |= 2;
								break;
							case 't':
							case 'T':
								kmer <<= 2;
								kmer |= 3;
								break;
							default:
								printf("hasher: unexpected nucleotide %c\n", ringbuffer[i].read2[j]);
								kmer <<= 2;
						} // switch
					} // for
					kmertemp = ~kmer;
					for(j=0; j<k2; j++){
						kmerrev <<= 2;
						kmerrev |= (kmertemp & 3);
						kmertemp >>= 2;
					}
					if(verbosity > 8){
						printf("hasher: kmer = %x\t kmerrev = %x\n", kmer, kmerrev);
					}
					if(kmer > kmerrev){
						do_hash(kmer, rnum, 2, &ringbuffer[i]);
					} else {
						do_hash(kmerrev, rnum, 2, &ringbuffer[i]);
					}
					rnum++;
					for(j=k; j<ringbuffer[i].rlen2; j++){
						switch(ringbuffer[i].read2[j]){
                                        		case 'a':
                                        		case 'A':
                                        		case 'n':
                                        		case 'N':
                                                		kmer <<= 2;
								kmerrev >>= 2;
								kmerrev |= (3UL << shiftk);
                                               			break;
                                        		case 'c':
                                        		case 'C':
                                                		kmer <<= 2;
                                                		kmer |= 1;
								kmerrev >>= 2;
								kmerrev |= (2UL << shiftk);
                                                		break;
                                        		case 'g':
                                        		case 'G':
								kmer <<= 2;
                                         		       	kmer |= 2;
								kmerrev >>= 2;
								kmerrev |= (1UL << shiftk);
                                                		break;
                                        		case 't':
                                        		case 'T':
                                                		kmer <<= 2;
                                                		kmer |= 3;
								kmerrev >>= 2;
                                                		break;
                                        		default:
                                                		printf("hasher: unexpected nucleotide %c\n", ringbuffer[i].read2[j]);
                                                		kmer <<= 2;
								kmerrev >>= 2;				
						} // switch
						kmer &= kmer_mask;
						kmerrev &= kmer_mask;	
						if(kmer > kmerrev){
							do_hash(kmer, rnum, 2, &ringbuffer[i]);
						} else {
							do_hash(kmerrev, rnum, 2, &ringbuffer[i]);
						}
						rnum++;
					} // for
					ringbuffer[i].hlen2 = rnum;
	

				} // paired

				omp_unset_lock(&ringlocks[i]);
			} else { // s <= 0
				if(s == -1){ // done
					// save CM?
					if(dumpfile == 1){
						df = gzopen(dumpfilename, "w");
						if(df== NULL){
							printf("error opening dumpfile %s, skipping dump\n", dumpfilename);
						} else {
							gzprintf(df,"CMdump\n%u\n", t * m);
							for(j=0; j < t*m; j++){ // write 1 MB per call
								if(gzwrite(df, &countmin[j * 1024 * 1024], 1024 * 1024) == 0){
									perror("hasher: error writing dump");
									break;
								}
							}
							gzclose(df);	
						}

					}
					omp_unset_lock(&ringlocks[i]);
					return 0;
				} else {
					if((s & READ_REJ) == READ_REJ){
						// read already marked as trash, nothing to do
						omp_unset_lock(&ringlocks[i]);
					} else {
						printf("hasher: unknown ringbuffer state %i\n", s);
						omp_unset_lock(&ringlocks[i]);
					}
				}
			} // else		
#ifdef HAVE_RINGBUFFER_MASK
			i = (i+1) & RINGBUFFER_MASK;
#else
			i = (i+1)%RINGBUFFER_SIZE;
#endif
		}
	}
}	
