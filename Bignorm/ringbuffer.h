// ringbuffer.h

#ifndef RINGBUFFER
#define RINGBUFFER

#include <stdint.h>
#include <zlib.h>
#include <omp.h>

#define RINGBUFFER_SIZE 32
#define HAVE_RINGBUFFER_MASK

#ifdef HAVE_RINGBUFFER_MASK
#define RINGBUFFER_MASK (RINGBUFFER_SIZE - 1)
#endif

#define R_MAX 300
#define T_MAX 15
#define Name_MAX 50
#define Comment_MAX 200

#define STAT_READER 0
#define STAT_TRIM 2
#define STAT_HASH 2
#define STAT_DECIDE 3
#define STAT_WRITE 4

#define READ_SINGLE 1
#define READ_PAIRED 2
#define READ_SHUFF 4
#define READ_KEEP 8
#define READ_REJ 16
#define IS_FASTQ 32

#define COUNTMAX 253

#define WAITTIME 100

struct ringbuffer_t
{
	int status, flag, filenr, rlen1, rlen2, hlen1, hlen2;
	char read1[R_MAX],read2[R_MAX],name1[Name_MAX],name2[Name_MAX],comment1[Comment_MAX],comment2[Comment_MAX];
	char qual1[R_MAX],qual2[R_MAX];
	uint8_t cmin1[R_MAX], cmin2[R_MAX];
	uint8_t count1[R_MAX][T_MAX], count2[R_MAX][T_MAX];
	uint64_t hash1[R_MAX][T_MAX], hash2[R_MAX][T_MAX];	
};

extern struct ringbuffer_t ringbuffer[RINGBUFFER_SIZE];
extern omp_lock_t ringlocks[RINGBUFFER_SIZE];

// functions around ringbuffer
int reader();
int hasher();
int dummy();
int writer();
int decider();

// CM
extern uint8_t *countmin;
extern uint64_t memneeded;
#endif
