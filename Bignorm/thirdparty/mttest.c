#include <stdio.h>
#include <time.h>
#include "mt64.h"


int main(){
	uint64_t a,min;
	uint32_t i, count;
	struct timespec ts;

	count = 0;
	clock_gettime(0, &ts);
	init_genrand64(ts.tv_nsec);
	min = genrand64_int64();

	for(i=0; i < 300000000; i++){
		a = genrand64_int64();
		if(a < min) min = a;
		if((a >> 32) > 0) count++;
	}
	printf("# > 2^32: %i (%.2f %), min: %u %u\n", count, ((float)count / 3000000), min >> 32, (uint32_t)min);
}
