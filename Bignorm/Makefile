Bignorm: main.o options.o reader.o writer.o stat.o sysinfo.o hasher.o decider.o mt19937-64.o
	gcc -o Bignorm  main.o options.o reader.o writer.o stat.o sysinfo.o hasher.o decider.o mt19937-64.o -lz -fopenmp 

clean:
	rm -f *.o Bignorm
main.o: main.c options.h ringbuffer.h stat.h
	gcc -c main.c -fopenmp	-O2

options.o: options.c options.h stat.h
	gcc -c options.c -O2

reader.o: reader.c options.h ringbuffer.h stat.h thirdparty/kseq.h
	gcc -c reader.c -fopenmp -O3

stat.o: stat.c stat.h
	gcc -c stat.c -O2

writer.o: writer.c options.h ringbuffer.h stat.h
	gcc -c writer.c -fopenmp -O3

hasher.o: hasher.c options.h ringbuffer.h stat.h thirdparty/mt64.h thirdparty/sysinfo.h
	gcc -c hasher.c -fopenmp -O3

decider.o: decider.c options.h ringbuffer.h stat.h
	gcc -c decider.c  -fopenmp -O3

mt19937-64.o: thirdparty/mt19937-64.c thirdparty/mt64.h
	gcc -c thirdparty/mt19937-64.c -O3

sysinfo.o: thirdparty/sysinfo.c thirdparty/sysinfo.h thirdparty/procps.h thirdparty/alloc.h
	gcc -c thirdparty/sysinfo.c -O2
