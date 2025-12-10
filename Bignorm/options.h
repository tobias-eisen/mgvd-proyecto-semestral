#ifndef OPTIONS_H
#define OPTIONS_H

#define MAX_FILENAMELENGTH 150

int options(int argc, char **argv);

extern int verbosity;
extern int k,C,m,t,h,s,d,dumpfile, do_hist, qbase, opt_N, opt_n, opt_f, opt_b, opt_q, opt_A, opt_B;
extern float mem_frac;
extern int count_unpaired, count_paired;
extern char *dumpfilename;
extern int compressOut;
extern int outputdetail;
extern int writerej;

struct infiles
{
	int flag;
	char infile1[MAX_FILENAMELENGTH], infile2[MAX_FILENAMELENGTH];
	struct infiles *next;
};

extern struct infiles *input_start;
extern struct infiles *input_akt;	
#endif
