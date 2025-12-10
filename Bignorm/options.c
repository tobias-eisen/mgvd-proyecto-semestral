/* options.c: parse the command line and set options

 Thanks to 'Step-by-Step into Argp' by  Ben Asselstine

*/
#include <stdio.h>
#include <argp.h>
// #include <argz.h>  // Not available on macOS and not used in this file
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "options.h"

// default values
int verbosity = 3;
int k = 32; // k-mer size
int C = 20; // parameter 1 for decision (strategy 1: cut-off constant for median
int s = 1;  // which memory strategy to use
float mem_frac = 0.9; // use max. 90% of RAM / available memory
int d = 6;  // which decision strategy to use
int dumpfile = 0; // dump CM when done?
char *dumpfilename;
int compressOut = 0;
int outputdetail = 0;
int writerej = 0;
int m = 0;
int t = 0;
int h = 1;
int do_hist = 0;
int qbase = 33; // default for phred score base
int opt_q = 20; // defaul minimum phred score above base
int opt_N = -1;
int opt_n = 1;
int opt_f = 0;
int opt_b = 0;
int opt_A = 3;
int opt_B = 3;

// input file handling
int count_unpaired = 0;
int count_paired = 0;
int pairflag = 0;
struct infiles *input_start = NULL;
struct infiles *input_akt = NULL;

// Argp stuff

// const char *argp_program_bug_address = "someone@example.com";
const char *argp_program_version = "version 0.01";

static int
parse_opt (int key, char *arg,
		struct argp_state *state)

{
	switch (key)
	{
		case 'v': verbosity++; break;
		case 'q': verbosity--; break;
		case 'z': compressOut = 1; break;
		case 'r': writerej = 1; break;
		case 'H': do_hist = 1; break;
		case 'f': opt_f = 1; break;
		case 'b': opt_b = 1; break;
		case 'n': opt_n = 0; break;
		case 'Q': opt_q=atoi(arg); break;
		case '0': qbase =atoi(arg);
			if((qbase != 33) & (qbase != 59) & (qbase != 64)){
				argp_failure(state, 1, 0, "error: got unknown qbase");
			}
			break;
		case 'k': k=atoi(arg);
			  if((k<3) | (k>32)){
				  argp_failure(state, 1, 0, "error: given k out of range\n\tk has to be between 3 and 32");
				  // return -1;
			  }
			  if(k < 15) printf("warning: the value of k you have choosen (%i) is very low, usable results should not be expected\n\n",k);
			  break;
		case 'A': opt_A = atoi(arg);
			break;
		case 'B': opt_B = atoi(arg);
			break;
		case 'm':
			m=atoi(arg);
			break;
		case 'N':
			opt_N=atoi(arg);
			break;
		case 't':
			t=atoi(arg);
			if(t < 0){
				printf("t has to be positive! setting t=0\n");
				t = 0;
			}
			if(t > 25){
				printf("warning: the value of t you have choosen (%i) is to big, reverting it to 25\n",t);
				t = 25;
			}
			break;
		case 'h':
			h=atoi(arg);
			if((h < 0) | (h > 3)){
				argp_failure(state, 1, 0, "error: given h out of range\n\ts has to be between 0  and 3");
			}
			break; 
		case 's': s=atoi(arg);
			if((s<0) | (s>3)){
				argp_failure(state, 1, 0, "error: given s out of range\n\ts has to be between 0 and 3");
			}
			break;
		case 'M':
			mem_frac=atof(arg);
			if(mem_frac <= 0.0) argp_failure(state, 1, 0, "error: mem_frac must be > 0");
			printf("mem_frac set to %f\n", mem_frac);
			break;
		case 'C':
			C=atoi(arg);
			break;
		case 'o':
			outputdetail=atoi(arg);
			break;
		case 'd':
			d=atoi(arg);
			//if((d<1) | (d>3)){
			//	argp_failure(state, 1, 0, "error: given d out of range\n\td has to be between 1 and 3");
			//}
			break; 
		case 'D':
			dumpfile = 1;
			dumpfilename=malloc(30);
			if (arg == NULL)
				dumpfilename = "CM.dump";
			else {
				strncpy(dumpfilename, arg, 29);
				dumpfilename[29] = '\0';
				}
			break;
		case 'u':
			if(pairflag != 0){
				printf("error: paired input files have to be given together\n");
				if(input_akt != NULL){
					printf("\t got filename %s, expecting paired file (-2)\n", input_akt->infile1);
				}
				exit(-1);
			}
			if(input_start == NULL){ // first input file
				input_start = malloc(sizeof(struct infiles));
				if(input_start == NULL){
					printf("error allocating memory for input file %s\n", arg);
					perror("Message:");
					exit(-1);
				}
				input_akt = input_start;
			} else {
				input_akt->next = malloc(sizeof(struct infiles));
				if(input_akt->next == NULL){
					printf("error allocating memory for input file %s\n", arg);
					perror("Message:");
					exit(-1);
				}
				input_akt = input_akt->next;
			}
			strncpy(input_akt->infile1, arg, MAX_FILENAMELENGTH);
			if(input_akt->infile1[MAX_FILENAMELENGTH - 1] != '\0'){
				printf("filename %s too long, pleae rename file\n", arg);
				exit(-1);
			}
			// is file accesible?
			if(access(input_akt->infile1, R_OK) != 0){
				printf("file %s is not accesible\n", input_akt->infile1);
				exit(-1);
			}
			input_akt->flag = 1;
			input_akt->next = NULL;
			count_unpaired++;
			break;
		case '1':
			if(pairflag != 0){
				printf("error: paired input files have to be given together\n");
				if(input_akt != NULL){
					printf("\t got filename %s, expecting paired file (-2)\n", input_akt->infile1);
				}
				exit(-1);
			}
			pairflag = 1;
			if(input_start == NULL){ // first input file
				input_start = malloc(sizeof(struct infiles));
				if(input_start == NULL){
					printf("error allocating memory for input file %s\n", arg);
					perror("Message:");
					exit(-1);
				}
				input_akt = input_start;
			} else {
				input_akt->next = malloc(sizeof(struct infiles));
				if(input_akt->next == NULL){
					printf("error allocating memory for input file %s\n", arg);
					perror("Message:");
					exit(-1);
				}
				input_akt = input_akt->next;
			}
			strncpy(input_akt->infile1, arg, MAX_FILENAMELENGTH);
			if(input_akt->infile1[MAX_FILENAMELENGTH - 1] != '\0'){
				printf("filename %s too long, pleae rename file\n", arg);
				exit(-1);
			}
			// is file accesible?
			if(access(input_akt->infile1, R_OK) != 0){
				printf("file %s is not accesible\n", input_akt->infile1);
				exit(-1);
			}
			input_akt->flag = 2;
			input_akt->next = NULL;
			break;
		case '2':
			if(pairflag != 1){
				printf("error: got second file of a pair (%s) before first file\n", arg);
				exit(-1);
			}
			pairflag = 0;
			strncpy(input_akt->infile2, arg, MAX_FILENAMELENGTH);
			if(input_akt->infile2[MAX_FILENAMELENGTH - 1] != '\0'){
				printf("filename %s too long, pleae rename file\n", arg);
				exit(-1);
			}
			if(access(input_akt->infile2, R_OK) != 0){
				printf("file %s is not accesible\n", input_akt->infile2);
				exit(-1);
			}
			input_akt->flag = 3;
			count_paired++;
			break;
		case ARGP_KEY_ARG: argp_failure(state, 1, 0,"got unexpected arg %s\n",arg); break;
	}
	return 0;
}

int options(int argc, char **argv){
	// Flags Feld 4: OPTION_ARG_OPTIONAL, 	
	struct argp_option options[] =
		{
			// { "dot", 'd', "NUM", OPTION_ARG_OPTIONAL, "Show a dot on the screen"},
			{0, 0, 0, 0, "Level of verbosity:",7},
			{ "verbose", 'v', 0, 0, "increase verbosity (use multiple times for much more verbosity)"},
			{ "quit", 'q', 0, 0, "decrease verbosity (use multiple times for much less verbosity)"},
			{0, 0, 0, 0, "options for normalization:",2},
			{ "table_width", 'm', "NUM", 0, "width of count-min sketch table given as power of two (size=2^m MB), default 0 (choose automatically)"},
			{ "num_rows", 't', "NUM", 0, "number of of count-min sketch rows, default 0 (choose automatically)"},
			{ "hash", 'h', "NUM", 0, "hash set to use (default: set 1, max 3, use 0 for random set)"},
			{ "kmer_size", 'k', "NUM", 0, "k-mer Size (default and maximum: 32)"},
			{ "qbase", '0', "NUM", 0, "base for phred score, default: 33, alternative: 64 or 59"},
			{ "qmin", 'Q', "NUM", 0, "minimum phred score above base to be accepted for counting (default: 20)"},
			{ "mem_strat", 's', "NUM", 0, "strategy for memory allocation to use"},
			{ "mem_frac", 'M', "NUM", 0, "parameter for memory allocation"},
			{ "decision", 'd', "NUM", 0, "strategy for deciding which read to keep"},
			{ "decPar1", 'C', "NUM", 0, "parameter 1 for deciding function"},
			{ "decPar2", 'A', "NUM", 0, "parameter 2 for deciding function"},
			{ "decPar3", 'B', "NUM", 0, "parameter 3 for deciding function"},
			{ "both", 'b', 0, 0, "for paired reads and some decision functions only: both reads have to be accepted (default: just one)"},
			{ "dumpCM", 'D', "FILE", OPTION_ARG_OPTIONAL, "dump the data of the count min sketch to file when program ends"},
			{0, 0, 0, 0, "N-handling:",3},
			{ "maxN", 'N', "NUM", 0, "maximal count of N (wildcard for unknown nucletide) per read"},
			{ "countN", 'n', 0, 0, "don't ignore k-mers including N when counting"},
			{0, 0, 0, 0, "input / output:",4},
			{ "unpaired", 'u', "FILE", OPTION_NO_USAGE, "unpaired input file, use multiple times for multiple files (fasta / fastq, may be gzip-ed)"},
			{ "paired_one", '1', "FILE", OPTION_NO_USAGE, "first file of a pair of input files"},
			{ "paired_two", '2', "FILE", OPTION_NO_USAGE, "second file of a pair of input files"}, 
			{ "output_detail", 'o', "NUM", 0, "level of detail in output files"},
			{ "force_fasta", 'f', 0, 0, "always output in FASTA format (even if input is FASTQ)"},
			{ "print_histo", 'H', 0, 0, "print a histogram of CM sketch"},
			{ "write_rejected", 'r', 0, 0, "write rejected reads too *_rej files"},
			{ "compress_output", 'z', 0, 0, "compress output files (using zlib)"},
			// { "long-only", 300, 0, 0, "test"},
			{0, 0, 0, 0, "Informational options:", -1},
			{0}
	};
	struct argp argp = { options, parse_opt, "-1 reads_1.fa -2 reads_2.fa [...]\n-u reads.fa [...]", "\nBignorm\n\t a faster diginorm implementation\vSee manual for further information" };
	return argp_parse (&argp, argc, argv, 0, 0, 0); // Flags: ARGP_IN_ORDER,
}



