
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include "strucs.h"
#include "functions.h"



/* =========================================================================================== */
void error_msg(char* msg) {

  fprintf(stderr,"\n%s\n\n", msg);
  fprintf(stderr,"This is fomc (force field parameter optimization via MC/MD\n");
  fprintf(stderr,"Version 1.0\n\n");
  fprintf(stderr,"usage: fomc [command-line parameters]\n\n");
  fprintf(stderr,"  -h          show this msg\n");
  fprintf(stderr,"  -i  string  id, basename of input files: parameter inp, xtal struc pdb\n");
  fprintf(stderr,"              as in: id.inp, id-01.pdb to id-0x.pdb\n");
  fprintf(stderr,"  -s  int     seed for random generator\n");
  fprintf(stderr,"  -n  int     number of xtal polymorphs\n");
  fprintf(stderr,"  -r  string  run id, used in combination with random seed for unique output filenames\n");
  fprintf(stderr,"  -t  float   MC temperature\n");
  fprintf(stderr,"  -d  float   RMSD threshold [nm]\n");
  fprintf(stderr,"  -l  float   latticeparameter deviation threshold [%%]\n");
  fprintf(stderr,"  -z  int     timeout for a single MD run in seconds\n");
  fprintf(stderr,"  -y  int     total number of MC cycles\n");
  fprintf(stderr,"  -x  float   max dtot at which to stop (usually 2)\n");
  fprintf(stderr,"  -a  string  input atp filename\n");

  fprintf(stderr,"\n");

  exit(1);
}

/* =========================================================================================== */
int words(const char sentence[]) {

    int counted = 0;
    const char* it = sentence;
    int inword = 0;

    do switch(*it) {
        case '\0': 
        case ' ': case '\t': case '\n': case '\r':
            if (inword) { inword = 0; counted++; }
            break;
        default: inword = 1;
    } while(*it++);

    return counted;
}

/* =========================================================================================== */
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

/* =========================================================================================== */
void timeval_print(struct timeval *tv)
{
    char buffer[30];
    time_t curtime;

    printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
    curtime = tv->tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
    printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}
