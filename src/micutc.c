#ifndef STDIOH
#include <stdio.h>
#define STDIOH
#endif
#ifndef TERMIOSH
#include <termios.h>
#define TERMIOSH
#endif
#ifndef STDLIBH
#include <stdlib.h>
#define STDLIBH
#endif
#ifndef UNISTDH
#include <unistd.h>
#define UNISTDH
#endif
#ifndef TIMEH
#include <time.h>
#define TIMEH
#endif

#define     NLo         6.022045e23
#define     kB          1.380662e-23
#define     Rgas        8.3143
#define     amu         1.6605655e-27
#define     epsilon_0   8.8542e-12
#define     el_Charge   1.6021e-19
#define     PI          3.14159265358979323846



struct termios saved_attributes;

void
reset_input_mode(void) {
    tcsetattr (STDIN_FILENO,TCSANOW,&saved_attributes);
    }

void 
set_input_mode(void) {
    struct termios tattr;

    if(!isatty(STDIN_FILENO)) {
        fprintf(stderr,"\nNot a terminal ...");
        exit(EXIT_FAILURE);
    }

    tcgetattr (STDIN_FILENO, &saved_attributes);
    atexit(reset_input_mode);

    tcgetattr(STDIN_FILENO, &tattr);
    tattr.c_lflag &= ~(ICANON|ECHO);
    tattr.c_cc[VMIN] = 1;
    tattr.c_cc[VTIME] = 0;
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &tattr);
}

int ask_stop()
    {
    char c;
    ssize_t s;
    set_input_mode();

    while(1)
        {
        printf("proceed ? (y/n)\n");
        s=read(STDIN_FILENO,&c,1);
        switch(c)
            {
            case 'y': reset_input_mode(); return 0;
            case 'n': reset_input_mode(); return 1;
            default:  printf("y or n !\n"); break;
            }
        }

    return EXIT_FAILURE;
    }
    
void presskey(void)
    {
    ssize_t s;
    char c;
    set_input_mode();
    printf("press any key to proceed ...\n");
    fflush(NULL);
    s=read(STDIN_FILENO,&c,1);
    reset_input_mode();
    }    

void mic_msg(char* msg)
    {
    char c;
    ssize_t s;
    set_input_mode();
    printf("%s", msg);
    printf(" ... press any key to proceed\n");
    s=read(STDIN_FILENO,&c,1);
    reset_input_mode();
    }    

   
char mic_getchar(char* text) 
    {
    ssize_t s;
    char c;
    set_input_mode();
    printf("%s\n",text);
    s=read(STDIN_FILENO,&c,1);
    reset_input_mode();
    return(c);
    }    

void mabort(const char* msg) {

    printf("runtime-error ocurred!\n");
    printf("%s\n", msg);
    exit(1);
}

double **dmat(int nrh, int nch) {
        int i;
        double **m;

        m=(double **) malloc((size_t)((nrh+1)*sizeof(double*)));
        if (!m) mabort("allocation failure 1 in dmat()");

        m[1]=(double *) malloc((size_t)((nrh*nch+1)*sizeof(double)));
        if (!m[1]) mabort("allocation failure 2 in dmat()");

        for(i=1+1;i<=nrh;i++) m[i]=m[i-1]+nch;

        return m;
}

double *dvec(int nh) {

        double *v;
        v=(double *)malloc((size_t) ((nh+1)*sizeof(double)));
        if (!v) mabort("allocation failure in dvec()");
        return v;
}


float **fmat(int nrh, int nch) {
        int i;
        float **m;

        m=(float **) malloc((size_t)((nrh+1)*sizeof(float*)));
        if (!m) mabort("allocation failure 1 in fmat()");

        m[1]=(float *) malloc((size_t)((nrh*nch+1)*sizeof(float)));
        if (!m[1]) mabort("allocation failure 2 in fmat()");

        for(i=1+1;i<=nrh;i++) m[i]=m[i-1]+nch;

        return m;
}

float *fvec(int nh) {

        float *v;
        v=(float *)malloc((size_t) ((nh+1)*sizeof(float)));
        if (!v) mabort("allocation failure in fvec()");
        return v;
}

int **imat(int nrh, int nch) {
        int i;
        int **m;

        m=(int **) malloc((size_t)((nrh+1)*sizeof(int*)));
        if (!m) mabort("allocation failure 1 in imat()");

        m[1]=(int *) malloc((size_t)((nrh*nch+1)*sizeof(int)));
        if (!m[1]) mabort("allocation failure 2 in imat()");

        for(i=1+1;i<=nrh;i++) m[i]=m[i-1]+nch;

        return m;
}

char **cmat(int nrh, int nch) {
        int i;
        char **m;

        m=(char **) malloc((size_t)((nrh+1)*sizeof(char*)));
        if (!m) mabort("allocation failure 1 in imat()");

        m[1]=(char *) malloc((size_t)((nrh*nch+1)*sizeof(char)));
        if (!m[1]) mabort("allocation failure 2 in imat()");

        for(i=1+1;i<=nrh;i++) m[i]=m[i-1]+nch;

        return m;
}

int *ivec(int nh) {

        int *v;
        v = (int *)calloc(nh+1,sizeof(int));
        if (!v) mabort("allocation failure in ivec()");
        return v;
}

char *cvec(int nh) {

        char *v;
        v=(char *)malloc((size_t) ((nh+1)*sizeof(char)));
        if (!v) mabort("allocation failure in cvec()");
        return v;
}

double **d0mat(int nrh, int nch) {
        int i;
        double **m;

        m=(double **) malloc((size_t)((nrh)*sizeof(double*)));
        if (!m) mabort("allocation failure 1 in dmat()");

        for(i=0;i<nrh;i++) {
          m[i]=(double *) malloc((size_t)(nch*sizeof(double)));
          if (!m[i]) mabort("allocation failure 2 in dmat()");
	}

        return m;
}

double *d0vec(int nh) {

        double *v;
        v=(double *)malloc((size_t) (nh*sizeof(double)));
        if (!v) mabort("allocation failure in dvec()");
        return v;
}

float **f0mat(int nrh, int nch) {
        int i;
        float **m;

        m=(float **) malloc((size_t)((nrh)*sizeof(float*)));
        if (!m) mabort("allocation failure 1 in dmat()");

        for(i=0;i<nrh;i++) {
          m[i]=(float *) malloc((size_t)(nch*sizeof(float)));
          if (!m[i]) mabort("allocation failure 2 in dmat()");
	}

        return m;
}

float *f0vec(int nh) {

        float *v;
        v=(float *)malloc((size_t) (nh*sizeof(float)));
        if (!v) mabort("allocation failure in dvec()");
        return v;
}

int **i0mat(int nrh, int nch) {
        int i;
        int **m;

        m=(int **) malloc((size_t)((nrh)*sizeof(int*)));
        if (!m) mabort("allocation failure 1 in dmat()");

        for(i=0;i<nrh;i++) {
          m[i]=(int *) malloc((size_t)(nch*sizeof(int)));
          if (!m[i]) mabort("allocation failure 2 in dmat()");
	}

        return m;
}

int *i0vec(int nh) {

        int *v;
        v=(int *)malloc((size_t) (nh*sizeof(int)));
        if (!v) mabort("allocation failure in dvec()");
        return v;
}

char **c0mat(int nrh, int nch) {
        int i;
        char **m;

        m=(char **) malloc((size_t)((nrh)*sizeof(char*)));
        if (!m) mabort("allocation failure 1 in dmat()");

        for(i=0;i<nrh;i++) {
          m[i]=(char *) malloc((size_t)(nch*sizeof(char)));
          if (!m[i]) mabort("allocation failure 2 in dmat()");
	}

        return m;
}

char *c0vec(int nh) {

        char *v;
        v=(char *)malloc((size_t) (nh*sizeof(char)));
        if (!v) mabort("allocation failure in dvec()");
        return v;
}

void waitsec(int secs) {

    time_t t1;
    time_t t2;

    time(&t1);

    do {
        time(&t2);
    } while(difftime(t2,t1)<(double)secs);
}
