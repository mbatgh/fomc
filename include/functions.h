

extern int get_parameters(int,char*[],t_parameters*);
extern int read_atp(const char*,t_atomtypes*);
extern int initialize_sigeps(const t_atomtypes*,const t_atomtypes*,t_ljparameters*);
extern int initialize_morphs(const int, const t_parameters*, t_morphdata**);

extern int words(const char*);
extern void error_msg(char*);

extern void reset_input_mode(void);
extern void set_input_mode(void);
extern int ask_stop(void);
extern void presskey(void);
extern void mic_msg(char*);
extern char mic_getchar(char*);
extern void mabort(const char*);

extern double **dmat(int,int);
extern double *dvec(int);
extern float **fmat(int,int);
extern float *fvec(int);
extern int **imat(int,int);
extern char **cmat(int,int);
extern int *ivec(int);
extern char *cvec(int);
extern double **d0mat(int,int);
extern double *d0vec(int);
extern float **f0mat(int,int);
extern float *f0vec(int);
extern int **i0mat(int,int);
extern int *i0vec(int);
extern char **c0mat(int,int);
extern char *c0vec(int);
extern void waitsec(int);


extern const double FRACDEPSMAX;
extern const double FRACDEPSMIN;
extern const double FRACDSIGMAX;
extern const double FRACDSIGMIN;
extern const int TEPS;
extern const int TSIG;
