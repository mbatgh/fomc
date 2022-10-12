

  typedef struct {


    int n_xtalpdbs;
    int n_mccycles;
    int mdtimeout;

    double mctemperature;
    double rmsthreshold;
    double abcthreshold;
    double max_dtot;

    int rseed;

    char* atp_filename;
    char* inputid;
    char* runid;

  } t_parameters;


  typedef struct {

    char**  atomtype;
    double* sigma;
    double* epsilon;
    
    int* modflag;
    int* i_type_mod;
    int* i_type_fix;
    
    int  ntypes_mod;
    int  ntypes_fix;
    
    int  n_types;

  } t_atomtypes;
  
  typedef struct {

    double* sig0;
    double* sigmin;
    double* sigmax;
    double* eps0;
    double* epsmin;
    double* epsmax;
    
  } t_ljparameters;
  
  typedef struct {
  
    char* pdbfn;
    char* ndxfn;

    double lata0;
    double latb0;
    double latc0;

    char** resname;
    int    nrestot;

    double mcrmsd;
    double mcdlat;

  } t_morphdata;
  
