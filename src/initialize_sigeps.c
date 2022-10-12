
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


int initialize_sigeps(const t_atomtypes* input_atp, const t_atomtypes* gaff2_atp, t_ljparameters* ljpars) {

  int npars;
  int i, j;
  int foundone;
  int nfound;
  
  nfound=0;
  
  npars = input_atp->n_types;

  ljpars->sig0   = calloc(npars, sizeof(double));
  ljpars->sigmin = calloc(npars, sizeof(double));
  ljpars->sigmax = calloc(npars, sizeof(double));
  ljpars->eps0   = calloc(npars, sizeof(double));
  ljpars->epsmin = calloc(npars, sizeof(double));
  ljpars->epsmax = calloc(npars, sizeof(double));

  for(i=0;i<npars;i++) {
    
    ljpars->sig0[i] = input_atp->sigma[i];
    ljpars->eps0[i] = input_atp->epsilon[i];
    
    foundone=0;
    
    for(j=0;j<gaff2_atp->n_types;j++) {
      if(strcmp(input_atp->atomtype[i], gaff2_atp->atomtype[j])==0) {foundone=1; break;}
    }
    
    if(foundone==0) {
      printf("looking for parameter nr %d: %s\n", i, input_atp->atomtype[i]);
      error_msg("did not find appropriate atom GAFF2 atom type");
    }

    else {
      printf("FOMC  found GAFF2 parameters for atom type %s\n", input_atp->atomtype[i]);
      nfound++;
    }

    ljpars->sigmax[i] = gaff2_atp->sigma[j]*FRACDSIGMAX;
    ljpars->sigmin[i] = gaff2_atp->sigma[j]*FRACDSIGMIN;
    ljpars->epsmax[i] = gaff2_atp->epsilon[j]*FRACDEPSMAX;
    ljpars->epsmin[i] = gaff2_atp->epsilon[j]*FRACDEPSMIN;
  }
  
  return(nfound);
}
