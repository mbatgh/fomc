
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


int read_atp(const char* filename, t_atomtypes* atp) {

  FILE* fp;
  char buf[256];
  char tmpchar1[256];
  char tmpchar2[256];
  char tmpchar3;
  int linesread;
  int nitems;
  int i;
  int nread;
  double ptmplf1, ptmplf2;
  
  atp->n_types=0;
  linesread=0;
  
  if(  (fp=fopen(filename,"r")) == NULL) { error_msg("FOMC ... Cannot open atp file"); }

  while(fgets(buf, 256, (FILE*)fp)) {
    nitems=words(buf);
    if(nitems==9 || nitems==7) { linesread++;}
  }

  atp->atomtype   = c0mat(linesread,4);
  atp->sigma      = d0vec(linesread);
  atp->epsilon    = d0vec(linesread);
  atp->i_type_mod = i0vec(linesread);
  atp->i_type_fix = i0vec(linesread);
  atp->modflag    = i0vec(linesread);
  
  rewind(fp);
  
  while(fgets(buf, 256, (FILE*)fp)) {

    nitems=words(buf);
    
    if(nitems==9) {

      nread=sscanf(buf,"%s %s %lf %lf %s %lf %lf %c %d",
             atp->atomtype[atp->n_types],
             tmpchar1,
             &ptmplf1,
             &ptmplf2,
             tmpchar2,
             &(atp->sigma[atp->n_types]),
             &(atp->epsilon[atp->n_types]),
             &tmpchar3,
             &(atp->modflag[atp->n_types])
            );

      if(nread==9) atp->n_types++; else error_msg("FOMC ... problem reading atp file");
      
    } 
  }
  
  fclose(fp);

  atp->ntypes_mod=0;
  atp->ntypes_fix=0;
  
  for(i=0;i<atp->n_types;i++) {
    
    if(atp->modflag[i]==1) {
      atp->i_type_mod[atp->ntypes_mod]=i;
      atp->ntypes_mod++;    
    } else if(atp->modflag[i]==0) {
      atp->i_type_fix[atp->ntypes_fix]=i;
      atp->ntypes_fix++;
    } else {
      error_msg("FOMC ... format issue in atp file");
    }
  }
  
  /*##### report contents of top #####*/

  printf("FOMC ... Read initial atom type parameters from %s\n", filename);
  printf("FOMC ... number of types:           %8d\n", atp->n_types);

  fflush(NULL);

  return(linesread);
}

/*
    if((cpos = strchr(buf, ';'))!=NULL){
      ii = (int)(cpos - buf);
      buf[ii]='\0';
    }  
*/
