
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


int initialize_morphs(const int nmorphs, const t_parameters* pars, t_morphdata** morph) {

  int nm;
  int nitemsread;
  int lastresnumber;
  int tmpnres;
  int resnumber;
  char s_resnumber[5];
  size_t len0;
  FILE* fpdb;
  char* oneline;
  char sdummy[256];
  char mycommand[5000];
  double dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6;
  
  oneline=NULL;
  len0=(size_t)0;

  *morph=calloc(nmorphs, sizeof(t_morphdata));

  for(nm=0;nm<nmorphs;nm++) {
    
    (*morph)[nm].pdbfn=(char*)calloc(256,sizeof(char));
    (*morph)[nm].ndxfn=(char*)calloc(256,sizeof(char));

    sprintf((*morph)[nm].pdbfn,"%s-%02d.pdb",pars->inputid, nm+1);
    sprintf((*morph)[nm].ndxfn,"%s-%02d.ndx",pars->inputid, nm+1);
  }

  for(nm=0;nm<nmorphs;nm++) {
    sprintf(mycommand,"cat 4ndx | gmx make_ndx -f %s -o %s > /dev/null 2>&1\n", (*morph)[nm].pdbfn, (*morph)[nm].ndxfn);
    if(system(mycommand)!=0) error_msg("gmx make_ndx did not work");
  }


  for(nm=0;nm<nmorphs;nm++) {
    
    (*morph)[nm].nrestot=0;
    lastresnumber=0;
    

    if((fpdb=fopen((*morph)[nm].pdbfn,"r"))==NULL) error_msg("cannot open pdb file");
    else {

      
      while((nitemsread=getline(&oneline,&len0,fpdb))>0) {
        
        if(strstr(oneline,"CRYST1")!=NULL) {
          
          if(sscanf(oneline,"%s %lf %lf %lf %lf %lf %lf",sdummy,&dtmp1,&dtmp2,&dtmp3,&dtmp4,&dtmp5,&dtmp6)!=7) {
            error_msg("cannot read CRYST1 line in pdb file");    
          } else {
            printf("FOMC ... read CRYST1 line in %s\n", (*morph)[nm].pdbfn);
            printf("FOMC ... %9.3lf %9.3lf %9.3lf %9.3lf %9.3lf %9.3lf\n",dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6);
            (*morph)[nm].lata0=dtmp1;
            (*morph)[nm].latb0=dtmp2;
            (*morph)[nm].latc0=dtmp3;
          }
        }


        if(strstr(oneline,"ATOM")!=NULL) {
          
          strncpy(s_resnumber, oneline+22, 4);
          resnumber=atoi(s_resnumber);
          if(lastresnumber!=resnumber) (*morph)[nm].nrestot++;
          lastresnumber=resnumber;
        }         
      }
  
      rewind(fpdb);
      
/*  morph[nm].resname=(char**)calloc(morph[nm].nrestot,sizeof(char));
 for(i=0;i<morph[nm].nrestot;i++) morph[nm].resname[i]=(char*)calloc(4,sizeof(char));*/

      (*morph)[nm].resname = c0mat((*morph)[nm].nrestot,4); 
  
      lastresnumber=0;
      tmpnres=0;
      
      while((nitemsread=getline(&oneline,&len0,fpdb))>0) {
        
        if(strstr(oneline,"ATOM")!=NULL) {
          
          strncpy(s_resnumber, oneline+22, 4);
          resnumber=atoi(s_resnumber);
        
          if(lastresnumber!=resnumber) {
            
            strncpy((*morph)[nm].resname[tmpnres], oneline+17, 3);
            (*morph)[nm].resname[tmpnres][3]='\0';
  
            tmpnres++;

          }
          lastresnumber=resnumber;
        }
      }
      
      fclose(fpdb);
    }
  }
  
  fflush(NULL);
 
  return((*morph)[0].nrestot);
}
