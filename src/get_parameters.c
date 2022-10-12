#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


int get_parameters(int ARGC,char** ARGV, t_parameters* pars) {

/***** variables set in input file or on the command line *****/

  int i;
  int irf,irc;

  FILE* fp;
  char buf[256];
  char tmpstring[256];
  char inputfilename[256];

  int IRUNID               = 0;
  int IATPFILE             = 0;
  int IINPID               = 0;

  pars->n_xtalpdbs         = -1;
  pars->n_mccycles         = -1;
  pars->mctemperature      = -1.0;
  pars->mdtimeout          = -1;
  pars->max_dtot           = -1.0;

  pars->rseed              = -1;

  pars->rmsthreshold       = -1.0;
  pars->abcthreshold       = -1.0;

  pars->atp_filename       = calloc(256, sizeof(char));
  pars->runid              = calloc(256, sizeof(char));
  pars->inputid            = calloc(256, sizeof(char));

/***** ... analyse the command line arguments ... *****/

  irc=0;

  for(i=1;i<ARGC;i++) {

    if(ARGV[i][0]!='-') error_msg("wrong format in command line");
      else switch(ARGV[i][1]) {

        case 'i': IINPID = ++i;
                  irc++;
                  break;
        case 's': pars->rseed = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'n': pars->n_xtalpdbs = atoi(ARGV[++i]);
                  irc++;
                  break;  
        case 'r': IRUNID = ++i;
                  irc++;
                  break;
        case 't': pars->mctemperature = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'z': pars->mdtimeout = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'y': pars->n_mccycles = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'x': pars->max_dtot = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'a': IATPFILE = ++i;
                  irc++;
                  break;
        case 'l': pars->abcthreshold = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'd': pars->rmsthreshold = atof(ARGV[++i]);
                  irc++;
                  break;                  
        case 'h': error_msg("");
                  break;
        default:  error_msg("unrecognized command line parameter");
    }
  }

  if(IATPFILE>0) strncpy(pars->atp_filename,ARGV[IATPFILE],strlen(ARGV[IATPFILE])+1);
  if(IRUNID>0)   strncpy(pars->runid,ARGV[IRUNID],strlen(ARGV[IRUNID])+1);

  if(IINPID==0)   {
    error_msg("# ERROR ... ID string missing");
  } else {

    strncpy(pars->inputid,ARGV[IINPID],strlen(ARGV[IINPID])+1);
    strncpy(inputfilename,pars->inputid,strlen(pars->inputid)+1);
    strncat(inputfilename,".inp",4);

    if(  (fp=fopen(inputfilename,"r"))==NULL  ) error_msg("# cannot open parameter input file");

    else {

      printf("FOMC  opened %s\n", inputfilename);

      irf=0;

/***** remove comments *****/

      while(fgets(buf, 256, (FILE*)fp)) {

        i=0;
        do {
          if(buf[i]=='#' || buf[i]==';') buf[i]='\0';
          i++;
        } while (buf[i]!='\n' && buf[i-1]!='\0');


/***** how many tokens are left? *****/

        if(words(buf)==2) {

          if(strstr(buf,"rseed")   !=NULL && pars->rseed==-1)
            if(sscanf(buf,"%s %d",  tmpstring, &(pars->rseed))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"n_xtalpdbs") !=NULL && pars->n_xtalpdbs==-1)
            if(sscanf(buf,"%s %d",  tmpstring, &(pars->n_xtalpdbs))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"atpfile") !=NULL && IATPFILE==0)
            if(sscanf(buf,"%s %s", tmpstring, pars->atp_filename)==2) {IATPFILE=1; irf++;}
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"runid")   !=NULL && IRUNID==0)
            if(sscanf(buf,"%s %s", tmpstring, pars->runid)==2)             {IRUNID=1; irf++;}
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"MCtemp")  !=NULL && pars->mctemperature==-1.0)
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->mctemperature))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"MDtimeout") !=NULL && pars->mdtimeout==-1)
            if(sscanf(buf,"%s %d", tmpstring, &(pars->mdtimeout))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"MCcycles") !=NULL && pars->n_mccycles==-1)
            if(sscanf(buf,"%s %d", tmpstring, &(pars->n_mccycles))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"thrRMSD") !=NULL && pars->rmsthreshold==-1.0)
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->rmsthreshold))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"thrDLAT") !=NULL && pars->abcthreshold==-1.0)
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->abcthreshold))==2) irf++;
            else error_msg("Problem reading parameter file");
            
          if(strstr(buf,"maxdtot") !=NULL && pars->max_dtot==-1.0)
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->max_dtot))==2) irf++;
            else error_msg("Problem reading parameter file");
	}
      }
    }
  }


/***** write info *****/

  if(IRUNID  ==0) error_msg("No run-id given");
  
  if(IATPFILE==0) {
    printf("FOMC  WARNING ... using default of atomtypes.atp for atp filename\n");
    strncpy(pars->atp_filename,"atomtypes.atp",13);
    pars->atp_filename[13]='\0';
  }
  if(pars->rseed              == -1) {
    printf("FOMC  WARNING ... using default of 1 for rseed\n");
    pars->rseed=1;
  }
  if(pars->mctemperature      == -1.0) {
    printf("FOMC  WARNING ... using default of 1.0 for initial mctemperature\n");
    pars->mctemperature=1.0;
  }
  if(pars->mdtimeout          == -1) {
    printf("FOMC  WARNING ... using default of 300 secs for mdtimeout\n");
    pars->mdtimeout=300.0;
  }
  if(pars->max_dtot           == -1) {
    printf("FOMC  WARNING ... using default of 2 for max_dtot\n");
    pars->max_dtot=2.0;
  }

  printf("#\nFOMC  This is fomc version 1\n");
  printf("FOMC \n");
  printf("FOMC  read %d parameters from command line and %d from %s\n", irc, irf, inputfilename);
  printf("FOMC \n");
  printf("FOMC  run id:                        %s\n",      pars->runid);
  printf("FOMC  input id:                      %s\n",      pars->inputid);
  printf("FOMC  number of polymorphs:          %d\n",      pars->n_xtalpdbs);

  printf("FOMC  MC temperature:                %12.3le\n", pars->mctemperature);
  printf("FOMC  MD timeout:                    %d\n",      pars->mdtimeout);
  printf("FOMC  MC number of cycles:           %d\n",      pars->n_mccycles);
  printf("FOMC  MC rmsd threshold:             %12.3le\n", pars->rmsthreshold);
  printf("FOMC  MC dlat threshold:             %12.3le\n", pars->abcthreshold);
  printf("FOMC  MC max dtot:                   %12.3le\n", pars->max_dtot);

  printf("FOMC  Input atp filename:            %s\n",      pars->atp_filename);
  printf("FOMC  Random seed:                   %12d\n",    pars->rseed);
  printf("FOMC \n");

  fflush(stdout);

  return(irc+irf);
}
