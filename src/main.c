/* ================================================================= */
/* === force field optimization with genetic algoriem, author MB === */
/* ================================================================= */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include "strucs.h"
#include "functions.h"

  const int TEPS=1;
  const int TSIG=2;
  const int TEST_INTERVAL=20;

  const double FRACDEPSMAX=2.0;
  const double FRACDEPSMIN=0.5;
  const double FRACDSIGMAX=1.4;
  const double FRACDSIGMIN=0.6;

/* ===== function declarations ============================================================= */

int get_parameters(int,char*[],t_parameters*);
int read_atp(const char*,t_atomtypes*);
int initialize_sigeps(const t_atomtypes*,const t_atomtypes*,t_ljparameters*);
int initialize_morphs(const int, const t_parameters*, t_morphdata**);

void error_msg(char*);



/* ===== main ============================================================================== */

int main(int ARGC, char* ARGV[]) {

/***** local variables *****/

  int j, k, l, nm, tmpmodindex;
  int nread;
  int iret;
  int GMXSTATUS;
  int nreject;
  int modindex;
  int changetype;
  int cycle_of_best_dtot;
  int n_cycle;
  int n_accepted;
  int n_ttest;

  double pvar;
  double dtmp;
  double dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6;
  double frand1,frand2,frand3;
  
  double mcrmsdcur;
  double mcdlatcur;

  double mcdtot_new;
  double mcdtot_old;
  double mcdtot_t0;

  double dlata;
  double dlatb;
  double dlatc;

  double origeps;
  double origsig;
  double bf;
  double best_dtot;

  FILE* ftop;
  FILE* fxvg;
  FILE* fone;

  char mycommand[5000];
  char sdummy[256];

  char mdbasename[256];
  char mdpfn[256];
  char topfn[256];
  char xvgfn[256];
  char onefn[256];
  
  char besttopfn[256];

  char* oneline;

  ssize_t nitemsread;
  size_t len0;

  t_parameters   pars;
  t_atomtypes    input_atp;
  t_atomtypes    gaff2_atp;
  t_ljparameters ljpars;
  t_morphdata*   morph;

  struct timeval tvBegin, tvEnd, tvDiff;

  len0=(size_t)0;

/***** record start time *****/
  gettimeofday(&tvBegin, NULL);

/***** get parameters from command line or input file *****/
  if(get_parameters(ARGC, ARGV, &pars)<1) error_msg("Something went wrong reading input parameters");

/***** read parameters from atp file *****/
  if(read_atp(pars.atp_filename,&input_atp)<1) error_msg("Something went wrong reading atp file");

/***** read gaff2 LJ parameters from atp file *****/
  if(read_atp("sigma-epsilon-gaff2.atp",&gaff2_atp)<1) error_msg("Something went wrong reading atp0 file");

/***** establish boundaries for LJ parameters *****/
  if(initialize_sigeps(&input_atp,&gaff2_atp,&ljpars)<1) error_msg("Something went wrong initializing LJ parameters");

/***** initiate data and allocate variables for all mrophic forms *****/
  if(initialize_morphs(pars.n_xtalpdbs, &pars, &morph)<1) error_msg("Something went wrong initializing morphs");

  fflush(NULL);

  sprintf(besttopfn,"best-%s-%d.top", pars.runid, pars.rseed);

/*###################################################*/
/*##### determine performance of input topology #####*/
/*###################################################*/

  GMXSTATUS=0;

/* ##### loop over all morphs ##### */

  for(nm=0;nm<pars.n_xtalpdbs;nm++) {

    sprintf(mdbasename,"%s-%d-%06d-%02d",pars.runid, pars.rseed, 0, nm+1);
    sprintf(topfn,"%s-%d-%06d-%02d.top",pars.runid, pars.rseed, 0, nm+1);
    sprintf(mdpfn,"md%02d.mdp",nm+1);

    if((ftop=fopen(topfn,"w"))==NULL) error_msg("cannot write topology file");
    else {
      fprintf(ftop, "[ defaults ]\n");
      fprintf(ftop, "1   2   yes   0.5  0.8333\n\n");
      fprintf(ftop, "[ atomtypes ]\n");

      for(k=0;k<input_atp.n_types;k++) {
        fprintf(ftop," %-6s  %-6s   0.00000  0.00000   A     %14.5e  %14.5e\n",
          input_atp.atomtype[k],input_atp.atomtype[k],input_atp.sigma[k],input_atp.epsilon[k]);
      }

      fprintf(ftop, "\n");
      fprintf(ftop, "#include \"moleculetypes.itp\"\n");
      fprintf(ftop, "\n");
      fprintf(ftop, "[ system ]\n");
      fprintf(ftop, "system\n\n");
      fprintf(ftop, "[ molecules ]\n");

      for(j=0;j<morph[nm].nrestot;j++) {
        fprintf(ftop, "%s 1\n", morph[nm].resname[j]);
      }
      fclose(ftop);
    }

/* ##### do MD ##### */

    sprintf(mycommand,"gmx grompp -f %s -c %s -p %s -o %s -po %s.mdp", mdpfn, morph[nm].pdbfn, topfn, mdbasename, mdbasename);
    if(system(mycommand)!=0) error_msg("MC: grompp command failed");
    sprintf(mycommand,"timeout %d gmx mdrun -deffnm %s", pars.mdtimeout, mdbasename);
    iret=system(mycommand);

    if(iret==0) {

/* ##### calc rmsd ##### */

      sprintf(mycommand,"cat 0 | gmx trjconv -f %s.trr -s %s.tpr -pbc nojump -o nj-%s.trr", mdbasename,mdbasename,mdbasename);
      iret=system(mycommand);
      sprintf(mycommand,"gmx rms -b 100 -f nj-%s.trr -s %s.tpr -o rmsd-%s.xvg -xvg none -n %s", mdbasename,mdbasename,mdbasename,morph[nm].ndxfn);
      iret=system(mycommand);

      sprintf(xvgfn,"rmsd-%s.xvg",mdbasename);
      if((fxvg=fopen(xvgfn,"r"))==NULL) error_msg("cannot open xvg file");
      else {
        nread=0;
        morph[nm].mcrmsd=0.0;
        while(fscanf(fxvg,"%lf %lf",&dtmp1,&dtmp2)==2) {
          morph[nm].mcrmsd+=dtmp2;
          nread++;
        }
        if(nread>0) morph[nm].mcrmsd /= (double)nread;
      }
      fclose(fxvg);

/* ##### calc deviations of lattice parameters a b c ##### */

      dlata=0.0;
      dlatb=0.0;
      dlatc=0.0;

      sprintf(onefn,"p-%s.pdb",mdbasename);

      sprintf(mycommand,"gmx trjconv -b 100 -f %s.trr -s %s.tpr -o %s -xvg none -n 1.ndx",mdbasename,mdbasename,onefn);
      iret=system(mycommand);
      if((fone=fopen(onefn,"r"))==NULL) error_msg("cannot open 1.pdb file");

      nread=0;

      while((nitemsread=getline(&oneline,&len0,fone))>0) {
        if(strstr(oneline,"CRYST1")!=NULL) {
          if(sscanf(oneline,"%s %lf %lf %lf %lf %lf %lf",sdummy,&dtmp1,&dtmp2,&dtmp3,&dtmp4,&dtmp5,&dtmp6)!=7) {
            error_msg("cannot read CRYST1 line in pdb file");
          } else {
            dlata += (dtmp1 - morph[nm].lata0);
            dlatb += (dtmp2 - morph[nm].latb0);
            dlatc += (dtmp3 - morph[nm].latc0);
            nread++;
          }
        }
      }

      fclose(fone);

      if(nread>0) {
        dlata/=(double)nread;
        dlatb/=(double)nread;
        dlatc/=(double)nread;
      }

      dlata   = (dlata/morph[nm].lata0*100.0);
      dlatb   = (dlatb/morph[nm].latb0*100.0);
      dlatc   = (dlatc/morph[nm].latc0*100.0);

      printf("FOMC  cycle %6d morph %d rmsd %9.3lf [Angstrom] dlat %8.2lf%8.2lf%8.2lf [%%]\n", 0, nm, 10.0*morph[nm].mcrmsd, dlata,dlatb,dlatc);

      morph[nm].mcdlat=fabs(dlata);
      if(fabs(dlatb)>morph[nm].mcdlat) morph[nm].mcdlat=fabs(dlatb);
      if(fabs(dlatc)>morph[nm].mcdlat) morph[nm].mcdlat=fabs(dlatc);

      fflush(NULL);

    } else {
      printf("# FOMC timeout gmx mdrun at cycle zero of morph %d returned: %d\n", nm, iret);
      exit(1);
    }

    sprintf(mycommand,"rm -f %s.top %s.mdp %s.tpr %s.trr %s.cpt %s_prev.cpt %s.edr %s.log nj-%s.trr %s.mdp e-%s.xvg p-%s.pdb %s.gro rmsd-%s.xvg",
            mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,
            mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,mdbasename);
    iret=system(mycommand);
  }

  mcrmsdcur=morph[0].mcrmsd;
  for(j=1;j<pars.n_xtalpdbs;j++) if(morph[j].mcrmsd > mcrmsdcur) mcrmsdcur = morph[j].mcrmsd;

  mcdlatcur=morph[0].mcdlat;
  for(j=1;j<pars.n_xtalpdbs;j++) if(morph[j].mcdlat > mcdlatcur) mcdlatcur = morph[j].mcdlat;

  mcdtot_new = (mcrmsdcur/pars.rmsthreshold)*(mcrmsdcur/pars.rmsthreshold) + (mcdlatcur/pars.abcthreshold)*(mcdlatcur/pars.abcthreshold);

  printf("FOMC  cycle %6d current dtot: %12.4f\n", 0, mcdtot_new);

  mcdtot_t0=mcdtot_new;
  best_dtot=mcdtot_new;
  cycle_of_best_dtot=0;

  printf("FOMC  cycle %6d best dtot = %12.4f from cycle %d\n", 0, best_dtot, cycle_of_best_dtot);

  fflush(NULL);


/*#################################################*/
/*########## begin main loop ######################*/
/*#################################################*/

  nreject=0;
  mcdtot_old=mcdtot_new;
  n_cycle=0;
  n_ttest=0;
  n_accepted=0;
  
  do {

      n_cycle++;
      n_ttest++;

/*****  adjust MC temperature *****/

      if(n_ttest==TEST_INTERVAL) {

        if(n_accepted<TEST_INTERVAL*0.3) {
          pars.mctemperature *= 1.2;
          printf("FOMC  cycle %6d raised MC temperature to %12.4e\n", n_cycle, pars.mctemperature);
        } else if(n_accepted>TEST_INTERVAL*0.7) {
          pars.mctemperature *= 0.8;
          printf("FOMC  cycle %6d lowered MC temperature to %12.4e\n", n_cycle, pars.mctemperature);
        }

        n_ttest=0;
        n_accepted=0;
      }

/***** apply random change to one randomly chosen LJ parameter *****/
      
      frand1=((double)rand()/RAND_MAX);
      
      if(frand1<0.5) {

        changetype=TEPS;
        frand2=((double)rand()/RAND_MAX);
        tmpmodindex=(int)(input_atp.ntypes_mod*frand2);
        modindex=input_atp.i_type_mod[tmpmodindex];
        origeps=input_atp.epsilon[modindex];
        frand3=((double)rand()/RAND_MAX);
        pvar=(0.1-frand3*0.2)*(ljpars.epsmax[modindex]-ljpars.epsmin[modindex]);
        dtmp=input_atp.epsilon[modindex]+pvar;
        if(dtmp<ljpars.epsmin[modindex]) dtmp=ljpars.epsmin[modindex]+sqrt(pvar*pvar)/10.0;
        if(dtmp>ljpars.epsmax[modindex]) dtmp=ljpars.epsmax[modindex]-sqrt(pvar*pvar)/10.0;
        printf("FOMC  ########################################################################\n");
        printf("FOMC  cycle %6d PARAMETER CHANGE: eps %4s   %12.4e -> ", n_cycle, input_atp.atomtype[modindex], input_atp.epsilon[modindex]);
        input_atp.epsilon[modindex]=dtmp;
        printf("%12.4e\n", input_atp.epsilon[modindex]);

      } else {

        changetype=TSIG;
        frand2=((double)rand()/RAND_MAX);
        tmpmodindex=(int)(input_atp.ntypes_mod*frand2);
        modindex=input_atp.i_type_mod[tmpmodindex];
        origsig=input_atp.sigma[modindex];
        frand3=((double)rand()/RAND_MAX);
        pvar=(0.1-frand3*0.2)*(ljpars.sigmax[modindex]-ljpars.sigmin[modindex]);
        dtmp=input_atp.sigma[modindex]+pvar;
        if(dtmp<ljpars.sigmin[modindex]) dtmp=ljpars.sigmin[modindex]+sqrt(pvar*pvar)/10.0;
        if(dtmp>ljpars.sigmax[modindex]) dtmp=ljpars.sigmax[modindex]-sqrt(pvar*pvar)/10.0;
        printf("FOMC  ########################################################################\n");
        printf("FOMC  cycle %6d PARAMETER CHANGE: sig %4s   %12.4e -> ", n_cycle, input_atp.atomtype[modindex], input_atp.sigma[modindex]);
        input_atp.sigma[modindex]=dtmp;
        printf("%12.4e\n", input_atp.sigma[modindex]);

      }

/* ##### use the resulting FF parameters to do MD with each polymorph #####*/

      GMXSTATUS=0;

      for(nm=0;nm<pars.n_xtalpdbs;nm++) {

        sprintf(mdbasename,"%s-%d-%06d-%02d",pars.runid, pars.rseed, n_cycle, nm+1);
        sprintf(topfn,"%s-%d-%06d-%02d.top",pars.runid, pars.rseed, n_cycle, nm+1);
        sprintf(mdpfn,"md%02d.mdp",nm+1);
        
        if((ftop=fopen(topfn,"w"))==NULL) error_msg("cannot write topology file");
        else {
          fprintf(ftop, "[ defaults ]\n");
          fprintf(ftop, "1   2   yes   0.5  0.8333\n\n");
          fprintf(ftop, "[ atomtypes ]\n");
          for(k=0;k<input_atp.n_types;k++) {
            fprintf(ftop," %-6s  %-6s   0.00000  0.00000   A     %14.5e  %14.5e\n",
              input_atp.atomtype[k],input_atp.atomtype[k],input_atp.sigma[k],input_atp.epsilon[k]);
          }
          fprintf(ftop, "\n");
          fprintf(ftop, "#include \"moleculetypes.itp\"\n");
          fprintf(ftop, "\n");
          fprintf(ftop, "[ system ]\n");
          fprintf(ftop, "system\n\n");
          fprintf(ftop, "[ molecules ]\n");
          for(j=0;j<morph[nm].nrestot;j++) fprintf(ftop, "%s 1\n", morph[nm].resname[j]);
          fclose(ftop);
        }  

        sprintf(mycommand,"gmx grompp -f %s -c %s -p %s -o %s -po %s.mdp", mdpfn, morph[nm].pdbfn, topfn, mdbasename, mdbasename);
        if(system(mycommand)!=0) error_msg("MC: grompp command failed");
        sprintf(mycommand,"timeout %d gmx mdrun -deffnm %s", pars.mdtimeout, mdbasename);
        iret=system(mycommand);

        if(iret==0) {

/* ##### determine RMSD ##### */

          sprintf(mycommand,"echo \"0\" | gmx trjconv -f %s.trr -s %s -pbc nojump -o nj-%s.trr", mdbasename,morph[nm].pdbfn,mdbasename);
          iret=system(mycommand);
          sprintf(mycommand,"gmx rms -b 100 -f nj-%s.trr -s %s -o rmsd-%s.xvg -xvg none -n %s",
                  mdbasename,morph[nm].pdbfn,mdbasename,morph[nm].ndxfn);
          iret=system(mycommand);

          sprintf(xvgfn,"rmsd-%s.xvg",mdbasename);
          if((fxvg=fopen(xvgfn,"r"))==NULL) error_msg("cannot open xvg file");
          else {
            nread=0;
            morph[nm].mcrmsd=0.0;
            while(fscanf(fxvg,"%lf %lf",&dtmp1,&dtmp2)==2) {
              morph[nm].mcrmsd+=dtmp2;
              nread++;
            }
            if(nread>0) morph[nm].mcrmsd /= (double)nread;
          }
          fclose(fxvg);

/* ##### determine avg lattice parameters ##### */

          dlata=0.0;
          dlatb=0.0;
          dlatc=0.0;
          nread=0;
          
          sprintf(onefn,"p-%s.pdb", mdbasename);

          sprintf(mycommand,"gmx trjconv -b 100 -f %s.trr -s %s.tpr -o %s -xvg none -n 1.ndx",mdbasename,mdbasename,onefn);
          iret=system(mycommand);

          if((fone=fopen(onefn,"r"))==NULL) error_msg("cannot open 1.pdb file");

          while((nitemsread=getline(&oneline,&len0,fone))>0) {
            if(strstr(oneline,"CRYST1")!=NULL) {
              if(sscanf(oneline,"%s %lf %lf %lf %lf %lf %lf",sdummy,&dtmp1,&dtmp2,&dtmp3,&dtmp4,&dtmp5,&dtmp6)!=7) {
                error_msg("cannot read CRYST1 line in pdb file");
              } else {
                dlata += (dtmp1-morph[nm].lata0);
                dlatb += (dtmp2-morph[nm].latb0);
                dlatc += (dtmp3-morph[nm].latc0);
                nread++;
              }
            }
          }

          fclose(fone);

          if(nread>0) {
            dlata/=(double)nread;
            dlatb/=(double)nread;
            dlatc/=(double)nread;
          }

          dlata   = (dlata/morph[nm].lata0*100.0);
          dlatb   = (dlatb/morph[nm].latb0*100.0);
          dlatc   = (dlatc/morph[nm].latc0*100.0);

          morph[nm].mcdlat=fabs(dlata);
          if(fabs(dlatb)>morph[nm].mcdlat) morph[nm].mcdlat=fabs(dlatb);
          if(fabs(dlatc)>morph[nm].mcdlat) morph[nm].mcdlat=fabs(dlatc);
          
          printf("FOMC  cycle %6d morph %d  rmsd %9.3lf [Angstrom]  dlatmax %8.2lf [%%] (%.2lf %.2lf %.2lf)\n",
                 n_cycle, nm, 10.0*morph[nm].mcrmsd, morph[nm].mcdlat,dlata,dlatb,dlatc);

        } else {
          
          printf("FOMC timeout gmx mdrun returned: %d\n", iret);
          exit(1);
          GMXSTATUS -= 1;
        }

        sprintf(mycommand,"rm -f %s.mdp %s.tpr %s.trr %s.cpt %s_prev.cpt %s.edr %s.log nj-%s.trr %s.mdp e-%s.xvg %s.gro %s rmsd-%s.xvg",
          mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,mdbasename,
          mdbasename,mdbasename,mdbasename,onefn,mdbasename);
        iret=system(mycommand);
      }

      if(GMXSTATUS==0) {

        mcrmsdcur = morph[0].mcrmsd;
        for(j=1;j<pars.n_xtalpdbs;j++) if(morph[j].mcrmsd>mcrmsdcur) mcrmsdcur = morph[j].mcrmsd;

        mcdlatcur = morph[0].mcdlat;
        for(j=1;j<pars.n_xtalpdbs;j++) if(morph[j].mcdlat>mcdlatcur) mcdlatcur = morph[j].mcdlat;

        mcdtot_new = (mcrmsdcur/pars.rmsthreshold)*(mcrmsdcur/pars.rmsthreshold) + (mcdlatcur/pars.abcthreshold)*(mcdlatcur/pars.abcthreshold);

        if(mcdtot_new>=mcdtot_old) bf=exp(-(mcdtot_new-mcdtot_old)/pars.mctemperature);
        frand1=((double)rand()/RAND_MAX);

        if(mcdtot_new<mcdtot_old || bf>frand1) {

/* #### report ##### */

          if(mcdtot_new<mcdtot_old) printf("FOMC  cycle %6d result dtot drops %9.2e -> %9.2e", n_cycle, mcdtot_old, mcdtot_new);
          else                      printf("FOMC  cycle %6d result you lucky  %9.2e -> %9.2e", n_cycle, mcdtot_old, mcdtot_new);
          if(changetype==TEPS)      printf("  eps %4s -> %12.4e\n", input_atp.atomtype[modindex],input_atp.epsilon[modindex]);
          if(changetype==TSIG)      printf("  sig %4s -> %12.4e\n", input_atp.atomtype[modindex],input_atp.sigma[modindex]);

/* #### report ##### */

          n_accepted++;
          nreject=0;
          mcdtot_old=mcdtot_new;

          if(mcdtot_new<best_dtot) {
            
            best_dtot=mcdtot_new;
            cycle_of_best_dtot=n_cycle;
            
            if((ftop=fopen(besttopfn,"w"))==NULL) error_msg("cannot write to best topology file"); 
            else {
              fprintf(ftop, "[ defaults ]\n");
              fprintf(ftop, "1   2   yes   0.5  0.8333\n\n");
              fprintf(ftop, "[ atomtypes ]\n");
              for(k=0;k<input_atp.n_types;k++) {
                fprintf(ftop," %-6s  %-6s   0.00000  0.00000   A     %14.5e  %14.5e ; %d\n",
                  input_atp.atomtype[k],input_atp.atomtype[k],input_atp.sigma[k],input_atp.epsilon[k], input_atp.modflag[k]);
              }
              fprintf(ftop, "\n");
            }
            fclose(ftop);
          }

        } else {

          printf("FOMC  cycle %6d result dtot poor  %9.2e -> %9.2e", n_cycle, mcdtot_old, mcdtot_new);
          if(changetype==TEPS)      printf("  eps %4s -> %12.4e\n", input_atp.atomtype[modindex],input_atp.epsilon[modindex]);
          if(changetype==TSIG)      printf("  sig %4s -> %12.4e\n", input_atp.atomtype[modindex],input_atp.sigma[modindex]);

          nreject++;
          if(changetype==TEPS) input_atp.epsilon[modindex]=origeps;
          if(changetype==TSIG) input_atp.sigma[modindex]=origsig;
        }

      } else {

        printf("FOMC  cycle %6d GMXSTATUS %d is bad\n", n_cycle, GMXSTATUS);
        printf("FOMC  cycle %6d current dtot: %12.4f\n", n_cycle, mcdtot_old);
        nreject++;
        if(changetype==TEPS) input_atp.epsilon[modindex]=origeps;
        if(changetype==TSIG) input_atp.sigma[modindex]=origsig;
      }

      if(nreject==20) {
        printf("FOMC  cycle %6d WE GOT STUCK - reset all parameters to initial values with small random variations!\n", n_cycle);
        printf("FOMC  cycle %6d current dtot: %12.4f\n", n_cycle, mcdtot_old);

        for(l=0;l<input_atp.n_types;l++) {
          frand1 = 0.01 - 0.02*((double)rand()/RAND_MAX);
          input_atp.epsilon[l] = ljpars.eps0[l]+frand1*(ljpars.epsmax[l]-ljpars.epsmin[l]);
        }
        for(l=0;l<input_atp.n_types;l++) {
          frand1 = 0.01 - 0.02*((double)rand()/RAND_MAX);
          input_atp.sigma[l]   = ljpars.sig0[l]+frand1*(ljpars.sigmax[l]-ljpars.sigmin[l]);
        }
        nreject=0;
        n_ttest=0;
        n_accepted=0;
        mcdtot_old=mcdtot_t0;
      }
      
      for(nm=0;nm<pars.n_xtalpdbs;nm++) {
        sprintf(mycommand,"rm -r %s-%d-%06d-%02d.top",pars.runid, pars.rseed, n_cycle, nm+1);
        iret=system(mycommand);
      }
      
      printf("FOMC  cycle %6d best dtot = %12.4e from cycle %d  last accepted: %12.4e  current: %12.4e\n",
             n_cycle, best_dtot, cycle_of_best_dtot, mcdtot_old, mcdtot_new);

      fflush(NULL);

  } while(n_cycle<pars.n_mccycles && best_dtot>pars.max_dtot);


/*#################################################*/
/*########## end main loop ########################*/
/*#################################################*/

  sprintf(mycommand,"cat moleculetypes.itp >> %s", besttopfn);
  iret=system(mycommand);

  printf("FOMC  best topology copied to: %s\n", besttopfn);

  gettimeofday(&tvEnd, NULL);
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
  printf("FOMC  total real time elapsed: %12.2lf hours\n", (double)tvDiff.tv_sec/3600.0);
  
  exit(0);
}
