/* ----------------------------------------------------------------------
   UPDTAE by  Himangsu Bhaumik
   
   In this file the VSP-BKS Silica potential developed keeping the parameter
   of the Paper I. Saika-Voivod, F. Sciortino, and P. H. Poole, Phys. Rev. E 69, 041503 (2004)
   I renamed the pair style as: 
   pair_style     vspbks
   
   

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_vspbks.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairVspbks::PairVspbks(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  writedata = 1;
  ftable = NULL;
}

/* ---------------------------------------------------------------------- */

PairVspbks::~PairVspbks()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(a);
    memory->destroy(rho);
    memory->destroy(c);
    memory->destroy(rhoinv);
    memory->destroy(buck1);
    memory->destroy(buck2);
    memory->destroy(offset);

    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);

    
    memory->destroy(eps);
    memory->destroy(sig);
    memory->destroy(Dp);
    memory->destroy(Ep);
    memory->destroy(Fp);
  }
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

void PairVspbks::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double rsq,r2inv,r6inv,forcecoul,forcebuck,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double r,rexp;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;   
	// (rsq < rsmsq[itype][jtype]) {
	//forcecoul = 0.0;
	//
        //else 
	//if (rsq < cut_coulsq) {
	if (rsq < rsmsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	   
	    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;


// 	if (rsq < rsmsq[itype][jtype]){
// 	  r = sqrt(rsq);
// 	  forcebuck = -(2.0*am[itype][jtype]*rsq+bm[itype][jtype]*r);
// 	}
// 
// 	
//        else if (rsq < cut_ljsq[itype][jtype]) {
//          r = sqrt(rsq);
//          r6inv = r2inv*r2inv*r2inv;
//          rexp = exp(-r*rhoinv[itype][jtype]);
//          forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
//        } else forcebuck = 0.0;

	if (rsq < rsmsq){
	  r = sqrt(rsq);
	  r6inv = r2inv*r2inv*r2inv;
	  rexp = exp(-r*rhoinv[itype][jtype]);
	  forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
	  forcelj = r6inv * (lj1[itype][jtype]*pow(r6inv,4.0) 
			     - lj2[itype][jtype]);
	  forcebuck = forcebuck + forcelj;
	}
 
	
	else if (rsq < cut_ljsq[itype][jtype]) {
	  r = sqrt(rsq);
	  rdiffR = r-Rc;
	  forcebuck = 5*Dp[itype][jtype]*pow(rdiffR,4.0) 
	    + 4*Ep[itype][jtype]*pow(rdiffR,3.0) 
	    + 3*Fp[itype][jtype]*pow(rdiffR,2.0);
	  forcebuck = -forcebuck*r;
	} else forcebuck = 0.0;

 
        fpair = (forcecoul + factor_lj*forcebuck) * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
	
        if (eflag) {
	  //if (rsq < rsmsq[itype][jtype]) ecoul = 0.0;

          //else 
	  //if (rsq < cut_coulsq) 
	      if(rsq <rsmsq)
	      {
            if (!ncoultablebits || rsq <= tabinnersq)
	      
		ecoul = prefactor*erfc;	
            else {
              table = etable[itable] + fraction*detable[itable];
              if(rsq <rsmsq)ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
	    } else ecoul = 0.0;

	  if (rsq < rsmsq){
	    //evdwl = am[itype][jtype]*rsq+bm[itype][jtype]*r+cm[itype][jtype];
	    
	    evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv 
	       -  offset[itype][jtype];
            
	    evdwl = evdwl + r6inv*(lj3[itype][jtype]*pow(r6inv,4.0)
				   -lj4[itype][jtype]);
	    
	    evdwl *= factor_lj;
	    
	  }
	  
          else if (rsq < cut_ljsq[itype][jtype]) {
	    rdiffR = sqrt(rsq)-Rc;
            evdwl = Dp[itype][jtype]*pow(rdiffR,5.0) 
	      + Ep[itype][jtype]*pow(rdiffR,4.0) 
	      + Fp[itype][jtype]*pow(rdiffR,3.0);
          } else evdwl = 0.0;
        }
	//if((itype==1 && jtype==2) || (itype==2 && jtype==1))printf("%e %e %e %e\n",sqrt(rsq),forcecoul,forcebuck,forcelj);
	//if((itype==1 && jtype==2) || (itype==2 && jtype==1))printf("%e %e %e %e %e\n",sqrt(rsq),fpair,evdwl,ecoul,(evdwl+ecoul));
	//printf("%e %e %e %e %e\n",sqrt(rsq),fpair,evdwl,ecoul,(evdwl+ecoul));
	//if(itype==1 && jtype==1)printf("%e %e %e %e %e\n",sqrt(rsq),fpair,evdwl,ecoul,(evdwl+ecoul));
	//if(itype==2 && jtype==2)printf("%e %e %e %e %e\n",sqrt(rsq),fpair,evdwl,ecoul,(evdwl+ecoul));
	//	 printf("prefactor=%e\n",prefactor);

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  } 

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairVspbks::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(a,n+1,n+1,"pair:a");
  memory->create(rho,n+1,n+1,"pair:rho");
  memory->create(c,n+1,n+1,"pair:c");
  memory->create(rhoinv,n+1,n+1,"pair:rhoinv");
  memory->create(buck1,n+1,n+1,"pair:buck1");
  memory->create(buck2,n+1,n+1,"pair:buck2");
  memory->create(offset,n+1,n+1,"pair:offset");
  
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  
  memory->create(eps,n+1,n+1,"pair:eps");
  memory->create(sig,n+1,n+1,"pair:sig");
  memory->create(Dp,n+1,n+1,"pair:Dp");
  memory->create(Ep,n+1,n+1,"pair:Ep");
  memory->create(Fp,n+1,n+1,"pair:Fp");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairVspbks::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);
  if (narg == 1) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairVspbks::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double a_one = force->numeric(FLERR,arg[2]);
  double rho_one = force->numeric(FLERR,arg[3]);
  if (rho_one <= 0) error->all(FLERR,"Incorrect args for pair coefficients");
  double c_one = force->numeric(FLERR,arg[4]);

  double cut_lj_one = cut_lj_global;
  if (narg == 6) cut_lj_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      rho[i][j] = rho_one;
      c[i][j] = c_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairVspbks::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j],cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  rhoinv[i][j] = 1.0/rho[i][j];
  buck1[i][j] = a[i][j]/rho[i][j];
  buck2[i][j] = 6.0*c[i][j];
  
  if (offset_flag) {
    double rexp = exp(-cut_lj[i][j]/rho[i][j]);
    offset[i][j] = a[i][j]*rexp - c[i][j]/pow(cut_lj[i][j],6.0);
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  a[j][i] = a[i][j];
  c[j][i] = c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];


  
  
  eps[1][1] = 0.0;
  eps[1][2] = 0.0714407145;
  eps[2][1] = 0.0714407145;
  eps[2][2] = 0.02423786381;
  
  sig[1][1] = 0.0;
  sig[1][2] = 1.313635;
  sig[2][1] = 1.313635;
  sig[2][2] = 1.779239;


  lj1[i][j] = 120.0 * eps[i][j] * pow(sig[i][j],30.0);
  lj2[i][j] = 24.0 * eps[i][j] * pow(sig[i][j],6.0);
  lj3[i][j] = 4.0 * eps[i][j] * pow(sig[i][j],30.0);
  lj4[i][j] = 4.0 * eps[i][j] * pow(sig[i][j],6.0);

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  
  
  Dp[1][1] = -0.0338751180316;
  Dp[1][2] = 0.0175621791329;
  Dp[2][1] = 0.0175621791329;
  Dp[2][2] = -0.007651894017;
  
  
  Ep[1][1] = -0.169552412209;
  Ep[1][2] = 0.088284968753;
  Ep[2][1] = 0.088284968753;
  Ep[2][2] = -0.0377950981;
  
  
  Fp[1][1] = -0.34310602604;
  Fp[1][2] = 0.17753394513;
  Fp[2][1] = 0.17753394513;
  Fp[2][2] = -0.07794280028;
  
    
  rsm = 7.7476;
  rsmsq = rsm*rsm;
  Rc=10.0;
  
  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double rho1 = rho[i][j];
    double rho2 = rho1*rho1;
    double rho3 = rho2*rho1;
    double rc = cut_lj[i][j];
    double rc2 = rc*rc;
    double rc3 = rc2*rc;
    etail_ij = 2.0*MY_PI*all[0]*all[1]*
      (a[i][j]*exp(-rc/rho1)*rho1*(rc2 + 2.0*rho1*rc + 2.0*rho2) -
       c[i][j]/(3.0*rc3));
    ptail_ij = (-1/3.0)*2.0*MY_PI*all[0]*all[1]*
      (-a[i][j]*exp(-rc/rho1)*
       (rc3 + 3.0*rho1*rc2 + 6.0*rho2*rc + 6.0*rho3) + 2.0*c[i][j]/rc3);
  }

  return cut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairVspbks::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style buck/coul/long/mod requires atom attribute q");

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  neighbor->request(this);
  
  // setup force tables

  if (ncoultablebits) init_tables(cut_coul,NULL);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairVspbks::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a[i][j],sizeof(double),1,fp);
        fwrite(&rho[i][j],sizeof(double),1,fp);
        fwrite(&c[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairVspbks::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&a[i][j],sizeof(double),1,fp);
          fread(&rho[i][j],sizeof(double),1,fp);
          fread(&c[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&rho[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairVspbks::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairVspbks::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
    fread(&ncoultablebits,sizeof(int),1,fp);
    fread(&tabinner,sizeof(double),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairVspbks::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,a[i][i],rho[i][i],c[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairVspbks::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,
              a[i][j],rho[i][j],c[i][j],cut_lj[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairVspbks::single(int i, int j, int itype, int jtype,
                                double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,r,rexp,grij,expm2,t,erfc,prefactor;
  double fraction,table,forcecoul,forcebuck,phicoul,phibuck;
  int itable;

  r2inv = 1.0/rsq;

  // if (rsq < rsmsq[itype][jtype]) {
  //  	  forcecoul = 0.0;
  // }
  // else

    if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq) {
      r = sqrt(rsq);
      grij = g_ewald * r;
      expm2 = exp(-grij*grij);
      t = 1.0 / (1.0 + EWALD_P*grij);
      erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
      prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
      
      if(rsq <rsmsq)forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
      else forcecoul = prefactor * (EWALD_F*grij*expm2);
      if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
    } else {
      union_int_float_t rsq_lookup;
      rsq_lookup.f = rsq;
      itable = rsq_lookup.i & ncoulmask;
      itable >>= ncoulshiftbits;
      fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
      table = ftable[itable] + fraction*dftable[itable];
      forcecoul = atom->q[i]*atom->q[j] * table;
      if (factor_coul < 1.0) {
        table = ctable[itable] + fraction*dctable[itable];
        prefactor = atom->q[i]*atom->q[j] * table;
        forcecoul -= (1.0-factor_coul)*prefactor;
      }
    }
  } else forcecoul = 0.0;


//  if (rsq < rsmsq[itype][jtype]){
// 	  r = sqrt(rsq);
// 	  forcebuck = -(2.0*am[itype][jtype]*rsq+bm[itype][jtype]*r);
// 	}
// 
//  
// else if (rsq < cut_ljsq[itype][jtype]) {
//    r6inv = r2inv*r2inv*r2inv;
//    r = sqrt(rsq);
//    rexp = exp(-r*rhoinv[itype][jtype]);
//    forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
//  } else forcebuck = 0.0;
 
    if (rsq < rsmsq){
	  r = sqrt(rsq);
	  r6inv = r2inv*r2inv*r2inv;
	  rexp = exp(-r*rhoinv[itype][jtype]);
	  forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
	  forcelj = r6inv * (lj1[itype][jtype]*pow(r6inv,4.0) 
			     - lj2[itype][jtype]);
	  forcebuck = forcebuck + forcelj;
	}
 
	
	else if (rsq < cut_ljsq[itype][jtype]) {
	  r = sqrt(rsq);
	  rdiffR = r-Rc;
	  forcebuck = 5*Dp[itype][jtype]*pow(rdiffR,4.0) 
	    + 4*Ep[itype][jtype]*pow(rdiffR,3.0) 
	    + 3*Fp[itype][jtype]*pow(rdiffR,2.0);
	  forcebuck = -forcebuck*r;
	} else forcebuck = 0.0;
    
    
    fforce = (forcecoul + factor_lj*forcebuck) * r2inv;
  
  double eng = 0.0;

  //if (rsq < rsmsq[itype][jtype]) phicoul = 0.0;

  //else 
  //if (rsq < cut_coulsq) {
  if (rsq < rsmsq) {
    if (!ncoultablebits || rsq <= tabinnersq)
       phicoul = prefactor*erfc;
    else {
      table = etable[itable] + fraction*detable[itable];
      phicoul = atom->q[i]*atom->q[j] * table;
    }
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    eng += phicoul;
  }
  if (rsq < rsmsq){
    //phibuck = am[itype][jtype]*rsq+bm[itype][jtype]*r+cm[itype][jtype];
    phibuck = a[itype][jtype]*rexp - c[itype][jtype]*r6inv 
	       -  offset[itype][jtype];
    phibuck = phibuck + r6inv*(lj3[itype][jtype]*pow(r6inv,4.0)
				   -lj4[itype][jtype]);

    eng += factor_lj*phibuck;
  }

  else if (rsq < cut_ljsq[itype][jtype]) {
    //phibuck = a[itype][jtype]*rexp - c[itype][jtype]*r6inv 
    //  - Ucij[itype][jtype] - offset[itype][jtype];

    rdiffR = sqrt(rsq)-Rc;
            phibuck = Dp[itype][jtype]*pow(rdiffR,5.0) 
	      + Ep[itype][jtype]*pow(rdiffR,4.0) 
	      + Fp[itype][jtype]*pow(rdiffR,3.0);
    eng += factor_lj*phibuck;
  }
  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairVspbks::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  return NULL;
}
