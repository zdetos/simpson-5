/*
    Hamiltonian calculation routines
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
				  2010 Zdenek Tosner

    This file is part of the SIMPSON General NMR Simulation Package

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
    Setup and calculations for the Hamiltonian for the spinsystem.
    Creates the spin and space tensors, rotates them and calculates
    fourier components.
    
    Called from readsys.c where the spin-system is set up,
    sim.c where the crystallite averaging takes place,
    and pulse.c where the pulse propagation is performed.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "matrix.h"
#include "cm.h"
#include "blockdiag.h"
#include "ham.h"
#include "wigner.h"
#include "spinsys.h"
#include "defs.h"

#ifdef INTEL_MKL
#include "mkl.h"
#include "mkl_spblas.h"
#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#elif defined(GSL)
#include <gsl/gsl_cblas.h>
#else
#include "cblas.h"
#endif

	/* for acurate timings on windows */
//#include <windows.h>
//#include <winbase.h>


void ham_set_offset(Sim_info *s, Sim_wsp *wsp, double* offset,int used)
{
  if (used) {
	  if (LEN(offset) != s->ss->nchan) {
		  fprintf(stderr,"Error: channels mismatch when setting offset\n");
		  exit(1);
	  }
	  if (wsp->offset) free_double_vector(wsp->offset);
	  wsp->offset = offset;
  } else {
	  if (wsp->offset) free_double_vector(wsp->offset);
	  wsp->offset = NULL;
  }
}

int shift_exist(Sim_info* s,int n)
{
	int i, res=-1;
	for (i=0;i<s->nCS;i++) {
		if (s->CS[i]->nuc == n) {
			res = i;
			break;
		}
	}
	return res;
}

int dipole_exist(Sim_info* s,int n1, int n2)
{
	int dum, i, res=-1;
	if (n1>n2) {
		dum = n1; n1 = n2; n2 = dum;
	}
	for (i=0; i<s->nDD; i++) {
		if ( (s->DD[i]->nuc[0] == n1) && (s->DD[i]->nuc[1] == n2) ) {
			res = i;
			break;
		}
	}
	return res;
}

int jcoupling_exist(Sim_info* s,int n1, int n2)
{
	int dum, i, res=-1;
	if (n1>n2) {
		dum = n1; n1 = n2; n2 = dum;
	}
	for (i=0; i<s->nJ; i++) {
		if ( (s->J[i]->nuc[0] == n1) && (s->J[i]->nuc[1] == n2) ) {
			res = i;
			break;
		}
	}
	return res;
}

int quadrupole_exist(Sim_info* s,int n)
{
	int i, res=-1;
	for (i=0;i<s->nQ;i++) {
		if (s->Q[i]->nuc == n) {
			res = i;
			break;
		}
	}
	return res;
}

int gtensor_exist(Sim_info* s,int n)
{
	int i, res=-1;
	for (i=0;i<s->nG;i++) {
		if (s->G[i]->nuc == n) {
			res = i;
			break;
		}
	}
	return res;
}

int hyperfine_exist(Sim_info* s,int n1, int n2)
{
	int i, res=-1;
	for (i=0; i<s->nDD; i++) {
		if ( (s->DD[i]->nuc[0] == n1) && (s->DD[i]->nuc[1] == n2) ) {
			res = i;
			break;
		}
		if ( (s->DD[i]->nuc[0] == n2) && (s->DD[i]->nuc[1] == n1) ) {
			res = i;
			break;
		}
	}
	return res;
}

int heisenberg_exist(Sim_info* s,int n1, int n2)
{
	int i, res=-1;
	for (i=0;i<s->nHEX;i++) {
		if ( (s->HEx[i]->electron[0] == n1) && (s->HEx[i]->electron[1] == n2) ) {
			res = i;
			break;
		}
		if ( (s->HEx[i]->electron[0] == n2) && (s->HEx[i]->electron[1] == n1) ) {
			res = i;
			break;
		}
	}
	return res;
}

int edipole_exist(Sim_info* s,int n1, int n2)
{
	int i, res=-1;
	for (i=0;i<s->nEDD;i++) {
		if ( (s->EDD[i]->electron[0] == n1) && (s->EDD[i]->electron[1] == n2) ) {
			res = i;
			break;
		}
		if ( (s->EDD[i]->electron[0] == n2) && (s->EDD[i]->electron[1] == n1) ) {
			res = i;
			break;
		}
	}
	return res;
}

/* this can be just wsp variable, set by turnon/turnoff commands */
int ham_ischanged(Sim_info *s, Sim_wsp *wsp)
{
  int i;

  /* return true if offset is set*/
  if (wsp->offset) return 1;

  /* return true if any interactions are disabled */
  for (i=0;i<s->nCS;i++) {
    if (!wsp->CS_used[i]) return 1;
  }
  for (i=0;i<s->nDD;i++) {
    if (!wsp->DD_used[i]) return 1;
  }
  for (i=0;i<s->nQ;i++) {
    if (!wsp->Q_used[i]) return 1;
  }
  for (i=0;i<s->nJ;i++) {
    if (!wsp->J_used[i]) return 1;
  }
  for (i=0;i<s->nG;i++) {
    if (!wsp->G_used[i]) return 1;
  }
  for (i=0;i<s->nHF;i++) {
    if (!wsp->HF_used[i]) return 1;
  }
  return 0;
}

/* create matrix with all interactions to be subtracted from total interactions */
void ham_turnoff(Sim_info *s, Sim_wsp *wsp,char* name)
{
  int n1,n2,i;
  Shift *csptr;
  Dipole *ddptr;
  Jcoupling *jptr;
  Quadrupole *qptr;
  mat_double *mx;
  blk_mat_double *blk_mx;

  n1 = s->matdim;
  if (wsp->Hiso_off == NULL) {
	  wsp->Hiso_off = create_blk_mat_double_copy(wsp->Hiso);
	  blk_dm_zero(wsp->Hiso_off);
  }
  if (s->Hassembly) {
	  for (i=0; i<5; i++) {
		  if (wsp->HQ_off[i] == NULL) {
			  wsp->HQ_off[i] = create_blk_mat_double_copy(wsp->HQ[i]);
			  blk_dm_zero(wsp->HQ_off[i]);
		  }
	  }
  }

  if (!strcmp(name,"all")) {
	  memset(wsp->CS_used,0,sizeof(int)*s->nCS);
	  memset(wsp->DD_used,0,sizeof(int)*s->nDD);
	  memset(wsp->Q_used,0,sizeof(int)*s->nQ);
	  memset(wsp->J_used,0,sizeof(int)*s->nJ);
	  memset(wsp->G_used,0,sizeof(int)*s->nG);
	  memset(wsp->HF_used,0,sizeof(int)*s->nHF);
	  memset(wsp->EDD_used,0,sizeof(int)*s->nEDD);
	  memset(wsp->HEx_used,0,sizeof(int)*s->nHEX);
	  blk_dm_copy(wsp->Hiso_off,wsp->Hiso);
	  if (s->Hassembly) {
		  for (i=0; i<5; i++) {
			  blk_dm_copy(wsp->HQ_off[i],wsp->HQ[i]);
		  }
	  }
	  wsp->Nint_off = s->nCS + s->nDD + s->nQ + s->nJ + s->nG + s->nHF + s->nHEX + s->nEDD;
     return;
  }

  if (sscanf(name,"shift_%d",&n1) == 1) {
	  for (i=0; i<s->nCS;i++) {
		  if (wsp->CS[i]->nuc == n1) {
			  if (wsp->CS_used[i]) {
				  csptr = wsp->CS[i];
				  if (csptr->T == NULL) {
					  /* MUTEX LOCK */
					  csptr->T = Iz_ham(s,n1);
					  /* MUTEX UNLOCK */
				  }
				  wsp->Nint_off++;
				  mx = dm_change_basis_2(csptr->T,wsp->Hiso->basis,s);
				  if (fabs(csptr->iso)>TINY) blk_dm_multod_diag(wsp->Hiso_off,mx,csptr->iso);
				  if (s->Hassembly && (fabs(csptr->delta)>TINY)) {
					  blk_dm_multod_diag(wsp->HQ_off[0],mx,csptr->Rmol[3].re);
					  blk_dm_multod_diag(wsp->HQ_off[1],mx,csptr->Rmol[4].re);
					  blk_dm_multod_diag(wsp->HQ_off[2],mx,csptr->Rmol[4].im);
					  blk_dm_multod_diag(wsp->HQ_off[3],mx,csptr->Rmol[5].re);
					  blk_dm_multod_diag(wsp->HQ_off[4],mx,csptr->Rmol[5].im);
				  }
				  free_double_matrix(mx);
			  }
			  wsp->CS_used[i] = 0;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnoff - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"quadrupole_%d",&n1) == 1) {
	  for (i=0; i<s->nQ;i++) {
		  if (wsp->Q[i]->nuc == n1) {
			  if (wsp->Q_used[i]) {
				  wsp->Nint_off++;
				  if (s->Hassembly) {
					  qptr = wsp->Q[i];
					  if (qptr->T == NULL) {
						  /* MUTEX LOCK */
						  qptr->T = T20II(s,n1);
						  /* MUTEX UNLOCK */
					  }
					  mx = dm_change_basis_2(qptr->T,wsp->HQ[0]->basis,s);
					  blk_dm_multod_diag(wsp->HQ_off[0],mx,qptr->Rmol[3].re);
					  blk_dm_multod_diag(wsp->HQ_off[1],mx,qptr->Rmol[4].re);
					  blk_dm_multod_diag(wsp->HQ_off[2],mx,qptr->Rmol[4].im);
					  blk_dm_multod_diag(wsp->HQ_off[3],mx,qptr->Rmol[5].re);
					  blk_dm_multod_diag(wsp->HQ_off[4],mx,qptr->Rmol[5].im);
					  free_double_matrix(mx);
				  }
			  }
			  wsp->Q_used[i] = 0;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnoff - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"dipole_%d_%d",&n1, &n2) == 2) {
	  if (n1 > n2) {
		  int dum;
		  dum = n1; n1 = n2; n2 = dum;
	  }
	  for (i=0; i<s->nDD;i++) {
		  if ( (wsp->DD[i]->nuc[0] == n1) && (wsp->DD[i]->nuc[1] == n2) ) {
			  if (wsp->DD_used[i]) {
				  wsp->Nint_off++;
				  if (s->Hassembly) {
					  ddptr = wsp->DD[i];
					  if (ddptr->blk_T == NULL) {
						  /* MUTEX LOCK */
						  if (ss_issame(s->ss,n1,n2)) {
							  ddptr->blk_T = T20(s,n1,n2);
						  } else {
							  ddptr->blk_T = IzIz_sqrt2by3(s,n1,n2);
						  }
						  /* MUTEX UNLOCK */
					  }
					  blk_mx = create_blk_mat_double_copy(wsp->HQ_off[0]);
					  blk_dm_change_basis(blk_mx, ddptr->blk_T,s);
					  blk_dm_multod(wsp->HQ_off[0],blk_mx,ddptr->Rmol[3].re);
					  blk_dm_multod(wsp->HQ_off[1],blk_mx,ddptr->Rmol[4].re);
					  blk_dm_multod(wsp->HQ_off[2],blk_mx,ddptr->Rmol[4].im);
					  blk_dm_multod(wsp->HQ_off[3],blk_mx,ddptr->Rmol[5].re);
					  blk_dm_multod(wsp->HQ_off[4],blk_mx,ddptr->Rmol[5].im);
					  free_blk_mat_double(blk_mx);
				  }
			  }
			  wsp->DD_used[i] = 0;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnoff - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"jcoupling_%d_%d",&n1, &n2) == 2) {
	  if (n1 > n2) {
		  int dum;
		  dum = n1; n1 = n2; n2 = dum;
	  }
	  for (i=0; i<s->nJ;i++) {
		  if ( (wsp->J[i]->nuc[0] == n1) && (wsp->J[i]->nuc[1] == n2) ) {
			  if (wsp->J_used[i]) {
				  wsp->Nint_off++;
				  jptr = wsp->J[i];
				  if (fabs(jptr->iso) > TINY) {
					  /* MUTEX LOCK */
					  if (jptr->blk_Tiso == NULL) {
						  if (ss_issame(s->ss,n1,n2)) {
							  jptr->blk_Tiso = II(s,n1,n2);
						  } else {
							  jptr->blk_Tiso = IzIz_sqrt2by3(s,n1,n2);
							  blk_dm_muld(jptr->blk_Tiso,1.0/SQRT2BY3);
						  }
					  }
					  /* MUTEX UNLOCK */
					  blk_mx = create_blk_mat_double_copy(wsp->Hiso_off);
					  blk_dm_change_basis(blk_mx, jptr->blk_Tiso, s);
					  blk_dm_multod(wsp->Hiso_off,blk_mx,jptr->iso);
					  free_blk_mat_double(blk_mx);
				  }
				  if ( (fabs(jptr->delta) > TINY) && s->Hassembly) {
					  /* MUTEX LOCK */
					  if (jptr->blk_T == NULL) {
						  if (ss_issame(s->ss,n1,n2)) {
							  jptr->blk_T = T20(s,n1,n2);
						  } else {
							  jptr->blk_T = IzIz_sqrt2by3(s,n1,n2);
						  }
					  }
					  /* MUTEX UNLOCK */
					  blk_mx = create_blk_mat_double_copy(wsp->HQ_off[0]);
					  blk_dm_change_basis(blk_mx, jptr->blk_T, s);
					  blk_dm_multod(wsp->HQ_off[0],blk_mx,jptr->Rmol[3].re);
					  blk_dm_multod(wsp->HQ_off[1],blk_mx,jptr->Rmol[4].re);
					  blk_dm_multod(wsp->HQ_off[2],blk_mx,jptr->Rmol[4].im);
					  blk_dm_multod(wsp->HQ_off[3],blk_mx,jptr->Rmol[5].re);
					  blk_dm_multod(wsp->HQ_off[4],blk_mx,jptr->Rmol[5].im);
					  free_blk_mat_double(blk_mx);
				  }
			  }
			  wsp->J_used[i] = 0;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnoff - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"gtensor_%d",&n1) == 1) {
	  Gtensor *gptr;
	  for (i=0; i<s->nG;i++) {
		  if (wsp->G[i]->nuc == n1) {
			  if (wsp->G_used[i]) {
				  gptr = wsp->G[i];
				  if (gptr->T == NULL) {
					  /* MUTEX LOCK */
					  gptr->T = Iz_ham(s,n1);
					  /* MUTEX UNLOCK */
				  }
				  wsp->Nint_off++;
				  mx = dm_change_basis_2(gptr->T,wsp->Hiso->basis,s);
				  if (fabs(gptr->iso)>TINY) blk_dm_multod_diag(wsp->Hiso_off,mx,gptr->iso);
				  if (s->Hassembly && (fabs(gptr->delta)>TINY)) {
					  blk_dm_multod_diag(wsp->HQ_off[0],mx,gptr->Rmol[3].re);
					  blk_dm_multod_diag(wsp->HQ_off[1],mx,gptr->Rmol[4].re);
					  blk_dm_multod_diag(wsp->HQ_off[2],mx,gptr->Rmol[4].im);
					  blk_dm_multod_diag(wsp->HQ_off[3],mx,gptr->Rmol[5].re);
					  blk_dm_multod_diag(wsp->HQ_off[4],mx,gptr->Rmol[5].im);
				  }
				  free_double_matrix(mx);
			  }
			  wsp->G_used[i] = 0;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnoff - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"hyperfine_%d_%d",&n1, &n2) == 2) {
	  i = hyperfine_exist(s,n1,n2); // returns index, if not found then i=-1
	  if ( i >= 0 ) {
		  if (wsp->HF_used[i]) {
			  wsp->Nint_off++;
			  Hyperfine *hfptr = wsp->HF[i];
			  if (fabs(hfptr->iso) > TINY) {
				  /* MUTEX LOCK */
				  if (hfptr->blk_Tiso == NULL) {
					  hfptr->blk_Tiso = IzIz_sqrt2by3(s,n1,n2);
					  blk_dm_muld(hfptr->blk_Tiso,1.0/SQRT2BY3);
				  }
				  /* MUTEX UNLOCK */
				  blk_mx = create_blk_mat_double_copy(wsp->Hiso_off);
				  blk_dm_change_basis(blk_mx, hfptr->blk_Tiso, s);
				  blk_dm_multod(wsp->Hiso_off,blk_mx,hfptr->iso);
				  free_blk_mat_double(blk_mx);
			  }
			  if ( (fabs(hfptr->delta) > TINY) && s->Hassembly) {
				  /* MUTEX LOCK */
				  if (hfptr->blk_T == NULL) {
					  hfptr->blk_T = IzIz_sqrt2by3(s,n1,n2);
				  }
				  /* MUTEX UNLOCK */
				  blk_mx = create_blk_mat_double_copy(wsp->HQ_off[0]);
				  blk_dm_change_basis(blk_mx, hfptr->blk_T, s);
				  blk_dm_multod(wsp->HQ_off[0],blk_mx,hfptr->Rmol[3].re);
				  blk_dm_multod(wsp->HQ_off[1],blk_mx,hfptr->Rmol[4].re);
				  blk_dm_multod(wsp->HQ_off[2],blk_mx,hfptr->Rmol[4].im);
				  blk_dm_multod(wsp->HQ_off[3],blk_mx,hfptr->Rmol[5].re);
				  blk_dm_multod(wsp->HQ_off[4],blk_mx,hfptr->Rmol[5].im);
				  free_blk_mat_double(blk_mx);
			  }
		  }
		  wsp->HF_used[i] = 0;
		  return;
	  }
	  fprintf(stderr,"Error: turnoff - unknown name '%s'\n",name);
	  exit(1);
  }

  fprintf(stderr,"error: turnoff: unknown interaction name '%s'\n",name);
  exit(1);
}

void clear_int_off(Sim_info *sim, Sim_wsp *wsp)
{
	int j;

	if (wsp->Hiso_off) {
		free_blk_mat_double(wsp->Hiso_off);
		wsp->Hiso_off = NULL;
	}
	if (sim->Hassembly) {
	  for (j=0; j<5; j++) {
		  if (wsp->HQ_off[j]) {
			  free_blk_mat_double(wsp->HQ_off[j]);
			  wsp->HQ_off[j] = NULL;
		  }
	  }
	}
	wsp->Nint_off = 0;
}

void ham_turnon(Sim_info *s, Sim_wsp *wsp,char* name)
{
  int n1,n2,i,j;
  mat_double *mx;
  blk_mat_double *blk_mx;

  if (!strcmp(name,"all")) {
	  for (i=0; i<s->nCS; i++) wsp->CS_used[i] = 1;
	  for (i=0; i<s->nDD; i++) wsp->DD_used[i] = 1;
	  for (i=0; i<s->nQ; i++) wsp->Q_used[i] = 1;
	  for (i=0; i<s->nJ; i++) wsp->J_used[i] = 1;
	  for (i=0; i<s->nG; i++) wsp->G_used[i] = 1;
	  for (i=0; i<s->nHF; i++) wsp->HF_used[i] = 1;
	  for (i=0; i<s->nHEX; i++) wsp->HEx_used[i] = 1;
	  for (i=0; i<s->nEDD; i++) wsp->EDD_used[i] = 1;
	  clear_int_off(s,wsp);
     return;
  }
  if (sscanf(name,"shift_%d",&n1) == 1) {
	  for (i=0; i<s->nCS;i++) {
		  if (s->CS[i]->nuc == n1) {
			  if (!(wsp->CS_used[i])) {
				  Shift *csptr;
				  csptr = wsp->CS[i];
				  assert(csptr->T);
				  if (--(wsp->Nint_off) > 0) {
					//  odecist
					  mx = dm_change_basis_2(csptr->T,wsp->Hiso_off->basis,s);
					  if (fabs(csptr->iso)>TINY) {
						  assert(wsp->Hiso_off);
						  blk_dm_multod_diag(wsp->Hiso_off,mx,-csptr->iso);
					  }
					  if (s->Hassembly && fabs(csptr->delta)>TINY) {
						  for (j=0; j<5; j++) {
							  assert(wsp->HQ_off[j]);
						  }
						  blk_dm_multod_diag(wsp->HQ_off[0],mx,-csptr->Rmol[3].re);
						  blk_dm_multod_diag(wsp->HQ_off[1],mx,-csptr->Rmol[4].re);
						  blk_dm_multod_diag(wsp->HQ_off[2],mx,-csptr->Rmol[4].im);
						  blk_dm_multod_diag(wsp->HQ_off[3],mx,-csptr->Rmol[5].re);
						  blk_dm_multod_diag(wsp->HQ_off[4],mx,-csptr->Rmol[5].im);
					  }
					  free_double_matrix(mx);
				  } else {
					//  zrusit
					  clear_int_off(s,wsp);
				  }
			  }
			  wsp->CS_used[i] = 1;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnon - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"quadrupole_%d",&n1) == 1) {
	  for (i=0; i<s->nQ;i++) {
		  if (s->Q[i]->nuc == n1) {
			  if (!(wsp->Q_used[i])) {
				  Quadrupole *qptr = wsp->Q[i];
				  assert(qptr->T);
				  if (--(wsp->Nint_off) > 0) {
					  if (s->Hassembly) {
						  for (j=0; j<5; j++) {
							  assert(wsp->HQ_off[j]);
						  }
						  mx = dm_change_basis_2(qptr->T,wsp->HQ_off[0]->basis,s);
						  blk_dm_multod_diag(wsp->HQ_off[0],mx,-qptr->Rmol[3].re);
						  blk_dm_multod_diag(wsp->HQ_off[1],mx,-qptr->Rmol[4].re);
						  blk_dm_multod_diag(wsp->HQ_off[2],mx,-qptr->Rmol[4].im);
						  blk_dm_multod_diag(wsp->HQ_off[3],mx,-qptr->Rmol[5].re);
						  blk_dm_multod_diag(wsp->HQ_off[4],mx,-qptr->Rmol[5].im);
						  free_double_matrix(mx);
					  }
				  } else {
					  clear_int_off(s,wsp);
				  }
			  }

			  wsp->Q_used[i] = 1;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnon - unknown name '%s'\n",name);
	  exit(1);
  }

  if (sscanf(name,"dipole_%d_%d",&n1, &n2) == 2) {
	  if (n1 > n2) {
		  int dum;
		  dum = n1; n1 = n2; n2 = dum;
	  }
	  for (i=0; i<s->nDD;i++) {
		  if ( (s->DD[i]->nuc[0] == n1) && (s->DD[i]->nuc[1] == n2) ) {
			  if (!(wsp->DD_used[i])) {
				  Dipole *ddptr = wsp->DD[i];
				  assert(ddptr->blk_T);
				  if (--(wsp->Nint_off) > 0) {
					  if (s->Hassembly) {
						  for (j=0; j<5; j++) {
							  assert(wsp->HQ_off[j]);
						  }
						  blk_mx = create_blk_mat_double_copy(wsp->HQ_off[0]);
						  blk_dm_change_basis(blk_mx, ddptr->blk_T,s);
						  blk_dm_multod(wsp->HQ_off[0],blk_mx,-ddptr->Rmol[3].re);
						  blk_dm_multod(wsp->HQ_off[1],blk_mx,-ddptr->Rmol[4].re);
						  blk_dm_multod(wsp->HQ_off[2],blk_mx,-ddptr->Rmol[4].im);
						  blk_dm_multod(wsp->HQ_off[3],blk_mx,-ddptr->Rmol[5].re);
						  blk_dm_multod(wsp->HQ_off[4],blk_mx,-ddptr->Rmol[5].im);
						  free_blk_mat_double(blk_mx);
					  }
				  } else {
					  clear_int_off(s,wsp);
				  }
			  }

			  wsp->DD_used[i] = 1;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnon - unknown name '%s'\n",name);
	  exit(1);
  }

  if (sscanf(name,"jcoupling_%d_%d",&n1, &n2) == 2) {
	  if (n1 > n2) {
		  int dum;
		  dum = n1; n1 = n2; n2 = dum;
	  }
	  for (i=0; i<s->nJ;i++) {
		  if ( (s->J[i]->nuc[0] == n1) && (s->J[i]->nuc[1] == n2) ) {
			  Jcoupling *jptr;
			  if (!(wsp->J_used[i])) {
				  jptr = wsp->J[i];

				  if (--(wsp->Nint_off) > 0) {
					//  odecist
					  if (fabs(jptr->iso)>TINY) {
						  assert(jptr->blk_Tiso && wsp->Hiso_off);
						  blk_mx = create_blk_mat_double_copy(wsp->Hiso_off);
						  blk_dm_change_basis(blk_mx, jptr->blk_Tiso,s);
						  blk_dm_multod(wsp->Hiso_off,blk_mx,-jptr->iso);
						  free_blk_mat_double(blk_mx);
					  }
					  if (s->Hassembly && fabs(jptr->delta)>TINY) {
						  assert(jptr->blk_T);
						  for (j=0; j<5; j++) {
							  assert(wsp->HQ_off[j]);
						  }
						  blk_mx = create_blk_mat_double_copy(wsp->HQ_off[0]);
						  blk_dm_change_basis(blk_mx, jptr->blk_T,s);
						  blk_dm_multod(wsp->HQ_off[0],blk_mx,-jptr->Rmol[3].re);
						  blk_dm_multod(wsp->HQ_off[1],blk_mx,-jptr->Rmol[4].re);
						  blk_dm_multod(wsp->HQ_off[2],blk_mx,-jptr->Rmol[4].im);
						  blk_dm_multod(wsp->HQ_off[3],blk_mx,-jptr->Rmol[5].re);
						  blk_dm_multod(wsp->HQ_off[4],blk_mx,-jptr->Rmol[5].im);
						  free_blk_mat_double(blk_mx);
					  }
				  } else {
					//  zrusit
					  clear_int_off(s,wsp);
				  }
			  }

			  wsp->J_used[i] = 1;
			  return;
		  }
	  }
	  fprintf(stderr,"Error: turnon - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"gtensor_%d",&n1) == 1) {
	  i = gtensor_exist(s,n1);
	  if (i >= 0) {
		  if (!(wsp->CS_used[i])) {
			  Gtensor *gptr;
			  gptr = wsp->G[i];
			  assert(gptr->T);
			  if (--(wsp->Nint_off) > 0) {
				  //  odecist
				  mx = dm_change_basis_2(gptr->T,wsp->Hiso_off->basis,s);
				  if (fabs(gptr->iso)>TINY) {
					  assert(wsp->Hiso_off);
					  blk_dm_multod_diag(wsp->Hiso_off,mx,-gptr->iso);
				  }
				  if (s->Hassembly && fabs(gptr->delta)>TINY) {
					  for (j=0; j<5; j++) {
						  assert(wsp->HQ_off[j]);
					  }
					  blk_dm_multod_diag(wsp->HQ_off[0],mx,-gptr->Rmol[3].re);
					  blk_dm_multod_diag(wsp->HQ_off[1],mx,-gptr->Rmol[4].re);
					  blk_dm_multod_diag(wsp->HQ_off[2],mx,-gptr->Rmol[4].im);
					  blk_dm_multod_diag(wsp->HQ_off[3],mx,-gptr->Rmol[5].re);
					  blk_dm_multod_diag(wsp->HQ_off[4],mx,-gptr->Rmol[5].im);
				  }
				  free_double_matrix(mx);
			  } else {
				  //  zrusit
				  clear_int_off(s,wsp);
			  }
		  }
		  wsp->G_used[i] = 1;
		  return;
	  }
	  fprintf(stderr,"Error: turnon - unknown name '%s'\n",name);
	  exit(1);
  }
  if (sscanf(name,"hyperfine_%d_%d",&n1, &n2) == 2) {
	  i = hyperfine_exist(s,n1,n2);
	  if ( i >= 0 ) {
		  Hyperfine *hfptr;
		  if (!(wsp->HF_used[i])) {
			  hfptr = wsp->HF[i];
			  if (--(wsp->Nint_off) > 0) {
				  //  odecist
				  if (fabs(hfptr->iso)>TINY) {
					  assert(hfptr->blk_Tiso && wsp->Hiso_off);
					  blk_mx = create_blk_mat_double_copy(wsp->Hiso_off);
					  blk_dm_change_basis(blk_mx, hfptr->blk_Tiso,s);
					  blk_dm_multod(wsp->Hiso_off,blk_mx,-hfptr->iso);
					  free_blk_mat_double(blk_mx);
				  }
				  if (s->Hassembly && fabs(hfptr->delta)>TINY) {
					  assert(hfptr->blk_T);
					  for (j=0; j<5; j++) {
						  assert(wsp->HQ_off[j]);
					  }
					  blk_mx = create_blk_mat_double_copy(wsp->HQ_off[0]);
					  blk_dm_change_basis(blk_mx, hfptr->blk_T,s);
					  blk_dm_multod(wsp->HQ_off[0],blk_mx,-hfptr->Rmol[3].re);
					  blk_dm_multod(wsp->HQ_off[1],blk_mx,-hfptr->Rmol[4].re);
					  blk_dm_multod(wsp->HQ_off[2],blk_mx,-hfptr->Rmol[4].im);
					  blk_dm_multod(wsp->HQ_off[3],blk_mx,-hfptr->Rmol[5].re);
					  blk_dm_multod(wsp->HQ_off[4],blk_mx,-hfptr->Rmol[5].im);
					  free_blk_mat_double(blk_mx);
				  }
			  } else {
				  //  zrusit
				  clear_int_off(s,wsp);
			  }
		  }
		  wsp->HF_used[i] = 1;
		  return;
	  }
	  fprintf(stderr,"Error: turnon - unknown name '%s'\n",name);
	  exit(1);
  }


  fprintf(stderr,"error: turnon: unknown interaction name '%s'\n",name);
  exit(1);
}

/* rotation MOL (=CRY) --> ROT */
void ham_rotate(Sim_info *s, Sim_wsp *wsp)
{
  int i;
  mat_complx *d2;
  
  //DEBUGPRINT("ham_rotate point 1\n");

  /* This procedure is called outside the pulse sequence and
     must therefore be performed even if the interactions are turned off. */
  d2 = wigner2(wsp->cryst.alpha,wsp->cryst.beta,wsp->cryst.gamma);

  if (s->Hassembly) {
	  if (wsp->Dmol_rot) free_complx_matrix(wsp->Dmol_rot);
	  wsp->Dmol_rot = d2;
	  for (i=0; i<s->nQ; i++) wig2rot(wsp->Q_Rrot[i],wsp->Q[i]->Rmol,d2);
	  for (i=0; i<s->nMIX; i++) {
		  if (s->MIX[i]->type == 1) {
			  wig2rot(wsp->CS_Rrot[s->MIX[i]->idx],wsp->CS[s->MIX[i]->idx]->Rmol,d2);
		  } else {
			  wig2rot(wsp->DD_Rrot[s->MIX[i]->idx],wsp->DD[s->MIX[i]->idx]->Rmol,d2);
		  }
	  }
	  // NOTE 30.5.2021 - Hassembly is disabled in DNPframe where hyperfine is used, so these lines can be obsolete
	  // I do this because I need Rrot for terms T2+-1 that are not included in HQ[]
	  for (i=0; i<s->nHF; i++) {
		  if (wsp->HF_Rrot[i] != NULL) wig2rot(wsp->HF_Rrot[i],wsp->HF[i]->Rmol,d2);
	  }
  } else {
	  for (i=0; i<s->nCS; i++) {
		  //DEBUGPRINT("   i = %d, w_iso = %f\n",i,s->CS[i]->iso);
		  if (wsp->CS_Rrot[i] != NULL) wig2rot(wsp->CS_Rrot[i],wsp->CS[i]->Rmol,d2);
	  }
	  //DEBUGPRINT("ham_rotate point 3 nDD=%d\n",s->nDD);
	  for (i=0; i<s->nDD; i++) wig2rot(wsp->DD_Rrot[i],wsp->DD[i]->Rmol,d2);
	  //DEBUGPRINT("ham_rotate point 4 nQ=%d\n",s->nQ);
	  for (i=0; i<s->nQ; i++) wig2rot(wsp->Q_Rrot[i],wsp->Q[i]->Rmol,d2);
	  //DEBUGPRINT("ham_rotate point 5 nJ=%d\n",s->nJ);
	  for (i=0; i<s->nJ; i++) {
		  if (wsp->J_Rrot[i] != NULL) wig2rot(wsp->J_Rrot[i],wsp->J[i]->Rmol,d2);
	  }
	  //DEBUGPRINT("ham_rotate point 6\n");
	  for (i=0; i<s->nG; i++) {
		  if (wsp->G_Rrot[i] != NULL) wig2rot(wsp->G_Rrot[i],wsp->G[i]->Rmol,d2);
	  }
	  //DEBUGPRINT("ham_rotate point 7\n");
	  for (i=0; i<s->nHF; i++) {
		  if (wsp->HF_Rrot[i] != NULL) wig2rot(wsp->HF_Rrot[i],wsp->HF[i]->Rmol,d2);
	  }
	  //DEBUGPRINT("ham_rotate point 8\n");
	  for (i=0; i<s->nEDD; i++) wig2rot(wsp->EDD_Rrot[i],wsp->EDD[i]->Rmol,d2);
	  free_complx_matrix(d2);
  }
}

/* Spatial rotation ROT->LAB. In LABframe, include everything.
 * This is called also by ham_hamilton_dnpframe() to create nuclear laboratory frame Hamiltonians
   Populate with anisotropic parts of interactions.
 */
void ham_hamilton_labframe(Sim_info *s, Sim_wsp *wsp, mat_complx *d2)
{
	int i, q, sgn;
	mat_complx *ham;

	// wsp->ham_blk is populated with offsets and isotropic parts in ham_hamilton()
	if (wsp->Hcplx == NULL) {
		wsp->Hcplx = create_blk_mat_complx_copy2(wsp->ham_blk);
	}
	blk_cm_zero(wsp->Hcplx);
	assert(wsp->Hcplx->Nblocks == 1);
    ham =  wsp->Hcplx->m;

	for (i=0; i<s->nCS; i++) {
		if (!wsp->CS_used[i] || wsp->CS_Rrot[i] == NULL) continue;
		wig2rot(wsp->CS_Rlab[i],wsp->CS_Rrot[i],d2);
		sgn = -1;
		for (q=1; q<4; q++) { // T2-2 and T22 are zero for CSA
			complx z = wsp->CS_Rlab[i][5-q];
			z.re *= sgn; z.im *= sgn;
			cm_multocr(ham,s->CS[i]->T2q[q], z);
			sgn *= -1;
		}
	}
	for (i=0; i<s->nDD; i++) {
		if (!wsp->DD_used[i]) continue;
		wig2rot(wsp->DD_Rlab[i],wsp->DD_Rrot[i],d2);
		sgn = 1;
		for (q=0; q<5; q++) {
			complx z = wsp->DD_Rlab[i][5-q];
			z.re *= sgn; z.im *= sgn;
			cm_multocr(ham,s->DD[i]->T2q[q], z);
			sgn *= -1;
		}
	}
	for (i=0; i<s->nJ; i++) {
		if (!wsp->J_used[i] || wsp->J_Rrot[i] == NULL) continue;
		wig2rot(wsp->J_Rlab[i],wsp->J_Rrot[i],d2);
		sgn = 1;
		for (q=0; q<5; q++) {
			complx z = wsp->J_Rlab[i][5-q];
			z.re *= sgn; z.im *= sgn;
			cm_multocr(ham,s->J[i]->T2q[q], z);
			sgn *= -1;
		}
	}
	for (i=0; i<s->nQ; i++) {
		if (!wsp->Q_used[i]) continue;
		wig2rot(wsp->Q_Rlab[i], wsp->Q_Rrot[i], d2);
		sgn = 1;
		for (q=0; q<5; q++) {
			complx z = wsp->Q_Rlab[i][5-q];
			z.re *= sgn; z.im *= sgn;
			cm_multocr(ham,s->Q[i]->T2q[q], z);
			sgn *= -1;
		}
	}
	if (s->frame != DNPFRAME) { // this function is called by ham_hamilton_dnpframe to prepare nuclear interactions.
		// For LABFRAME we include full forms also for electrons
		for (i=0; i<s->nG; i++) {
			fprintf(stderr,"Error: Gtensor not yet implemented in LABframe");
			exit(1);
		}
		for (i=0; i<s->nHF; i++) {
			fprintf(stderr,"Error: Hyperfine not yet implemented in LABframe");
			exit(1);
		}
		for (i=0; i<s->nEDD; i++) { // electron dipole-dipole is clear
			if (!wsp->EDD_used[i]) continue;
			wig2rot(wsp->EDD_Rlab[i],wsp->EDD_Rrot[i],d2);
			sgn = 1;
			for (q=0; q<5; q++) {
				complx z = wsp->EDD_Rlab[i][5-q];
				z.re *= sgn; z.im *= sgn;
				cm_multocr(ham,s->EDD[i]->T2q[q], z);
				sgn *= -1;
			}
		}

	}
}

/* Spatial rotation ROT->LAB. In DNPframe, nuclear interactions are as in LABframe, electrons are in rotating frame
 * Having it separate, is should allow implementation of block-diagonalization with respect to electrons
   Populate with anisotropic parts of interactions.
 ***/
void ham_hamilton_dnpframe(Sim_info *sim, Sim_wsp *wsp, complx *d20, mat_complx *d2)
{
	int i;
	double dw;
	complx cdw, zz;
	// in current implementation, no block_diag (30.5.2021)

	// wsp->ham_blk is populated with offsets and isotropic parts in ham_hamilton()
	if (wsp->Hcplx == NULL) {
		wsp->Hcplx = create_blk_mat_complx_copy2(wsp->ham_blk);
	}
	blk_cm_zero(wsp->Hcplx);

	// call labframe for nuclei, result is in wsp->Hcplx first block = assumes block_diag is OFF
	ham_hamilton_labframe(sim, wsp, d2);

	// add gtensor  anisotropy (isotropic part was added to Hiso in readsys)
	for (i=0; i<sim->nG; i++) {
		if (!wsp->G_used[i] || wsp->G_Rrot[i] == NULL) continue;
		dw = wig20rot(wsp->G_Rrot[i],d20);
		//DEBUGPRINT("ham_hamilton: adding ANISO gtensor %d (dw=%f)\n",i,dw);
		if (fabs(dw) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->G[i]->T,dw);
	}
	// add hyperfine anisotropic part (isotropic part was added to Hiso in readsys)
	for (i=0; i<sim->nHF; i++) {
		if (!wsp->HF_used[i] || wsp->HF_Rrot[i] == NULL) continue;
		dw = wig20rot(wsp->HF_Rrot[i],d20);
		//DEBUGPRINT("ham_hamilton: adding ANISO hyperfine T20 %d (dw=%f)\n",i,dw);
		if (fabs(dw) > TINY) blk_dm_multod(wsp->ham_blk,wsp->HF[i]->blk_T,dw);
		// T2+-1 terms
		complx *zptr1 = wsp->HF_Rrot[i];
		complx *zptr2 = &(d2->data[5]); // here is located D-2-1
		//cv_print(zptr1,"Rrot");
		cblas_zdotu_sub(LEN(d20),zptr1+1,1,zptr2,1,&cdw);
		//printf("hyperfine -1: (%g, %g)\n",cdw.re, cdw.im);
		blk_cm_multocr(wsp->Hcplx,wsp->HF[i]->blk_Tb,cdw);
		zptr2 = &(d2->data[15]); // here is D-2+1
		cblas_zdotu_sub(LEN(d20),zptr1+1,1,zptr2,1,&cdw);
		//printf("hyperfine +1: (%g, %g)\n",cdw.re, cdw.im);
		blk_cm_multocr(wsp->Hcplx,wsp->HF[i]->blk_Ta,cdw);
	}
	// add electron-electron dipolar part
	for (i=0; i<sim->nEDD; i++) {
		if (!wsp->EDD_used[i]) continue;
		dw = wig20rot(wsp->EDD_Rrot[i],d20);
		//DEBUGPRINT("ham_hamilton: adding electron dipole %d\n",i);
		if (fabs(dw) > TINY) blk_dm_multod(wsp->ham_blk,wsp->EDD[i]->blk_T,dw);
	}

	// wsp->ham_blk will be added in ham_hamilton()
}

/* Spatial rotation ROT->LAB. In ROTframe, include only T00 and T20 terms of all interactions, and quadrupole higher orders */
void ham_hamilton_rotframe(Sim_info *s, Sim_wsp *wsp, complx *d20, mat_complx *d2)
{
	if (s->do_avg == 1) {
		ham_hamilton_integrate(s,wsp,wsp->dtmax);
		blk_dm_muld(wsp->ham_blk,1.0e6/wsp->dtmax);
		return;
	}

  int i;
  double dw;
  complx cdw, zz;

  // wsp->ham_blk is prepared in ham_hamilton()
  if (s->Hassembly) {
	  complx res[6];
	  wig2rot_t(res,d20,wsp->Dmol_rot);
	  blk_dm_multod(wsp->ham_blk,wsp->HQ[0],res[3].re);
	  blk_dm_multod(wsp->ham_blk,wsp->HQ[1],2.0*res[4].re);
	  blk_dm_multod(wsp->ham_blk,wsp->HQ[2],-2.0*res[4].im);
	  blk_dm_multod(wsp->ham_blk,wsp->HQ[3],2.0*res[5].re);
	  blk_dm_multod(wsp->ham_blk,wsp->HQ[4],-2.0*res[5].im);
	  if (wsp->Nint_off) {
		  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[0],-res[3].re);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[1],-2.0*res[4].re);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[2],2.0*res[4].im);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[3],-2.0*res[5].re);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[4],2.0*res[5].im);
	  }
	  // in ROTframe we omit the pseudosecular T2+-1 terms of Hyperfine interaction
  } else {
	  for (i=0; i<s->nCS; i++) {
		  if (!wsp->CS_used[i] || wsp->CS_Rrot[i] == NULL) continue;
		  dw = wig20rot(wsp->CS_Rrot[i],d20);
		  //DEBUGPRINT("ham_hamilton: adding ANISO shift %d (dw=%f)\n",i,dw);
		  if (fabs(dw) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->CS[i]->T,dw);
	  }
	  for (i=0; i<s->nDD; i++) {
		  if (!wsp->DD_used[i]) continue;
		  dw = wig20rot(wsp->DD_Rrot[i],d20);
		  //DEBUGPRINT("ham_hamilton: adding dipole %d\n",i);
		  if (fabs(dw) > TINY) blk_dm_multod(wsp->ham_blk,wsp->DD[i]->blk_T,dw);
	  }
	  for (i=0; i<s->nJ; i++) {
		  //if (!wsp->J_used[i] || s->J[i]->blk_T == NULL) continue;
		  if (!wsp->J_used[i] || wsp->J_Rrot[i] == NULL) continue;
		  dw = wig20rot(wsp->J_Rrot[i],d20);
		  //DEBUGPRINT("ham_hamilton: adding ANISO Jcoupling %d",i);
		  //DEBUGPRINT("(matdim=%d,",s->J[i]->T->row);
		  //DEBUGPRINT("Janiso=%f)\n",dw);
		  if (fabs(dw) > TINY) blk_dm_multod(wsp->ham_blk,wsp->J[i]->blk_T,dw);
	  }
	  for (i=0; i<s->nG; i++) {
		  if (!wsp->G_used[i] || wsp->G_Rrot[i] == NULL) continue;
		  dw = wig20rot(wsp->G_Rrot[i],d20);
		  //DEBUGPRINT("ham_hamilton: adding ANISO gtensor %d (dw=%f)\n",i,dw);
		  if (fabs(dw) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->G[i]->T,dw);
	  }
	  for (i=0; i<s->nHF; i++) { // adds the secular part but not the pseudosecular part, so no DNP transfers possible
		  if (!wsp->HF_used[i] || wsp->HF_Rrot[i] == NULL) continue;
		  dw = wig20rot(wsp->HF_Rrot[i],d20);
		  //DEBUGPRINT("ham_hamilton: adding ANISO hyperfine T20 %d (dw=%f)\n",i,dw);
		  if (fabs(dw) > TINY) blk_dm_multod(wsp->ham_blk,wsp->HF[i]->blk_T,dw);
	  }
	  for (i=0; i<s->nEDD; i++) {
		  if (!wsp->EDD_used[i]) continue;
		  dw = wig20rot(wsp->EDD_Rrot[i],d20);
		  //DEBUGPRINT("ham_hamilton: adding dipole %d\n",i);
		  if (fabs(dw) > TINY) blk_dm_multod(wsp->ham_blk,wsp->EDD[i]->blk_T,dw);
	  }
  }

  /* quadrupoles are special */
  for (i=0; i<s->nQ; i++) {
	  if (!wsp->Q_used[i]) continue;
	  wig2rot(wsp->Q_Rlab[i], wsp->Q_Rrot[i], d2);
	  if (!s->Hassembly) {
		  dw = wsp->Q_Rlab[i][3].re;
		  //DEBUGPRINT("ham_hamilton: adding quad %d order 1\n",i);
		  if (fabs(dw) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->Q[i]->T,dw);
	  }
	  if (s->Q[i]->order > 1) {
		  cdw = Cmul(wsp->Q_Rlab[i][1],wsp->Q_Rlab[i][5]);
		  cdw.re /= 2.0*wsp->Q[i]->w0; cdw.im /= 2.0*wsp->Q[i]->w0;
		  assert(fabs(cdw.im) < TINY);
		  if (fabs(cdw.re) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->QTa[i],cdw.re);
		  cdw = Cmul(wsp->Q_Rlab[i][2],wsp->Q_Rlab[i][4]);
		  cdw.re /= 2.0*wsp->Q[i]->w0; cdw.im /= 2.0*wsp->Q[i]->w0;
		  assert(fabs(cdw.im) < TINY);
		  if (fabs(cdw.re) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->QTb[i],cdw.re);
	  }
	  if (s->Q[i]->order == 3) {
		  dw = (wsp->Q[i]->w0)*(wsp->Q[i]->w0);
		  cdw = Cmul(wsp->Q_Rlab[i][1],wsp->Q_Rlab[i][5]);
		  cdw = Cmul(cdw,wsp->Q_Rlab[i][3]);
		  blk_dm_multod_diag(wsp->ham_blk,wsp->QT3b[i],cdw.re/(4*dw));
		  cdw = Cmul(wsp->Q_Rlab[i][4],wsp->Q_Rlab[i][4]);
		  cdw = Cmul(cdw,wsp->Q_Rlab[i][1]);
		  zz = Cmul(wsp->Q_Rlab[i][2],wsp->Q_Rlab[i][2]);
		  zz = Cmul(zz,wsp->Q_Rlab[i][5]);
		  blk_dm_multod_diag(wsp->ham_blk,wsp->QT3c[i],(cdw.re+zz.re)/dw);
		  cdw = Cmul(wsp->Q_Rlab[i][2],wsp->Q_Rlab[i][4]);
		  cdw = Cmul(cdw,wsp->Q_Rlab[i][3]);
		  blk_dm_multod_diag(wsp->ham_blk,wsp->QT3a[i],cdw.re/dw);
	  }
  }

  /* mixing terms are special as well */
  for (i=0; i<s->nMIX; i++) {
	  assert(d2 != NULL);
	  int j1 = s->MIX[i]->qidx;
	  if (!wsp->Q_used[j1]) continue;
	  int j2 = s->MIX[i]->idx;
	  if (s->MIX[i]->type == 1) {
		  /* Q-CSA */
		  wig2rot(wsp->CS_Rlab[j2], wsp->CS_Rrot[j2], d2);
		  cdw = Cmul(wsp->Q_Rlab[j1][2],wsp->CS_Rlab[j2][4]);
		  cdw = Cadd(cdw, Cmul(wsp->Q_Rlab[j1][4],wsp->CS_Rlab[j2][2]));
		  blk_dm_multod_diag(wsp->ham_blk, wsp->MT[i], cdw.re/2.0/wsp->Q[j1]->w0);
	  } else {
		  /* Q-DD */
		  wig2rot(wsp->DD_Rlab[j2], wsp->DD_Rrot[j2], d2);
		  cdw = Cmul(wsp->Q_Rlab[j1][2],wsp->DD_Rlab[j2][4]);
		  cdw = Cadd(cdw, Cmul(wsp->Q_Rlab[j1][4],wsp->DD_Rlab[j2][2]));
		  //printf("vystup_hetero: (%g, %g)\n",cdw.re, cdw.im);
		  //assert( fabs(cdw.im) < TINY );
		  blk_dm_multod_diag(wsp->ham_blk, wsp->MT[i], -cdw.re/2.0/wsp->Q[j1]->w0);
		  if (wsp->MTa[i] != NULL) {
			  cdw = Cmul(wsp->Q_Rlab[j1][2],wsp->DD_Rlab[j2][4]);
			  cdw = Cadd(cdw, Cmul(wsp->Q_Rlab[j1][1],wsp->DD_Rlab[j2][5]));
			  printf("vystup_homo1: (%g, %g)\n",cdw.re, cdw.im);
			  assert( fabs(cdw.im) < TINY);
			  blk_dm_multod(wsp->ham_blk, wsp->MTa[i], cdw.re/4.0/wsp->Q[j1]->w0);
			  cdw = Cmul(wsp->Q_Rlab[j1][4],wsp->DD_Rlab[j2][2]);
			  cdw = Cadd(cdw, Cmul(wsp->Q_Rlab[j1][5],wsp->DD_Rlab[j2][1]));
			  printf("vystup_homo2: (%g, %g)\n",cdw.re, cdw.im);
			  assert( fabs(cdw.im) < TINY);
			  blk_dm_multod(wsp->ham_blk, wsp->MTb[i], cdw.re/4.0/wsp->Q[j1]->w0);
		  }
	  }
  }
  // wsp->ham_blk is DONE.
}

/* this makes final rotation ROT->LAB and creates interaction Hamiltonian */
void ham_hamilton(Sim_info *sim, Sim_wsp *wsp)
{
	int i;
	complx *d20=NULL;
	mat_complx *d2=NULL;

	// transformation matrices
	if (sim->dor == 1) {
		mat_complx *md1, *md2;
		md1 = wigner2(sim->wr1*wsp->t*1.0e-6*RAD2DEG+wsp->gamma_add,sim->brl1,0.0);
		md2 = wigner2(sim->wr2*wsp->t*1.0e-6*RAD2DEG,sim->brl2,0.0);
		d2 = cm_mul(md1,md2);
		d20 = complx_vector(5);
		for (i=1;i<=5;i++) d20[i] = cm_getelem(d2,i,3);
		free_complx_matrix(md1);
		free_complx_matrix(md2);
	} else {
		d20 = wigner20(sim->wr*wsp->t*1.0e-6*RAD2DEG+wsp->gamma_add,wsp->brl);
		d2 = wigner2(sim->wr*wsp->t*1.0e-6*RAD2DEG+wsp->gamma_add,wsp->brl,0.0);
	}

	// adding interactions
	assert(wsp->ham_blk != NULL);
	blk_dm_zero(wsp->ham_blk);

	/* offset terms */
	if (wsp->inhom_offset != NULL) {
		fprintf(stderr,"Error: checkpoint in ham.c, function ham_hamiltion() - inhomogeneity offsets are disabled\n");
		exit(1);
	}
	if (wsp->offset != NULL) {
		for (i=1; i<=sim->ss->nchan; i++) {
			//printf("ham_hamilton: adding offsets %g to channel %d\n",ovals[i],i);
			switch (sim->frame) {
			case ROTFRAME:
				blk_dm_multod_diag(wsp->ham_blk,wsp->chan_Iz[i],wsp->offset[i]);
				break;
			case DNPFRAME:
				if (sim->ss->iso[sim->ss->chan[i][1]]->number == 0) { // it is electron channel
					//printf("ham_hamilton is adding electron offset %f\n",wsp->offset[i]);
					blk_dm_multod_diag(wsp->ham_blk,wsp->chan_Iz[i],wsp->offset[i]);
					//dm_print(wsp->chan_Iz[i],"chan Iz");
				}
				break;
			case LABFRAME:
				// do nothing
				break;
			default :
				fprintf(stderr,"Error in ham_hamiltion() - invalid sim->frame %d\n",sim->frame);
				exit(1);
			}
		}
	}

	// add isotropic part of Hamiltonian
	blk_dm_multod(wsp->ham_blk,wsp->Hiso,1.0);
	if (wsp->Nint_off) {
		assert(wsp->Hiso_off);
		blk_dm_multod(wsp->ham_blk,wsp->Hiso_off,-1.0);
	}

	// add anisotropic parts, depending on calculation frame
	switch (sim->frame) {
	case ROTFRAME :
		ham_hamilton_rotframe(sim, wsp, d20, d2);
		// Hamiltonian is in wsp->ham_blk and it is block-diagonal and real
		break;
	case DNPFRAME:
		ham_hamilton_dnpframe(sim, wsp, d20, d2);
		// Hamiltonian is in wsp->Hcplx and it is block-diagonal and complex; wsp->ham_blk was also modified
		blk_cm_multocr(wsp->Hcplx,wsp->ham_blk,Cunit);
		break;
	case LABFRAME:
		ham_hamilton_labframe(sim, wsp, d2);
		// Hamiltonian is in wsp->Hcplx first block and it is complex
		// now add isotropic part to Hlab
	    assert(wsp->ham_blk->Nblocks == 1);
	    blk_cm_multocr(wsp->Hcplx,wsp->ham_blk,Cunit);
		break;
	default:
		fprintf(stderr,"Error: ham_hamilton - simulation frame (%d) not recognized\n",sim->frame);
		exit(1);
	}

	// clean up
	if (d20 != NULL) free_complx_vector(d20);
	if (d2 != NULL) free_complx_matrix(d2);
}


/*  int_t1^t2 exp(-i m wr t) = -i/(m wr) (exp(-i m wr t2) - exp(-i m wr t1)) */
void integ(complx* R,double t1,double t2,double wr)
{
   double wrm;

   wrm= wr * (-2.0); R[1] = Cmul(Complx(0,1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm)));
   wrm= wr * (-1.0); R[2] = Cmul(Complx(0,1.0/wrm),Csub(Cexpi(-t2*wrm),Cexpi(-t1*wrm)));
                     R[3] = Complx(t2-t1,0);
  R[4] = Conj(R[2]);
  R[5] = Conj(R[1]);
}

typedef struct {
	double times[2], beta[2], wr[2];
	mat_complx *d2mas, *d2dor1, *d2dor2;
	int type, Istate;
	complx *Imas, *Idor;
	complx *q2a, *q2b, *q3a, *q3b, *q3c, *qcsadd, *qddhomoa, *qddhomob;
} wigints;

void wigints_init(wigints *p, double t1, double t2, Sim_info *sim, Sim_wsp *wsp)
{
	if (p == NULL) {
		fprintf(stderr,"Error: wigints_init - got NULL pointer\n");
		exit(1);
	}
	p->times[0] = t1;
	p->times[1] = t2;
	if (sim->dor == 1) {  // we have DOR
		p->type = 1;
		p->beta[0] = sim->brl1;
		p->beta[1] = sim->brl2;
		p->wr[0] = sim->wr1;
		p->wr[1] = sim->wr2;
	} else {  // we have MAS
		p->type = 0;
		p->beta[0] = wsp->brl;
		p->wr[0] = sim->wr;
	}
	p->d2mas = NULL;
	p->d2dor1 = NULL;
	p->d2dor2 = NULL;
	p->Istate = 0;
	p->Imas = NULL;
	p->Idor = NULL;
	p->q2a = NULL;
	p->q2b = NULL;
	p->q3a = NULL;
	p->q3b = NULL;
	p->q3c = NULL;
	p->qcsadd = NULL;
	p->qddhomoa = NULL;
	p->qddhomob = NULL;
}

void wigints_destroy(wigints *p)
{
	assert(p != NULL);
	if (p->d2mas != NULL) { free_complx_matrix(p->d2mas); p->d2mas = NULL;}
	if (p->d2dor1 != NULL) {free_complx_matrix(p->d2dor1); p->d2dor1 = NULL;}
	if (p->d2dor2 != NULL) {free_complx_matrix(p->d2dor2); p->d2dor2 = NULL;}
	p->Istate = 0;
	if (p->Imas != NULL) {free(p->Imas); p->Imas = NULL;}
	if (p->Idor != NULL) {free(p->Idor); p->Idor = NULL;}
	if (p->q2a != NULL) {free(p->q2a); p->q2a = NULL;}
	if (p->q2b != NULL) {free(p->q2b); p->q2b = NULL;}
	if (p->q3a != NULL) {free(p->q3a); p->q3a = NULL;}
	if (p->q3b != NULL) {free(p->q3b); p->q3b = NULL;}
	if (p->q3c != NULL) {free(p->q3c); p->q3c = NULL;}
	if (p->qcsadd != NULL) {free(p->qcsadd); p->qcsadd = NULL;}
	if (p->qddhomoa != NULL) {free(p->qddhomoa); p->qddhomoa = NULL;}
	if (p->qddhomob != NULL) {free(p->qddhomob); p->qddhomob = NULL;}
}

void wigints_prepare(wigints *p, int code)
{
	int i;

	assert(p != NULL);
	if (p->type == 0) { // we have MAS
		if (p->d2mas == NULL) {
			p->d2mas = wigner2(0,p->beta[0],0);
		}
		if (p->Istate == 0) {
			p->Imas = (complx*)malloc(13*sizeof(complx));
			if (p->Imas == NULL) {
				fprintf(stderr,"Error: wigints_prepare - Imas alloc failure\n");
				exit(1);
			}
			for (i=2;i<6;i++) {
				double dw = p->wr[0] * (double)(i-6);
				p->Imas[i] = Cmul(Complx(0,1.0/dw),Csub(Cexpi(-dw*p->times[1]),Cexpi(-dw*p->times[0])));
			}
			p->Imas[6] = Complx(p->times[1]-p->times[0],0.0);
			p->Imas[7] = Conj(p->Imas[5]);
			p->Imas[8] = Conj(p->Imas[4]);
			p->Imas[9] = Conj(p->Imas[3]);
			p->Imas[10] = Conj(p->Imas[2]);
			p->Istate = 2; // good for up to second order
		}
	} else {            // we have DOR
		if (p->d2dor1 == NULL) {
			p->d2dor1 = wigner2(0,p->beta[0],0);
			p->d2dor2 = wigner2(0,p->beta[1],0);
		}
		if (p->Istate == 0) {
			p->Idor = (complx*)malloc(13*13*sizeof(complx));
			if (p->Idor == NULL) {
				fprintf(stderr,"Error: wigints_prepare - Idor alloc failure\n");
				exit(1);
			}
			int j;
			for (i=2; i<11; i++) {
				double dw1 = p->wr[0]*(double)(i-6);
				for (j=2; j<11; j++) {
					double dw2 = p->wr[1]*(double)(j-6);
					p->Idor[i*13+j] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));
				}
			}
			p->Istate = 2; // good for up to second order
		}
	}
	switch (code) {
	case 2:  // quadrupole second order
		if (p->q2a != NULL) break; // work is already done
		p->q2a = (complx *)malloc(5*5*sizeof(complx));  // d(m,-1)*d(n,+1)*I(m+n)
		p->q2b = (complx *)malloc(5*5*sizeof(complx));  // d(m,-2)*d(n,+2)*I(m+n)
		if (p->q2a == NULL || p->q2b == NULL) {
			fprintf(stderr,"Error: wigints_prepare can not allocate q2a, q2b\n");
			exit(1);
		}
		if (p->type == 2) { //we have MAS
			int j;
			for (i=0; i<5; i++) {
				for (j=0; j<5; j++) {
					p->q2a[i+j*5] = Cmul(p->d2mas->data[i+1*5],Cmul(p->d2mas->data[j+3*5],p->Imas[i+j+2]));
					p->q2b[i+j*5] = Cmul(p->d2mas->data[i+0*5],Cmul(p->d2mas->data[j+4*5],p->Imas[i+j+2]));
				}
			}
		} else {    // we have DOR
			int j, k, l;
			complx z1, z2, z3;
			for (i=0; i<5; i++) {
				for (j=0; j<5; j++) {
					for (k=0; k<5; k++) {
						for (l=0; l<5; l++) {
							z1 = Cmul(p->d2dor1->data[i+k*5],p->d2dor2->data[k+1*5]);
							z2 = Cmul(p->d2dor1->data[j+l*5],p->d2dor2->data[l+3*5]);
							z3.re =z1.re*z2.re - z1.im*z2.im; z3.im = z1.re*z2.im+z1.im*z2.re;
							p->q2a[i+j*5] = Cmul(z3,p->Idor[(i+j)*13+(k+l)]);
							z1 = Cmul(p->d2dor1->data[i+k*5],p->d2dor2->data[k+0*5]);
							z2 = Cmul(p->d2dor1->data[j+l*5],p->d2dor2->data[l+4*5]);
							z3.re =z1.re*z2.re - z1.im*z2.im; z3.im = z1.re*z2.im+z1.im*z2.re;
							p->q2b[i+j*5] = Cmul(z3,p->Idor[(i+j)*13+(k+l)]);
						}
					}
				}
			}
		}
		break;
	case 3: // quadrupole third order
		if (p->Istate != 3) {
			assert(p->Istate == 2);
			if (p->type == 0) {  // we have MAS
				double dw = -6*p->wr[0];
				p->Imas[0] = Cmul(Complx(0,1.0/dw),Csub(Cexpi(-dw*p->times[1]),Cexpi(-dw*p->times[0])));
				dw = -5*p->wr[0];
				p->Imas[1] = Cmul(Complx(0,1.0/dw),Csub(Cexpi(-dw*p->times[1]),Cexpi(-dw*p->times[0])));
				p->Imas[11] = Conj(p->Imas[1]);
				p->Imas[12] = Conj(p->Imas[0]);
			} else {             // we have DOR
				double dw1 = -6*p->wr[0];
				double dw2 = -6*p->wr[1];
				p->Idor[0*13+0] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));
				dw2 = -dw2; // = +6
				p->Idor[0*13+12] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));
				dw2 = -5*p->wr[1];
				p->Idor[0*13+1] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));
				dw2 = -dw2; // = +5
				p->Idor[0*13+11] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));

				dw1 = -5*p->wr[0];  // dw2 remains +5
				p->Idor[1*13+11] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));
				dw2 = -dw2; // = -5
				p->Idor[1*13+1] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));
				dw2 = -6*p->wr[1];
				p->Idor[1*13+0] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));
				dw2 = -dw2; // = +6
				p->Idor[1*13+12] = Cmul(Complx(0,1.0/(dw1+dw2)),Csub(Cexpi(-(dw1+dw2)*p->times[1]),Cexpi(-(dw1+dw2)*p->times[0])));

				p->Idor[12*13+0] = Conj(p->Idor[0*13+12]);
				p->Idor[12*13+1] = Conj(p->Idor[0*13+11]);
				p->Idor[12*13+12] = Conj(p->Idor[0*13+0]);
				p->Idor[12*13+11] = Conj(p->Idor[0*13+1]);
				p->Idor[11*13+0] = Conj(p->Idor[1*13+12]);
				p->Idor[11*13+1] = Conj(p->Idor[1*13+11]);
				p->Idor[11*13+12] = Conj(p->Idor[1*13+0]);
				p->Idor[11*13+11] = Conj(p->Idor[1*13+1]);
			}
			p->Istate = 3;
		}
		break;
	case 4:  // mixing Q-CSA or Q-DDhetero
		break;
	case 5:  // mixing Q-DDhomo
		break;
	default:
		fprintf(stderr,"Error: wigints_prepare - invalid code\n");
		exit(1);
	}
}

/* this makes final rotation ROT->LAB, integrates and creates interaction Hamiltonian */
void ham_hamilton_integrate(Sim_info *s, Sim_wsp *wsp, double dur)
{
	int i, j;
	double dw, t;
	complx *d20, R[6], RR[6], II[13], cdw;
	mat_complx *d2=NULL;

	d20 = wigner20(0.0,wsp->brl);
	t = wsp->t * 1.0e-6;
	dur *= 1.0e-6;
	integ(R, t, t+dur, s->wr);
	if (s->nQ) {
		d2 = wigner2(0,wsp->brl,0);
		for (i=0;i<6;i++) {
			j = i-6;
			dw = s->wr * (double)j;
			II[i] = Cmul(Complx(0,1.0/dw),Csub(Cexpi(-(t+dur)*dw),Cexpi(-t*dw)));
		}
		II[6] = Complx(dur,0.0);
		II[7] = Conj(II[5]);
		II[8] = Conj(II[4]);
		II[9] = Conj(II[3]);
		II[10] = Conj(II[2]);
		II[11] = Conj(II[1]);
		II[12] = Conj(II[0]);
	}

	blk_dm_zero(wsp->ham_blk);
	blk_dm_multod(wsp->ham_blk,wsp->Hiso,dur);
	if (wsp->Nint_off) {
		assert(wsp->Hiso_off);
		blk_dm_multod(wsp->ham_blk,wsp->Hiso_off,-dur);
	}

	if (s->Hassembly) {
		  for (i=1; i<6; i++) {
			  RR[i].re = d20[i].re*R[i].re - d20[i].im*R[i].im;
			  RR[i].im = d20[i].im*R[i].re + d20[i].re*R[i].im;
		  }
		  wig2rot_t(R,RR,wsp->Dmol_rot);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ[0],R[3].re);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ[1],2.0*R[4].re);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ[2],-2.0*R[4].im);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ[3],2.0*R[5].re);
		  blk_dm_multod(wsp->ham_blk,wsp->HQ[4],-2.0*R[5].im);
		  if (wsp->Nint_off) {
			  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[0],-R[3].re);
			  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[1],-2.0*R[4].re);
			  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[2],2.0*R[4].im);
			  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[3],-2.0*R[5].re);
			  blk_dm_multod(wsp->ham_blk,wsp->HQ_off[4],2.0*R[5].im);
		  }
		  // but we miss Hyperfine terms T2+-1
	} else {
		for (i=0; i<s->nCS; i++) {
			if (!wsp->CS_used[i] || wsp->CS_Rrot[i] == NULL) continue;
			cdw = Cnull;
			for (j=1; j<6; j++) {
				cdw = Cadd(cdw,Cmul(Cmul(wsp->CS_Rrot[i][j],d20[j]),R[j]));
			}
			assert(fabs(cdw.im) < TINY);
			if (fabs(cdw.re) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->CS[i]->T,cdw.re);
		}
		for (i=0; i<s->nDD; i++) {
			if (!wsp->DD_used[i]) continue;
			cdw = Cnull;
			for (j=1; j<6; j++) {
				cdw = Cadd(cdw,Cmul(Cmul(wsp->DD_Rrot[i][j],d20[j]),R[j]));
			}
			assert(fabs(cdw.im) < TINY);
			if (fabs(cdw.re) > TINY) blk_dm_multod(wsp->ham_blk,wsp->DD[i]->blk_T,cdw.re);
		}
		for (i=0; i<s->nJ; i++) {
			if (!wsp->J_used[i] || wsp->J_Rrot[i]==NULL) continue;
			cdw = Cnull;
			for (j=1; j<6; j++) {
				cdw = Cadd(cdw,Cmul(Cmul(wsp->J_Rrot[i][j],d20[j]),R[j]));
			}
			assert(fabs(cdw.im) < TINY);
			if (fabs(cdw.re) > TINY) blk_dm_multod(wsp->ham_blk,wsp->J[i]->blk_T,cdw.re);
		}
		for (i=0; i<s->nG; i++) {
			if (!wsp->G_used[i] || wsp->G_Rrot[i] == NULL) continue;
			cdw = Cnull;
			for (j=1; j<6; j++) {
				cdw = Cadd(cdw,Cmul(Cmul(wsp->G_Rrot[i][j],d20[j]),R[j]));
			}
			assert(fabs(cdw.im) < TINY);
			if (fabs(cdw.re) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->G[i]->T,cdw.re);
		}
	}

	/* quadrupoles are special */
	for (i=0; i<s->nQ; i++) {
		if (!wsp->Q_used[i]) continue;
		if (!s->Hassembly) {
			cdw = Cnull;
			for (j=1; j<6; j++) {
				cdw = Cadd(cdw,Cmul(Cmul(wsp->Q_Rrot[i][j],d20[j]),R[j]));
			}
			assert(fabs(cdw.im) < TINY);
			//printf("vypis: (%g,%g)\n",cdw.re,cdw.im);
			if (fabs(cdw.re) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->Q[i]->T,cdw.re);
		}
		if (s->Q[i]->order > 1) { // second order
			complx cc1, cc2, c1, c2, c4, c5;
			int l, m;
			cc1 = cc2 = Cnull;
			for (l=0;l<5;l++) {
				c2 = Cmul(wsp->Q_Rrot[i][l+1],d2->data[l+1*5]);
				c1 = Cmul(wsp->Q_Rrot[i][l+1],d2->data[l+0*5]);
				for (m=0;m<5;m++) {
					c4 = Cmul(wsp->Q_Rrot[i][m+1],d2->data[m+3*5]);
					c5 = Cmul(wsp->Q_Rrot[i][m+1],d2->data[m+4*5]);
					cc2 = Cadd(cc2,Cmul(Cmul(c2,c4),II[l+m+2]));
					cc1 = Cadd(cc1,Cmul(Cmul(c1,c5),II[l+m+2]));
				}
			}
			//printf("vypis: (%g,%g) a (%g,%g)\n",cc1.re,cc1.im,cc2.re,cc2.im);
			//assert(fabs(cc1.im) < TINY);
			//assert(fabs(cc2.im) < TINY);
			if (fabs(cc1.re) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->QTa[i],cc1.re/2.0/wsp->Q[i]->w0);
			if (fabs(cc2.re) > TINY) blk_dm_multod_diag(wsp->ham_blk,wsp->QTb[i],cc2.re/2.0/wsp->Q[i]->w0);
		}
		if (s->Q[i]->order == 3) {
			complx cc1, cc2, cc3, z1a, z2a, z3a, z1b1, z2b1, z3b1, z1b2, z2b2, z3b2, z1c, z2c, z3c;
			int l, m, n;
			cc1 = cc2 = cc3 = Cnull;
			for (l=0; l<5; l++) {
				z1a = Cmul(wsp->Q_Rrot[i][l+1],d2->data[l+4*5]);
				z1b1 = Cmul(wsp->Q_Rrot[i][l+1],d2->data[l+3*5]);
				z1b2 = Cmul(wsp->Q_Rrot[i][l+1],d2->data[l+1*5]);
				z1c = Cmul(wsp->Q_Rrot[i][l+1],d2->data[l+1*5]);
				for (m=0; m<5; m++) {
					z2a = Cmul(z1a,Cmul(wsp->Q_Rrot[i][m+1],d2->data[m+2*5]));
					z2b1 = Cmul(z1b1,Cmul(wsp->Q_Rrot[i][m+1],d2->data[m+3*5]));
					z2b2 = Cmul(z1b2,Cmul(wsp->Q_Rrot[i][m+1],d2->data[m+1*5]));
					z2c = Cmul(z1c,Cmul(wsp->Q_Rrot[i][m+1],d2->data[m+2*5]));
					for (n=0; n<5; n++) {
						z3a = Cmul(z2a,Cmul(wsp->Q_Rrot[i][n+1],d2->data[n+0*5]));
						cc1 = Cadd(cc1,Cmul(z3a,II[l+m+n]));
						z3b1 = Cmul(z2b1,Cmul(wsp->Q_Rrot[i][n+1],d2->data[n+0*5]));
						z3b2 = Cmul(z2b2,Cmul(wsp->Q_Rrot[i][n+1],d2->data[n+4*5]));
						cc2 = Cadd(cc2,Cmul(Cadd(z3b1,z3b2),II[l+m+n]));
						z3c = Cmul(z2c,Cmul(wsp->Q_Rrot[i][n+1],d2->data[n+3*5]));
						cc3 = Cadd(cc3,Cmul(z3c,II[l+m+n]));
					}
				}
			}
			double w2 = (wsp->Q[i]->w0)*(wsp->Q[i]->w0);
			blk_dm_multod_diag(wsp->ham_blk,wsp->QT3b[i],cc1.re/(4.0*w2));
			blk_dm_multod_diag(wsp->ham_blk,wsp->QT3c[i],cc2.re/w2);
			blk_dm_multod_diag(wsp->ham_blk,wsp->QT3a[i],cc3.re/w2);
		}
	}
	/* mixing terms */
	for (i=0; i<s->nMIX; i++) {
		int j1 = s->MIX[i]->qidx;
		if (!wsp->Q_used[j1]) continue;
		int j2 = s->MIX[i]->idx;
		complx cc1, c1, c2, c3, c4;
		int l, m;
		if (s->MIX[i]->type == 1) {
			/* Q-CSA */
			cc1 = Cnull;
			for (l=0;l<5;l++) {
				c2 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+3*5]);
				c1 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+1*5]);
				for (m=0;m<5;m++) {
					c3 = Cmul(wsp->CS_Rrot[j2][m+1],d2->data[m+1*5]);
					c4 = Cmul(wsp->CS_Rrot[j2][m+1],d2->data[m+3*5]);
					cc1 = Cadd(cc1,Cmul(Cmul(c2,c3),II[l+m+2]));
					cc1 = Cadd(cc1,Cmul(Cmul(c1,c4),II[l+m+2]));
				}
			}
			blk_dm_multod_diag(wsp->ham_blk,wsp->MT[i],cc1.re/2.0/wsp->Q[j1]->w0);
		} else {
			/* Q-DD */
			cc1 = Cnull;
			for (l=0;l<5;l++) {
				c2 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+3*5]);
				c1 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+1*5]);
				for (m=0;m<5;m++) {
					c3 = Cmul(wsp->DD_Rrot[j2][m+1],d2->data[m+1*5]);
					c4 = Cmul(wsp->DD_Rrot[j2][m+1],d2->data[m+3*5]);
					cc1 = Cadd(cc1,Cmul(Cmul(c2,c3),II[l+m+2]));
					cc1 = Cadd(cc1,Cmul(Cmul(c1,c4),II[l+m+2]));
				}
			}
			//printf("uistup: (%g, %g)\n",cc1.re, cc1.im);
			//assert(fabs(cc1.im) < TINY);
			// following is according to Brinkmann
			blk_dm_multod_diag(wsp->ham_blk,wsp->MT[i],-cc1.re/2.0/wsp->Q[j1]->w0);
			// following is according Larsen
			//blk_dm_multod_diag(wsp->ham_blk,wsp->MT[i],cc1.re/2.0/wsp->Q[j1]->w0);
			if (wsp->MTa[i] != NULL) {
				complx cc2, z1, z2, z3, z4;
				cc1 = cc2 = Cnull;
				for (l=0;l<5;l++) {
					c1 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+1*5]);
					c2 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+0*5]);
					z1 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+3*5]);
					z2 = Cmul(wsp->Q_Rrot[j1][l+1],d2->data[l+4*5]);
					for (m=0;m<5;m++) {
						c3 = Cmul(wsp->DD_Rrot[j2][m+1],d2->data[m+3*5]);
						c4 = Cmul(wsp->DD_Rrot[j2][m+1],d2->data[m+4*5]);
						z3 = Cmul(wsp->DD_Rrot[j2][m+1],d2->data[m+1*5]);
						z4 = Cmul(wsp->DD_Rrot[j2][m+1],d2->data[m+0*5]);
						cc1 = Cadd(cc1,Cmul(Cmul(c1,c3),II[l+m+2]));
						cc1 = Cadd(cc1,Cmul(Cmul(c2,c4),II[l+m+2]));
						cc2 = Cadd(cc2,Cmul(Cmul(z1,z3),II[l+m+2]));
						cc2 = Cadd(cc2,Cmul(Cmul(z2,z4),II[l+m+2]));
					}
				}
				printf("uistup_homo1: (%g, %g)\n",cc1.re, cc1.im);
				printf("uistup_homo2: (%g, %g)\n",cc2.re, cc2.im);
				assert(fabs(cc1.im) < TINY);
				assert(fabs(cc2.im) < TINY);
				blk_dm_multod(wsp->ham_blk,wsp->MTa[i],cc1.re/4.0/wsp->Q[j1]->w0);
				blk_dm_multod(wsp->ham_blk,wsp->MTb[i],cc2.re/4.0/wsp->Q[j1]->w0);
			}
		}
	}
	// Hyperfine is special
	for (i=0; i<s->nHF; i++) {
		if (!wsp->HF_used[i] || wsp->HF_Rrot[i] == NULL) continue;
		if (wsp->Hcplx == NULL) {
			wsp->Hcplx = create_blk_mat_complx_copy2(wsp->ham_blk);
			blk_cm_zero(wsp->Hcplx);
		}
		if (!s->Hassembly) { // T20 term
			cdw = Cnull;
			for (j=1; j<6; j++) {
				cdw = Cadd(cdw,Cmul(Cmul(wsp->HF_Rrot[i][j],d20[j]),R[j]));
			}
			assert(fabs(cdw.im) < TINY);
			if (fabs(cdw.re) > TINY) blk_dm_multod(wsp->ham_blk,wsp->HF[i]->blk_T,cdw.re);
		}
		// now T2+-1 terms
		if (d2 == NULL) d2 = wigner2(0,wsp->brl,0);
		cdw = Cnull; // D2-1
		for (j=1; j<6; j++) {
			cdw = Cadd(cdw,Cmul(Cmul(wsp->HF_Rrot[i][j],d2->data[j+4]),R[j]));
		}
		blk_cm_multocr(wsp->Hcplx,wsp->HF[i]->blk_Tb,cdw);
		cdw = Cnull; // D2+1
		for (j=1; j<6; j++) {
			cdw = Cadd(cdw,Cmul(Cmul(wsp->HF_Rrot[i][j],d2->data[j+14]),R[j]));
		}
		blk_cm_multocr(wsp->Hcplx,wsp->HF[i]->blk_Ta,cdw);
	}

	/* offset terms */
	if (wsp->offset != NULL || wsp->inhom_offset != NULL) {
		int nch = s->ss->nchan;
		double *ovals = double_vector(nch);
		dv_zero(ovals);
		if (wsp->offset != NULL) dv_multod(ovals,wsp->offset,dur);
		if (wsp->inhom_offset != NULL) dv_multod(ovals,wsp->inhom_offset,dur);
		for (i=1; i<=nch; i++) {
			blk_dm_multod_diag(wsp->ham_blk,wsp->chan_Iz[i],ovals[i]);
		}
		free_double_vector(ovals);
	}


	free_complx_vector(d20);
	if (d2 != NULL) free_complx_matrix(d2);

	if (wsp->Hcplx != NULL) blk_cm_multocr(wsp->Hcplx,wsp->ham_blk,Cunit);

}

/* Caution: assumes _setrfprop has been already called */
blk_mat_complx * ham_rf(Sim_wsp *wsp)
{
	int i;
	blk_mat_complx *HRF = create_blk_mat_complx_copy2(wsp->sumHrf);

	for (i=0; i<HRF->Nblocks; i++) {
		if (HRF->blk_dims[i] == 1) {
			HRF->m[i].data[0] = Complx(wsp->sumHrf->m[i].data[0],0.0);
		} else {
			dm_copy2cm(wsp->sumHrf->m + i, HRF->m + i);
		}
	}
	blk_simtrans_zrot2(HRF, wsp->sumUph);

	return HRF;
}




