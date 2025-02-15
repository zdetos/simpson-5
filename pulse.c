/*
    Pulse propagation
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
     Zdenek Tosner & Niels Chr. Nielsen: changes for optimal control
     
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
    
    Performs propagation with the pulses given in the 'pulseq' procedure
    given in the input file. Makes available several commands that
    can be read about in the manual.
    
    Called from sim.c that setup and perform the simulation.    
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "matrix.h"
#include "cm.h"
#include "blockdiag.h"
#include "pulse.h"
#include "tclutil.h"
#include "defs.h"
#include "rfshapes.h"
#include "iodata.h"
#include "OCroutines.h"
#include "B0inhom.h"
#include "relax.h"
#include "sim.h"
#include "ham.h"
#include "spinsys.h"
#include "spinach.h"
#include "fidcalc.h"
	/* for acurate timings on windows */
//#define TIMING
#include "timing.h"

extern glob_info_type glob_info;

/*
  tpropstart_usec : the time when the propagator was reset last time
  tpulsestart_usec : the time when the pulse sequence was started
  t_usec : the current time.
  tstartincr_usec : the time that will be added to the current time after resetting.
  t_proplength : the length of a propagator
 
  In this file t*_usec is in microseconds, the other time variables
  in seconds. t_usec alias 't' can be changed from the tcl script and must therefore be in microseconds

  ZT: all times are in microseconds...
*/

void direct_propagate(Sim_info *sim, Sim_wsp *wsp)
{
	int i;

	DEBUGPRINT("direct_propagate start\n");
  /* enable all interactions, turn off offset and undo other possible remainders from previous pulseq calls) */
  ham_turnon(sim,wsp,"all");
  if (wsp->offset) {
	  free_double_vector(wsp->offset);
	  wsp->offset = NULL;
  }
  wsp->brl = sim->brl;
  wsp->gamma_add = 0.0;
  // optim. tip: fstart and fdetect are not usually changed in pulseq
  // use sim pointers and do not copy
  // commands changing them will create new matrices if necessary
  for (i=0; i<sim->Nfstart; i++) {
	  if (wsp->fstart[i] != sim->fstart[i]) {
		  free_complx_matrix(wsp->fstart[i]);
		  wsp->fstart[i] = sim->fstart[i];
	  }
	  if (wsp->sigma[i] == NULL)
		  wsp->sigma[i] = cm_dup(wsp->fstart[i]);
	  else
		  cm_copy(wsp->sigma[i],wsp->fstart[i]);
  }
  for (i=0; i<sim->Nfdetect; i++) {
	  if (wsp->fdetect[i] != sim->fdetect[i]) {
		  free_complx_matrix(wsp->fdetect[i]);
		  wsp->fdetect[i] = sim->fdetect[i];
	  }
  }

  wsp->Uisunit=1;
  wsp->isselectivepulse=0;
  wsp->cannotbestored=0;
  wsp->waspulse=0;
  wsp->tpropstart = wsp->t = wsp->tstart;
  wsp->curr_nsig = 0;

  //printf("direct_propagate wsp->tstart = %g\n",wsp->tstart);
  DEBUGPRINT("direct_propagate -> calling pulseq '%s'\n",sim->pulsename);
  //printf("direct_propagate - interp pointer %p\n",wsp->interp);
  //Tcl_Eval(wsp->interp,"puts \"NAZDAREK\"");
  if (Tcl_Eval(wsp->interp,sim->pulsename) != TCL_OK) {
    fprintf(stderr,"error: in evaluation of pulse sequence: %s\n",Tcl_GetStringResult(wsp->interp));
    exit(1);
  }
  
  /* ZT: when calculating gradients, reset mx_pos when crystallite calculation finishes */
  if (OCpar.gradmode) wsp->OC_mxpos=0;

}


/* Internal routines to be called from the TCL routines */

mat_complx * get_matrix_1(Sim_info *sim, Sim_wsp *wsp,char* name)
{
  int num, dupl=1;
  mat_complx *m, *res;
  Tcl_Interp* interp = wsp->interp;
  
  if (Tcl_GetInt(interp,name,&num) == TCL_OK) {    
    if (num < 1 || num > MAXSTO) {
       TclError(interp,"matrix: called with number larger than %d or smaller than 1\n",MAXSTO);
       return NULL;
    }
    if (!wsp->matrix[num]) {
       TclError(interp,"matrix number %d is undefined\n",num);
       return NULL;
    }
    m = wsp->matrix[num];
  } else if (!strcmp(name,"start")) {
	  m = wsp->fstart[0];
  } else if (!strcmp(name,"detect")) {
	  m = wsp->fdetect[0];
  } else if (!strcmp(name,"density")) {
	  m = wsp->sigma[0];
  } else if (!strcmp(name,"propagator")) {
      //m = wsp->U;
	  //blk_cm_print(wsp->U,"get_matrix PROP");
	  if (wsp->Uisunit == 1) {
		  m = cm_creatediag('d',sim->matdim,Cunit,wsp->U->basis);
	  } else {
		  m = complx_matrix(sim->matdim,sim->matdim,MAT_DENSE,1,wsp->U->basis);
		  blk_cm_copy_2(m,wsp->U);
	  }
	  dupl = 0;
  } else if (!strcmp(name,"hamiltonian")) {
    ham_hamilton(sim,wsp);
    if (sim->frame == DNPFRAME) {
    	m = complx_matrix(sim->matdim,sim->matdim,MAT_DENSE,1,wsp->Hcplx->basis);
    	blk_cm_copy_2(m,wsp->Hcplx);
    	dupl = 0;
    } else {
    	mat_double *fham = double_matrix(sim->matdim,sim->matdim,MAT_DENSE,0,wsp->ham_blk->basis);
    	blk_dm_copy_2(fham,wsp->ham_blk);
    	m = dm_complx(fham);
    	free_double_matrix(fham);
    	dupl = 0;
    }
  } else if (!strcmp(name,"avgham")) {
    double dt = wsp->t - wsp->tpropstart;
    blk_mat_complx *ah = blk_cm_ln(wsp->U);
    m = complx_matrix(sim->matdim, sim->matdim, MAT_DENSE, 0, wsp->U->basis);
    blk_cm_copy_2(m,ah);
    cm_mulc(m, Complx(0.0, 1.0e6/dt));
    dupl = 0;
    free_blk_mat_complx(ah);
  } else {
    TclError(interp,"error: matrix: argument must be 'start', 'detect', 'density', 'propagator', 'avgham', or "
                    "'hamiltonian'\n");
    return NULL;
  }
  
  if (dupl) {
    res = cm_dup(m);
  } else {
    res =  m;
  }
  cm_dense(res);
  //cm_print(res,"get_matrix_1 output");
  return res;
}

mat_complx* get_matrix_2(Sim_info *sim, Sim_wsp *wsp,char* name,char* list)
{
  mat_complx * m;
  int N;
  Tcl_Interp* interp = wsp->interp;
  
  if (!strcmp(name,"operator")) {
    return ss_readoper(sim, list);
  }
  N = sim->matdim;
  m = complx_matrix(N,N,MAT_DENSE,0,0);
  cm_zero(m);

  if (!strcmp(name,"notelements") || !strcmp(name,"elements")) {
      const char **list1,**list2;
      int i,nlist1,nlist2,val1,val2;
      
      if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK)
        return NULL;

      for (i=0;i<nlist1;i++) {
        if (Tcl_SplitList(interp,list1[i],&nlist2,&list2) != TCL_OK) {
          Tcl_Free((char *) list1);
          return NULL;
        }
        if (nlist2 != 2) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2);
          TclError(interp,"matrix: expecting two elements (matrix row and column) in list");
          return NULL;
        }
        if (Tcl_GetInt(interp,list2[0],&val1) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2); 
          return NULL;
        }
        if (Tcl_GetInt(interp,list2[1],&val2) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2); 
          return NULL;
        }
        if ( val1 < 1 || val1 > N || val2 < 1 || val2 > N) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2); 
          TclError(interp,"matrix: value out of range in list");
          return NULL;
        }
        m->data[val1-1+N*(val2-1)].re=1.0;
        Tcl_Free((char *) list2);        
      }
      if (!strcmp(name,"notelements")) {
	for (i=0; i<N*N; i++) {
	   if (m->data[i].re != 0.0)
	      m->data[i].re = 0.0;
	   else
	      m->data[i].re = 1.0;
	}
      }
      Tcl_Free((char *) list1);
  } else if (!strcmp(name,"totalcoherence")) {
    int coh[MAXSPINS+1];
    const char **list1;
    int nlist1,i,j,k;    
    mat_complx * Q;
    
    if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK)
      return NULL;
    if (nlist1 > MAXSPINS) {      
      TclError(interp,"matrix: internal array overflow\n");
      return NULL;
    }

    for (i=0;i<nlist1;i++) {
      if (Tcl_GetInt(interp,list1[i],&coh[i]) != TCL_OK) {
        Tcl_Free((char *) list1);
        return NULL;
      }
    }
    Q=Iqdelta(sim); cm_dense(Q);

    for (i=0;i<N;i++) {
      for (j=0;j<N;j++) {
        int qdelta = (int)( Q->data[i+j*N].re );
        for (k=0;k<nlist1;k++) {
          if (coh[k] == qdelta) {
             m->data[i+j*N].re = 1.0;
             break;
          }
        }
      }
    }
    Tcl_Free((char *) list1);
    free_complx_matrix(Q);
  } else if (!strcmp(name,"coherence")) {
    const char **list1,**list2;
    int i,j,nlist1,nlist2;      
    double coh[MAXSPINS+1];
    mat_complx *tmp;


    if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK)
      return NULL;

    for (i=0;i<nlist1;i++) {
      if (Tcl_SplitList(interp,list1[i],&nlist2,&list2) != TCL_OK) {
        Tcl_Free((char *) list1);
        return NULL;
      }
      if (nlist2 != sim->ss->nspins) {
        Tcl_Free((char *) list1);
        Tcl_Free((char *) list2);
        TclError(interp,"matrix: expecting %d elements (number of spins)",sim->ss->nspins);
        return NULL;
      }
      for (j=0;j<nlist2;j++) {
        if (Tcl_GetDouble(interp,list2[j],&coh[j+1]) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2);
          return NULL;
        }
      }
      tmp=Icoherence(sim,coh);
      /* here cm_or(m,m, tmp) was used. It did only that nonzero real elements of tmp were
       * replaced by Complx(1,0) and stored in m. New code does the same...
       */
      //cm_print(tmp,"get_matrix_2 Icoherence");
      cm_nnz2one(tmp);
      cm_addto(m,tmp);
      //cm_print(m,"get_matrix_2 -> m<-");
      free_complx_matrix(tmp);
      Tcl_Free((char *) list2);
    }
    Tcl_Free((char *) list1);
    
  } else if (!strcmp(name,"list")) {
    const char **list1,**list2,**list3;
    int i,j,nlist1,nlist2,nlist3;

    if (Tcl_SplitList(interp,list,&nlist1,&list1) != TCL_OK) {
      TclError(interp,"matrix: specify list");
      return NULL;
    }
    if (nlist1 != N) {
      TclError(interp,"matrix: expecting %d elements (matrix dimension)",
	  N);
      return NULL;
    }

    for (i=0;i<nlist1;i++) {
      if (Tcl_SplitList(interp,list1[i],&nlist2,&list2) != TCL_OK) {
        Tcl_Free((char *) list1);
        TclError(interp,"matrix: specify list (in list)");
        return NULL;
      }
      if (nlist2 != N) {
        Tcl_Free((char *) list1);
        Tcl_Free((char *) list2);
        TclError(interp,"matrix: expecting %d elements (matrix dimension)",
	  N);
        return NULL;
      }
      for (j=0;j<nlist2;j++) {
        if (Tcl_SplitList(interp,list2[j],&nlist3,&list3) != TCL_OK) {
	  Tcl_Free((char *) list1);
	  Tcl_Free((char *) list2);
          TclError(interp,"matrix: specify {re im} for each element");
	  return NULL;
	}
	if (nlist3 != 1 && nlist3 != 2) {
	  Tcl_Free((char *) list1);
	  Tcl_Free((char *) list2);
	  Tcl_Free((char *) list3);
          TclError(interp,"matrix: specify {re im} for each element");
	  return NULL;
	}
        if (Tcl_GetDouble(interp,list3[0],&(m->data[i+j*N].re)) != TCL_OK) {
          Tcl_Free((char *) list1);
          Tcl_Free((char *) list2);
	  Tcl_Free((char *) list3);
          TclError(interp,"matrix: real element (%d,%d) does not appear to be a double",
	    i,j);
          return NULL;
        }
        if (nlist3 == 2) {
	  if (Tcl_GetDouble(interp,list3[1],&(m->data[i+j*N].im)) != TCL_OK) {
            Tcl_Free((char *) list1);
            Tcl_Free((char *) list2);
  	    Tcl_Free((char *) list3);
            TclError(interp,"matrix: imaginary element (%d,%d) does not appear to be a double",
	      i,j);
            return NULL;
          }
	}
	Tcl_Free((char *) list3);
      }
      Tcl_Free((char *) list2);
    }
    Tcl_Free((char *) list1);
  } else {
    TclError(interp, "Usage: matrix: arguments must be"
                            "  elements {{i j} {i j} ..}\n"
                            "  notelements {{i j} {i j} ..}\n"                            
                            "  totalcoherence {dm1 dm2 ..}\n"
                            "  coherence {dm1 dm2 dmN} {dm1 dm2 dmN} ..., where N is number of nuclei\n"
                            "  list {{c11 c12 ..} {c21 c22 ..} ..}, where c is a real number or {re im}\n"
			    );
    return NULL;
  }
  return m;
}


/****
 * ZT: modification in scaling rf power by rf inhomogeneity factor
 ****/
void _rf(Sim_wsp *wsp, int channel,double rffield)
{
  if (channel < 1 || channel > LEN(wsp->rfv)) {
    fprintf(stderr,"error: channel %d out of range in _rf()\n",channel);
    exit(1);
  }
  wsp->rfv[channel] = rffield*M_PI*2.0*(wsp->rffac[channel]);
}

void _ph(Sim_wsp *wsp, int channel, double phase)
{
  if (channel < 1 || channel > LEN(wsp->phv)) {
    fprintf(stderr,"error: channel %d out of range in _ph()\n",channel);
    exit(1);
  }
  wsp->phv[channel] = phase*DEG2RAD;
}

void _reset_prop(Sim_info *sim, Sim_wsp *wsp)
{
  wsp->tpropstart = wsp->t;
  wsp->Uisunit = 1;
  wsp->waspulse = 0;
  wsp->isselectivepulse = 0;
  wsp->cannotbestored = 0;
}

void _evolve_with_prop(Sim_info *sim, Sim_wsp *wsp)
{
  if (!wsp->Uisunit) {
	  //printf("\n=== _evolve_with_prop ===\n");
	  //cm_print(wsp->sigma,"initial sigma");
	  //blk_cm_print(wsp->U,"propagator");
	  int i;
	  for (i=0; i<sim->Nfstart; i++) {
		  // basis compatibility check
		  if (wsp->sigma[i]->basis != wsp->U->basis) {
			  mat_complx *dum = cm_change_basis_2(wsp->sigma[i],wsp->U->basis,sim);
			  free_complx_matrix(wsp->sigma[i]);
			  wsp->sigma[i] = dum;
		  }
		  blk_simtrans(wsp->sigma[i],wsp->U,sim);
	  }
	  //cm_print(wsp->sigma,"resulting sigma");
	  //printf("    ======\n\n");
  }
}

void _reset(Sim_info *sim, Sim_wsp *wsp, double tstartincr_usec)
{
  int i;
  for (i=0; i<sim->Nfstart; i++)
	  cm_copy(wsp->sigma[i],wsp->fstart[i]);
  wsp->t = wsp->tstart + tstartincr_usec;
  _reset_prop(sim,wsp);
  if (wsp->offset) {
	  free_double_vector(wsp->offset);
	  wsp->offset = NULL;
  }
  ham_turnon(sim,wsp,"all");
}


void _filter(Sim_info *sim, Sim_wsp *wsp, int num, int proj)
{
  mat_complx *mask;
  int i;

  mask = wsp->matrix[num];
  if (!mask) {
    fprintf(stderr,"error: illegal filter number %d\n",num);
    exit(1);
  }
  if ( mask->row != sim->matdim || mask->col != sim->matdim) {
    fprintf(stderr,"error: filter: matrix %d must be a %dx%d matrix\n",num,sim->matdim,sim->matdim);
    exit(1);
  }

  //DEBUGPRINT("_filter 1: ham dim %d, type '%s'\n",wsp->ham->row, matrix_type(wsp->ham->type));
  _evolve_with_prop(sim,wsp);
  //DEBUGPRINT("_filter 2: ham dim %d, type '%s'\n",wsp->ham->row, matrix_type(wsp->ham->type));
  _reset_prop(sim,wsp);
  //DEBUGPRINT("_filter 3: ham dim %d, type '%s'\n",wsp->ham->row, matrix_type(wsp->ham->type));
  wsp->cannotbestored=1;

  for (i=0; i<sim->Nfstart; i++) {
	  // basis compatibility check
	  if (wsp->sigma[i]->basis != mask->basis) {
		  mask = cm_change_basis_2(wsp->matrix[num], wsp->sigma[i]->basis,sim);
		  free_complx_matrix(wsp->matrix[num]);
		  wsp->matrix[num] = mask;
	  }
	  if (proj == 1) {
		  cm_filter_proj(wsp->sigma[i],mask);
	  } else {
		  cm_filter(wsp->sigma[i], mask);
	  }
  }
  //DEBUGPRINT("_filter 4: ham dim %d, type '%s'\n",wsp->ham->row, matrix_type(wsp->ham->type));
}

void _store(Sim_info *sim, Sim_wsp *wsp, int num, int adj)
{
  if (wsp->cannotbestored) {
    fprintf(stderr,"error: the command 'filter' cannot be stored\n");
    exit(-1);
  }

  if (wsp->STO[num] != NULL) {
	  free_blk_mat_complx(wsp->STO[num]);
  }
  if (wsp->Uisunit) {
	  wsp->STO[num] = create_blk_mat_complx(sim->matdim,1,NULL,MAT_DENSE_DIAG,sim->basis);
	  blk_cm_unit(wsp->STO[num]);
  } else {
	  if (adj) {
		  wsp->STO[num] = blk_cm_adjoint(wsp->U);
	  } else {
		  wsp->STO[num] = blk_cm_dup(wsp->U);
	  }
  }
  wsp->STO_tproplength_usec[num] = wsp->t - wsp->tpropstart;
  wsp->STO_tpropstart_usec[num] = wsp->tpropstart;
  /* save whether a pulse was applied */
  wsp->STO_waspulse[num] = wsp->waspulse;
  wsp->STO_brl[num] = wsp->brl;
}

/******
 *  check if the time the propagator is applied is the
 *  time it was calculated. Time is relative to the
 *  rotor period because the Hamiltonian is periodic
 *  with the rotor period.
 *  ZT: added condition about non-negativity of propagator start time
 *     negative start time is used with propagators created from matrices,lists...
 *     for which I don't want to check start time
 *****/
void check_prop_reuse(Sim_info *sim, Sim_wsp *wsp, int num)
{
	double tr,tdiff,t,tpropstart;

	//if ( (sim->wr != 0.0) && (wsp->STO_tpropstart_usec[num] >= 0.0) && !(VARIOUS_NOCHECKPROPTIME & various)) {
	if ( (sim->taur > 0) && (wsp->STO_tpropstart_usec[num] >= 0.0) && !(VARIOUS_NOCHECKPROPTIME & various)) {

	     //tr = M_PI*2.0e6/sim->wr;
		tr = sim->taur;
	     tpropstart = wsp->STO_tpropstart_usec[num];
	     t = wsp->t;
	     tdiff = t-tpropstart;
	     //DEBUGPRINT("check_prop_reuse: tpropstart=%f, t=%f, tdiff=%.10g\n",tpropstart,t,tdiff);
	     if (fabs(tdiff-floor(tdiff/tr+0.5)*tr) > 1e-5) {
	       fprintf(stderr,"error: a propagator was calculated at time %.10gusec\nrelative to the"
	       " start of the rotor period but reused at time %.10gusec.\n",
	        tpropstart-floor(tpropstart/tr)*tr,
	        t-floor(t/tr)*tr );
	        exit(1);
	     }
	     /* check for dynamic rotor angle */
	     if (fabs(wsp->STO_brl[num] - wsp->brl) > TINY) {
	       fprintf(stderr,"error: a propagator was calculated using a rotor angle of %g\nand reused with a rotor angle of %g\n",
	         wsp->STO_brl[num],wsp->brl);
	       exit(1);
	     }
	  }
}

void _prop(Sim_info *sim, Sim_wsp *wsp, int num,int ntimes)
{
  int i;

  if (!wsp->STO[num]) {
    fprintf(stderr,"error: illegal prop number %d\n",num);
    exit(1);
  }
  if (ham_ischanged(sim,wsp)) {
    fprintf(stderr,"error: 'store', 'turnon' and 'turnoff' will not have any"
                   " effect\n on a stored propagator\n"
                   "Call 'reset' before using stored propagators\n");
    exit(1);
  }

  for (i=0; i<ntimes; i++) {
	  check_prop_reuse(sim, wsp, num);
	  update_propagator(wsp->U,wsp->STO[num],sim,wsp);
	  wsp->t += wsp->STO_tproplength_usec[num];
	  if (wsp->evalmode == EM_ACQBLOCK) {
		  if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY ) {
			  _store(sim,wsp,wsp->acqblock_sto,0);
			  wsp->acqblock_t0 = wsp->t;
			  acqblock_sto_incr(wsp);
		  } else if ( wsp->t-wsp->acqblock_t0 > wsp->dw) {
			  fprintf(stderr,"Error: acq_block - 'prop' event not synchronized with dwelltime\n");
			  exit(1);
		  }
	  }
  }
  /* does the current propagator contain any pulses ? */  
  //wsp->waspulse = wsp->STO_waspulse[num] || wsp->waspulse;
  wsp->waspulse += wsp->STO_waspulse[num];
}


void _acq(Sim_info *sim, Sim_wsp *wsp, double phase, int do_reset)
{
  complx z, phfac, *ptr;
  int i, is, id, NN;
  
  NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);
  if (wsp->curr_nsig + 1 > sim->ntot) {
    fprintf(stderr,"error: acq overflow in fid points\n");
    exit(1);
  }
  (wsp->curr_nsig)++;
  if (fabs(phase)>TINY) {
	  phfac.re = cos(phase*DEG2RAD);
	  phfac.im = sin(phase*DEG2RAD);
  } else {
	  phfac = Cnull;
  }
  /* in relaxation mode rho is already evolved */
  if (!(sim->relax)) {
     _evolve_with_prop(sim,wsp);
  }
  if (do_reset) _reset_prop(sim,wsp);

  /*** would this work??? ***/
  for (i=0; i<NN; i++) {
	  is = i % sim->Nfstart;
	  id = i % sim->Nfdetect;
	  // basis compatibility check
	  if (wsp->fdetect[id]->basis != wsp->sigma[is]->basis) {
		  mat_complx *dum = cm_change_basis_2(wsp->fdetect[i],wsp->sigma[is]->basis,sim);
		  if (wsp->fdetect[id] != sim->fdetect[id]) free_complx_matrix(wsp->fdetect[id]);
		  wsp->fdetect[id] = dum;
	  }
	  if (sim->acq_adjoint == 0) {
		  z = cm_trace(wsp->fdetect[id],wsp->sigma[is]);
	  } else {
		  z = cm_trace_adjoint(wsp->fdetect[id],wsp->sigma[is]);
	  }
	  ptr = &(wsp->fid[wsp->curr_nsig+i*sim->ntot]);
	  if (fabs(phase) > TINY) {
	      ptr->re += phfac.re*z.re+phfac.im*z.im;
	      ptr->im += -phfac.im*z.re+phfac.re*z.im;
	  } else {
		  ptr->re += z.re;
		  ptr->im += z.im;
	  }
  }

/****
  if (sim->Nfstart == sim->Nfdetect ) {
	  for (i=0; i<sim->Nfstart; i++) {
		  if (sim->acq_adjoint == 0) {
			  z = cm_trace(wsp->fdetect[i],wsp->sigma[i]);
		  } else {
			  z = cm_trace_adjoint(wsp->fdetect[i],wsp->sigma[i]);
		  }
		  ptr = &(wsp->fid[wsp->curr_nsig+i*sim->ntot]);
		  if (fabs(phase) > TINY) {
		      ptr->re += phfac.re*z.re+phfac.im*z.im;
		      ptr->im += -phfac.im*z.re+phfac.re*z.im;
		  } else {
			  ptr->re += z.re;
			  ptr->im += z.im;
		  }
	  }
  } else if (sim->Nfstart == 1) {
	  for (i=0; i<sim->Nfdetect; i++) {
		  if (sim->acq_adjoint == 0) {
			  z = cm_trace(wsp->fdetect[i],wsp->sigma[0]);
		  } else {
			  z = cm_trace_adjoint(wsp->fdetect[i],wsp->sigma[0]);
		  }
		  ptr = &(wsp->fid[wsp->curr_nsig+i*sim->ntot]);
		  if (fabs(phase) > TINY) {
		      ptr->re += phfac.re*z.re+phfac.im*z.im;
		      ptr->im += -phfac.im*z.re+phfac.re*z.im;
		  } else {
			  ptr->re += z.re;
			  ptr->im += z.im;
		  }
	  }
  } else {
	  assert(sim->Nfdetect == 1);
	  for (i=0; i<sim->Nfstart; i++) {
		  if (sim->acq_adjoint == 0) {
			  z = cm_trace(wsp->fdetect[0],wsp->sigma[i]);
		  } else {
			  z = cm_trace_adjoint(wsp->fdetect[0],wsp->sigma[i]);
		  }
		  ptr = &(wsp->fid[wsp->curr_nsig+i*sim->ntot]);
		  if (fabs(phase) > TINY) {
		      ptr->re += phfac.re*z.re+phfac.im*z.im;
		      ptr->im += -phfac.im*z.re+phfac.re*z.im;
		  } else {
			  ptr->re += z.re;
			  ptr->im += z.im;
		  }
	  }
  }
***/
  //printf("_acq: fid(%d) = (%f, %f) ; z = (%f, %f)\n",wsp->curr_nsig,ptr->re,ptr->im,z.re,z.im);
}

void _nacq(Sim_info *sim, Sim_wsp *wsp, int np,int num,double phase)
{
  int k;
//  double cosph=0,sinph=0;
//  complx *ptr, z;

  if (!wsp->STO[num]) {
    fprintf(stderr,"error: illegal prop number %d\n",num);
    exit(1);
  }
//  if (wsp->curr_nsig + np > LEN(wsp->fid)) {
  if (wsp->curr_nsig + np > sim->ntot) {
	  fprintf(stderr,"error: acq overflow in fid points\n");
	  exit(1);
  }

  if (ham_ischanged(sim,wsp)) {
    fprintf(stderr,"error: 'store', 'turnon' and 'turnoff' will not have any"
                   " effect\n on a stored propagator (indirectly used"
                   " by 'acq <prop> <ntimes>' in this case)\n"
                   "Call 'reset' before using stored propagators\n");
    exit(1);
  }
  _acq(sim,wsp,phase,1);

  wsp->waspulse = wsp->STO_waspulse[num];
  blk_cm_copy(wsp->U,wsp->STO[num]);
  wsp->Uisunit = 0;
 // if (wsp->sigma->basis != wsp->U->basis) {
//	  DEBUGPRINT("_nacq : changing basis of sigma %d --> %d\n",wsp->sigma->basis,wsp->U->basis);
//	  mat_complx *dum = cm_change_basis_2(wsp->sigma,wsp->U->basis,sim);
//	  free_complx_matrix(wsp->sigma);
//	  wsp->sigma = dum;
// }
//  if (wsp->fdetect->basis != wsp->U->basis) {
//	  DEBUGPRINT("_nacq : changing basis of fdetect %d --> %d\n",wsp->fdetect->basis,wsp->U->basis);
//	  mat_complx *dum = cm_change_basis_2(wsp->fdetect,wsp->U->basis,sim);
//	  if (wsp->fdetect != sim->fdetect) free_complx_matrix(wsp->sigma);
//	  wsp->fdetect = dum;
//  }
//  if (phase != 0.0) {
//    phase *= DEG2RAD;
//    cosph=cos(phase);
//    sinph=sin(phase);
//  }
//  ptr = &(wsp->fid[(wsp->curr_nsig)]);
  for (k=2;k<=np;k++) {
	  check_prop_reuse(sim,wsp,num);
	  //DEBUGPRINT("_nacq: wsp->t = %f, proplength=%f\n",wsp->t,wsp->STO_tproplength_usec[num]);
      wsp->t += wsp->STO_tproplength_usec[num];
      //DEBUGPRINT("_nacq: k=%d, wsp->t=%f\n",k,wsp->t);
      //_evolve_with_prop(sim,wsp);
      _acq(sim,wsp,phase,0);
      //      ptr++; (wsp->curr_nsig)++;
//      if (sim->acq_adjoint == 0) {
//        z = cm_trace(wsp->fdetect,wsp->sigma);
//      } else {
//        z = cm_trace_adjoint(wsp->fdetect,wsp->sigma);
//      }
  
//      if (phase != 0.0) {
//          ptr->re += cosph*z.re+sinph*z.im;
//          ptr->im += -sinph*z.re+cosph*z.im;
//      } else {
//    	  ptr->re += z.re;
//    	  ptr->im += z.im;
//      }
      //DEBUGPRINT("_nacq: fid(%d) = (%f, %f)\n",wsp->curr_nsig,ptr->re,ptr->im);
      //if (k==10) exit(1);
  }
  _reset_prop(sim,wsp);
}

void _delay_simple(Sim_info *sim, Sim_wsp *wsp, double duration)
{
	int i,n;
	double dt, dt_us;

	DEBUGPRINT("_delay_simple with duration %f\n",duration);

	assert(wsp->dU != NULL);
	blk_cm_unit(wsp->dU);
	if (sim->wr < TINY) {
		/* static case */
		ham_hamilton(sim,wsp);
		if ( (sim->frame == LABFRAME) || (sim->frame == DNPFRAME) ) {
			blk_prop_complx_3(wsp->dU,wsp->Hcplx,duration*1e-6,sim);
		} else {
			blk_prop_real(wsp->dU,wsp->ham_blk,duration*1e-6,sim);
		}
		wsp->t += duration;
	} else {
		/* MAS */
		if (sim->Hint_isdiag) {
			DEBUGPRINT("_delay_simple integration of diagonal Hamiltonian\n");
			/* delay in pulse sequence */
			assert(sim->frame == ROTFRAME);
			ham_hamilton_integrate(sim,wsp,duration);
			blk_prop_real(wsp->dU,wsp->ham_blk,1.0, sim); // unit duration since time was integrated
			wsp->t += duration;
		} else {
			//assert(wsp->tmpU != NULL);
			n=(int)ceil(duration/wsp->dtmax);
			if (n < 1) n=1;
			dt_us = duration/(double)n;
			dt = dt_us*1.0e-6;
			DEBUGPRINT("_delay_simple split into %d steps of %f us length\n",n,dt_us);
			for (i=1;i<=n;i++) {
				ham_hamilton(sim,wsp);
				//printf("STEP %i, time %f, ",i,wsp->t);
				//blk_dm_print(wsp->ham_blk," Ham");
				//if (i==1) {
				if ( (sim->frame == LABFRAME) || (sim->frame == DNPFRAME) ) {
					blk_prop_complx_3(wsp->dU,wsp->Hcplx,dt,sim);
				} else {
					blk_prop_real(wsp->dU, wsp->ham_blk, dt, sim);
				}
				//} else {
				//	blk_prop_real(wsp->tmpU, wsp->ham_blk, dt, sim->propmethod);
				//	update_propagator(wsp->dU, wsp->tmpU, sim, NULL);
				//}
				wsp->t += dt_us;
			}
		}
	}
	//blk_cm_print(wsp->dU,"_delay_simple step propagator");
}


void _delay(Sim_info *sim, Sim_wsp *wsp, double duration)
{
	double dt, t_tmp;
	int i, n;

	if (duration > wsp->dt_gcompute) {
		duration = wsp->dt_gcompute;
	}

	if (wsp->evalmode == EM_ACQBLOCK) {
		if ( fabs(wsp->t - wsp->acqblock_t0) > TINY ) {
			/* there was some previous event not synchronized with dw */
			dt = wsp->dw - (wsp->t - wsp->acqblock_t0);
			if ( duration < dt) {
				_delay_simple(sim,wsp,duration);
				t_tmp = 0;
			} else {
				_delay_simple(sim,wsp,dt);
				t_tmp = duration - dt;
			}
			update_propagator(wsp->U,wsp->dU,sim,wsp);
			/* did it fill time up to dw? */
			if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY ) {
				_store(sim,wsp,wsp->acqblock_sto,0);
				wsp->acqblock_t0 = wsp->t;
				acqblock_sto_incr(wsp);
			}
		} else {
			t_tmp = duration;
		}
		if ( fabs(t_tmp) > TINY) {
			n = (int)floor(t_tmp/wsp->dw+1e-6);
			dt = t_tmp - wsp->dw*(double)n;
			for (i=1; i<=n; i++) {
				_delay_simple(sim,wsp,wsp->dw);
				update_propagator(wsp->U, wsp->dU, sim, wsp);
				_store(sim,wsp,wsp->acqblock_sto,0);
				wsp->acqblock_t0 = wsp->t;
				acqblock_sto_incr(wsp);
			}
			if (dt > TINY) {
				/* remaining time not synchronized with dw */
				_delay_simple(sim,wsp,dt);
				update_propagator(wsp->U, wsp->dU, sim, wsp);
			}
		}
	} else {
		_delay_simple(sim,wsp,duration);  // this just creates wsp->dU, prop of event
		update_propagator(wsp->U, wsp->dU, sim, wsp);
	}
}


/* creates the propagators used for the pulses */
int _setrfprop(Sim_info *sim, Sim_wsp *wsp)
{
  int i,j,spin,isany;

  isany=0;
  blk_dm_zero(wsp->sumHrf);
  dm_zero(wsp->sumUph);

  for (i=1; i<=sim->ss->nchan; i++) {
    if (fabs(wsp->rfv[i]) > TINY) {
      isany++;
      if (wsp->isselectivepulse) {
        for (j=1; j<=sim->ss->nchanelem[i]; j++) {
          spin = sim->ss->chan[i][j];
          if (wsp->spinused[spin]) {
        	  dm_multod(wsp->sumUph,wsp->Iz[spin],wsp->phv[i]);
        	  blk_dm_multod_extract(wsp->sumHrf, wsp->Ix[spin], wsp->rfv[i]);
          }
        }
      } else {
    	  DEBUGPRINT("\n\n --- setrfprop ---\n");
    	  //dm_print(wsp->sumUph,"sumUph");
    	  //dm_print(wsp->chan_Iz[i],"chan_Iz");
    	  if (fabs(wsp->phv[i]) > TINY) dm_multod(wsp->sumUph, wsp->chan_Iz[i], wsp->phv[i]);
    	  //dm_print(wsp->chan_Ix[i],"_setrfprop chan_Ix");
    	  blk_dm_multod_extract(wsp->sumHrf, wsp->chan_Ix[i], wsp->rfv[i]);
      }
    }
  }
  /* only use the selection of spins once */
  wsp->isselectivepulse=0;
  if (wsp->spinused) free_int_vector(wsp->spinused);
  wsp->spinused = NULL;

  //blk_dm_print(wsp->sumHrf,"_setrfprop sumHrf");
  return isany;
}

/* creates the irradiation Hamiltonian for LABframe and DNPframe every maxdt
 * result is in wsp->Hrf_blk
 ******/
int _setrfprop_labdnp(Sim_info *sim, Sim_wsp *wsp)
{
  int i,j,spin,isany1, isany2;
  double w1, wrf;

  assert(wsp->Hrf_blk != NULL);
  assert(wsp->Hrf_blk->Nblocks == 1); // assume no blockdiagonalization in the current implementation

  isany1 = 0;
  isany2 = 0;
  blk_cm_zero(wsp->Hrf_blk);
  blk_dm_zero(wsp->sumHrf);
  dm_zero(wsp->sumUph);

  for (i=1; i<=sim->ss->nchan; i++) {
    if (fabs(wsp->rfv[i]) > TINY) {
      if (wsp->isselectivepulse) {
    	if (sim->frame == DNPFRAME) {
    	  if (sim->ss->iso[spin]->number == 0) { // it is electron
    	    for (j=1; j<=sim->ss->nchanelem[i]; j++) {
    	      spin = sim->ss->chan[i][j];
              if (wsp->spinused[spin]) {
    			cm_multocr(wsp->Hrf_blk->m, wsp->Ix[spin], Complx(wsp->rfv[i],0));
    			if (fabs(wsp->phv[i]) > TINY) {
    			  dm_multod(wsp->sumUph, wsp->Iz[spin], wsp->phv[i]);
    			}
    			isany1++;
    		  }
    		}
  		    //make the rotation for all elements in this channel
    	    if (fabs(wsp->phv[i]) > TINY) blk_simtrans_zrot2(wsp->Hrf_blk,wsp->sumUph);
    		continue; // continue with other channel
    	  }
        }
    	// do other channels as lab frame. Electron in DNPframe should not get here
	    for (j=1; j<=sim->ss->nchanelem[i]; j++) {
	      spin = sim->ss->chan[i][j];
          if (wsp->spinused[spin]) {
        	wrf = sim->ss->chanfreq[i]*2*M_PI;
        	if (wsp->offset != NULL) wrf += wsp->offset[i];  // rf frequency with offset
        	w1 = 2.0*wsp->rfv[i]*cos(wrf*wsp->t + wsp->phv[i]); // amplitude in front of Ix
        	blk_dm_multod_extract(wsp->sumHrf, wsp->Ix[spin], w1);
			isany2++;
		  }
	    }
    	// end of selective
      } else {
        if (sim->frame == DNPFRAME) {
          if (sim->ss->iso[sim->ss->chan[i][1]]->number == 0) { // it is electron channel
    	    cm_multocr(wsp->Hrf_blk->m, wsp->chan_Ix[i], Complx(wsp->rfv[i],0));
    	    if (fabs(wsp->phv[i]) > TINY) {
    		  dm_multod(wsp->sumUph, wsp->chan_Iz[i], wsp->phv[i]);
    		  //make the rotation
    		  blk_simtrans_zrot2(wsp->Hrf_blk,wsp->sumUph);
    	    }
    	    isany1++;
    	    continue; // continue with other channel
          }
        }
        // do labframe stuff for all channels. Electron in DNPframe should not get here
        wrf = sim->ss->chanfreq[i]*2*M_PI;
        if (wsp->offset != NULL) wrf += wsp->offset[i];  // rf frequency with offset
        w1 = 2.0*wsp->rfv[i]*cos(wrf*wsp->t + wsp->phv[i]); // amplitude in front of Ix
        blk_dm_multod_extract(wsp->sumHrf, wsp->chan_Ix[i], w1);
        isany2++;
      }
    }
  }
  if (isany2) { // add real block matrix wsp->sumHrf to complex block matrix wsp->Hrf_blk
	  blk_cm_multocr(wsp->Hrf_blk,wsp->sumHrf,Cunit);
  }

  return isany1+isany2;
}

void _pulse_simple(Sim_info *sim, Sim_wsp *wsp, double duration)
{
	/* this assumes that wsp->sumHrf and wsp->sumUph are already made (at _pulse) */
	int i, n;
	double dt, dt_us;

	assert(wsp->dU != NULL);
	blk_cm_unit(wsp->dU);
	//assert(wsp->tmpU != NULL);
	n = (int)ceil(duration/wsp->dtmax);
	if (n < 1) n = 1;
	dt_us = duration/(double)n;
	dt = dt_us*1.0e-6;
	//printf("_pulse_simple duration of %f us split into %d steps of %f us\n",duration,n,dt*1.0e+6);
	for (i=1;i<=n;i++) {
		//blk_dm_print(wsp->Hiso,"_pulse Hiso PRED");
		ham_hamilton(sim,wsp);
		//blk_dm_print(wsp->ham_blk,"_pulse HAM_int");
		//blk_dm_print(wsp->sumHrf,"_pulse HAM_rf");
		if ( (sim->frame == LABFRAME) || (sim->frame == DNPFRAME) ) {
			//blk_cm_print(wsp->Hcplx,"Hamilton interactions");
			_setrfprop_labdnp(sim, wsp);  // it needs to be called here as the pulse parameters are changing every maxdt (wsp->dtmax)
			blk_cm_multod(wsp->Hcplx,wsp->Hrf_blk,1.0);
			//blk_cm_print(wsp->Hrf_blk,"Hamilton pulse");
			//blk_cm_print(wsp->Hcplx,"Ham total");
			//blk_cm_print(wsp->dU,"dU clean propagator");
			blk_prop_complx_3(wsp->dU,wsp->Hcplx,dt,sim);
			//blk_cm_print(wsp->dU,"propagator");
		} else {
			//blk_dm_print(wsp->ham_blk,"_pulse HAM ints");
			blk_dm_multod(wsp->ham_blk,wsp->sumHrf,1.0);
		//blk_dm_print(wsp->ham_blk,"_pulse HAM_TOT");
		//printf("thread %d, cryst %d: i %d: Ham(1,1) = %10.5f\n",wsp->thread_id,wsp->cryst_idx,i,blk_dm_getelem(wsp->ham_blk,1,1));
		//if (i == 1) {
			blk_prop_real(wsp->dU,wsp->ham_blk,dt,sim);
		}
		//} else {
		//	blk_prop_real(wsp->tmpU,wsp->ham_blk,dt,sim->propmethod);
		//	update_propagator(wsp->dU, wsp->tmpU, sim, NULL);
		//}
		wsp->t += dt_us;
		//complx cprt = blk_cm_getelem(wsp->dU,1,1);
		//printf("thread %d, cryst %d: i %d: dU(1,1) = (%10.5f %10.5f)\n",wsp->thread_id,wsp->cryst_idx,i,cprt.re,cprt.im);

	}
	//exit(1);
	//blk_cm_print(wsp->dU,"_pulse pre-phase propagator");
	if (sim->frame == ROTFRAME) {
		blk_simtrans_zrot2(wsp->dU,wsp->sumUph);
	}
	//dm_print(wsp->sumUph,"_pulse sumUph");
	//blk_cm_print(wsp->dU,"_pulse step propagator");
}

void _pulse(Sim_info *sim, Sim_wsp *wsp, double duration)
{
	int i,n;
	double dt, t_tmp;

	if (duration > wsp->dt_gcompute) {
		duration = wsp->dt_gcompute;
	}
	if (sim->frame == ROTFRAME) { // evaluate setrfprop here, pulse parameters are constant, Hamiltonian can be time dependent due to MAS
		if (!_setrfprop(sim,wsp)) {
			/* ZT: relaxation? */
			if (sim->relax) {
				_delay_relax(sim,wsp,duration);
			} else {
				_delay(sim,wsp,duration);
			}
			return;
		}
	}

	if (wsp->evalmode == EM_ACQBLOCK) {
		if ( fabs(wsp->t - wsp->acqblock_t0) > TINY ) {
			/* there was some previous event not synchronized with dw */
			dt = wsp->dw - (wsp->t - wsp->acqblock_t0);
			if ( duration < dt) {
				_pulse_simple(sim,wsp,duration);
				t_tmp = 0;
			} else {
				_pulse_simple(sim,wsp,dt);
				t_tmp = duration - dt;
			}
			update_propagator(wsp->U, wsp->dU, sim, wsp);
			/* did it fill time up to dw? */
			if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY*2 ) {
				_store(sim,wsp,wsp->acqblock_sto,0);
				wsp->acqblock_t0 = wsp->t;
				acqblock_sto_incr(wsp);
			}
		} else {
			t_tmp = duration;
		}
		if ( fabs(t_tmp) > TINY) {
			n = (int)floor(t_tmp/wsp->dw+1e-6);
			dt = t_tmp - wsp->dw*(double)n;
			for (i=1; i<=n; i++) {
				_pulse_simple(sim,wsp,wsp->dw);
				update_propagator(wsp->U, wsp->dU, sim, wsp);
				_store(sim,wsp,wsp->acqblock_sto,0);
				wsp->acqblock_t0 = wsp->t;
				acqblock_sto_incr(wsp);
			}
			if (dt > TINY) {
				/* remaining time not synchronized with dw */
				_pulse_simple(sim,wsp,dt);
				update_propagator(wsp->U, wsp->dU, sim, wsp);
			}
		}
	} else {
		_pulse_simple(sim,wsp,duration);
		update_propagator(wsp->U, wsp->dU, sim, wsp);
	}
	wsp->waspulse += 1;
}


void _pulseid(Sim_info *sim, Sim_wsp *wsp, double duration)
{
	if (sim->frame != ROTFRAME) {
		fprintf(stderr,"Error: pulseid not implemented for frame %d (checkpoint in _pulseid)\n",sim->frame);
		exit(1);
	}

	if (!_setrfprop(sim,wsp)) {
		/* no delay if there is no rf-pulses */
		return;
	}

  /****** replaced, Uisunit is tested inside update_propagator *************
  if (wsp->Uisunit) {
	  blk_prop_real(wsp->U,wsp->sumHrf,duration*1.0e-6,sim->propmethod);
	  blk_simtrans_zrot2(wsp->U,wsp->sumUph);
	  wsp->Uisunit = 0;
  } else {
	  assert(wsp->dU);
	  blk_prop_real(wsp->dU,wsp->sumHrf,duration*1.0e-6,sim->propmethod);
	  blk_simtrans_zrot2(wsp->dU,wsp->sumUph);
	  update_propagator(wsp->U, wsp->dU, sim, wsp);
  }
  wsp->Uisunit=0;
  ************************** end replaced *********************************/

  assert(wsp->dU != NULL);
  blk_cm_unit(wsp->dU);
  blk_prop_real(wsp->dU,wsp->sumHrf,duration*1.0e-6,sim);
  blk_simtrans_zrot2(wsp->dU,wsp->sumUph);
  update_propagator(wsp->U, wsp->dU, sim, wsp);


  wsp->waspulse += 1;
}

/****
 *  ZT: propagation with simple, static, average hamiltonian - 'avgham_static'
 *      This implementation is as I think it should look like.
 ****/
void _avgham_static(Sim_info *sim, Sim_wsp *wsp, double duration, char* expr)
{
	mat_complx *mtx;
	int i, n;
	double dt, t_tmp;

	if (duration > wsp->dt_gcompute) {
		duration = wsp->dt_gcompute;
	}
	/***
	 * in this function, duration is scaled by 2*PI in all prop calls
	 * in order to take care of Hz->rad.s-1 conversion in mtx
	 ***/

	mtx = ss_readoper(sim,expr);
	if (!cm_ishermit(mtx)) {
		fprintf(stderr,"avgham_static error: expression doesn't give hermitian operator\n");
		exit(1);
	}

	//cm_print(mtx,"_avgham_static matrix");

	assert(wsp->dU != NULL);
	if (wsp->evalmode == EM_ACQBLOCK) {
		if ( fabs(wsp->t - wsp->acqblock_t0) > TINY ) {
			/* there was some previous event not synchronized with dw */
			dt = wsp->dw - (wsp->t - wsp->acqblock_t0);
			if ( duration < dt) {
				blk_prop_complx_2(wsp->dU,mtx,duration*2.0e-6*M_PI,sim->propmethod);
				wsp->t += duration;
				t_tmp = 0;
			} else {
				blk_prop_complx_2(wsp->dU,mtx,dt*2.0e-6*M_PI,sim->propmethod);
				wsp->t += dt;
				t_tmp = duration - dt;
			}
			update_propagator(wsp->U,wsp->dU,sim,wsp);
			/* did it fill time up to dw? */
			if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY ) {
				_store(sim,wsp,wsp->acqblock_sto,0);
				wsp->acqblock_t0 = wsp->t;
				acqblock_sto_incr(wsp);
			}
		} else {
			t_tmp = duration;
		}
		if ( fabs(t_tmp) > TINY) {
			n = (int)floor(t_tmp/wsp->dw+1e-6);
			dt = t_tmp - wsp->dw*(double)n;
			for (i=1; i<=n; i++) {
				blk_prop_complx_2(wsp->dU,mtx,wsp->dw*2.0e-6*M_PI,sim->propmethod);
				update_propagator(wsp->U, wsp->dU, sim, wsp);
				wsp->t += wsp->dw;
				_store(sim,wsp,wsp->acqblock_sto,0);
				wsp->acqblock_t0 = wsp->t;
				acqblock_sto_incr(wsp);
			}
			if (dt > TINY) {
				/* remaining time not synchronized with dw */
				blk_prop_complx_2(wsp->dU,mtx,dt*2.0e-6*M_PI,sim->propmethod);
				update_propagator(wsp->U, wsp->dU, sim, wsp);
				wsp->t += dt;
			}
		}
	} else {
		blk_prop_complx_2(wsp->dU,mtx,duration*2.0e-6*M_PI,sim->propmethod);
		update_propagator(wsp->U, wsp->dU, sim, wsp);
		wsp->t += duration;
	}

	free_complx_matrix(mtx);
	wsp->Uisunit=0;
	wsp->waspulse += 1;
}


int TclGetPhase(Tcl_Interp* interp,Tcl_Obj* obj,double* phase)
{
  char *buf;

  *phase=0;
  if (Tcl_GetDoubleFromObj(interp,obj,phase) != TCL_OK) {
	  buf = Tcl_GetString(obj);
	  if (!strcasecmp(buf,"-X")) {
	    *phase=180;
	  } else if (!strcasecmp(buf,"-Y")) {
	    *phase=270;
	  } else if (!strcasecmp(buf,"X")) {
	    *phase=0;
	  } else if (!strcasecmp(buf,"Y")) {
	    *phase=90;
	  } else {
		  fprintf(stderr,"Error: can not get phase from '%s'",buf);
		  return TCL_ERROR;
	  }
  }
  return TCL_OK;
}

int tclAcq(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double phase;
  int np,prop;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"acq");
  acqblock_disable_command(wsp, "acq");

  if (argc > 4)
    return TclError(interp,"Usage: acq ?phase?\n   or: acq <np> <prop> ?phase?");

//  if (sim->imethod != M_DIRECT) {
//    fprintf(stderr,"error: the 'acq' command can only be used when the computing method is 'direct'\n");
//    exit(1);
//  }
  phase=0.0;
  if (argc == 2 || argc == 4) {
    if (TclGetPhase(interp,argv[argc-1],&phase) != TCL_OK)
       return TCL_ERROR;
  }
  if (argc > 2) {
    /* disable in relaxation mode */
    if (sim->relax)
       return TclError(interp,"acq: using propagators is prohibited in combination with relaxation");

    if (Tcl_GetIntFromObj(interp,argv[1],&np) != TCL_OK)
       return TCL_ERROR;
    if (np < 2)
      return TclError(interp,"acq: number of data points must be larger than 1 if a "
                              "propagator is specified");
    if (Tcl_GetIntFromObj(interp,argv[2],&prop) != TCL_OK) return TCL_ERROR;
    _nacq(sim,wsp,np,prop,phase);
  } else {
    _acq(sim,wsp,phase,1);
  }
  return TCL_OK;
}

void acq_modulus(Sim_info *sim, Sim_wsp *wsp, int do_reset)
{
  complx *ptr;
  double val;
  int i, is, id, NN, k, l, dim;

  NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);
  if (wsp->curr_nsig + 1 > sim->ntot) {
    fprintf(stderr,"error: acq_modulus overflow in fid points\n");
    exit(1);
  }
  (wsp->curr_nsig)++;
  /* in relaxation mode rho is already evolved */
  if (!(sim->relax)) {
     _evolve_with_prop(sim,wsp);
  }
  if (do_reset) _reset_prop(sim,wsp);
  dim = sim->ss->matdim;

  /*** would this work??? ***/
  for (i=0; i<NN; i++) {
	  is = i % sim->Nfstart;
	  id = i % sim->Nfdetect;
	  // basis compatibility check
	  if (wsp->fdetect[id]->basis != wsp->sigma[is]->basis) {
		  mat_complx *dum = cm_change_basis_2(wsp->fdetect[i],wsp->sigma[is]->basis,sim);
		  if (wsp->fdetect[id] != sim->fdetect[id]) free_complx_matrix(wsp->fdetect[id]);
		  wsp->fdetect[id] = dum;
	  }
	  val = 0;
	  for (k=1; k<=dim; k++) {
		  for (l=1; l<=dim; l++) {
			  complx zd = cm_getelem(wsp->fdetect[id],k,l);
			  complx zs = cm_getelem(wsp->sigma[is],k,l);
			  val += zd.re*(sqrt(zs.re*zs.re+zs.im*zs.im));
		  }
	  }
	  ptr = &(wsp->fid[wsp->curr_nsig+i*sim->ntot]);
      ptr->re += val;
  }

}

int tclAcqModulus(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double phase;
  int np,prop;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"acq_modulus");
  acqblock_disable_command(wsp, "acq_modulus");

  if (argc != 1 && argc != 3)
    return TclError(interp,"Usage: acq_modulus\n   or: acq_modulus <np> <prop>");

//  if (sim->imethod != M_DIRECT) {
//	    return TclError(interp,"Error: acq_modulus can be used only for direct method");

  if (argc == 3) {
    /* disable in relaxation mode */
    if (sim->relax)
       return TclError(interp,"acq_modulus: using propagators is prohibited in combination with relaxation");

    if (Tcl_GetIntFromObj(interp,argv[1],&np) != TCL_OK)
       return TCL_ERROR;
    if (np < 2)
      return TclError(interp,"acq_modulus: number of data points must be larger than 1 if a "
                              "propagator is specified");
    if (Tcl_GetIntFromObj(interp,argv[2],&prop) != TCL_OK) return TCL_ERROR;
    return TclError(interp,"acq_modulus: not implemented when using propagators");
  } else {
    acq_modulus(sim,wsp,1);
  }
  return TCL_OK;
}


int tclEvolve(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"evolve");
  acqblock_disable_command(wsp, "evolve");

  /* using within pulseq optimizing propagator is prohibited */
  if (OCpar.gradmodeprop) {
     return TclError(interp,"evolve can not be used when optimizing propagator");
  }

  if (argc != 1)
    return TclError(interp,"Usage: evolve");

  /* check if called in gradient mode */
  if (OCpar.gradmode) {
     /* yes, take care of this situation */
     incr_OCmx_pos(wsp);
     store_OCprop(wsp);
     _evolve_with_prop(sim,wsp);
     _reset_prop(sim,wsp);
     store_OCdens(sim,wsp);
     set_OCmx_code(wsp,"P");
  } else {
     /* no, do usual things */  
     _evolve_with_prop(sim,wsp);
     _reset_prop(sim,wsp);
  }
  
  return TCL_OK;
}


int tclOffset(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  double *offset,off=0.0;
  int i,nchan;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  if (argc < 2)
    return TclError(interp,"Usage: offset <offset(1)/Hz> ?<offset(2)/Hz>? ...");

  nchan=argc-1;
  if (nchan != sim->ss->nchan)
    return TclError(interp,"\nerror: offset: arguments must match the number of channels\n");

  offset = double_vector(nchan);
  for (i=1;i<=nchan;i++) {
    if (Tcl_GetDoubleFromObj(interp,argv[i],&offset[i]) != TCL_OK) return TCL_ERROR;
    offset[i] *= 2.0*M_PI;
    off += fabs(offset[i]);
  }
  if (off < TINY) {
	  free_double_vector(offset);
	  if (wsp->offset) free_double_vector(wsp->offset);
	  wsp->offset =NULL;
  } else {
	  if (wsp->offset) free_double_vector(wsp->offset);
	  wsp->offset = offset;
  }
  return TCL_OK;
}

int tclMaxdt(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double duration;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  if (argc != 2)
     return TclError(interp,"Usage: maxdt <duration/usec>");

  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;

  if (duration < 1e-6)
     return TclError(interp,"maxdt: argument must be larger than 1e-6\n");
  /* ZT: modification to allow maxdt also in relaxation mode */
  if ( (sim->wr == 0.0) && (sim->relax == 0) )
     return TclError(interp,"maxdt: cannot be called if spin-rate is zero\n");
     
  wsp->dtmax = duration;
  return TCL_OK;
}

int tclReset(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double tim_usec;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  acqblock_disable_command(wsp, "reset");

  if (argc != 1 && argc != 2)
      return TclError(interp,"Usage: reset [time increment/usec]");

  tim_usec=0.0;
  if (argc == 2) {
    if (Tcl_GetDoubleFromObj(interp,argv[1],&tim_usec) != TCL_OK) return TCL_ERROR;
  }
  _reset(sim,wsp,tim_usec);
  
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 0;
  }

  return TCL_OK;
}


int tclPulseid(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double duration,rf,ph;
  int i, chan, basis = 0;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  if (argc < 2 || ((argc % 2) != 0))
    return TclError(interp,"Usage: pulseid <duration/usec> ?<rf(1)/Hz> <phase(1)/degrees>? ?<rf(2)> <phase(2)>? ...");

  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
    return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"pulseid: duration must be zero or positive");

  if (argc > 2) {
    if ((argc/2-1) != sim->ss->nchan)
      return TclError(interp,"pulseid: arguments must match number of channels\n");
 
    for (chan=1,i=2;i<argc;i += 2,chan++) {
      if (Tcl_GetDoubleFromObj(interp,argv[i],&rf) != TCL_OK) return TCL_ERROR;

      if (rf < 0.0)
        return TclError(interp,"pulseid: rf field must be positive (change pulse phases if you"
                               " want to rotate the other way.");
      _rf(wsp,chan,rf/wsp->rffac[chan]); // to undo scaling by rf inhomogeneity
      if (rf < TINY) {
    	  // use channels with no rf for blocking
    	  basis += 1 << (chan-1);
      }

      if (TclGetPhase(interp,argv[i+1],&ph) != TCL_OK) return TCL_ERROR;
      _ph(wsp,chan,ph);
    }
  }

  if (sim->Nmz > sim->ss->nchan) {
	  // use dummy channels for blocking
	  for (i=sim->ss->nchan; i<sim->Nmz; i++) basis += 1 << i;
  }

  if (sim->imethod == M_SPINACH) {
	  spinach_pulseid(sim,wsp,duration);
  } else {
	  if (sim->block_diag) {
		  change_basis_pulse(sim, wsp, basis);
		  // IMPROVE: this changes nasis also for H_int which is not necessary!!!
		  if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
		  if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
	  }
    _pulseid(sim,wsp,duration);
  }

  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 1;
  }

  /* do immediate evolution in relaxation mode */
  if (sim->relax) {
     _evolve_with_prop(sim,wsp);
     _reset_prop(sim,wsp);
  }
  
  return TCL_OK;
}

int tclPulse(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double duration,rf,ph;
  int i, chan, basis = 0;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  TIMING_INIT;
  TIMING_INIT_VAR(tv1);
  TIMING_INIT_VAR(tv2);
  TIMING_TIC(tv1);

  read_sim_pointers(interp, &sim, &wsp);

  if (argc < 2 || ((argc % 2) != 0))
    return TclError(interp,"Usage: pulse <duration/usec> ?<rf(1)/Hz> <phase(1)/degrees>? ?<rf(2)> <phase(2)>? ...");

  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;

  if (argc > 2) {
    if ((argc/2-1) != sim->ss->nchan)
      return TclError(interp,"\nerror: pulse: arguments must match number of channels\n");

    for (chan=1,i=2;i<argc;i += 2,chan++) {
      if (Tcl_GetDoubleFromObj(interp,argv[i],&rf) != TCL_OK) return TCL_ERROR;

      if (rf < 0.0)
        return TclError(interp,"pulse: rf field must be positive (change pulse phases if you"
                               " want to rotate the other way.");
      _rf(wsp,chan,rf);

      if (rf < TINY) basis += 1 << (chan-1); // use channels with no rf for blocking

      if (TclGetPhase(interp,argv[i+1],&ph) != TCL_OK) return TCL_ERROR;
      _ph(wsp,chan,ph);
    }
  }

  if (sim->Nmz > sim->ss->nchan) {
	  // use dummy channels for blocking
	  for (i=sim->ss->nchan; i<sim->Nmz; i++) basis += 1 << i;
  }

  if (duration < 0.0)
      return TclError(interp,"pulse: duration must be zero or positive");
  DEBUGPRINT("...pulse duration = %g, maxdt = %g\n",duration,wsp->dtmax);

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  return TCL_OK;
  }

  if (sim->imethod == M_SPINACH) {
	  spinach_pulse(sim,wsp,duration);
  } else {
	  if (sim->block_diag) {
	  	  change_basis_pulse(sim, wsp, basis);
	  	  if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
	  	  if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
	  	  //if (wsp->tmpU == NULL) wsp->tmpU = create_blk_mat_complx_copy2(wsp->ham_blk);
	  }
     _pulse(sim,wsp,duration);
  }

  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 1;
  }

  /* do immediate evolution in relaxation mode */
  if (sim->relax) {
     _evolve_with_prop(sim,wsp);
     _reset_prop(sim,wsp);
  }

  TIMING_TOC(tv1,tv2,"tcl pulse");
  return TCL_OK;
}

/****
 * ZT: helper function for normal pulse_shaped
 ****/
 void _pulse_shaped(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *mask, double steptime)
{
  int i, j;
  
  // tip: in case of normal evolution and blocking it is possible to
  //      create propagator of the whole shape and then update wsp->U:
  //      store wsp->U on safe place, set wsp->U to unit and it will
  //      accumulate the propagator of the rf shape; update wsp->U properly
  for (j=1; j<=Nelem; j++) {
	  for (i=1; i<=sim->ss->nchan; i++) {
		  if (mask[i] == -1) {
			  _rf(wsp,i,0.0);
			  _ph(wsp,i,0.0);
		  } else {
			  _rf(wsp,i,RFshapes[mask[i]][j].ampl);
			  _ph(wsp,i,RFshapes[mask[i]][j].phase);
		  }
	  }
	  _pulse(sim,wsp,steptime);
  }
  
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 1;
  }
  
}

/****
 * ZT: implementation of shaped pulse, using global variable RFshapes[]
 ****/
int tclPulseShaped(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i, j, slot, Nelem=-1, Nch=0, basis=0;
  double duration, steptime;
  int *mask, *OCchanmap = NULL;
  char cd[128], buf2[64];;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"pulse_shaped");

  /* disable when relaxation is ON */
  if (sim->relax) {
     fprintf(stderr,"pulse_shaped error: not supported in combination with relaxation.\n");
     exit(1);
  }
  
  mask = int_vector(sim->ss->nchan);
  
  /* READING ARGUMENTS  */
  if (argc <= 2)
    return TclError(interp,"Usage: pulse_shaped <duration/usec> <RFshape on channel 1> ?<RFshape on channel 2>? ...");
  if (argc-2 != sim->ss->nchan)
    return TclError(interp,"\nerror: pulse_shaped: arguments must match number of channels\n"
                            "        use 'nothing' as place holder when no rf is applied\n");
  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"pulse: duration must be zero or positive");
      
  for (i=1; i<=sim->ss->nchan; i++) {
     if (!strcmp(Tcl_GetString(argv[i+1]),"nothing")) {
        mask[i] = -1;
        basis += 1 << (i-1); // use channels with no rf for blocking
	    DEBUGPRINT("Channel %d: no rf applied\n",i);
	    continue;
     } else {
        /* read RFshape and check it */
        if (Tcl_GetIntFromObj(interp,argv[i+1],&slot) != TCL_OK) {
           return TclError(interp,"error in pulse_shaped: argument %d must be interger <RFshape>",i+1 );
        }
	    if (!RFshapes[slot]) {
           return TclError(interp,"pulse_shaped: argument %d points to non-existing RFshape",i+1);
	    }
	    if (Nelem == -1) {
	    	/* set variable for length of RFshapes */
	    	Nelem=RFshapes_len(slot);
	    } else {
	    	/* check consistency of RFshape lengths */
	    	if ( RFshapes_len(slot) != Nelem )
	    		return TclError(interp,"pulse_shaped: inconsistent number of elements in RFshapes!");
	    }
        DEBUGPRINT("Channel %d: Number of elements in the RFshape = %d\n",i-1,Nelem);
        mask[i] = slot;
     }
   }
   
  if (sim->Nmz > sim->ss->nchan) {
	  // use dummy channels for blocking
	  for (i=sim->ss->nchan; i<sim->Nmz; i++) basis += 1 << i;
  }

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  free_int_vector(mask);
	  return TCL_OK;
  }

   steptime = duration/Nelem;
   
   /* setting propagator flag is done only within _pulse_shaped function */

   if (sim->block_diag) {
	   change_basis_pulse(sim, wsp, basis);
	   if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
	   if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
	   //if (wsp->tmpU == NULL) wsp->tmpU = create_blk_mat_complx_copy2(wsp->ham_blk);
   }

   if (OCpar.gradmode) {
      DEBUGPRINT("pulse_shaped does gradients\n");
      /* determine which channels are active for gradients */
      if (!OCpar.grad_shapes) 
         return TclError(interp,"error when calculating propagators for gradients: grad_shapes not defined");
      int Nsh=LEN(OCpar.grad_shapes);
      if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
    	  cd[0]='\0';
    	  for (i=1; i<=sim->ss->nchan; i++) {
    		  for (j=1; j<=Nsh; j++) {
    			  if (OCpar.grad_shapes[j] == mask[i]) {
    				  sprintf(buf2," I%dC%d",j,i);
    				  strcat(cd,buf2);
    				  Nch++;
    				  break;
    			  }
    		  }
    	  }
    	  /* check if there is any gradient to calculate */
    	  if (Nch==0) {
    		  /* no, then do usual things */
    		  DEBUGPRINT("pulse_shaped - no variable shapes, just creates propagator\n");
    		  _pulse_shaped(sim,wsp,Nelem, mask, steptime);
    	  } else {
    		  /* yes, check for type of optimization */
    		  DEBUGPRINT("pulse_shaped created code '%s'\n",cd);
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  _pulse_shapedOCprops(sim,wsp,cd,Nch,Nelem,mask,steptime);
    		  } else {
    			  /* do stuff for state to state optimization */
    			  _pulse_shapedOC(sim,wsp,cd,Nch,Nelem,mask,steptime);
    		  }
    	  }
      } else { /* GRAPE with advanced gradients */
    	  OCchanmap = int_vector(Nsh);
    	  for (i=1; i<=Nsh; i++) {
    		  OCchanmap[i] = -1;
    		  //printf("grad_shapes[%d] = %d \n",i,OCpar.grad_shapes[i]);
    		  for (j=1; j<=sim->ss->nchan; j++) {
        		  //printf("\t mask[%d] = %d \n",j,mask[j]);
    			  if (OCpar.grad_shapes[i] == mask[j]) {
    				  OCchanmap[i] = j;
    				  Nch++;
        			  break;
    			  }
    		  }
    		  //printf("\t OCchanmap[%d] = %d \n",i,OCchanmap[i]);
    	  }
    	  if (Nch == 0) {
    		  _pulse_shaped(sim,wsp,Nelem, mask, steptime);
    	  } else {
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  switch (sim->frame) {
    			  case ROTFRAME:
    				  _pulse_shapedOCprops_2(sim,wsp,Nelem,OCchanmap,mask,steptime);
    				  break;
    			  case DNPFRAME:
    			  case LABFRAME:
    				  fprintf(stderr,"Error: optimal control on propagators not supported in DNP/LAB frames (checkpoint tclPulseShaped)\n");
    				  exit(1);
    				  break;
    			  default:
    				  fprintf(stderr,"Error: unknown calculation frame %d (checkpoint tclPulseShaped 1)\n",sim->frame);
    				  exit(1);
    			  }

    		  } else {
    			  /* do stuff for state to state optimization */
    			  switch (sim->frame) {
    			  case ROTFRAME:
    				  _pulse_shapedOC_2(sim,wsp,Nelem,OCchanmap,mask,steptime);
    				  break;
    			  case DNPFRAME:
    				  _pulse_shapedOC_2_dnpframe(sim,wsp,Nelem,OCchanmap,mask,steptime);
    				  break;
    			  case LABFRAME:
    				  fprintf(stderr,"Error: optimal control not supported in LABframe (checkpoint tclPulseShaped)\n");
    				  exit(1);
    				  break;
    			  default:
    				  fprintf(stderr,"Error: unknown calculation frame %d (checkpoint tclPulseShaped 1)\n",sim->frame);
    				  exit(1);
    			  }
    		  }
    	  }
    	  free_int_vector(OCchanmap);
      }
   } else {
      /* do just actual pulsing */
      _pulse_shaped(sim,wsp,Nelem, mask, steptime);
   }
   
  free_int_vector(mask);
    
  return TCL_OK;
}

/****
 * ZT: implementation of shaped pulse whose shape gets modulated by rfmap
 ****/
int tclPulseShapedRotorModulated(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i, j, slot, iz, iphi, Nelem=-1, Nch=0, basis=0;
  double duration, steptime, am, ph, bx, by;
  int *mask, *OCchanmap = NULL, *mask2;
  char cd[128], buf2[64];;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"pulse_shaped_rotormodulated");

  /* disable when relaxation is ON */
  if (sim->relax) {
     fprintf(stderr,"pulse_shaped_rotormodulated error: not supported in combination with relaxation.\n");
     exit(1);
  }

  mask = int_vector(sim->ss->nchan);

  /* READING ARGUMENTS  */
  if (argc <= 2)
    return TclError(interp,"Usage: pulse_shaped_rotormodulated <duration/usec> <RFshape on channel 1> ?<RFshape on channel 2>? ... ");
  if (argc-2 != sim->ss->nchan)
    return TclError(interp,"\nerror: pulse_shaped_rotormodulated: arguments must match number of channels\n"
                            "        use 'nothing' as place holder when no rf is applied\n");
  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"pulse_shaped_rotormodulated: duration must be zero or positive");

  for (i=1; i<=sim->ss->nchan; i++) {
     if (!strcmp(Tcl_GetString(argv[i+1]),"nothing")) {
        mask[i] = -1;
        basis += 1 << (i-1); // use channels with no rf for blocking
	    DEBUGPRINT("Channel %d: no rf applied\n",i);
	    continue;
     } else {
        /* read RFshape and check it */
        if (Tcl_GetIntFromObj(interp,argv[i+1],&slot) != TCL_OK) {
           return TclError(interp,"error in pulse_shaped_rotormodulated: argument %d must be interger <RFshape>",i+1 );
        }
	    if (!RFshapes[slot]) {
           return TclError(interp,"pulse_shaped_rotormodulated: argument %d points to non-existing RFshape",i+1);
	    }
	    if (Nelem == -1) {
	    	/* set variable for length of RFshapes */
	    	Nelem=RFshapes_len(slot);
	    } else {
	    	/* check consistency of RFshape lengths */
	    	if ( RFshapes_len(slot) != Nelem )
	    		return TclError(interp,"pulse_shaped_rotormodulated: inconsistent number of elements in RFshapes!");
	    }
        DEBUGPRINT("Channel %d: slot %d Number of elements in the RFshape = %d\n",i-1,slot,Nelem);
        mask[i] = slot;
     }
   }

  if (sim->Nmz > sim->ss->nchan) {
	  // use dummy channels for blocking
	  for (i=sim->ss->nchan; i<sim->Nmz; i++) basis += 1 << i;
  }

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  free_int_vector(mask);
	  return TCL_OK;
  }

   steptime = duration/Nelem;

   /* setting propagator flag is done only within _pulse_shaped function */

   if (sim->block_diag) {
	   change_basis_pulse(sim, wsp, basis);
	   if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
	   if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
	   //if (wsp->tmpU == NULL) wsp->tmpU = create_blk_mat_complx_copy2(wsp->ham_blk);
   }

   // tady by asi mela bejt modulace shape, vytvori se novej mask, kterej se podstrci jen do _pulse*
   // nakonec se  modulovanej shape i mask zahodi
   // !!! vytvareni a ruseni RFshape neni thread safe, je to globalni promenna
   //     proto kazdy vlakno musi mit vlastni RFslot
   mask2 = int_vector(sim->ss->nchan);
   for (i=1; i<=sim->ss->nchan; i++) {
	   if (mask[i] == -1) {
		   mask2[i] = mask[i];
//		   continue;
	   } else {
		   slot = MAXRFSHAPES-1-(i-1)*glob_info.num_threads-wsp->thread_id;  // get new slot index
		   //printf("thread %d uses slot %d for shape %d\n",wsp->thread_id,slot,mask[i]);
		   if (RFshapes[slot] != NULL) {
			   fprintf(stderr,"pulse rotormodulated error: RF slot %d seems occupied\n",slot);
			   exit(1);
		   }
		   RFshapes[slot] = RFshapes_alloc(Nelem); // allocate new shape there
		   mask2[i] = slot;
		   iz = sim->rfmap->loop[(wsp->rf_idx)*2+0];   // z coil coordinate index
		   iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1]; // coil initial "phase" index
		   int Nphi = sim->rfmap->z[iz*(sim->rfmap->Nch+1)+sim->rfmap->Nch]; // number of phi elements for current z position
		   double steptimephi = sim->taur/Nphi;
		   //printf("Timings: shape steptime - rfmap steptimephi = %g - %g = %g\n",steptime,steptimephi,steptime-steptimephi);
		   assert(steptimephi - steptime >= -TINY); // limitation for fine-digitized shapes
		   int iphi_shift = 0;
		   double phitime = 0.0;
		   //printf("Timings: shape steptime = %g, rfmap steptimephi = %g\n",steptime,steptimephi);
		   for (j=1; j<=Nelem; j++) {
//			   rfmap_get(sim->rfmap,iz,iphi+j-1,i-1,&bx,&by); // get distortion coefs
			   rfmap_get(sim->rfmap,iz,iphi+iphi_shift,i-1,&bx,&by); // get distortion coefs
			   am = RFshapes[mask[i]][j].ampl;  // original ampl/phase element
			   ph = RFshapes[mask[i]][j].phase;
			   RFshapes[mask2[i]][j].ampl = am*sqrt(bx*bx+by*by);
			   RFshapes[mask2[i]][j].phase = ph+RAD2DEG*atan2(by,bx);
			   //printf("Thread %d: chan %d, slot %3d: phitime = %g ... %d %g %g\n",wsp->thread_id,i,slot,phitime,j,bx,by);
			   phitime += steptime;
			   //if (phitime >= steptimephi) {
			   if (steptimephi-phitime < 0.5*steptime) {
				   phitime -= steptimephi;
				   iphi_shift++;
			   }
		   }
	   }
   }

   if (OCpar.gradmode) {
      DEBUGPRINT("pulse_shaped_rotormodulated does gradients\n");
      /* determine which channels are active for gradients */
      if (!OCpar.grad_shapes)
         return TclError(interp,"error when calculating propagators for gradients: grad_shapes not defined");
      int Nsh=LEN(OCpar.grad_shapes);
      if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
    	  cd[0]='\0';
    	  for (i=1; i<=sim->ss->nchan; i++) {
    		  for (j=1; j<=Nsh; j++) {
    			  if (OCpar.grad_shapes[j] == mask[i]) {
    				  sprintf(buf2," I%dC%d",j,i);
    				  strcat(cd,buf2);
    				  Nch++;
    				  break;
    			  }
    		  }
    	  }
    	  /* check if there is any gradient to calculate */
    	  if (Nch==0) {
    		  /* no, then do usual things */
    		  DEBUGPRINT("pulse_shaped_rotormodulated - no variable shapes, just creates propagator\n");
    		  _pulse_shaped(sim,wsp,Nelem, mask2, steptime); // propagate with distorted shapes
    	  } else {
    		  /* yes, check for type of optimization */
    		  DEBUGPRINT("pulse_shaped_rotormodulated created code '%s'\n",cd);
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  _pulse_shapedOCprops(sim,wsp,cd,Nch,Nelem,mask2,steptime); // propagate with distorted shapes
    		  } else {
    			  /* do stuff for state to state optimization */
    			  _pulse_shapedOC(sim,wsp,cd,Nch,Nelem,mask2,steptime); // propagate with distorted shapes
    		  }
    	  }
      } else { /* GRAPE with advanced gradients */
    	  OCchanmap = int_vector(Nsh);
    	  for (i=1; i<=Nsh; i++) {
    		  OCchanmap[i] = -1;
    		  //printf("grad_shapes[%d] = %d \n",i,OCpar.grad_shapes[i]);
    		  for (j=1; j<=sim->ss->nchan; j++) {
        		  //printf("\t mask[%d] = %d \n",j,mask[j]);
    			  if (OCpar.grad_shapes[i] == mask[j]) {
    				  OCchanmap[i] = j;
    				  Nch++;
        			  break;
    			  }
    		  }
    		  //printf("\t OCchanmap[%d] = %d \n",i,OCchanmap[i]);
    	  }
    	  if (Nch == 0) {
    		  _pulse_shaped(sim,wsp,Nelem, mask2, steptime); // propagate with distorted shapes
    	  } else {
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  _pulse_shapedOCprops_2(sim,wsp,Nelem,OCchanmap,mask2,steptime); // propagate with distorted shapes
    		  } else {
    			  /* do stuff for state to state optimization */
    			  _pulse_shapedOC_2(sim,wsp,Nelem,OCchanmap,mask2,steptime); // propagate with distorted shapes
    		  }

    	  }
    	  free_int_vector(OCchanmap);
      }
   } else {
      /* do just actual pulsing */
      _pulse_shaped(sim,wsp,Nelem, mask2, steptime); // propagate with distorted shapes
   }

  free_int_vector(mask);
  for (i=1; i<=LEN(mask2); i++) {
	  if (mask2[i] != -1) free_RFshapes(mask2[i]);
  }
  free_int_vector(mask2);

  return TCL_OK;
}

/****
 * ZT: helper function for normal pulse_shaped_linear
 ****/
 void _pulse_shaped_linear(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *mask, double steptime)
{
	 int i, n, Ntot, irf;
	 //double dt;
	 int Nchan = sim->ss->nchan;
	 int chan;
	 double *wx_left, *wy_left, *wx_right, *wy_right, *w_swap;
	 double wam, wph, wxA, wxB, wyA, wyB;
	 blk_mat_complx *HRF_left, *HRF_right;

	 /* current limitation */
	 if (wsp->evalmode == EM_ACQBLOCK) {
		 fprintf(stderr,"pulse_shaped_linear cannot be used inside acq_block - not implemented\n");
		 exit(1);
	 }

	 /* initialize */
	 assert(wsp->dU != NULL);
	 blk_cm_unit(wsp->dU);
	 wx_left = double_vector(Nchan);
	 wy_left = double_vector(Nchan);
	 wx_right = double_vector(Nchan);
	 wy_right = double_vector(Nchan);
	 /* time steps */
	 n = (int)ceil(steptime/wsp->dtmax);
	 Ntot = (Nelem-1)*n;
	 printf("pulse_shaped_linear: steptime %g split into %d steps of %g\n",steptime,n,wsp->dtmax);
	 if (fabs(n*wsp->dtmax - steptime) > TINY) {
		 fprintf(stderr,"pulse_shaped_linear error: steptime is not multiple of maxdt\n");
		 exit(1);
	 }
	 printf("Total pulse_shaped_linear duration is (%d - 1) x %g = %g us, \n",Nelem, steptime, (Nelem-1)*steptime);
	 printf("   --> will execute (%d - 1) x %d = %d steps\n",Nelem, n, Ntot);
	 /* calculate left Hamiltonian at time zero */
	 irf = 1;
	 printf("pulse_shaped_linear starts at time %g us\n",wsp->t);
	 ham_hamilton(sim,wsp); // H_int is in wsp->ham_blk, real matrix
	 for (chan=1; chan<=Nchan; chan++) { // prepares (wx,wy) rf parameters
		 if (mask[chan] == -1) {
			 /* no rf applied on this channel */
			 wam = 0.0;
			 wph = 0.0;
		 } else {
			 wam = RFshapes[mask[chan]][irf].ampl;
			 wph = RFshapes[mask[chan]][irf].phase*DEG2RAD;
		 }
		 /* set Ix and Iy components of the rf Hamiltonian */
		 wx_left[chan] = wam*cos(wph);
		 wy_left[chan] = wam*sin(wph);
		 //printf("   chan %d : %g and %g\n",chan,wx_left[chan], wy_left[chan]);
	 }
	 // _setrfprop_linear(sim, wsp, wx_left, wy_left); // create RF Hamiltonian
	 //
	 for (i=1; i<=Ntot; i++) {
		 printf("Executing step %d\n",i);
		 /* calculate right Hamiltonian */
		 //printf("Ham-R at time %g (irf = %d)\n",i*wsp->dtmax, irf);
		 wsp->t += wsp->dtmax;
		 printf("Ham-R calculated at time %g us\n",wsp->t);
		 for (chan=1; chan<=Nchan; chan++) {
			 if (mask[chan] == -1) {
				 /* no rf applied on this channel */
				 wxA = 0.0;
				 wyA = 0.0;
				 wxB = 0.0;
				 wyB = 0.0;
			 } else {
				 wam = RFshapes[mask[chan]][irf].ampl;
				 wph = RFshapes[mask[chan]][irf].phase*DEG2RAD;
				 wxA = wam*cos(wph);
				 wyA = wam*sin(wph);
				 wam = RFshapes[mask[chan]][irf+1].ampl;
				 wph = RFshapes[mask[chan]][irf+1].phase*DEG2RAD;
				 wxB = wam*cos(wph);
				 wyB = wam*sin(wph);
			 }
			 /* set Ix and Iy components of the rf Hamiltonian */
			 wx_right[chan] = wx_left[chan]+(wxB-wxA)/steptime*wsp->dtmax;
			 wy_right[chan] = wy_left[chan]+(wyB-wyA)/steptime*wsp->dtmax;
			 printf("   chan %d : L = (%g, %g); R = (%g, %g)\n",chan,wx_left[chan], wy_left[chan],wx_right[chan], wy_right[chan]);
		 }

		 /* assemble exponent for propagator calculation */
		 printf("calculating exponent\n");
		 /* create propagator and update wsp-dU */
		 printf("prop and update of wsp->dU\n");
		 /* right Hamiltonian becomes left Hamiltonian */
		 printf("Ham-R becomes Ham-L at relative-pulse-time %g\n",i*wsp->dtmax);
		 w_swap = wx_left; wx_left = wx_right; wx_right = w_swap;
		 w_swap = wy_left; wy_left = wy_right; wy_right = w_swap;
		 if (i % n == 0) {
			 irf++;
		 }
	 }

	 // clear memory usage
	 free_double_vector(wx_left);
	 free_double_vector(wy_left);
	 free_double_vector(wx_right);
	 free_double_vector(wy_right);
	 printf("finishing with pulse_shaped_linear at time %g us\n",wsp->t);

	 /* add total pulse propagator to total evolution propagator */
	 update_propagator(wsp->U, wsp->dU, sim, wsp);
	 /* book keeping at the end */
	 wsp->waspulse += 1;
	 /* if used within pulseq optimizing propagator set this flag */
	if (OCpar.gradmodeprop) {
		wsp->OC_propstatus = 1;
	}

}

/****
 * ZT: implementation of shaped pulse whose shape gets modulated by rfmap
 *     implementing PIECEWISE-LINEAR rf shape
 ****/
int tclPulseShapedRotorModulatedLinear(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i, j, slot, iz, iphi, Nelem=-1, Nch=0, basis=0;
  double duration, steptime, am, ph, bx, by;
  int *mask, *OCchanmap = NULL, *mask2;
  char cd[128], buf2[64];;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"pulse_shaped_linear");

  /* disable when relaxation is ON */
  if (sim->relax) {
     fprintf(stderr,"pulse_shaped_linear error: not supported in combination with relaxation.\n");
     exit(1);
  }

  mask = int_vector(sim->ss->nchan);

  /* READING ARGUMENTS  */
  if (argc <= 2)
    return TclError(interp,"Usage: pulse_shaped_linear <duration/usec> <RFshape on channel 1> ?<RFshape on channel 2>? ... ");
  if (argc-2 != sim->ss->nchan)
    return TclError(interp,"\nerror: pulse_shaped_linear: arguments must match number of channels\n"
                            "        use 'nothing' as place holder when no rf is applied\n");
  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"pulse_shaped_linear: duration must be zero or positive");

  for (i=1; i<=sim->ss->nchan; i++) {
     if (!strcmp(Tcl_GetString(argv[i+1]),"nothing")) {
        mask[i] = -1;
        basis += 1 << (i-1); // use channels with no rf for blocking
	    DEBUGPRINT("Channel %d: no rf applied\n",i);
	    continue;
     } else {
        /* read RFshape and check it */
        if (Tcl_GetIntFromObj(interp,argv[i+1],&slot) != TCL_OK) {
           return TclError(interp,"error in pulse_shaped_linear: argument %d must be interger <RFshape>",i+1 );
        }
	    if (!RFshapes[slot]) {
           return TclError(interp,"pulse_shaped_linear: argument %d points to non-existing RFshape",i+1);
	    }
	    if (Nelem == -1) {
	    	/* set variable for length of RFshapes */
	    	Nelem=RFshapes_len(slot);
	    } else {
	    	/* check consistency of RFshape lengths */
	    	if ( RFshapes_len(slot) != Nelem )
	    		return TclError(interp,"pulse_shaped_linear: inconsistent number of elements in RFshapes!");
	    }
        DEBUGPRINT("Channel %d: slot %d Number of elements in the RFshape = %d\n",i-1,slot,Nelem);
        mask[i] = slot;
     }
   }

  if (sim->Nmz > sim->ss->nchan) {
	  // use dummy channels for blocking
	  for (i=sim->ss->nchan; i<sim->Nmz; i++) basis += 1 << i;
  }

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  free_int_vector(mask);
	  return TCL_OK;
  }

  /* piecewise-linear rfshape contains (n+1) values and n intervals */
  steptime = duration/(Nelem-1);

   /* setting propagator flag is done only within _pulse_shaped function */

   if (sim->block_diag) {
	   change_basis_pulse(sim, wsp, basis);
	   if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
	   if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
   }

   // tady by asi mela bejt modulace shape, vytvori se novej mask, kterej se podstrci jen do _pulse*
   // nakonec se  modulovanej shape i mask zahodi
   // !!! vytvareni a ruseni RFshape neni thread safe, je to globalni promenna
   //     proto kazdy vlakno musi mit vlastni RFslot
   mask2 = int_vector(sim->ss->nchan);
   for (i=1; i<=sim->ss->nchan; i++) {
	   if (mask[i] == -1) {
		   mask2[i] = mask[i];
//		   continue;
	   } else {
		   slot = MAXRFSHAPES-1-(i-1)*glob_info.num_threads-wsp->thread_id;  // get new slot index
		   //printf("thread %d uses slot %d for shape %d\n",wsp->thread_id,slot,mask[i]);
		   if (RFshapes[slot] != NULL) {
			   fprintf(stderr,"pulse_shaped_linear error: RF slot %d seems occupied\n",slot);
			   exit(1);
		   }
		   RFshapes[slot] = RFshapes_alloc(Nelem); // allocate new shape there
		   mask2[i] = slot;
		   iz = sim->rfmap->loop[(wsp->rf_idx)*2+0];   // z coil coordinate index
		   iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1]; // coil initial "phase" index
		   int Nphi = sim->rfmap->z[iz*(sim->rfmap->Nch+1)+sim->rfmap->Nch]; // number of phi elements for current z position
		   double steptimephi = sim->taur/Nphi;
		   //printf("Timings: shape steptime - rfmap steptimephi = %g - %g = %g\n",steptime,steptimephi,steptime-steptimephi);
		   assert(steptimephi - steptime >= -TINY); // limitation for fine-digitized shapes
		   int iphi_shift = 0;
		   double phitime = 0.0;
		   //printf("Timings: shape steptime = %g, rfmap steptimephi = %g\n",steptime,steptimephi);
		   for (j=1; j<=Nelem; j++) {
//			   rfmap_get(sim->rfmap,iz,iphi+j-1,i-1,&bx,&by); // get distortion coefs
			   rfmap_get(sim->rfmap,iz,iphi+iphi_shift,i-1,&bx,&by); // get distortion coefs
			   am = RFshapes[mask[i]][j].ampl;  // original ampl/phase element
			   ph = RFshapes[mask[i]][j].phase;
			   RFshapes[mask2[i]][j].ampl = am*sqrt(bx*bx+by*by);
			   RFshapes[mask2[i]][j].phase = ph+RAD2DEG*atan2(by,bx);
			   //printf("Thread %d: chan %d, slot %3d: phitime = %g ... %d %g %g\n",wsp->thread_id,i,slot,phitime,j,bx,by);
			   phitime += steptime;
			   //if (phitime >= steptimephi) {
			   if (steptimephi-phitime < 0.5*steptime) {
				   phitime -= steptimephi;
				   iphi_shift++;
			   }
		   }
	   }
   }

   if (OCpar.gradmode) {
      printf("pulse_shaped_linear does not know how to do gradients...");
      exit(1);
      DEBUGPRINT("pulse_shaped_linear does gradients\n");
      /* determine which channels are active for gradients */
      if (!OCpar.grad_shapes)
         return TclError(interp,"error when calculating propagators for gradients: grad_shapes not defined");
      int Nsh=LEN(OCpar.grad_shapes);
      if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
    	  cd[0]='\0';
    	  for (i=1; i<=sim->ss->nchan; i++) {
    		  for (j=1; j<=Nsh; j++) {
    			  if (OCpar.grad_shapes[j] == mask[i]) {
    				  sprintf(buf2," I%dC%d",j,i);
    				  strcat(cd,buf2);
    				  Nch++;
    				  break;
    			  }
    		  }
    	  }
    	  /* check if there is any gradient to calculate */
    	  if (Nch==0) {
    		  /* no, then do usual things */
    		  DEBUGPRINT("pulse_shaped_linear - no variable shapes, just creates propagator\n");
    		  _pulse_shaped_linear(sim,wsp,Nelem, mask2, steptime); // propagate with distorted shapes
    	  } else {
    		  /* yes, check for type of optimization */
    		  DEBUGPRINT("pulse_shaped_linear created code '%s'\n",cd);
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  _pulse_shapedOCprops(sim,wsp,cd,Nch,Nelem,mask2,steptime); // propagate with distorted shapes
    		  } else {
    			  /* do stuff for state to state optimization */
    			  _pulse_shapedOC(sim,wsp,cd,Nch,Nelem,mask2,steptime); // propagate with distorted shapes
    		  }
    	  }
      } else { /* GRAPE with advanced gradients */
    	  OCchanmap = int_vector(Nsh);
    	  for (i=1; i<=Nsh; i++) {
    		  OCchanmap[i] = -1;
    		  //printf("grad_shapes[%d] = %d \n",i,OCpar.grad_shapes[i]);
    		  for (j=1; j<=sim->ss->nchan; j++) {
        		  //printf("\t mask[%d] = %d \n",j,mask[j]);
    			  if (OCpar.grad_shapes[i] == mask[j]) {
    				  OCchanmap[i] = j;
    				  Nch++;
        			  break;
    			  }
    		  }
    		  //printf("\t OCchanmap[%d] = %d \n",i,OCchanmap[i]);
    	  }
    	  if (Nch == 0) {
    		  _pulse_shaped_linear(sim,wsp,Nelem, mask2, steptime); // propagate with distorted shapes
    	  } else {
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  _pulse_shapedOCprops_2(sim,wsp,Nelem,OCchanmap,mask2,steptime); // propagate with distorted shapes
    		  } else {
    			  /* do stuff for state to state optimization */
    			  _pulse_shapedOC_2(sim,wsp,Nelem,OCchanmap,mask2,steptime); // propagate with distorted shapes
    		  }

    	  }
    	  free_int_vector(OCchanmap);
      }
   } else {
      /* do just actual pulsing */
      _pulse_shaped_linear(sim,wsp,Nelem, mask2, steptime); // propagate with distorted shapes
   }

  free_int_vector(mask);
  for (i=1; i<=LEN(mask2); i++) {
	  if (mask2[i] != -1) free_RFshapes(mask2[i]);
  }
  free_int_vector(mask2);

  return TCL_OK;
}



int tclDelay(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double duration;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  //DEBUGPRINT("--- sim=%p, wsp=%p ---\n",sim,wsp);
  read_sim_pointers(interp, &sim, &wsp);
  //DEBUGPRINT("--- sim=%p, wsp=%p ---\n",sim,wsp);

  if (argc != 2)
      return TclError(interp,"Usage: delay <duration/usec>");

  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
      return TCL_ERROR;
  //DEBUGPRINT("tclDelay duration = %f\n",duration);
  //DEBUGPRINT("tclDelay - accessing sim->pulsename: \n");
  //DEBUGPRINT("'%s'\n",sim->pulsename);
  if (duration < 0.0)
      return TclError(interp,"delay: duration must be zero or positive");
  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  return TCL_OK;
  }

  if (sim->block_diag) {
  	  change_basis_delay(sim, wsp);
  	  if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
  	  if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
  	  //if (wsp->tmpU == NULL) wsp->tmpU = create_blk_mat_complx_copy2(wsp->ham_blk);
  }

  /* ZT: relaxation? */
  if (sim->relax) {
       _delay_relax(sim,wsp,duration);
  } else if (sim->imethod == M_SPINACH) {
	  spinach_delay(sim,wsp,duration);
  } else {
       _delay(sim,wsp,duration);
  }
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 1;
  }
  
  return TCL_OK;
}

int tclSelect(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int spin,i;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  if (argc < 2)
    return TclError(interp,"Usage: select <n> <n> ...");

  wsp->isselectivepulse=1;
  if (wsp->spinused == NULL) {
	  wsp->spinused = int_vector(sim->ss->nspins);
  }
  iv_zero(wsp->spinused);


  for (i=1;i<argc;i++) {
    if (Tcl_GetIntFromObj(interp,argv[i],&spin) != TCL_OK)
      return TCL_ERROR;
    if (spin < 1 || spin > sim->ss->nspins)
       return TclError(interp,"\nerror: select: spin out of range (%d not between 1 and %d)\n",spin,sim->ss->nspins);
    wsp->spinused[spin]=1;
    if (wsp->Ix[spin] == NULL) {
    	wsp->Ix[spin] = Ix_real(sim,spin);
    	if (sim->basis != 0) {
    		dm_permute(wsp->Ix[spin],sim->perm_table[0+sim->basis*sim->Nbasis]);
    		wsp->Ix[spin]->basis = sim->basis;
    	}
    }
    if (wsp->Iz[spin] == NULL) {
    	wsp->Iz[spin] = Iz_ham(sim,spin); // this already permutes...
    }
  }
  return TCL_OK;
}


/****
 * ZT: added more options so the 'store' command can be used to store propagator from a matrix
 ****/
int tclStore(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int num, adj=0;
  char *buf;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"store");
  acqblock_disable_command(wsp, "store");

  /* disable when relaxation is ON */
  if (sim->relax) {
     fprintf(stderr,"store error: cannot be used in combination with relaxation.\n");
     exit(1);
  }
  
  if ( (argc<2) || (argc>5) ) {
     return TclError(interp,"Usage: store <number>  [adjoint]             \n"
                            "       store <number> <matrix number>  [adjoint]\n"
			    "       store <number> list {...}  [adjoint]     \n"
			    "       store <number> elements {...}  [adjoint] \n"
			    "       store <number> notelements {...} [adjoint]\n");
  }
  
  if (Tcl_GetIntFromObj(interp,argv[1],&num) != TCL_OK)
     return TclError(interp,"store first argument must be integer");
  buf = Tcl_GetString(argv[argc-1]);
  if (!strcmp(buf,"adjoint")) {
	  adj = 1;
	  argc--;
  }
  if (argc==2) {
     /* usual, original way */ 
     _store(sim,wsp,num,adj);
     return TCL_OK;
  } else if (argc==3) {
     int mxidx;
     
     if (Tcl_GetIntFromObj(interp,argv[2],&mxidx) != TCL_OK)
        return TclError(interp,"store second argument must be integer, usage: store <num> <matrix num>");
     if (mxidx<1 || mxidx>MAXSTO)
        return TclError(interp,"store: matrix argument out of range <1,%d>",MAXSTO);
     if (!wsp->matrix[mxidx])
        return TclError(interp,"store: matrix number %d is undefined",mxidx);
     if ( (wsp->matrix[mxidx]->row != sim->matdim) || (wsp->matrix[mxidx]->col != sim->matdim) )
        return TclError(interp,"store: mismatch in matrix dimensions");
	
     if (wsp->STO[num] != NULL) free_blk_mat_complx(wsp->STO[num]);
     wsp->STO[num] = create_blk_mat_complx(sim->matdim,1,NULL,wsp->matrix[mxidx]->type,wsp->matrix[mxidx]->basis);
     cm_copy(wsp->STO[num]->m,wsp->matrix[mxidx]);
     if (adj) cm_adjointi(wsp->STO[num]->m);
     
  } else {
     /* case with argc==4 */
     mat_complx *m;
     buf = Tcl_GetString(argv[2]);
     if ( !strcmp(buf,"list") || !strcmp(buf,"elements") || !strcmp(buf,"notelements") ) {
        m = get_matrix_2(sim,wsp,buf,Tcl_GetString(argv[3]));
	if (!m) return TclError(interp,"store could not set propagator from %s",buf);
     } else {
        return TclError(interp,"store called with wrong argument %s",buf);
     }
     if ( (m->row != sim->matdim) || (m->col != sim->matdim) )
        return TclError(interp,"store: mismatch in matrix dimensions");
     if (adj) cm_adjointi(m);
     /* take care if STO[num] was already alocated */
     if (wsp->STO[num] != NULL) free_blk_mat_complx(wsp->STO[num]);
     wsp->STO[num] = create_blk_mat_complx(sim->matdim,1,NULL,m->type,m->basis);
     cm_copy(wsp->STO[num]->m,m);
     free_complx_matrix(m);
  }
  /* this is done only for cases argc==3 or 4. As a consequence, time checking is not done when reusing this propagator */
  wsp->STO_tproplength_usec[num]=0.0;
  wsp->STO_tpropstart_usec[num]=-1.0; /*note negative time here */
  wsp->STO_waspulse[num]=wsp->waspulse;
     
  return TCL_OK;
}


int tclProp(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int num,ntimes = 1;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;
  const char usage[] = "Usage: prop <number> [<times>] [-phase <degrees>]";
  double phase = 0;
  char *buf;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"prop");

  if (argc < 2 || argc > 5 )
	  return TclError(interp,"Error in prop: wrong number of arguments\n%s\n",usage);
  if (Tcl_GetIntFromObj(interp,argv[1],&num) != TCL_OK)
      return TclError(interp,"Error in prop: first argument must be integer\n%s\n",usage);

  switch (argc) {
  case 2:
	  ntimes = 1;
	  break;
  case 3:
	  if (Tcl_GetIntFromObj(interp,argv[2],&ntimes) != TCL_OK)
		  return TclError(interp,"Error in prop: can not get integer from second argument\n%s\n",usage);
	  break;
  case 4:
	  buf = Tcl_GetString(argv[2]);
	  if (!strcmp(buf,"-phase")) {
		  return TclError(interp,"Error in prop: option -phase not implemented yet\n");
		  //if (Tcl_GetDoubleFromObj(interp,argv[3],&phase) != TCL_OK)
		  //	  return TclError(interp,"Error in prop: can not get double from 3rd argument\n%s\n",usage);
	  } else {
		  return TclError(interp,"Error in prop: wrong option '%s'\n%s\n",buf,usage);
	  }
	  break;
  case 5:
	  if (Tcl_GetIntFromObj(interp,argv[2],&ntimes) != TCL_OK)
		  return TclError(interp,"Error in prop: can not get integer from second argument\n%s\n",usage);
	  buf = Tcl_GetString(argv[3]);
	  if (!strcmp(buf,"-phase")) {
		  return TclError(interp,"Error in prop: option -phase not implemented yet\n");
		  //if (Tcl_GetDoubleFromObj(interp,argv[4],&phase) != TCL_OK)
		  //	  return TclError(interp,"Error in prop: can not get double from 4th argument\n%s\n",usage);
	  } else {
		  return TclError(interp,"Error in prop: wrong option '%s'\n%s\n",buf,usage);
	  }
  }

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += ntimes*wsp->STO_tproplength_usec[num];
	  return TCL_OK;
  }

  if (sim->block_diag) {
  	  if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
  }

  if (ntimes > 0) _prop(sim,wsp,num,ntimes);
  
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 1;
  }
  
  return TCL_OK;
}

 
/****
 * ZT: modification to properly deal OC situations
 ****/
int tclFilter(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int num, proj = 0;
  char *buf;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"filter");
  acqblock_disable_command(wsp, "filter");

  /* usage within pulseq optimizing propagator is prohibited */
  if (OCpar.gradmodeprop) {
     return TclError(interp,"filter can not be used when optimizing propagator!\n");     
  }
  
  if (argc != 2 && argc != 3)
      return TclError(interp,"Usage: filter <number> [-projection]");
  if (Tcl_GetIntFromObj(interp,argv[1],&num) != TCL_OK)
     return TclError(interp,"Error in filter: argument must be integer\n");
  if (argc == 3) {
	  buf = Tcl_GetString(argv[2]);
	  if (!strcmp(buf,"-projection")) {
		  proj = 1;
	  } else {
		  return TclError(interp,"Error in filter: wrong option '%s'\n",buf);
	  }
  }

  /* check if called in gradient mode */
  if (OCpar.gradmode) {
     /* yes, take care of this situation */
	 if (proj == 1) {
		//return TclError(interp,"Error: filter -projection not allowed in OC\n");
		_filter(sim,wsp,num,proj);
	 } else {
        _filterOC(sim,wsp,num);
	 }
  } else {
     /* no, do usual things */ 
     _filter(sim,wsp,num,proj);
  }
  return TCL_OK;
}



int tclMatrixoper(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i, j[3];
  mat_complx *m;
  complx val1, val2;
  char *buf;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  
  if (argc != 4 && argc != 7) 
    return TclError(interp,
    "Usage:\n"
    "matrixoper <whichmatrix> <ir> <ic> [(+|-|*|/) <jr> <jc>]\n");

  if (!(m = get_matrix_1(sim, wsp, Tcl_GetString(argv[1]))))
    return TclError(interp, "Error: matrixoper - Specify valid matrix\n");
  
  for (i=1; i<=2; i++) {
    if (Tcl_GetIntFromObj(interp, argv[i+1], &j[i]) != TCL_OK)
      return TclError(interp, "Error: matrixoper - Specify valid number\n");
    if (j[i] < 1 || j[i] > sim->matdim)
      return TclError(interp, "Error: matrixoper - Invalid range %d (should be 1 <= n <= %d)\n", j[i], sim->matdim);
  }
  val1 = cm_getelem(m,j[1],j[2]);

  if (argc == 7) {
    buf = Tcl_GetString(argv[4]);
    for (i=1; i<=2; i++) {
      if (Tcl_GetIntFromObj(interp, argv[i+4], &j[i]) != TCL_OK)
        return TclError(interp, "Error: matrixoper - Specify valid number\n");
      if (j[i] < 1 || j[i] > sim->matdim)
        return TclError(interp, "Error: matrixoper - Invalid range %d (should be 1 <= n <= %d)\n", j[i],sim->matdim);
    }
    val2 = cm_getelem(m,j[1],j[2]);
    switch (*buf) {
    case '+': val1 = Cadd(val1, val2); break;
    case '-': val1 = Csub(val1, val2); break;
    case '*': val1 = Cmul(val1, val2); break;
    case '/': val1 = Cdiv(val1, val2); break;
    default:
      return TclError(interp, "Error: matrixoper - Specify valid operator (+|-|*|/)");
    }
  }
  
  TclAppendResult(interp, "%.8e %.8e", val1.re, val1.im);
  free_complx_matrix(m);
  return TCL_OK;
}

int tclMatrix(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int num;
  mat_complx * m;
  int arg;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  if (argc < 3 || argc > 5)
      return TclError(interp,
      "Usage:\n"
      "  matrix set <to> <from>\n"
      "  <result> matrix get <from>\n"
      "\n"
      "  <to> can be:\n"
      "       <number>\n"
      "       start\n"
      "       detect\n"
      "\n"
      "  <from> can be:\n"
      "       <number>\n"
      "       start\n"
      "       detect\n"
      "       density\n"
      "       propagator\n"
      "       hamiltonian\n"
      "       avgham\n"
      "       operator {I2x+I3y}\n"
      "       totalcoherence {dm1 dm2 ..}\n"
      "       coherence {{dm1 .. dmN} {dm1 .. dmN}} ..., where N is number of nuclei\n"
      "       list {{c11 c12 ..} {c21 c22 ..} ..}, where c is <re> or {<re> <im>}\n"
      "       elements {{i j} {i j} ..}\n"
      "       notelements {{i j} {i j} ..}\n"
      );
  
  arg=1;
  if (!strcmp(argv[arg],"get")) {
    if (arg + 2 == argc) {
      if ( !(m = get_matrix_1(sim,wsp,argv[arg+1])))
        return TCL_ERROR;
        
    } else if (arg + 3 == argc) {
      if ( !(m = get_matrix_2(sim,wsp,argv[arg+1],argv[arg+2])))
        return TCL_ERROR;
    } else
      return TclError(interp,"error: 'matrix get' must have one or two additional arguments\n");

    Tcl_ResetResult(interp);
    if (TclAppendMatrix(interp,m) == TCL_ERROR)
      return TCL_ERROR;
    free_complx_matrix(m);
    
  } else if (!strcmp(argv[arg],"set")) {
    if (arg + 3 == argc) {
      if ( !(m = get_matrix_1(sim,wsp,argv[arg+2])))
        return TCL_ERROR;
    } else if (arg + 4 == argc) {
      if ( !(m = get_matrix_2(sim,wsp,argv[arg+2],argv[arg+3])))
        return TCL_ERROR;
    } else
      return TclError(interp,"error: 'matrix set' must have two or three additional arguments\n");
    arg++;
    if (Tcl_GetInt(interp,argv[arg],&num) == TCL_OK) {    
      if (num < 1 || num > MAXSTO)
        return TclError(interp,"matrix: called with number larger than %d\n",MAXSTO);
      if (wsp->matrix[num]) {
        free_complx_matrix(wsp->matrix[num]);
      }
      wsp->matrix[num] = m;
      //DEBUGPRINT("tclMatrix sets slot %d\n",num);
      //cm_print(m,"tclMatrix");
    } else if (!strcmp(argv[arg],"start")) {
      if (m->row != sim->matdim || m->col != sim->matdim ) {
        fprintf(stderr,"error: 'matrix set start': matrix %d must be a %dx%d matrix\n",num,sim->matdim,sim->matdim);
        exit(1);
      }
      if (wsp->fstart[0] != sim->fstart[0]) free_complx_matrix(wsp->fstart[0]);
      wsp->fstart[0] = m;
    } else if (!strcmp(argv[arg],"detect")) {
      if (m->row != sim->matdim  || m->col != sim->matdim ) {
        fprintf(stderr,"error: 'matrix set detect': matrix %d must be a %dx%d matrix\n",num,sim->matdim,sim->matdim);
        exit(1);
      }
      if (wsp->fdetect[0] != sim->fdetect[0]) free_complx_matrix(wsp->fdetect[0]);
      wsp->fdetect[0] = m;
    } else {
      return TclError(interp,"error: first argument to 'matrix set' must be <number>, 'start' or 'detect', "
                           "but not '%s'\n",argv[arg]);
    }
  } else
    return TclError(interp,"error: first argument to 'matrix' must be 'get' or 'set' "
                           "but not '%s'\n",argv[arg]);

  return TCL_OK;
}

int tclTurnon(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  if (argc < 2)
    return TclError(interp,"Usage: turnon <int> <int> ...");

  for (i=1;i<argc;i++) {
    ham_turnon(sim,wsp,argv[i]);
  }
  return TCL_OK;
}

int tclTurnoff(ClientData data,Tcl_Interp* interp,int argc, char *argv[])
{
  int i;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  
  if (argc < 2)
    return TclError(interp,"Usage: turnoff <int> <int> ...");

  for (i=1;i<argc;i++) {
    ham_turnoff(sim,wsp,argv[i]);
  }
  return TCL_OK;
}

int tclRotorangle(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double brl,gamma;
  char *buf;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  if (argc != 2 && argc != 3)
    return TclError(interp,"Usage: rotorangle (<angle>|reset|ma) ?gamma_add?");

  buf = Tcl_GetString(argv[1]);
  if (!strcmp(buf, "reset")) {
    wsp->brl = sim->brl;
  } else if (!strcmp(buf, "ma")) {
    wsp->brl = MAGIC_ANGLE;
  } else {
    if (Tcl_GetDoubleFromObj(interp,argv[1],&brl) != TCL_OK) return TCL_ERROR;
    wsp->brl = brl;
  }
  if (argc == 3) {
    if (Tcl_GetDoubleFromObj(interp,argv[2],&gamma) != TCL_OK) return TCL_ERROR;
    wsp->gamma_add = gamma;
  }
  return TCL_OK;
}

int tclGetInteractions(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  
  if (argc != 1)
    return TclError(interp,"Usage: {list} getinteractions");

  Tcl_ResetResult(interp);

  for (i=0;i<sim->nCS;i++)
	  TclAppendResult(interp,"shift_%d %d",sim->CS[i]->nuc, wsp->CS_used[i]);
  for (i=0;i<sim->nJ;i++)
	  TclAppendResult(interp,"jcoupling_%d_%d %d",sim->J[i]->nuc[0],sim->J[i]->nuc[1], wsp->J_used[i]);
  for (i=0;i<sim->nDD;i++)
	  TclAppendResult(interp,"dipole_%d_%d %d",sim->DD[i]->nuc[0],sim->DD[i]->nuc[1], wsp->DD_used[i]);
  for (i=0;i<sim->nQ;i++)
	  TclAppendResult(interp,"quadrupole_%d %d",sim->Q[i]->nuc, wsp->Q_used[i]);
  
  return TCL_OK;
}


int tclCurrenttime(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  
  return TclSetResult(interp,"%g",wsp->t);
}


/****
 * ZT: Zdenek Tosner - implementation of avgham_static command
 ****/
int tclAvgHam_static(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double duration;
  char* expr;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"avgham_static");

  if (argc != 3)
    return TclError(interp,"Usage: avgham_static <duration/usec> <expr>");

  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;

  expr=Tcl_GetString(argv[2]);
  
  if (duration < 0.0)
      return TclError(interp,"avgham_static: duration must be zero or positive");

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  return TCL_OK;
  }

  if (sim->block_diag) {
  	  if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
  	  if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
  }

  if (duration != 0.0) 
     _avgham_static(sim,wsp,duration,expr);
     
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 1;
  }
  
  return TCL_OK;
}

/******
 * TV: Routine to return Euler Angles
 ******/
int tclEulerAngles(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  char buf[256];
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  if (argc != 1)
    return TclError(interp,"Usage: {alpha beta gamma weight} eulerangles");
  
  Tcl_ResetResult(interp);

  sprintf(buf, "%g", wsp->cryst.alpha);
  Tcl_AppendElement(interp, buf);
  sprintf(buf, "%g", wsp->cryst.beta);
  Tcl_AppendElement(interp, buf);
  sprintf(buf, "%g", wsp->cryst.gamma);
  Tcl_AppendElement(interp, buf);
  sprintf(buf, "%g", wsp->cryst.weight);
  Tcl_AppendElement(interp, buf);
  return TCL_OK;
}

/****
 * ZT: implementation of Tcl zgrad_pulse_shape command 
 ****/
int tclZgrPulse(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  double duration, zcoor;
  int slot, Nelem, i, Nchan;
  double *ovals;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"zgrad_pulse_shape");

  if (argc != 3)
      return TclError(interp,"Usage: zgrad_pulse_shape <duration/usec> <zgrad shape>");

  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
      return TclError(interp,"zgrad_pulse_shape: cannot convert argument 1 to duration");
  if (duration < 0.0)
      return TclError(interp,"zgrad_pulse_shape: duration must be zero or positive");
  if (Tcl_GetIntFromObj(interp,argv[2],&slot) != TCL_OK)
      return TclError(interp,"zgrad_pulse_shape: cannot convert argument 2 to zgrad_shape");
  if (!ZgradShapes[slot])
      return TclError(interp,"zgrad_pulse_shape: zgrad_shape was not activated/does not exist");

  if (sim->block_diag) {
  	  change_basis_delay(sim, wsp);
  	  if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
  	  if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
  	  //if (wsp->tmpU == NULL) wsp->tmpU = create_blk_mat_complx_copy2(wsp->ham_blk);
  }

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  return TCL_OK;
  }

  if (duration != 0.0) {
    Nelem = ZgradShapes_len(slot);
    duration /= (double)Nelem;
    Nchan = sim->ss->nchan;
    ovals = double_vector(Nchan);
    wsp->inhom_offset = ovals;
    zcoor = wsp->zcoor;
    for (i=1; i<=Nelem; i++) {
       get_chan_offset_ratios(sim->ss, 2.0*M_PI*ZgradShapes[slot][i], ovals);
       dv_muld(ovals,zcoor);
       /* ZT: relaxation? */
       if (sim->relax) {
          _delay_relax(sim,wsp,duration);
       } else {
          _delay(sim,wsp,duration);
       }
    }
    free_double_vector(ovals);
    wsp->inhom_offset = NULL;
    /* if used within pulseq optimizing propagator set this flag */
    if (OCpar.gradmodeprop) {
       wsp->OC_propstatus = 1;
    }
  }  
  return TCL_OK;
}

/****
 * ZT: helper function for normal pulse_and_zgrad_shaped
 ****/
 void _pulse_and_zgrad_shaped(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *mask, int zgrslot, double steptime)
{
  int i, j, Nchan;
  double *ovals, zcoor;

  /* prepare for offset term */
  Nchan = sim->ss->nchan;
  ovals = double_vector(Nchan);
  zcoor = wsp->zcoor;
  
  for (j=1; j<=Nelem; j++) {
     /* offset term */
	 get_chan_offset_ratios(sim->ss, 2.0*M_PI*ZgradShapes[zgrslot][j], ovals);
     dv_muld(ovals,zcoor);
     wsp->inhom_offset = ovals;
     for (i=1; i<=Nchan; i++) {
        /* prepare rf term */
        if (mask[i] == -1) {
	     _rf(wsp,i,0.0);
	     _ph(wsp,i,0.0);
        } else {
	     _rf(wsp,i,RFshapes[mask[i]][j].ampl);
	     _ph(wsp,i,RFshapes[mask[i]][j].phase);
        }
     }
     /* do step pulse */
     _pulse(sim, wsp, steptime);
  }
  /* clean up after z grad offsets */
  free_double_vector(ovals);
  wsp->inhom_offset = NULL;
  /* if used within pulseq optimizing propagator set this flag */
  if (OCpar.gradmodeprop) {
     wsp->OC_propstatus = 1;
  }
  
}

/****
 * ZT: implementation of shaped pulse with z gradient on
 *      - using global variables RFshapes[], ZgradShapes[]
 ****/
int tclPulseZgrShaped(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i, j, slot, zgrslot, Nelem, Nch=0;
  double duration, steptime;
  int *mask, basis = 0;
  char buf[256], cd[128], buf2[32];
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"pulse_and_zgrad_shaped");

  /* disable when relaxation is ON */
  if (sim->relax) {
     fprintf(stderr,"pulse_and_zgrad_shaped error: not supported in combination with relaxation.\n");
     exit(1);
  }
  
  mask = int_vector(sim->ss->nchan);

  /* READING ARGUMENTS  */
  if (argc <= 3)
    return TclError(interp,"Usage: pulse_and_zgrad_shaped <duration/usec> <zgrad shape> <RFshape on channel 1> ?<RFshape on channel 2>? ...");
  if (argc-3 != sim->ss->nchan)
    return TclError(interp,"\nerror: pulse_and_zgrad_shaped: arguments must match number of channels\n"
                            "        use 'nothing' as place holder when no rf is applied\n");
  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"pulse_and_zgrad_shaped: duration must be zero or positive");
  if (Tcl_GetIntFromObj(interp,argv[2],&zgrslot) != TCL_OK)
     return TCL_ERROR;
  if (!ZgradShapes[zgrslot]) 
     return TclError(interp,"pulse_and_zgrad_shaped: second argument points to non-existing zgrad shape");

  Nelem = ZgradShapes_len(zgrslot);
      
  for (i=1; i<=sim->ss->nchan; i++) {
	  if (!strcmp(Tcl_GetString(argv[i+2]),"nothing")) {
		  mask[i] = -1;
		  basis += 1 << (i-1);
		  /* Debug output */
		  /*
		printf("Channel %d: no rf applied\n",i);
		   */
		  continue;
	  } else {
		  /* read RFshape and check it */
		  if (Tcl_GetIntFromObj(interp,argv[i+2],&slot) != TCL_OK) {
			  sprintf(buf,"error in pulse_and_zgrad_shaped: argument %d must be interger <RFshape>",i+2);
			  return TclError(interp, buf);
		  }
		  if (!RFshapes[slot]) {
			  sprintf(buf,"pulse_and_zgrad_shaped: argument %d points to non-existing RFshape",i+2);
			  return TclError(interp,buf);
		  }
		  /* check consistency of RFshape lengths */
		  if ( RFshapes_len(slot) != Nelem )
			  return TclError(interp,"pulse_and_zgrad_shaped: inconsistend number of elements in RFshapes and zgrad shape!");
		  /* Debug output */
		  /*
		  printf("Channel %d: Number of elements in the RFshape = %d\n",i-1,Nelem);
		  */
		  mask[i] = slot;
	  }
  }
   
  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  free_int_vector(mask);
	  return TCL_OK;
  }

   steptime = duration/Nelem;
   
   /* setting propagator flag is done only within _pulse_shaped function */

   if (sim->block_diag) {
	   change_basis_pulse(sim, wsp, basis);
	   if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
	   if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
	   //if (wsp->tmpU == NULL) wsp->tmpU = create_blk_mat_complx_copy2(wsp->ham_blk);
   }

   if (OCpar.gradmode) {
      /* printf("pulse_and_zgrad_shaped does gradients\n"); */
      /* determine which channels are active for gradients */
      if (!OCpar.grad_shapes) 
         return TclError(interp,"error when calculating propagators for gradients: grad_shapes not defined");
      int Nsh=LEN(OCpar.grad_shapes);      
      cd[0]='\0';    
      for (i=1; i<=sim->ss->nchan; i++) {
         for (j=1; j<=Nsh; j++) {
	    if (OCpar.grad_shapes[j] == mask[i]) {
	       sprintf(buf2," I%dC%d",j,i);
	       strcat(cd,buf2);
	       Nch++;
	       break;
	     }
	  }
      }  
      /* check if there is any gradient to claculate */
      if (Nch==0) {
         /* no, then do usual things */
	 /* printf("pulse_shaped - no variable shapes, just creates propagator\n"); */
	 _pulse_and_zgrad_shaped(sim, wsp, Nelem, mask, zgrslot, steptime);
      } else {
         /* yes, check for type of optimization */
	 /* printf("pulse_shaped created code '%s'\n",cd); */
	 if ( OCpar.gradmodeprop == 1 ) {
	    /* do stuff for propagator optimization */
	    _pulse_and_zgrad_shapedOCprops(sim, wsp, cd,Nch,Nelem,mask,zgrslot,steptime);
	 } else {
	    /* do stuff for state to state optimization */
   	    _pulse_and_zgrad_shapedOC(sim, wsp, cd,Nch,Nelem,mask,zgrslot,steptime);
	 }
      }

   } else {
      /* do just actual pulsing */
      _pulse_and_zgrad_shaped(sim, wsp, Nelem, mask, zgrslot, steptime);
   }
   
  free_int_vector(mask);
    
  return TCL_OK;
}


int tclAcqBlock(ClientData data,Tcl_Interp* interp,int objc, Tcl_Obj *objv[])
{
	int i, np=-1;
	double phase=0.0;
	char *s;
	Sim_info *sim = NULL;
	Sim_wsp *wsp = NULL;

	read_sim_pointers(interp, &sim, &wsp);

	if ( objc%2 != 0 )
		return TclError(interp,"acq_block argument count error. usage: \n\t"
								"acq_block -np X -ph x { pulse code }");
	for (i=1; i<objc-1; i++) {
		s = Tcl_GetString(objv[i]);
		/* printf("   analyzing '%s'\n",s); */
		i++;
		if (!strncmp(s,"-np",3)) {
			if ( Tcl_GetIntFromObj(interp,objv[i],&np) != TCL_OK )
				return TclError(interp,"acq_block error: can not get int for -np parameter");
		} else if (!strncmp(s,"-ph",3)) {
			if ( Tcl_GetDoubleFromObj(interp,objv[i],&phase) != TCL_OK ) {
				/* perhaps not a number... */
				s = Tcl_GetString(objv[i]);
				if (!strcasecmp(s,"-X")) {
					phase=180;
				} else if (!strcasecmp(s,"-Y")) {
					phase=270;
				} else if (!strcasecmp(s,"X")) {
					phase=0;
				} else if (!strcasecmp(s,"Y")) {
					phase=90;
				} else {
					return TclError(interp,"acq_block error: can not get phase for -ph parameter");
				}
			}
		} else {
			return TclError(interp,"acq_block error: unknown option %s\n",s);
		}
	}

	if (np<0) np = sim->np;
//	if (wsp->curr_nsig + np > LEN(wsp->fid))
	if (wsp->curr_nsig + np > sim->ntot)
		return TclError(interp,"Error: acq_block - acq overflow in fid points\n");

	wsp->Nacq = np;
	wsp->acqphase = phase;

	switch (sim->imethod) {
	case M_DIRECT_TIME:
		if (sim->interpolation == 1) {
			direct_acqblock_time_FWT(interp,objv[objc-1],sim,wsp);
		} else {
			direct_acqblock(interp,objv[objc-1],sim,wsp);
		}
		break;
	case M_DIRECT_FREQ:
		if (sim->interpolation == INTERPOL_NOT_USED) {
			direct_acqblock_freq(interp,objv[objc-1],sim,wsp);
		} else if (sim->interpolation == INTERPOL_ASG) {
			direct_acqblock_ASGdata(interp,objv[objc-1],sim,wsp);
		} else { // any variant of FWT
			direct_acqblock_FWTdata(interp,objv[objc-1],sim,wsp);
		}
		break;
	case M_GCOMPUTE_TIME:
	case M_GCOMPUTE_FREQ:
		gcompute_acqblock(interp,objv[objc-1],sim,wsp);
		break;
	case M_SPINACH:
		wsp->evalmode = EM_ACQBLOCK;
		spinach_acqblock(interp,objv[objc-1],sim,wsp);
		break;
	default:
		return TclError(interp,"Error: acq_block has not recognized the current method\n");
	}

	wsp->evalmode = EM_NORMAL;
	return TCL_OK;
}

#include <time.h>
int tclTestFunc(ClientData data,Tcl_Interp* interp,int objc, Tcl_Obj *objv[])
{
	Sim_info *sim = NULL;
	Sim_wsp *wsp = NULL;
	int i, j, n, nnz;
	//LARGE_INTEGER tv1, tv2, tickpsec;
	extern void prop_pade_real(mat_complx *prop, mat_double *ham,double dt);
	extern void prop_cheby_real(mat_complx *prop, mat_double *ham,double dt);
	extern void prop_diag1_real(mat_complx *prop, mat_double *ham, double dt);
	extern void prop_diag2_real(mat_complx *prop, mat_double *ham, double dt);
	// Taylor pres prop_real(), method = 3
	//QueryPerformanceFrequency(&tickpsec);

	//if (Tcl_GetIntFromObj(interp,objv[1],&n) != TCL_OK)
	//     return TCL_ERROR;
	//if (Tcl_GetDoubleFromObj(interp,objv[2],&dens) != TCL_OK)
	//     return TCL_ERROR;

	read_sim_pointers(interp, &sim, &wsp);

	const double nrm[] = {0.5, 1.2, 2.5, 5, 7.5, 10};
	//const int Nrep = 6;
	int Nrep = 1;
	const double dens[] = {0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0};
	const int Ndens = 7;
	int idens;

	for (n=3; n<11; n++) {
	for (idens=0; idens<Ndens; idens++) {
	n=4;
	idens = 2;
		fprintf(stderr,"Doing n = %d, dens = %g\n",n,dens[idens]);

	int dim = 2<<(n-1);
	nnz = (int)(dens[idens]*dim*dim);
	if (nnz < 3) continue;
	assert(nnz <= dim*dim);
	int ond = nnz / dim;
	int offd = ( nnz * (dim-1) / dim ) / 2;
	i = nnz - ond - 2*offd;
	if (i%2) {
		ond++;
		i--;
	}
	offd += i/2;
	assert(ond<=dim);
	assert(offd<=dim*dim-dim);
	assert(nnz-ond-2*offd==0);
	mat_double *ham = double_matrix(dim,dim,MAT_DENSE,0,0);
	dm_zero(ham);
	mat_complx *prop = cm_creatediag('d',dim,Complx(1,0),0);
	mat_double *hh;

	printf("ond=%d, offd=%d\n",ond,offd);
	srand(time(NULL));
	if (nnz < dim*dim) {
		int step;
		if (ond == 0) {
			step = 0;
		} else {
			step = dim/ond;
		}
		for (i=0;i<ond;i++) {
			int r = i*step % dim;
			//printf("diag: %d\n",r);
			ham->data[r+r*dim] = rand();
		}
		if (offd == 0) {
			step = 0;
		} else {
			step = (dim*dim-dim)/2 / offd;
		}
		for (i=0; i<offd; i++) {
			int r = 0;
			int c = (i*step) % ( (dim*dim-dim)/2 );
			//printf("\tc=%d   ",c);
			while ( c - (dim-r-1) > 0 ) {
				c -= dim-r-1;
				r++;
			}
			c += r;
			//printf("OFF-diag: r=%d, c=%d\n",r,c);
			double dum = rand();
			ham->data[r+c*dim] = dum;
			ham->data[c+r*dim] = dum;
		}
		/*
		while (i < nnz) {
			int idx = (int)(rand()/(double)RAND_MAX*dim*dim);
			int r = idx / dim;
			int c = idx % dim;
			if (fabs(ham->data[r+c*dim]) > 0) continue;
			double dum = rand();
			if (r != c) {
				ham->data[r+c*dim] = dum;
				i++;
				ham->data[c+r*dim] = dum;
				i++;
			} else {
				ham->data[r+c*dim] = dum;
				i++;
			}
		}
		*/
	} else {
		for (i=0; i<dim*dim; i++) ham->data[i] = rand();
	}
	nnz = dm_nnz(ham);
	printf("\nGenerated matrix dim = %d, nnz = %d\n",dim,nnz);
	dm_print(ham,"HAM");
	void slanpro(mat_double *A, double *dg0, double *dg1, mat_double *Q);
	double *dg0 = (double*)malloc(dim*sizeof(double));
	double *dg1 = (double*)malloc((dim-1)*sizeof(double));
	mat_double *Q = double_matrix(dim,dim,MAT_DENSE,0,0);
	slanpro(ham,dg0,dg1,Q);
	printf("diagonal:\n");
	for (j=0; j<dim; j++) printf("%g ",dg0[j]);
	printf("\nsubdiagonal:\n");
	for (j=0; j<dim-1; j++) printf("%g ",dg1[j]);
	printf("\n");
	dm_print(Q,"transf. matrix Q");
	exit(1);
	printf("Norm\tDIAG_1\tDIAG_2\tPade_D\tTaylor_D\tCheby_D\tPade_S\tTaylor_S\tCheby_S\n");

	for (j=0; j<Nrep; j++) {
		double norm = dm_normest(ham);
		dm_muld(ham,nrm[j]/norm);
		norm = dm_normest(ham);
		printf("%.9f\t",norm);
		hh = dm_dup(ham);
		dm_dense(hh);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_diag1_real(prop,hh,1.0);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		hh = dm_dup(ham);
		dm_dense(hh);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_diag2_real(prop,hh,1.0);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		hh = dm_dup(ham);
		dm_dense(hh);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_pade_real(prop,hh,1.0);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		hh = dm_dup(ham);
		dm_dense(hh);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_real(prop,hh,1.0,3);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		hh = dm_dup(ham);
		dm_dense(hh);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_cheby_real(prop,hh,1.0);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		if (nnz > floor(dim*dim*(1.0-SPARSITY)) ) {
			printf("NaN\tNaN\tNaN\t\n");
			continue;
		}
		hh = dm_dup(ham);
		dm_sparse(hh,TINY);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_pade_real(prop,hh,1.0);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		hh = dm_dup(ham);
		dm_sparse(hh,TINY);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_real(prop,hh,1.0,3);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		hh = dm_dup(ham);
		dm_sparse(hh,TINY);
		cm_unit(prop);
		//QueryPerformanceCounter(&tv1);
		prop_cheby_real(prop,hh,1.0);
		//QueryPerformanceCounter(&tv2);
		//printf("%.9f\t",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		free_double_matrix(hh);

		printf("\n");
	}

	free_double_matrix(ham);
	free_complx_matrix(prop);
	} // for over idens
	} // end for over n

	return TCL_OK;
}
/* ZT: single purpose shaped pulse for testing 3-point rule of propagation
 *     Assumes shaped pulse elements are one maxdt long and constant
 *     Evaluates interaction Hamiltonian in left, middle, and right points
 *     Accepting possible trouble arising from discontinuity of rf field between shape elements
 */
void local_prop_complx(mat_complx *prop, mat_complx *ham, double dt)
{
	int i;
	int dim = ham->row;
	double *eigs=double_vector(dim);
	mat_complx *T=complx_matrix(dim,dim,MAT_DENSE,0,ham->basis);

	cm_muld(ham,dt);
	//cm_print(ham,"------> H*dt");
	cm_diag_hermit(ham,eigs,T);
	//for (i=1; i<=dim; i++) printf("   eig[%i] = (%g,%g)\n",i,eigs[i].re,eigs[i].im);
	//cm_print(T,"====> T");
	if (prop->type != MAT_DENSE_DIAG) {
		mat_complx *dum = complx_matrix(dim,dim,MAT_DENSE_DIAG, dim, ham->basis);
		cm_swap_innards_and_destroy(prop,dum);
	}
	prop->basis = ham->basis;
	for (i=0; i<dim; i++) {
		prop->data[i] = Cexpi(-eigs[i+1]);
	}
	//cm_print(prop,"+++++>prop");
	simtrans(prop,T);
	free_double_vector(eigs);
	free_complx_matrix(T);
	//cm_print(prop,"~~~~~> prop");
}
void _pulse_simple_3pointruleconstantrf(Sim_info *sim, Sim_wsp *wsp, double duration)
{
	/* this assumes that wsp->sumHrf and wsp->sumUph are already made */
	int i, n;
	double dt, dt_us;
	blk_mat_double *Hleft, *Hmiddle, *Hright, *Htemp, *temp_ham_blk;

	// initialization and allocation
	assert(wsp->dU != NULL);
	assert(wsp->dU->Nblocks == 1); // assume no blockdiagonalization in the current implementation
	assert(wsp->ham_blk->Nblocks == 1); // assume no blockdiagonalization in the current implementation
	assert(wsp->sumHrf->Nblocks == 1); // assume no blockdiagonalization in the current implementation
	blk_cm_unit(wsp->dU);
	Hleft = create_blk_mat_double_copy(wsp->ham_blk);
	Hmiddle = create_blk_mat_double_copy(wsp->ham_blk);
	Hright = create_blk_mat_double_copy(wsp->ham_blk);
	if (wsp->Hcplx == NULL) {
		wsp->Hcplx = create_blk_mat_complx_copy2(wsp->ham_blk);
	}

	// split duration in maxdt elements
	n = (int)ceil(duration/wsp->dtmax);
	if (n < 1) n = 1; // use 3point rule also for short durations (<maxdt)
	dt_us = duration/(double)n;
	dt = dt_us*1.0e-6;
	//printf("_pulse_simple_3pointruleconstantrf duration of %f us split into %d steps of %f us\n",duration,n,dt*1.0e+6);
	// create LEFT hamiltonian
	temp_ham_blk = wsp->ham_blk;
	wsp->ham_blk = Hleft;
	ham_hamilton(sim,wsp);
	blk_dm_multod(Hleft,wsp->sumHrf,1.0);
	for (i=1;i<=n;i++) {
		// create MIDDLE hamiltonian
		wsp->t += dt_us/2.0;
		//Htemp = wsp->ham_blk;
		wsp->ham_blk = Hmiddle;
		ham_hamilton(sim,wsp);
		blk_dm_multod(Hmiddle,wsp->sumHrf,1.0);
		// create RIGHT hamiltonian
		wsp->t += dt_us/2.0;
		wsp->ham_blk = Hright;
		ham_hamilton(sim,wsp);
		blk_dm_multod(Hright,wsp->sumHrf,1.0);
		// create exponent for propagation
		blk_cm_zero(wsp->Hcplx);
		cm_multocr(wsp->Hcplx->m, Hleft->m, Complx(1.0/6.0,0));
		cm_multocr(wsp->Hcplx->m, Hmiddle->m, Complx(4.0/6.0,0));
		cm_multocr(wsp->Hcplx->m, Hright->m, Complx(1.0/6.0,0));
		mat_double *dum = dm_commutator(Hleft->m,Hright->m);
		cm_multocr(wsp->Hcplx->m, dum, Complx(0.0,1.0/12.0*dt));
		free_double_matrix(dum);
		// calculate propagator
		//blk_cm_print(wsp->Hcplx,"Hamiltonian");
		//blk_cm_print(wsp->dU,"propagator");
		//blk_prop_complx_3(wsp->dU,wsp->Hcplx,dt,sim);
			mat_complx *cum = complx_matrix(wsp->Hcplx->blk_dims[0],wsp->Hcplx->blk_dims[0],MAT_DENSE_DIAG, 0, wsp->Hcplx->basis);
			local_prop_complx(cum,wsp->Hcplx->m,dt);
			cm_multo_rev(wsp->dU->m,cum);
			free_complx_matrix(cum);
		//blk_cm_print(wsp->dU,"updated propagator");
		//printf("time step is %g\n",dt);
		//exit(1);
		// swap Hleft and Hright
		//wsp->ham_blk = Htemp;
		Htemp = Hleft;
		Hleft = Hright;
		Hright = Htemp;
	}
	// add phase of the pulse element
	blk_simtrans_zrot2(wsp->dU,wsp->sumUph);
	// put ham_blk back
	wsp->ham_blk = temp_ham_blk;
	// cleaning memory
	free_blk_mat_double(Hleft);
	free_blk_mat_double(Hmiddle);
	free_blk_mat_double(Hright);
}
void _pulse_shaped_3pointruleconstantrf(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *mask, double duration)
{
	int i, j, ii, n;
	double dt, t_tmp;

	for (j=1; j<=Nelem; j++) {
		// this prepares RF Hamiltonian data for using the trick of real interaction Hamiltonian
		//  commuting with Iz operator, so the rf phase can be applied after the propagator
		//  is calculated from a real matrix
		for (i=1; i<=sim->ss->nchan; i++) {
			if (mask[i] == -1) {
				_rf(wsp,i,0.0);
				_ph(wsp,i,0.0);
			} else {
				_rf(wsp,i,RFshapes[mask[i]][j].ampl);
				_ph(wsp,i,RFshapes[mask[i]][j].phase);
			}
		}
		// replace _pulse with new calculations
		//_pulse(sim,wsp,steptime);
		// new calculations start here; just replacing _pulse_simple for _pulse_simple_3pointruleconstantrf
		_setrfprop(sim,wsp); //use this and do not care about possible zero rf amplitudes

		if (wsp->evalmode == EM_ACQBLOCK) {
			if ( fabs(wsp->t - wsp->acqblock_t0) > TINY ) {
				/* there was some previous event not synchronized with dw */
				dt = wsp->dw - (wsp->t - wsp->acqblock_t0);
				if ( duration < dt) {
					_pulse_simple_3pointruleconstantrf(sim,wsp,duration);
					t_tmp = 0;
				} else {
					_pulse_simple_3pointruleconstantrf(sim,wsp,dt);
					t_tmp = duration - dt;
				}
				update_propagator(wsp->U, wsp->dU, sim, wsp);
				/* did it fill time up to dw? */
				if (fabs(wsp->t - wsp->acqblock_t0 - wsp->dw) < TINY*2 ) {
					_store(sim,wsp,wsp->acqblock_sto,0);
					wsp->acqblock_t0 = wsp->t;
					acqblock_sto_incr(wsp);
				}
			} else {
				t_tmp = duration;
			}
			if ( fabs(t_tmp) > TINY) {
				n = (int)floor(t_tmp/wsp->dw+1e-6);
				dt = t_tmp - wsp->dw*(double)n;
				for (ii=1; ii<=n; ii++) {
					_pulse_simple_3pointruleconstantrf(sim,wsp,wsp->dw);
					update_propagator(wsp->U, wsp->dU, sim, wsp);
					_store(sim,wsp,wsp->acqblock_sto,0);
					wsp->acqblock_t0 = wsp->t;
					acqblock_sto_incr(wsp);
				}
				if (dt > TINY) {
					/* remaining time not synchronized with dw */
					_pulse_simple_3pointruleconstantrf(sim,wsp,dt);
					update_propagator(wsp->U, wsp->dU, sim, wsp);
				}
			}
		} else {
			_pulse_simple_3pointruleconstantrf(sim,wsp,duration);
			update_propagator(wsp->U, wsp->dU, sim, wsp);
		}
		wsp->waspulse += 1;
		// new calculations end here
	}
	/* if used within pulseq optimizing propagator set this flag */
	if (OCpar.gradmodeprop) {
		wsp->OC_propstatus = 1;
	}
 /*** END of _pulse***/
}
int tclPulseShaped3PointRuleConstantRF(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
  int i, j, slot, Nelem=-1, Nch=0, basis=0;
  double duration, steptime;
  int *mask, *OCchanmap = NULL;
  char cd[128], buf2[64];;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;
  const char *cmdname = "pulse_shaped_3pointruleconstantrf";

  read_sim_pointers(interp, &sim, &wsp);
  spinach_disable_command(sim,"pulse_shaped");

  /* disable when relaxation is ON */
  if (sim->relax) {
	  return TclError(interp,"%s error: not supported in combination with relaxation.\n",cmdname);
  }
  if (sim->frame != ROTFRAME) {
	  return TclError(interp,"%s error: works only in ROTFRAME.\n",cmdname);
  }
  if (sim->imethod != M_DIRECT_TIME) {
	  return TclError(interp,"%s error: works only with direct method in time domain.\n",cmdname);
  }

  mask = int_vector(sim->ss->nchan);
  /* READING ARGUMENTS  */
  if (argc <= 2)
    return TclError(interp,"Usage: %s <duration/usec> <RFshape on channel 1> ?<RFshape on channel 2>? ...",cmdname);
  if (argc-2 != sim->ss->nchan)
    return TclError(interp,"\nerror: %s: arguments must match number of channels\n"
                            "        use 'nothing' as place holder when no rf is applied\n",cmdname);
  if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
     return TCL_ERROR;
  if (duration < 0.0)
      return TclError(interp,"%s: duration must be zero or positive",cmdname);

  for (i=1; i<=sim->ss->nchan; i++) {
     if (!strcmp(Tcl_GetString(argv[i+1]),"nothing")) {
        mask[i] = -1;
        basis += 1 << (i-1); // use channels with no rf for blocking
	    DEBUGPRINT("Channel %d: no rf applied\n",i);
	    continue;
     } else {
        /* read RFshape and check it */
        if (Tcl_GetIntFromObj(interp,argv[i+1],&slot) != TCL_OK) {
           return TclError(interp,"error in %s: argument %d must be integer <RFshape>",cmdname,i+1 );
        }
	    if (!RFshapes[slot]) {
           return TclError(interp,"%s: argument %d points to non-existing RFshape",cmdname,i+1);
	    }
	    if (Nelem == -1) {
	    	/* set variable for length of RFshapes */
	    	Nelem=RFshapes_len(slot);
	    } else {
	    	/* check consistency of RFshape lengths */
	    	if ( RFshapes_len(slot) != Nelem )
	    		return TclError(interp,"%s: inconsistent number of elements in RFshapes!",cmdname);
	    }
        DEBUGPRINT("Channel %d: Number of elements in the RFshape = %d\n",i-1,Nelem);
        mask[i] = slot;
     }
   }

  if (sim->Nmz > sim->ss->nchan) {
	  // use dummy channels for blocking
	  for (i=sim->ss->nchan; i<sim->Nmz; i++) basis += 1 << i;
  }

  if (wsp->evalmode == EM_MEASURE) {
	  wsp->t += duration;
	  free_int_vector(mask);
	  return TCL_OK;
  }

   steptime = duration/Nelem;

   if (sim->block_diag) {
	   change_basis_pulse(sim, wsp, basis);
	   if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
	   if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
   }

   if (OCpar.gradmode) {
      DEBUGPRINT("%s does gradients\n",cmdname);
      /* determine which channels are active for gradients */
      if (!OCpar.grad_shapes)
         return TclError(interp,"error when calculating propagators for gradients: grad_shapes not defined");
      int Nsh=LEN(OCpar.grad_shapes);
      if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
    	  cd[0]='\0';
    	  for (i=1; i<=sim->ss->nchan; i++) {
    		  for (j=1; j<=Nsh; j++) {
    			  if (OCpar.grad_shapes[j] == mask[i]) {
    				  sprintf(buf2," I%dC%d",j,i);
    				  strcat(cd,buf2);
    				  Nch++;
    				  break;
    			  }
    		  }
    	  }
    	  /* check if there is any gradient to calculate */
    	  if (Nch==0) {
    		  /* no, then do usual things */
    		  DEBUGPRINT("%s - no variable shapes, just creates propagator\n",cmdname);
    		  _pulse_shaped_3pointruleconstantrf(sim,wsp,Nelem, mask, steptime);
    	  } else {
    		  /* yes, check for type of optimization */
    		  DEBUGPRINT("%s created code '%s'\n",cmdname,cd);
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  return TclError(interp,"%s: OC of prop ala Khaneja NOT implemented!",cmdname);
    			  //_pulse_shapedOCprops(sim,wsp,cd,Nch,Nelem,mask,steptime);
    		  } else {
    			  /* do stuff for state to state optimization */
    			  return TclError(interp,"%s: OC state-state ala Khaneja NOT implemented!",cmdname);
    			  //_pulse_shapedOC(sim,wsp,cd,Nch,Nelem,mask,steptime);
    		  }
    	  }
      } else { /* GRAPE with advanced gradients */
    	  OCchanmap = int_vector(Nsh);
    	  for (i=1; i<=Nsh; i++) {
    		  OCchanmap[i] = -1;
    		  //printf("grad_shapes[%d] = %d \n",i,OCpar.grad_shapes[i]);
    		  for (j=1; j<=sim->ss->nchan; j++) {
        		  //printf("\t mask[%d] = %d \n",j,mask[j]);
    			  if (OCpar.grad_shapes[i] == mask[j]) {
    				  OCchanmap[i] = j;
    				  Nch++;
        			  break;
    			  }
    		  }
    		  //printf("\t OCchanmap[%d] = %d \n",i,OCchanmap[i]);
    	  }
    	  if (Nch == 0) {
    		  _pulse_shaped_3pointruleconstantrf(sim,wsp,Nelem, mask, steptime);
    	  } else {
    		  if ( OCpar.gradmodeprop == 1 ) {
    			  /* do stuff for propagator optimization */
    			  return TclError(interp,"%s: OC of prop with advanced grads NOT implemented!",cmdname);
    			  /***
    			  switch (sim->frame) {
    			  case ROTFRAME:
    				  _pulse_shapedOCprops_2(sim,wsp,Nelem,OCchanmap,mask,steptime);
    				  break;
    			  case DNPFRAME:
    			  case LABFRAME:
    				  fprintf(stderr,"Error: optimal control on propagators not supported in DNP/LAB frames (checkpoint tclPulseShaped)\n");
    				  exit(1);
    				  break;
    			  default:
    				  fprintf(stderr,"Error: unknown calculation frame %d (checkpoint tclPulseShaped 1)\n",sim->frame);
    				  exit(1);
    			  }
    			  ***/
    		  } else {
    			  /* do stuff for state to state optimization */
    			  return TclError(interp,"%s: OC of state-state with advanced grads NOT implemented!",cmdname);
    			  /***
    			  switch (sim->frame) {
    			  case ROTFRAME:
    				  _pulse_shapedOC_2(sim,wsp,Nelem,OCchanmap,mask,steptime);
    				  break;
    			  case DNPFRAME:
    				  _pulse_shapedOC_2_dnpframe(sim,wsp,Nelem,OCchanmap,mask,steptime);
    				  break;
    			  case LABFRAME:
    				  fprintf(stderr,"Error: optimal control not supported in LABframe (checkpoint tclPulseShaped)\n");
    				  exit(1);
    				  break;
    			  default:
    				  fprintf(stderr,"Error: unknown calculation frame %d (checkpoint tclPulseShaped 1)\n",sim->frame);
    				  exit(1);
    			  }
    			  ***/
    		  }
    	  }
    	  free_int_vector(OCchanmap);
      }
   } else {
      /* do just actual pulsing */
      _pulse_shaped_3pointruleconstantrf(sim,wsp,Nelem, mask, steptime);
   }

  free_int_vector(mask);

  return TCL_OK;
}


/*****
 * Implementation of special purpose pulse_shaped_3pointrulerfmap command
 * Intention: test 3-point rule propagation over one single element of the
 * shaped pulse. And, optimal control.
 * Restrictions: when rf shape has N elements in one rotor period,
 *   rfmap MUST have 2*N elements over one rotor period
 * RF parameters are assumed constant over those 3 points, accepting discontinuity
 *   between pulse elements.
 */
void local_setrfprop_complex(Sim_info *sim, Sim_wsp *wsp, blk_mat_complx *HRF)
{  // creates RF Hamiltonian in wsp->Hrf_blk, ignores selective pulses
	int i;

	assert(HRF->Nblocks == 1); // assume no blockdiagonalization in the current implementation
	blk_cm_zero(HRF);
	dm_zero(wsp->sumUph);

	for (i=1; i<=sim->ss->nchan; i++) {
		cm_multocr(HRF->m, wsp->chan_Ix[i], Complx(wsp->rfv[i],0));
		dm_multod(wsp->sumUph, wsp->chan_Iz[i], wsp->phv[i]);
	}
	// apply phases as z-rotation
	blk_simtrans_zrot2(HRF,wsp->sumUph);
}
void _pulse_shaped_3pointrulerfmap(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *mask, double duration)
{
	int j, iz, iphi, iphi_shift, chan;
	double steptimephi, bx, by, am, ph;
	blk_mat_complx *Hleft, *Hmiddle, *Hright;
	mat_complx *cmat;
	int Nchan = sim->ss->nchan;
	double dt = duration*1e-6; // element duration in seconds

	if (wsp->evalmode == EM_ACQBLOCK) {
		fprintf(stderr,"pulse_shaped_3pointrulerfmap cannot be used inside acq_block\n");
		exit(1);
	}

	// initialization and allocation
	assert(wsp->dU != NULL);
	assert(wsp->dU->Nblocks == 1); // assume no blockdiagonalization in the current implementation
	assert(wsp->ham_blk->Nblocks == 1); // assume no blockdiagonalization in the current implementation
	blk_cm_unit(wsp->dU);
	Hleft = create_blk_mat_complx_copy2(wsp->ham_blk);
	Hmiddle = create_blk_mat_complx_copy2(wsp->ham_blk);
	Hright = create_blk_mat_complx_copy2(wsp->ham_blk);
	if (wsp->Hcplx == NULL) {
		wsp->Hcplx = create_blk_mat_complx_copy2(wsp->ham_blk);
	}

	/*** rfmap data coordinates and book keeping; check of timings done before! ***/
	iz = sim->rfmap->loop[(wsp->rf_idx)*2+0];   // z coil coordinate index
	iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1]; // coil initial "phase" index
	iphi_shift = 0; // index of rotor "phase", can run to infinity, periodicity assumed in rfmap_get()
	steptimephi = duration*0.5;  // assuming duration = 2*steptimephi !!!

	for (j=1; j<=Nelem; j++) { // main loop over all shape elements
		//printf("Calculating element %d\n",j);
		// preparing LEFT Hamiltonian
		for (chan=1; chan<=Nchan; chan++) {
			if (mask[chan] == -1) { // no rf applied on this channel
				_rf(wsp,chan,0.0);
				_ph(wsp,chan,0.0);
			} else {
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,chan-1,&bx,&by); // get distortion coefs
				am = RFshapes[mask[chan]][j].ampl;  // original ampl/phase element
				ph = RFshapes[mask[chan]][j].phase;
				am *= sqrt(bx*bx+by*by);  // modify rf parameters at this element
				ph += RAD2DEG*atan2(by,bx);
				_rf(wsp,chan,am);  // set them on
				_ph(wsp,chan,ph);
			}
		}
		local_setrfprop_complex(sim, wsp, Hleft); // creates complex RF Hamiltonian in Hleft
		ham_hamilton(sim,wsp);  // creates interaction Hamiltonian in wsp->ham_blk
		blk_cm_multocr(Hleft,wsp->ham_blk,Complx(1.0,0.0)); // add the two, this is final Hleft
		// preparing MIDDLE Hamiltonian
		wsp->t += steptimephi; // move time forward
		iphi_shift++;         // move rotor forward
		for (chan=1; chan<=Nchan; chan++) {
			if (mask[chan] == -1) { // no rf applied on this channel
				_rf(wsp,chan,0.0);
				_ph(wsp,chan,0.0);
			} else {
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,chan-1,&bx,&by); // get distortion coefs
				am = RFshapes[mask[chan]][j].ampl;  // original ampl/phase element
				ph = RFshapes[mask[chan]][j].phase;
				am *= sqrt(bx*bx+by*by);  // modify rf parameters at this element
				ph += RAD2DEG*atan2(by,bx);
				_rf(wsp,chan,am);  // set them on
				_ph(wsp,chan,ph);
			}
		}
		local_setrfprop_complex(sim, wsp, Hmiddle); // creates complex RF Hamiltonian in Hmiddle
		ham_hamilton(sim,wsp);  // creates interaction Hamiltonian in wsp->ham_blk
		blk_cm_multocr(Hmiddle,wsp->ham_blk,Complx(1.0,0.0)); // add the two, this is final Hmiddle
		// preparing RIGHT Hamiltonian
		wsp->t += steptimephi; // move time forward
		iphi_shift++;         // move rotor forward
		for (chan=1; chan<=Nchan; chan++) {
			if (mask[chan] == -1) { // no rf applied on this channel
				_rf(wsp,chan,0.0);
				_ph(wsp,chan,0.0);
			} else {
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,chan-1,&bx,&by); // get distortion coefs
				am = RFshapes[mask[chan]][j].ampl;  // original ampl/phase element
				ph = RFshapes[mask[chan]][j].phase;
				am *= sqrt(bx*bx+by*by);  // modify rf parameters at this element
				ph += RAD2DEG*atan2(by,bx);
				_rf(wsp,chan,am);  // set them on
				_ph(wsp,chan,ph);
			}
		}
		local_setrfprop_complex(sim, wsp, Hright); // creates complex RF Hamiltonian in Hright
		ham_hamilton(sim,wsp);  // creates interaction Hamiltonian in wsp->ham_blk
		blk_cm_multocr(Hright,wsp->ham_blk,Complx(1.0,0.0)); // add the two, this is final Hright
		// combine all Hamiltonians into the complex exponent
		blk_cm_zero(wsp->Hcplx);
		cm_multoc(wsp->Hcplx->m, Hleft->m, Complx(1.0/6.0,0));
		cm_multoc(wsp->Hcplx->m, Hmiddle->m, Complx(4.0/6.0,0));
		cm_multoc(wsp->Hcplx->m, Hright->m, Complx(1.0/6.0,0));
		cmat = cm_commutator(Hleft->m,Hright->m); // NOTE: cmat is aalocated here!!!
		cm_multoc(wsp->Hcplx->m, cmat, Complx(0.0,1.0/12.0*dt));
		// calculate propagator of rf element and update the total propagator
		local_prop_complx(cmat,wsp->Hcplx->m,dt);
		cm_multo_rev(wsp->dU->m,cmat);
		free_complx_matrix(cmat); // NOTE: need to free this as it is not used anymore
		/*** timing for the next shape element starts at RIGHT Hamiltonian, no adjustment necessary ***/
	} // end of loop over shape elements
	update_propagator(wsp->U, wsp->dU, sim, wsp);

	// cleaning memory
	free_blk_mat_complx(Hleft);
	free_blk_mat_complx(Hmiddle);
	free_blk_mat_complx(Hright);
}
int tclPulseShaped3PointRuleRFmap(ClientData data,Tcl_Interp* interp,int argc, Tcl_Obj *argv[])
{
	int i, j, slot, iz, Nphi, Nelem=-1, Nch=0, basis=0;
	double duration, steptime, steptimephi;
	int *mask, *OCchanmap = NULL;
	Sim_info *sim = NULL;
	Sim_wsp *wsp = NULL;
	const char *cmdname = "pulse_shaped_3pointrulerfmap";

	read_sim_pointers(interp, &sim, &wsp);
	spinach_disable_command(sim,cmdname);

	/* disable when relaxation is ON */
	if (sim->relax) {
		return TclError(interp,"%s error: not supported in combination with relaxation.\n",cmdname);
	}
	if (sim->frame != ROTFRAME) {
		return TclError(interp,"%s error: works only in ROTFRAME.\n",cmdname);
	}
	if (sim->imethod != M_DIRECT_TIME) {
		return TclError(interp,"%s error: works only with direct method in time domain.\n",cmdname);
	}

	mask = int_vector(sim->ss->nchan);
	/***** READING ARGUMENTS  start ****/
	if (argc <= 2)
		return TclError(interp,"Usage: %s <duration/usec> <RFshape on channel 1> ?<RFshape on channel 2>? ...",cmdname);
	if (argc-2 != sim->ss->nchan)
		return TclError(interp,"\nerror: %s: arguments must match number of channels\n"
	                            "        use 'nothing' as place holder when no rf is applied\n",cmdname);
	if (Tcl_GetDoubleFromObj(interp,argv[1],&duration) != TCL_OK)
		return TclError(interp,"%s: cannot convert duration to double\n",cmdname);
	if (duration < 0.0)
		return TclError(interp,"%s: duration must be zero or positive",cmdname);
	// read shpaes for all channels
	for (i=1; i<=sim->ss->nchan; i++) {
		if (!strcmp(Tcl_GetString(argv[i+1]),"nothing")) {
			mask[i] = -1;
	        basis += 1 << (i-1); // use channels with no rf for blocking
		    DEBUGPRINT("Channel %d: no rf applied\n",i);
		    continue;
		} else {
			/* read RFshape and check it */
			if (Tcl_GetIntFromObj(interp,argv[i+1],&slot) != TCL_OK) {
				return TclError(interp,"error in %s: argument %d must be integer <RFshape>",cmdname,i+1 );
	        }
		    if (!RFshapes[slot]) {
		    	return TclError(interp,"%s: argument %d points to non-existing RFshape",cmdname,i+1);
		    }
		    if (Nelem == -1) {
		    	/* set variable for length of RFshapes */
		    	Nelem=RFshapes_len(slot);
		    } else {
		    	/* check consistency of RFshape lengths */
		    	if ( RFshapes_len(slot) != Nelem )
		    		return TclError(interp,"%s: inconsistent number of elements in RFshapes!",cmdname);
		    }
		    DEBUGPRINT("Channel %d: Number of elements in the RFshape = %d\n",i-1,Nelem);
		    mask[i] = slot;
		}
	}
	/****** reading parameters DONE here ***/

	if (wsp->evalmode == EM_MEASURE) {
		wsp->t += duration;
		free_int_vector(mask);
		return TCL_OK;
	}

	/* disable when block diagonalization is active */
	if (sim->block_diag) {
		return TclError(interp,"%s: not supported in combination with block diagonalization feature", cmdname);
		/* remove blockdiag related code
		change_basis_pulse(sim, wsp, basis);
		if (wsp->U == NULL) wsp->U = create_blk_mat_complx_copy2(wsp->ham_blk);
		if (wsp->dU == NULL) wsp->dU = create_blk_mat_complx_copy2(wsp->ham_blk);
		END: blockdiag related code */
	}
	/* remove blockdiag related code
	if (sim->Nmz > sim->ss->nchan) {
		// use dummy channels for blocking
		for (i=sim->ss->nchan; i<sim->Nmz; i++) basis += 1 << i;
	}
    END: blockdiag related code */

	/*** check timings of shaped pulse and rfmap ***/
	steptime = duration/Nelem;
	iz = sim->rfmap->loop[(wsp->rf_idx)*2+0];   // z coil coordinate index
	Nphi = sim->rfmap->z[iz*(sim->rfmap->Nch+1)+sim->rfmap->Nch]; // number of phi elements for current z position
	steptimephi = sim->taur/Nphi;
	if (Nphi > 1) { // allow for no modulations with rotation
		if (fabs(steptime - 2*steptimephi) > TINY) {
			return TclError(interp,"%s: bad timing of shape and rfmap. Pulse %g, rfmap %g\n", cmdname, steptime, steptimephi);
		}
	}

	if (OCpar.gradmode) {
		if (!OCpar.grad_shapes) {
			return TclError(interp,"%s: error when calculating propagators for gradients: grad_shapes not defined", cmdname);
		}
		if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
			return TclError(interp,"%s: original GRAPE ala Khaneja not implemented", cmdname);
		} else { /* GRAPE with advanced gradients */
			int Nsh=LEN(OCpar.grad_shapes);
			OCchanmap = int_vector(Nsh);
			for (i=1; i<=Nsh; i++) { // assign rfshapes their respective rf channels
				OCchanmap[i] = -1;
				//printf("grad_shapes[%d] = %d \n",i,OCpar.grad_shapes[i]);
				for (j=1; j<=sim->ss->nchan; j++) {
					//printf("\t mask[%d] = %d \n",j,mask[j]);
					if (OCpar.grad_shapes[i] == mask[j]) {
						OCchanmap[i] = j;
						Nch++;
						break;
					}
				}
				//printf("\t OCchanmap[%d] = %d \n",i,OCchanmap[i]);
			}
			if (Nch == 0) {  // no rf shape to take gradients
				_pulse_shaped_3pointrulerfmap(sim, wsp, Nelem, mask, steptime); // propagate with distorted shapes
			} else {  // lets do gradients
				if ( OCpar.gradmodeprop == 1 ) {
					/* do stuff for propagator optimization */
					return TclError(interp,"%s: does not support OC of propagators\n",cmdname);
					//_pulse_shapedOCprops_2(sim,wsp,Nelem,OCchanmap,mask2,steptime); // propagate with distorted shapes
				} else {
					/* do stuff for state to state optimization */
					//return TclError(interp,"%s: not yet supported in combination with OC", cmdname);
					_pulse_shapedOC_3pointrulerfmap(sim,wsp,Nelem,OCchanmap,mask,steptime);
				}
			}
			free_int_vector(OCchanmap);
		}
	} else {
		/* do just actual pulsing */
		_pulse_shaped_3pointrulerfmap(sim, wsp, Nelem, mask, steptime);
	}

	free_int_vector(mask);
	return TCL_OK;
}
/**** END of pulse_shaped_3pointrulerfmap implementation *****/

void tclcmd_pulse(Tcl_Interp* interp) {

  Tcl_CreateObjCommand(interp,"filter",(Tcl_ObjCmdProc *)tclFilter,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"pulse",(Tcl_ObjCmdProc *)tclPulse,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"reset",(Tcl_ObjCmdProc *)tclReset,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"pulseid",(Tcl_ObjCmdProc *)tclPulseid,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"delay",(Tcl_ObjCmdProc *)tclDelay,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"maxdt",(Tcl_ObjCmdProc *)tclMaxdt,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"acq",(Tcl_ObjCmdProc *)tclAcq,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"evolve",(Tcl_ObjCmdProc *)tclEvolve,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"prop",(Tcl_ObjCmdProc *)tclProp,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateObjCommand(interp,"store",(Tcl_ObjCmdProc *)tclStore,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"offset",(Tcl_ObjCmdProc *)tclOffset,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"select",(Tcl_ObjCmdProc *)tclSelect,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"turnon",(Tcl_CmdProc *)tclTurnon,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"turnoff",(Tcl_CmdProc *)tclTurnoff,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateCommand(interp,"matrix",(Tcl_CmdProc *)tclMatrix,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"matrixoper",(Tcl_ObjCmdProc *)tclMatrixoper,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"getinteractions",(Tcl_ObjCmdProc *)tclGetInteractions,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

  Tcl_CreateObjCommand(interp,"currenttime",(Tcl_ObjCmdProc *)tclCurrenttime,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"rotorangle",(Tcl_ObjCmdProc *)tclRotorangle,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  
/************ ZT: new commands: ***************/
  Tcl_CreateObjCommand(interp,"avgham_static",(Tcl_ObjCmdProc *)tclAvgHam_static,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(interp,"pulse_shaped",(Tcl_ObjCmdProc *)tclPulseShaped,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

/************ TV: new command: ****************/
  Tcl_CreateObjCommand(interp,"eulerangles",(Tcl_ObjCmdProc *)tclEulerAngles,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

/**** ZT: commands for zgrad pulses ****/
Tcl_CreateObjCommand(interp,"zgrad_pulse_shaped",(Tcl_ObjCmdProc *)tclZgrPulse,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"pulse_and_zgrad_shaped",(Tcl_ObjCmdProc *)tclPulseZgrShaped,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

/**** ZT: new acq commands ****/
Tcl_CreateObjCommand(interp,"acq_block",(Tcl_ObjCmdProc *)tclAcqBlock,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"acq_modulus",(Tcl_ObjCmdProc *)tclAcqModulus,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"pulse_shaped_rotormodulated",(Tcl_ObjCmdProc *)tclPulseShapedRotorModulated,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

Tcl_CreateObjCommand(interp,"test_function",(Tcl_ObjCmdProc *)tclTestFunc,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"pulse_shaped_linear",(Tcl_ObjCmdProc *)tclPulseShapedRotorModulatedLinear,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"pulse_shaped_3pointruleconstantrf",(Tcl_ObjCmdProc *)tclPulseShaped3PointRuleConstantRF,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"pulse_shaped_3pointrulerfmap",(Tcl_ObjCmdProc *)tclPulseShaped3PointRuleRFmap,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}
