#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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


#include "defs.h"
#include "cm.h"
#include "blockdiag.h"
#include "tclutil.h"
#include "cryst.h"
#include "OCroutines.h"
#include "rfshapes.h"
#include "pulse.h"
#include "iodata.h"
#include "tclutil.h"
#include "B0inhom.h"
#include "ham.h"
#include "lbfgs.h"

/* added for testing speed improvements */
//#include <windows.h>
//#include <winbase.h>
//LARGE_INTEGER _tickpsec_, _tv1_, _tv2_;

/* global variable holding all OC parameters */
OCoptPars OCpar;

/****
 * initialize global variable OCpar
 ****/ 
void OCpar_initialize(void) 
{
	OCpar.gradmode=0;
	OCpar.gradmodeprop=0;
	OCpar.var_shapes=NULL;
	OCpar.var_shapes_min=NULL;
	OCpar.var_shapes_max=NULL;
	OCpar.var_shapes_rmsmax=NULL;
	OCpar.grad_shapes=NULL;
	OCpar.var_shapes_penalty_order = NULL;
	OCpar.var_shapes_penalty_factor = NULL;
	OCpar.fvals[0] = 0.0;
	OCpar.fvals[1] = 0.0;
	OCpar.ispenalty = 0;
	OCpar.isinit=1;
}

/****
 * clear global variable OCpar
 ****/
void OCpar_destroy(void)
{
  OCpar.gradmode=0;
  OCpar.gradmodeprop=0;
  
  if (OCpar.var_shapes) {
    free_int_vector(OCpar.var_shapes);
    OCpar.var_shapes=NULL;
  }
  if (OCpar.var_shapes_min) {
    free_double_vector(OCpar.var_shapes_min);
    OCpar.var_shapes_min=NULL;
  }
  if (OCpar.var_shapes_max) {
    free_double_vector(OCpar.var_shapes_max);
    OCpar.var_shapes_max=NULL;
  }
  if (OCpar.var_shapes_rmsmax) {
    free_double_vector(OCpar.var_shapes_rmsmax);
    OCpar.var_shapes_rmsmax=NULL;
  }
  if (OCpar.grad_shapes) {
    free_int_vector(OCpar.grad_shapes);
    OCpar.grad_shapes=NULL;
  }
  if (OCpar.var_shapes_penalty_order) {
    free_int_vector(OCpar.var_shapes_penalty_order);
    OCpar.var_shapes_penalty_order=NULL;
  }
  if (OCpar.var_shapes_penalty_factor) {
    free_double_vector(OCpar.var_shapes_penalty_factor);
    OCpar.var_shapes_penalty_factor=NULL;
  }
  OCpar.isinit=0;
}

/****
 * increase matrix pointer
 ****/
void incr_OCmx_pos(Sim_wsp *wsp)
{
  if (wsp->OC_mxpos == MAXOCPROPS-1) {
     fprintf(stderr,"error: there is no more space for OC propagators\n");
     exit(1);
  }
  
  (wsp->OC_mxpos)++;
}

/****
 * store current propagator to wsp->OC_props[i]
 ****/
void store_OCprop(Sim_wsp *wsp)
{
  int i = wsp->OC_mxpos;
  
  /* debug print */
  //printf("Storing OC propagator to slot %d\n",i);
   
  if (wsp->OC_props[i] != NULL) free_blk_mat_complx(wsp->OC_props[i]);
  wsp->OC_props[i] = blk_cm_dup(wsp->U);

}

/****
 * store current density matrix to OCpar->dens[i]
 ****/
void store_OCdens(Sim_info *sim, Sim_wsp *wsp)
{
  int i = wsp->OC_mxpos;
  int ii;

  /*  printf("Storing OC matrix density to slot %d\n",i); */
  for (ii=0; ii<sim->Nfstart; ii++) {
	  if (wsp->OC_dens[i+MAXOCPROPS*ii] == NULL)
		  wsp->OC_dens[i+MAXOCPROPS*ii] = cm_dup(wsp->sigma[ii]);
	  else
		  cm_copy(wsp->OC_dens[i+MAXOCPROPS*ii],wsp->sigma[ii]);
  }

}

/****
 * this serves debugging purposes only
 ****/
void test_print_codes(Sim_wsp *wsp)
{
  int i;
  for (i=1; i<=wsp->OC_mxpos; i++) {
     if (wsp->OC_mxcode[i]) {
        printf("/%d/ code is '%s'\n",i,wsp->OC_mxcode[i]);
     }
  }
}

/****
 * store code to matrix to OCpar.mx_code[i]
 ****/
void set_OCmx_code(Sim_wsp *wsp, char *c)
{
  int i = wsp->OC_mxpos;
  
  if (wsp->OC_mxcode[i] == NULL) {
     wsp->OC_mxcode[i] = malloc(strlen(c)+1);
  } else {
     free( (char*)(wsp->OC_mxcode[i]) );
     wsp->OC_mxcode[i] = malloc(strlen(c)+1);
  }
  strcpy(wsp->OC_mxcode[i],c);
  //printf("OC_mxcode: %3d) -%s-\n",wsp->OC_mxpos,wsp->OC_mxcode[wsp->OC_mxpos]);
  
}

/****
 * store filter matrix as OCprop
 ****/
void store_OCfilter(Sim_wsp *wsp,int num)
{
  int i = wsp->OC_mxpos;
  blk_mat_complx *obj;
  mat_complx *cm = wsp->matrix[num];
  
  /* debug print */
  /*   printf("Storing OC filter matrix to slot %d\n",i); */

  if (wsp->OC_props[i] != NULL) free_blk_mat_complx(wsp->OC_props[i]);

  obj = (blk_mat_complx*)malloc(sizeof(blk_mat_complx));
  obj->dim = cm->row;
  obj->Nblocks = 1;
  obj->basis = cm->basis;
  obj->blk_dims = (int*)malloc(sizeof(int));
  obj->blk_dims[0] = obj->dim;
  obj->m = cm_dup(cm);

  wsp->OC_props[i] = obj;
}

char* my_strtok_r(char *str, const char *delim, char **nextp)
{
    char *ret;

    //if (str == NULL) str = *nextp;
    assert(str != NULL);
    // clear initial sequence of delimiters
    str += strspn(str, delim);
    if (*str == '\0') {
        return NULL;
    }
    ret = str;
    str += strcspn(str, delim);
    // if not at the end of string put there end of string
    if (*str) {
        *str++ = '\0';
    }
    *nextp = str;

    return ret;
}

/****
 * check and set switch for optimization of propagators
 *  (in that mode cummulative propagators are calculated)
 ****/
void test_pulseq_for_acqOC_prop(Tcl_Interp *interp)
{
  char *res, *lin, *dum, *chptr;
  char buf[256],pulseq[128];
  int l=0;
  
  TclGetString(interp,pulseq,"par","pulse_sequence",0,"pulseq");
  strcpy(buf,"info body ");
  strcat(buf,pulseq);
  if ( Tcl_Eval(interp,buf) != TCL_OK ) {
     fprintf(stderr,"test_pulseq_for_acqOC_prop can not execute %s :\n",buf);
     fprintf(stderr,"%s",Tcl_GetStringResult(interp));
     exit(1);
  }
     
  //printf("test_pulseq_for_acqOC_prop:\n%s\n---\n",buf);

  res = (char *)Tcl_GetStringResult(interp);
  //printf("%s\n///\n",res);

  lin = my_strtok_r(res,"\n",&chptr);
  //printf("First line: '%s'\n///\n",lin);
  while ( lin != NULL) {
     l++;
     dum = strstr(lin,"oc_acq_prop");
     if ( dum ) {
        /* found the command, now check if it is commented out */
    	DEBUGPRINT("found oc_acq_prop on line %d ",l);
    	char *comt;
    	comt = strchr(lin,'#');
    	if ( comt ) {
    		/* comment sign found, check if it is before the command */
    		if ( comt < dum ) {
    			DEBUGPRINT("but it is within a comment\n");
    		} else {
    			DEBUGPRINT(" and it seems to be active command\n");
    			OCpar.gradmodeprop++;
    		}
        } else {
        	DEBUGPRINT(" and it seems to be active command\n");
        	OCpar.gradmodeprop++;
        }
     }  
     lin = my_strtok_r(chptr,"\n",&chptr);
     //printf("next line: '%s'\n///\n",lin);
  }
  
  if (OCpar.gradmodeprop > 1) {
     fprintf(stderr,"test_pulseq_for_acqOC_prop error: the command 'oc_acq_prop' can be used\n   only once in the pulse sequence!!! Please do not use reserved word\n   'oc_acq_prop' neither as a text argument for other commands.\n");
     exit(1);
  }

}



/****
 * filter command in gradient mode 
 ****/
void _filterOC(Sim_info *sim, Sim_wsp *wsp, int num)
{
  int ii;
  mat_complx *mask;

  mask = wsp->matrix[num];

  if (!mask) {
    fprintf(stderr,"error: illegal filter number %d\n",num);
    exit(1);
  }
  if (mask->row != sim->matdim || mask->col != sim->matdim) {
    fprintf(stderr,"error: filter: matrix %d must be a %dx%d matrix\n",num,sim->matdim,sim->matdim);
    exit(1);
  }

  if (wsp->Uisunit != 1) {
     incr_OCmx_pos(wsp);
     store_OCprop(wsp);
     _evolve_with_prop(sim,wsp);
     _reset_prop(sim,wsp);
     store_OCdens(sim,wsp);
     set_OCmx_code(wsp,"P");
  }
  wsp->cannotbestored=1;
  for (ii=0; ii<sim->Nfstart; ii++)
	  cm_filter(wsp->sigma[ii],mask);

  incr_OCmx_pos(wsp);
  store_OCfilter(wsp,num);
  store_OCdens(sim,wsp);
  set_OCmx_code(wsp,"F");
}

/****
 * filter command in backevolution when calculating gradients 
 ****/
void lambda_filterOC( Sim_wsp *wsp, int num, mat_complx *lam)
{
  assert(wsp->OC_props[num]->Nblocks == 1);
  cm_filter(lam, wsp->OC_props[num]->m);
}

/****
 * creates complex Ham including interactions and rf (with phase)
 ****/
void ham_hamrf_complex(mat_complx *cH, Sim_info *sim, Sim_wsp *wsp)
{
	assert(cH != NULL);
	assert(cH->row == wsp->ham_blk->dim);
	assert(cH->col == cH->row);
	assert(cH->type == MAT_DENSE);
	assert(wsp->ham_blk->Nblocks == 1); /* no block_diag */
	assert(wsp->chan_Ix[1]->type == MAT_DENSE);
	assert(wsp->chan_Iy[1]->type == MAT_DENSE);

	int i;
	int dim = cH->row;
	int dim2 = dim*dim;
	cm_zero(cH);
	/* add interactions */
	mat_double *smx = wsp->ham_blk->m;
	switch (smx->type) {
		case MAT_DENSE:
			cblas_dcopy(dim2,smx->data,1,(double*)(cH->data),2);
			break;
		case MAT_DENSE_DIAG:
			cblas_dcopy(dim,smx->data,1,(double*)(cH->data),2*(dim+1));
			//printf("H_1_1 = %g (%g, %g)\n",smx->data[0],cH->data[0].re,cH->data[0].im);
			break;
		case MAT_SPARSE:
		default:
			fprintf(stderr,"Error: ham_hamrf_complex - unsupported matrix type (%d)\n",smx->type);
			exit(1);
	}
	/* add rf fields */
	double c, s;
	for (i=1; i<=sim->ss->nchan; i++) {
  	  if (fabs(wsp->rfv[i]) > TINY) {
  		  c = wsp->rfv[i]*cos(wsp->phv[i]);
  		  s = wsp->rfv[i]*sin(wsp->phv[i]);
  		  cblas_daxpy(dim2,c,wsp->chan_Ix[i]->data,1,(double*)(cH->data),2);
  		  cblas_daxpy(dim2,s,wsp->chan_Iy[i]->data,1,(double*)(cH->data)+1,2);
  	  }
	}
}

/****
 * fills auxiliary complex matrix with complex Ham including interactions
 *     and rf (with phase)
 *     For exact derivatives via exponential
 ****/
void ham_hamrf_complex_2(mat_complx *cH, int chan, int ph, Sim_info *sim, Sim_wsp *wsp)
{
	assert(cH != NULL);
	assert(cH->row == 2*wsp->ham_blk->dim);
	assert(cH->col == cH->row);
	assert(cH->type == MAT_DENSE);
	assert(wsp->ham_blk->Nblocks == 1); /* no block_diag */
	assert(wsp->chan_Ix[chan]->type == MAT_DENSE);
	assert(wsp->chan_Iy[chan]->type == MAT_DENSE);

	int i,j;
	int dim = cH->row;
	int dim2 = dim*dim;
	cm_zero(cH);
	/* add interactions to first block */
	mat_double *smx = wsp->ham_blk->m;
	switch (smx->type) {
		case MAT_DENSE:
			for (i=0; i<dim; i++) {
				cblas_dcopy(dim,smx->data+i*dim,1,(double*)(cH->data)+2*2*dim*i,2);
			}
			break;
		case MAT_DENSE_DIAG:
			cblas_dcopy(dim,smx->data,1,(double*)(cH->data),2*2*dim+2);
			break;
		case MAT_SPARSE:
		default:
			fprintf(stderr,"Error: ham_hamrf_complex_2 - unsupported matrix type (%d)\n",smx->type);
			exit(1);
	}
	/* add rf fields to first block */
	double c, s;
	for (i=1; i<=sim->ss->nchan; i++) {
  	  if (fabs(wsp->rfv[i]) > TINY) {
  		  c = wsp->rfv[i]*cos(wsp->phv[i]);
  		  s = wsp->rfv[i]*sin(wsp->phv[i]);
  		  for (j=0; j<dim; j++) {
  			  cblas_daxpy(dim,c,wsp->chan_Ix[i]->data+j*dim,1,(double*)(cH->data)+2*2*dim*j,2);
  			  cblas_daxpy(dim,s,wsp->chan_Iy[i]->data+j*dim,1,(double*)(cH->data)+2*2*dim*j+1,2);
  		  }
  	  }
	}
	/* copy first block to the second, add channel for derivative */
	for (i=0; i<dim; i++) {
		cblas_zcopy(dim, cH->data+2*dim*i, 1, cH->data+2*dim2+dim+2*dim*i,1);
		switch (ph) {
		case 1:
			cblas_dcopy(dim,wsp->chan_Ix[chan]->data,1,(double*)(cH->data)+2*2*dim2+i*2*2*dim,2);
			break;
		case 2:
			cblas_dcopy(dim,wsp->chan_Iy[chan]->data,1,(double*)(cH->data)+2*2*dim2+i*2*2*dim+1,2);
			break;
		default:
			fprintf(stderr,"Error: ham_hamrf_complex_2 - invalid argument ph (%d)\n",ph);
			exit(1);
		}
	}
}

void OC_commutator(mat_complx *cH,mat_complx *kom,int iter,mat_complx *cdum, double *tol)
{
	complx *p1, *p2, *p3, *p4, *pstop;
	int dim = cH->row;
	int dim2 = dim*dim;
	int sign = -1;

	cblas_zhemm(CblasColMajor,CblasLeft,CblasUpper,dim,dim,&Cunit,cH->data,dim,kom->data,dim,&Cnull,cdum->data,dim);
	//cm_print(cdum,"cH*kom");

	p1 = kom->data; /* result */
	p2 = p3 = cdum->data;
	pstop = kom->data + dim2;
	p4 = cdum->data + dim2;
	while (--iter > 0) sign *= -1;
	*tol = -1000;
	while (p1<pstop) {
		p1->re = p2->re + sign*p3->re;
		p1->im = p2->im - sign*p3->im;
		//printf("\t p2 = %g, %g, p3 = %g, %g\n",p2->re,p2->im,p3->re,p3->im);
		if (fabs(p1->re)+fabs(p1->im) > *tol) *tol = fabs(p1->re)+fabs(p1->im);
		p1++;
		p2++;
		p3 += dim;
		if (p3>=p4) p3 = p3 - dim2 + 1;
	}
}

/*****
 * for DNPframe and electron channel only
 * grads via expansion in nested commutators. Mostly copy-paste from _pulse_shapedOC_2
 */
void _pulse_shapedOC_2_dnpframe(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *OCchanmap, int *mask, double duration)
{
	int i, j, k, n;
	double dt, dt_us;
	int dim = wsp->ham_blk->dim;
	int dim2 = dim*dim;
	int Nsh = LEN(OCchanmap);
	mat_complx *cH = NULL, *kom = NULL, *cdum = NULL;
	//complx *cvec = NULL;
	mat_complx *gr = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	// variables for rfmap hack
	int iz = 0, iphi = 0, Nphi = 0, iphi_shift = 0;
	double steptimephi = 0.0, phitime = 0.0;
	// duration is steptime

	assert(wsp->dU != NULL);
	assert(wsp->dU->Nblocks == 1); /* no block_diag */
	assert(wsp->Hcplx != NULL);
	assert(wsp->Hcplx->Nblocks == 1);

	if (wsp->Uisunit != 1) {
		/* there is some pending propagator, store it but make no gradient */
		incr_OCmx_pos(wsp);
		store_OCprop(wsp);
		store_OCdens(sim,wsp); /* store sigma from previous step */
		_evolve_with_prop(sim,wsp); /* get sigma ready for next store */
		_reset_prop(sim,wsp);
	}

	int iter;
	int maxiter = -1; /* exact gradients calculated via exponential */
	double tol;
	double maxtol = 1e-6;
	if (fabs(OCpar.grad_level)<0.5) {
		/* gradient with higher order corrections up to grad_level accuracy */
		maxiter = 10;
		maxtol = fabs(OCpar.grad_level);
	} else if (OCpar.grad_level > 1.5){
		/* gradient with higher order corrections up to grad_level order */
		maxiter = (int)round(OCpar.grad_level);
		maxtol = 1e-10;
	}
	//printf("pulse_shaped_OC: maxiter = %d, tol = %g\n",maxiter,maxtol);

	/* do pulsing, gradients, evolving, storing the results */
	n = (int)ceil(duration/wsp->dtmax);
	if (n < 1) n = 1;
	dt_us = duration/(double)n;
	dt = dt_us*1.0e-6;
	//printf("_pulse_shapedOC_2 duration of %f us split into %d steps of %f us\n",duration,n,dt*1.0e+6);
	// pre-hack block of rfmap and RotorModulated
	if (sim->rfmap != NULL) { // we have RotorModulated, initialize variables
		iz = sim->rfmap->loop[(wsp->rf_idx)*2+0]; // z coil coordinate index
		iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1];  // coil initial "phase" index
		Nphi = sim->rfmap->z[iz*(sim->rfmap->Nch+1)+sim->rfmap->Nch]; // number of phi elements for current z position
		steptimephi = sim->taur/Nphi;
		assert(steptimephi - duration >= -TINY); // limitation for fine-digitized shapes, steptime <=> duration
		iphi_shift = 0;
		phitime = 0.0;
	} // end of pre-hack block of rfmap and RotorModulated
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
		incr_OCmx_pos(wsp);
		/* like _pulse_simple */
		for (i=1;i<=n;i++) { // loop over maxdt steps
			ham_hamilton(sim,wsp);
			_setrfprop_labdnp(sim,wsp); // it needs to be called here as the pulse parameters are changing every maxdt (wsp->dtmax)
			//blk_dm_print(wsp->ham_blk,"Ham int real");
			if (maxiter<0) {
				/* exact gradients calculated via exponential */
				fprintf(stderr,"Error: _pulse_shapedOC_2_dnpframe - not implemented exact grads via expm\n");
				exit(1);
			}else {
				double scl;
				complx cscl;
				/* gradients improved with higher order terms */
				if (kom == NULL) {
					kom = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
					cdum = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
				}
				//ham_hamrf_complex(cH, sim, wsp);
				blk_cm_multod(wsp->Hcplx,wsp->Hrf_blk,1.0);
				cH = wsp->Hcplx->m;
				blk_cm_unit(wsp->dU);
				blk_prop_complx_3(wsp->dU,wsp->Hcplx,dt,sim);
				for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
					//printf("OCchanmap[%d] = %d\n",k,OCchanmap[k]);
					if (OCchanmap[k] < 0) { /* this shape is not active */
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = NULL;
							//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
						}
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = NULL;
							//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
						}
						continue;
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
						//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
						//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
					}
					/* x channel */
					cm_zero(gr);
					cblas_daxpy(dim2,-1.0,wsp->chan_Ix[OCchanmap[k]]->data,1,(double*)(gr->data)+1,2); /* first order */
					cm_zero(kom);
					cblas_dcopy(dim2,wsp->chan_Ix[OCchanmap[k]]->data,1,(double*)(kom->data),2);
					//cm_print(kom,"channel Ix");
					iter=1;
					tol = 1000;
					scl = dt/2.0;
					while (iter < maxiter && tol > maxtol) {
						//printf("Iter %d, tol %g\n============\n",iter,tol);
						//cm_print(gr,"gr matrix");
						OC_commutator(cH,kom,iter,cdum, &tol);
						//cm_print(kom,"komutator");
						switch (iter % 4) {
						case 0: cscl.re = 0; cscl.im = -scl;
							break;
						case 1: cscl.re = scl; cscl.im = 0;
							break;
						case 2: cscl.re = 0; cscl.im = scl;
							break;
						case 3: cscl.re = -scl; cscl.im = 0;
							break;
						}
						cblas_zaxpy(dim2,&cscl, kom->data,1,gr->data,1);
						tol *= scl;
						scl *= (dt/((double)(iter)+2.0));
						iter++;
					}
					cm_multo_rev(gr,wsp->dU->m); /* gr is done here */
					//cm_print(gr,"READY grad");
					//exit(1);
					if (i==1) {
						cm_copy(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)],gr);
					} else {
						cm_multo_rev(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)],wsp->dU->m);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cunit,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]->data,dim);
					}
					/* y channel */
					cm_zero(gr);
					cblas_daxpy(dim2,1.0,wsp->chan_Iy[OCchanmap[k]]->data,1,(double*)(gr->data),2); /* first order */
					cm_zero(kom);
					cblas_dcopy(dim2,wsp->chan_Iy[OCchanmap[k]]->data,1,(double*)(kom->data)+1,2);
					//cm_print(kom,"channel Iy");
					iter=1;
					tol = 1000;
					scl = dt/2.0;
					while (iter < maxiter && tol > maxtol) {
						//printf("Iter %d, tol %g\n============\n",iter,tol);
						//cm_print(gr,"gr matrix");
						OC_commutator(cH,kom,iter,cdum, &tol);
						//cm_print(kom,"komutator");
						switch (iter % 4) {
						case 0: cscl.re = 0; cscl.im = -scl;
							break;
						case 1: cscl.re = scl; cscl.im = 0;
							break;
						case 2: cscl.re = 0; cscl.im = scl;
							break;
						case 3: cscl.re = -scl; cscl.im = 0;
							break;
						}
						cblas_zaxpy(dim2,&cscl, kom->data,1,gr->data,1);
						tol *= scl;
						scl *= (dt/((double)(iter)+2.0));
						iter++;
					}
					cm_multo_rev(gr,wsp->dU->m); /* gr is done here */
					//cm_print(gr,"READY grad");
					//exit(1);
					if (i==1) {
						cm_copy(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1],gr);
					} else {
						cm_multo_rev(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1],wsp->dU->m);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cunit,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]->data,dim);
					}
				}
			}
			/* update kumulativni propagator ve wsp->U */
			update_propagator(wsp->U, wsp->dU, sim, wsp);
			wsp->t += dt_us;
		} /* end for i over maxdt steps inside pulse step */
		store_OCprop(wsp);
		store_OCdens(sim,wsp); /* sigma of previous step */
		_evolve_with_prop(sim,wsp); /* get sigma ready for the next step */
		_reset_prop(sim,wsp);
		// hack for RotorModulated
		if (sim->rfmap != NULL) { // we have RotorModulated
			for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
				if (OCchanmap[k] < 0) { /* this shape is not active */
					continue;
				}
			//-	int iz = sim->rfmap->loop[(wsp->rf_idx)*2+0];
			//-	int iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1];
				double bx, by;
				//rfmap_get(sim->rfmap,iz,iphi+j-1,OCchanmap[k]-1,&bx,&by);
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,OCchanmap[k]-1,&bx,&by); // get distortion coefs
				mat_complx *dUx = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)];   // dU/dwx
				mat_complx *dUy = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]; // dU/dwy
				mat_complx *dUx2 = cm_dup(dUx);
				cm_muld(dUx,bx);
				cm_multod(dUx,dUy,by);
				cm_muld(dUy,bx);
				cm_multod(dUy,dUx2,-by);
				free_complx_matrix(dUx2);
				//printf("GRAD: Thread %d: chan %d, slot %3d: phitime = %g ... %d %g %g\n",wsp->thread_id,OCchanmap[k],-1,phitime,j,bx,by);
			}
			// this needs to be last at the end of loop over j (shape elements)
			phitime += duration;  // steptime <=> duration
  		    if (steptimephi-phitime < 0.5*duration) {
			    phitime -= steptimephi;
				iphi_shift++;
			}
		} // end of hack
	} /* end for j over Nelem */

	// cH is a pointer to existing complex Hamiltonian wsp->Hcplx->m, do nto free it here
	if (kom != NULL) free_complx_matrix(kom);
	if (cdum != NULL) free_complx_matrix(cdum);
	free_complx_matrix(gr);
	//printf("_pulse_shapedOC_2 done\n");

}

/****
 * NOT FINISHED exact gradients via exponentialization.
 * Finished: grads via expansion in nested commutators
 ****/
void _pulse_shapedOC_2(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *OCchanmap, int *mask, double duration)
{
	int i, j, k, l, n;
	double dt, dt_us;
	int dim = wsp->ham_blk->dim;
	int dim2 = dim*dim;
	int Nsh = LEN(OCchanmap);
	mat_complx *cH = NULL, *kom = NULL, *cdum = NULL;
	complx *cvec = NULL;
	mat_complx *gr = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	// variables for rfmap hack
	int iz = 0, iphi = 0, Nphi = 0, iphi_shift = 0;
	double steptimephi = 0.0, phitime = 0.0;
	// duration is steptime

	assert(wsp->dU != NULL);
	assert(wsp->dU->Nblocks == 1); /* no block_diag */

	if (wsp->Uisunit != 1) {
		/* there is some pending propagator, store it but make no gradient */
		incr_OCmx_pos(wsp);
		store_OCprop(wsp);
		store_OCdens(sim,wsp); /* store sigma from previous step */
		_evolve_with_prop(sim,wsp); /* get sigma ready for next store */
		_reset_prop(sim,wsp);
	}

	int iter;
	int maxiter = -1; /* exact gradients calculated via exponential */
	double tol;
	double maxtol = 1e-6;
	if (fabs(OCpar.grad_level)<0.5) {
		/* gradient with higher order corrections up to grad_level accuracy */
		maxiter = 10;
		maxtol = fabs(OCpar.grad_level);
	} else if (OCpar.grad_level > 1.5){
		/* gradient with higher order corrections up to grad_level order */
		maxiter = (int)round(OCpar.grad_level);
		maxtol = 1e-10;
	}
	//printf("pulse_shaped_OC: maxiter = %d, tol = %g\n",maxiter,maxtol);

	/* do pulsing, gradients, evolving, storing the results */
	n = (int)ceil(duration/wsp->dtmax);
	if (n < 1) n = 1;
	dt_us = duration/(double)n;
	dt = dt_us*1.0e-6;
	//printf("_pulse_shapedOC_2 duration of %f us split into %d steps of %f us\n",duration,n,dt*1.0e+6);
	// pre-hack block of rfmap and RotorModulated
	if (sim->rfmap != NULL) { // we have RotorModulated, initialize variables
		iz = sim->rfmap->loop[(wsp->rf_idx)*2+0]; // z coil coordinate index
		iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1];  // coil initial "phase" index
		Nphi = sim->rfmap->z[iz*(sim->rfmap->Nch+1)+sim->rfmap->Nch]; // number of phi elements for current z position
		steptimephi = sim->taur/Nphi;
		assert(steptimephi - duration >= -TINY); // limitation for fine-digitized shapes, steptime <=> duration
		iphi_shift = 0;
		phitime = 0.0;
	} // end of pre-hack block of rfmap and RotorModulated
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
		_setrfprop(sim,wsp);
		incr_OCmx_pos(wsp);
		/* like _pulse_simple */
		for (i=1;i<=n;i++) {
			ham_hamilton(sim,wsp);
			//blk_dm_print(wsp->ham_blk,"Ham int real");
			if (maxiter<0) {
				/* exact gradients calculated via exponential */
				if (cH == NULL) {
					cH = complx_matrix(wsp->ham_blk->dim*2,wsp->ham_blk->dim*2,MAT_DENSE,0,wsp->ham_blk->basis);
					cdum = complx_matrix(wsp->ham_blk->dim*2,wsp->ham_blk->dim*2,MAT_DENSE,0,wsp->ham_blk->basis);
					cvec = complx_vector(2*dim);
				}
				for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
					if (OCchanmap[k] < 0) { /* this shape is not active */
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = NULL;
						}
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = NULL;
						}
						continue;
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
					}
					/* x channel */
					ham_hamrf_complex_2(cH,k,1,sim,wsp);
					fprintf(stderr,"Error: gradient via exponential not finished.\n");
					exit(1);
					cm_diag(cH,cvec,cdum); /* problem jak expm obecne komplexni matice... */
					cm_zero(cH);
					for (l=0; l<2*dim; l++) {
						double f = exp(cvec[l+1].im*dt);
						cH->data[l].re = f*cos(-cvec[l+1].re*dt);
						cH->data[l].im = f*sin(-cvec[l+1].re*dt);
					}
					simtrans(cH,cdum); /* tihle nemusi platit, obecne potrebuju inverzi matice, ne adjoint */
					fprintf(stderr,"Error: _pulse_shapedOC_2 - not implemented exact grads via expm\n");
					exit(1);
					/* y channel */
					/* buffery s derivacema vynasobit  wsp->dU, krome i==1 */
					/* do bufferu pricist grad_i * wsp->U */
				}
			} else {
				double scl;
				complx cscl;
				/* gradients improved with higher order terms */
				if (cH == NULL) {
					cH = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
					kom = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
					cdum = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
				}
				ham_hamrf_complex(cH, sim, wsp);
				//cm_print(cH, "Hamiltonian");
				blk_dm_multod(wsp->ham_blk,wsp->sumHrf,1.0);
				//blk_dm_print(wsp->ham_blk,"Ham total real");
				blk_cm_unit(wsp->dU);
				blk_prop_real(wsp->dU,wsp->ham_blk,dt,sim);
				blk_simtrans_zrot2(wsp->dU,wsp->sumUph);
				//blk_cm_print(wsp->dU,"Propagator");
				for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
					//printf("OCchanmap[%d] = %d\n",k,OCchanmap[k]);
					if (OCchanmap[k] < 0) { /* this shape is not active */
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = NULL;
							//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
						}
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = NULL;
							//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
						}
						continue;
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
						//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
						//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
					}
					/* x channel */
					cm_zero(gr);
					cblas_daxpy(dim2,-1.0,wsp->chan_Ix[OCchanmap[k]]->data,1,(double*)(gr->data)+1,2); /* first order */
					cm_zero(kom);
					cblas_dcopy(dim2,wsp->chan_Ix[OCchanmap[k]]->data,1,(double*)(kom->data),2);
					//cm_print(kom,"channel Ix");
					iter=1;
					tol = 1000;
					scl = dt/2.0;
					while (iter < maxiter && tol > maxtol) {
						//printf("Iter %d, tol %g\n============\n",iter,tol);
						//cm_print(gr,"gr matrix");
						OC_commutator(cH,kom,iter,cdum, &tol);
						//cm_print(kom,"komutator");
						switch (iter % 4) {
						case 0: cscl.re = 0; cscl.im = -scl;
							break;
						case 1: cscl.re = scl; cscl.im = 0;
							break;
						case 2: cscl.re = 0; cscl.im = scl;
							break;
						case 3: cscl.re = -scl; cscl.im = 0;
							break;
						}
						cblas_zaxpy(dim2,&cscl, kom->data,1,gr->data,1);
						tol *= scl;
						scl *= (dt/((double)(iter)+2.0));
						iter++;
					}
					cm_multo_rev(gr,wsp->dU->m); /* gr is done here */
					//cm_print(gr,"READY grad");
					//exit(1);
					if (i==1) {
						cm_copy(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)],gr);
					} else {
						cm_multo_rev(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)],wsp->dU->m);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cunit,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]->data,dim);
					}
					/* y channel */
					cm_zero(gr);
					cblas_daxpy(dim2,1.0,wsp->chan_Iy[OCchanmap[k]]->data,1,(double*)(gr->data),2); /* first order */
					cm_zero(kom);
					cblas_dcopy(dim2,wsp->chan_Iy[OCchanmap[k]]->data,1,(double*)(kom->data)+1,2);
					//cm_print(kom,"channel Iy");
					iter=1;
					tol = 1000;
					scl = dt/2.0;
					while (iter < maxiter && tol > maxtol) {
						//printf("Iter %d, tol %g\n============\n",iter,tol);
						//cm_print(gr,"gr matrix");
						OC_commutator(cH,kom,iter,cdum, &tol);
						//cm_print(kom,"komutator");
						switch (iter % 4) {
						case 0: cscl.re = 0; cscl.im = -scl;
							break;
						case 1: cscl.re = scl; cscl.im = 0;
							break;
						case 2: cscl.re = 0; cscl.im = scl;
							break;
						case 3: cscl.re = -scl; cscl.im = 0;
							break;
						}
						cblas_zaxpy(dim2,&cscl, kom->data,1,gr->data,1);
						tol *= scl;
						scl *= (dt/((double)(iter)+2.0));
						iter++;
					}
					cm_multo_rev(gr,wsp->dU->m); /* gr is done here */
					//cm_print(gr,"READY grad");
					//exit(1);
					if (i==1) {
						cm_copy(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1],gr);
					} else {
						cm_multo_rev(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1],wsp->dU->m);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cunit,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]->data,dim);
					}
				}
			}
			/* update kumulativni propagator ve wsp->U */
			update_propagator(wsp->U, wsp->dU, sim, wsp);
			wsp->t += dt_us;
		} /* end for i over maxdt steps inside pulse step */
		store_OCprop(wsp);
		store_OCdens(sim,wsp); /* sigma of previous step */
		_evolve_with_prop(sim,wsp); /* get sigma ready for the next step */
		_reset_prop(sim,wsp);
		// hack for RotorModulated
		if (sim->rfmap != NULL) { // we have RotorModulated
			for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
				if (OCchanmap[k] < 0) { /* this shape is not active */
					continue;
				}
			//-	int iz = sim->rfmap->loop[(wsp->rf_idx)*2+0];
			//-	int iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1];
				double bx, by;
				//rfmap_get(sim->rfmap,iz,iphi+j-1,OCchanmap[k]-1,&bx,&by);
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,OCchanmap[k]-1,&bx,&by); // get distortion coefs
				mat_complx *dUx = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)];   // dU/dwx
				mat_complx *dUy = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]; // dU/dwy
				mat_complx *dUx2 = cm_dup(dUx);
				cm_muld(dUx,bx);
				cm_multod(dUx,dUy,by);
				cm_muld(dUy,bx);
				cm_multod(dUy,dUx2,-by);
				free_complx_matrix(dUx2);
				//printf("GRAD: Thread %d: chan %d, slot %3d: phitime = %g ... %d %g %g\n",wsp->thread_id,OCchanmap[k],-1,phitime,j,bx,by);
			}
			// this needs to be last at the end of loop over j (shape elements)
			phitime += duration;  // steptime <=> duration
  		    if (steptimephi-phitime < 0.5*duration) {
			    phitime -= steptimephi;
				iphi_shift++;
			}
		} // end of hack
	} /* end for j over Nelem */

	if (cH != NULL) free_complx_matrix(cH);
	if (kom != NULL) free_complx_matrix(kom);
	if (cdum != NULL) free_complx_matrix(cdum);
	free_complx_matrix(gr);
	//printf("_pulse_shapedOC_2 done\n");
}

/****
 * helper function for pulse_shaped in gradient mode and state to state optimization
 *        (here step by step propagators are saved together with step-evolved
 *         density matrices)
 * Original GRAPE ala Khaneja
 ****/
 void _pulse_shapedOC(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, double steptime)
{
  int i, j;
  char cd[128];

   if (wsp->Uisunit != 1) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      _evolve_with_prop(sim,wsp);
      store_OCdens(sim,wsp);
      _reset_prop(sim,wsp);
      set_OCmx_code(wsp,"P");
   }
   /* do pulsing and evolving, storing the results */
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
      incr_OCmx_pos(wsp);
      _pulse(sim,wsp,steptime);
      store_OCprop(wsp);
      _evolve_with_prop(sim,wsp);
      store_OCdens(sim,wsp);
      _reset_prop(sim,wsp);
      if (j==1)
         set_OCmx_code(wsp,"G");
   }
   /* set the code only for the last memory slot */ 
   sprintf(cd,"G%dE%d",Nelem,Nch); 
   strcat(cd,code);
   set_OCmx_code(wsp,cd);

}

void _pulse_shapedOCprops_2(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *OCchanmap, int *mask, double duration)
{
	int i, j, k, n;
	double dt, dt_us;
	int dim = wsp->ham_blk->dim;
	int dim2 = dim*dim;
	int Nsh = LEN(OCchanmap);
	mat_complx *cH = NULL, *kom = NULL, *cdum = NULL;
	mat_complx *gr = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	// variables for rfmap hack
	int iz = 0, iphi = 0, Nphi = 0, iphi_shift = 0;
	double steptimephi = 0.0, phitime = 0.0;
	// duration is steptime

	assert(wsp->dU != NULL);
	assert(wsp->dU->Nblocks == 1); /* no block_diag */

	if ( (wsp->Uisunit != 1) && (wsp->OC_propstatus == 1) ) {
		/* there is some pending propagator, store it but make no gradient */
		incr_OCmx_pos(wsp);
		store_OCprop(wsp);
	}

	int iter;
	int maxiter = -1; /* exact gradients calculated via exponential */
	double tol;
	double maxtol = 1e-6;
	if (fabs(OCpar.grad_level)<0.5) {
		/* gradient with higher order corrections up to grad_level accuracy */
		maxiter = 10;
		maxtol = fabs(OCpar.grad_level);
	} else if (OCpar.grad_level > 1.5){
		/* gradient with higher order corrections up to grad_level order */
		maxiter = (int)round(OCpar.grad_level);
		maxtol = 1e-10;
	}

	/* do pulsing, element gradients, storing cumulative props */
	n = (int)ceil(duration/wsp->dtmax);
	if (n < 1) n = 1;
	dt_us = duration/(double)n;
	dt = dt_us*1.0e-6;
	//printf("_pulse_shapedOC_2 duration of %f us split into %d steps of %f us\n",duration,n,dt*1.0e+6);

	// pre-hack block of rfmap and RotorModulated
	if (sim->rfmap != NULL) { // we have RotorModulated, initialize variables
		iz = sim->rfmap->loop[(wsp->rf_idx)*2+0]; // z coil coordinate index
		iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1];  // coil initial "phase" index
		Nphi = sim->rfmap->z[iz*(sim->rfmap->Nch+1)+sim->rfmap->Nch]; // number of phi elements for current z position
		steptimephi = sim->taur/Nphi;
		assert(steptimephi - duration >= -TINY); // limitation for fine-digitized shapes, steptime <=> duration
		iphi_shift = 0;
		phitime = 0.0;
	} // end of pre-hack block of rfmap and RotorModulated

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
		_setrfprop(sim,wsp);
		incr_OCmx_pos(wsp);
		/* like _pulse_simple */
		for (i=1;i<=n;i++) {
			ham_hamilton(sim,wsp);
			if (maxiter<0) {
				/* exact gradients calculated via exponential */
				fprintf(stderr,"Error: _pulse_shapedOCprops_2 - not implemented exact grads via expm\n");
				exit(1);
			} else {
				double scl;
				complx cscl;
				/* gradients improved with higher order terms */
				if (cH == NULL) {
					cH = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
					kom = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
					cdum = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
				}
				ham_hamrf_complex(cH, sim, wsp);
				//cm_print(cH, "Hamiltonian");
				blk_dm_multod(wsp->ham_blk,wsp->sumHrf,1.0);
				//blk_dm_print(wsp->ham_blk,"Ham total real");
				blk_cm_unit(wsp->dU);
				blk_prop_real(wsp->dU,wsp->ham_blk,dt,sim);
				blk_simtrans_zrot2(wsp->dU,wsp->sumUph);
				//blk_cm_print(wsp->dU,"Propagator");
				for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
					//printf("OCchanmap[%d] = %d\n",k,OCchanmap[k]);
					if (OCchanmap[k] < 0) { /* this shape is not active */
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = NULL;
							//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
						}
						if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] != NULL) {
							free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]);
							wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = NULL;
							//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
						}
						continue;
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
						//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
					}
					if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] == NULL) {
						wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
						//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
					}
					/* x channel */
					cm_zero(gr);
					cblas_daxpy(dim2,-1.0,wsp->chan_Ix[OCchanmap[k]]->data,1,(double*)(gr->data)+1,2); /* first order */
					cm_zero(kom);
					cblas_dcopy(dim2,wsp->chan_Ix[OCchanmap[k]]->data,1,(double*)(kom->data),2);
					//cm_print(kom,"channel Ix");
					iter=1;
					tol = 1000;
					scl = dt/2.0;
					while (iter < maxiter && tol > maxtol) {
						//printf("Iter %d, tol %g\n============\n",iter,tol);
						//cm_print(gr,"gr matrix");
						OC_commutator(cH,kom,iter,cdum, &tol);
						//cm_print(kom,"komutator");
						switch (iter % 4) {
						case 0: cscl.re = 0; cscl.im = -scl;
							break;
						case 1: cscl.re = scl; cscl.im = 0;
							break;
						case 2: cscl.re = 0; cscl.im = scl;
							break;
						case 3: cscl.re = -scl; cscl.im = 0;
							break;
						}
						cblas_zaxpy(dim2,&cscl, kom->data,1,gr->data,1);
						tol *= scl;
						scl *= (dt/((double)(iter)+2.0));
						iter++;
					}
					cm_multo_rev(gr,wsp->dU->m); /* gr is done here */
					//cm_print(gr,"READY grad");
					//exit(1);
					if (i == 1) {
						if (wsp->Uisunit) {
							cm_copy(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)],gr);
						} else {
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cnull,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]->data,dim);
						}
					} else {
						cm_multo_rev(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)],wsp->dU->m);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cunit,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]->data,dim);
					}
					/* y channel */
					cm_zero(gr);
					cblas_daxpy(dim2,1.0,wsp->chan_Iy[OCchanmap[k]]->data,1,(double*)(gr->data),2); /* first order */
					cm_zero(kom);
					cblas_dcopy(dim2,wsp->chan_Iy[OCchanmap[k]]->data,1,(double*)(kom->data)+1,2);
					//cm_print(kom,"channel Iy");
					iter=1;
					tol = 1000;
					scl = dt/2.0;
					while (iter < maxiter && tol > maxtol) {
						//printf("Iter %d, tol %g\n============\n",iter,tol);
						//cm_print(gr,"gr matrix");
						OC_commutator(cH,kom,iter,cdum, &tol);
						//cm_print(kom,"komutator");
						switch (iter % 4) {
						case 0: cscl.re = 0; cscl.im = -scl;
							break;
						case 1: cscl.re = scl; cscl.im = 0;
							break;
						case 2: cscl.re = 0; cscl.im = scl;
							break;
						case 3: cscl.re = -scl; cscl.im = 0;
							break;
						}
						cblas_zaxpy(dim2,&cscl, kom->data,1,gr->data,1);
						tol *= scl;
						scl *= (dt/((double)(iter)+2.0));
						iter++;
					}
					cm_multo_rev(gr,wsp->dU->m); /* gr is done here */
					//cm_print(gr,"READY grad");
					//exit(1);
					if (i == 1) {
						if (wsp->Uisunit) {
							cm_copy(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1],gr);
						} else {
							cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cnull,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]->data,dim);
						}
					} else {
						cm_multo_rev(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1],wsp->dU->m);
						cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,gr->data,dim,wsp->U->m->data,dim,&Cunit,wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]->data,dim);
					}
				}
			}
			/* update kumulativni propagator ve wsp->U */
			update_propagator(wsp->U, wsp->dU, sim, wsp);
			wsp->t += dt_us;
		} /* end for i over maxdt steps inside pulse step */
		store_OCprop(wsp);
		// !!!
		//blk_cm_print(wsp->OC_props[wsp->OC_mxpos],"U(1..k)");
		//cm_print(wsp->OC_deriv[wsp->OC_mxpos][0],"dU U(1..k)");

		// hack for RotorModulated
		if (sim->rfmap != NULL) { // we have RotorModulated
			for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
				if (OCchanmap[k] < 0) { /* this shape is not active */
					continue;
				}
			//-	int iz = sim->rfmap->loop[(wsp->rf_idx)*2+0];
			//-	int iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1];
				double bx, by;
				//rfmap_get(sim->rfmap,iz,iphi+j-1,OCchanmap[k]-1,&bx,&by);
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,OCchanmap[k]-1,&bx,&by);
				mat_complx *dUx = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)];   // dU/dwx
				mat_complx *dUy = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]; // dU/dwy
				mat_complx *dUx2 = cm_dup(dUx);
				cm_muld(dUx,bx);
				cm_multod(dUx,dUy,by);
				cm_muld(dUy,bx);
				cm_multod(dUy,dUx2,-by);
				free_complx_matrix(dUx2);
				//printf("GRAD: Thread %d: chan %d, slot %3d: phitime = %g ... %d %g %g\n",wsp->thread_id,OCchanmap[k],-1,phitime,j,bx,by);
			}
			// this needs to be last at the end of loop over j (shape elements)
			phitime += duration;  // steptime <=> duration
  		    if (steptimephi-phitime < 0.5*duration) {
			    phitime -= steptimephi;
				iphi_shift++;
			}
		} // end of hack

	} /* end for j over Nelem */

	if (cH != NULL) free_complx_matrix(cH);
	if (kom != NULL) free_complx_matrix(kom);
	if (cdum != NULL) free_complx_matrix(cdum);
	free_complx_matrix(gr);
	//printf("_pulse_shapedOC_2 done\n");

	/* mark current propagator as saved */
	wsp->OC_propstatus = 0;
}

/****
 * helper function for pulse_shaped in gradient mode and propagator optimization
 *        (here cumulative propagators are stored and no density matrices)
 ****/
void _pulse_shapedOCprops(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, double steptime)
{
  int i, j;
  char cd[128];

   if ( (wsp->Uisunit != 1) && (wsp->OC_propstatus == 1) ) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      set_OCmx_code(wsp,"P");
   }
   /* do pulsing and storing cumulative propagators */
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
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      if (j==1) {
         /* set the code only for the first memory slot */
         sprintf(cd,"G%dE%d",Nelem,Nch);
         strcat(cd,code);
         set_OCmx_code(wsp,cd);
      }

   }
   /* this is just to mark last slot with variable shape propagator*/
   set_OCmx_code(wsp,"G");
   /* mark current propagator as saved */
   wsp->OC_propstatus = 0;
}





/****
 * helper function for pulse_and_zgrad_shaped in gradient mode and state to state optimization
 *        (here step by step propagators are saved together with step-evolved
 *         density matrices)
 ****/
 void _pulse_and_zgrad_shapedOC(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, int zgrslot, double steptime)
{
  int i, j, Nchan;
  char cd[128];
  double zcoor;

  Nchan = sim->ss->nchan;
  if (wsp->inhom_offset) wsp->inhom_offset = double_vector(Nchan);
  zcoor = wsp->zcoor;

   if (wsp->Uisunit != 1) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      _evolve_with_prop(sim,wsp);
      store_OCdens(sim,wsp);
      _reset_prop(sim,wsp);
      set_OCmx_code(wsp,"P");
   }
   /* do pulsing and evolving, storing the results */
   for (j=1; j<=Nelem; j++) {
      get_chan_offset_ratios(sim->ss, zcoor*2.0*M_PI*ZgradShapes[zgrslot][j], wsp->inhom_offset);
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
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      _evolve_with_prop(sim,wsp);
      store_OCdens(sim,wsp);
      _reset_prop(sim,wsp);
      if (j==1)
         set_OCmx_code(wsp,"G");
   }
   /* set the code only for the last memory slot */
   sprintf(cd,"G%dE%d",Nelem,Nch);
   strcat(cd,code);
   set_OCmx_code(wsp,cd);
   /* clean up after z grad offsets */
   free_double_vector(wsp->inhom_offset);
   wsp->inhom_offset = NULL;
}


/****
 * helper function for pulse_and_zgrad_shaped in gradient mode and propagator optimization
 *        (here cummulative propagators are stored and no density matrices)
 ****/
void _pulse_and_zgrad_shapedOCprops(Sim_info *sim, Sim_wsp *wsp, char *code, int Nch, int Nelem, int *mask, int zgrslot, double steptime)
{
  int i, j, Nchan;
  char cd[128];
  double zcoor;

  Nchan = sim->ss->nchan;
  if (wsp->inhom_offset == NULL) wsp->inhom_offset = double_vector(Nchan);
  zcoor = wsp->zcoor;

   if ( (wsp->Uisunit != 1) && (wsp->OC_propstatus == 1) ) {
      /* there is some pending propagator, store it but set code for not take gradient */
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      set_OCmx_code(wsp,"P");
   }
   /* do pulsing and storing cumulative propagators */
   for (j=1; j<=Nelem; j++) {
      get_chan_offset_ratios(sim->ss, zcoor*2.0*M_PI*ZgradShapes[zgrslot][j], wsp->inhom_offset);
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
      _pulse(sim,wsp,steptime);
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      if (j==1) {
         /* set the code only for the first memory slot */
         sprintf(cd,"G%dE%d",Nelem,Nch);
         strcat(cd,code);
         set_OCmx_code(wsp,cd);
      }

   }
   /* this is just to mark last slot with variable shape propagator*/
   set_OCmx_code(wsp,"G");
   /* mark current propagator as saved */
   wsp->OC_propstatus = 0;
   /* clean up after z grad offsets */
   free_double_vector(wsp->inhom_offset);
   wsp->inhom_offset = NULL;
}

/*****
 * special purpose OC with 3-point rule propagation
 * called from tclPulseShaped3PointRuleRFmap()
 * Mostly copy-paste of _pulse_shapedOC_2,
 * grads via expansion in nested commutators, assumes rfmap modulations,
 * removed un-necessary options and tests, not fool-proofed!!
 */
extern void local_prop_complx(mat_complx *prop, mat_complx *ham, double dt); // defined in pulse.c
void _pulse_shapedOC_3pointrulerfmap(Sim_info *sim, Sim_wsp *wsp, int Nelem, int *OCchanmap, int *mask, double duration)
{
	int j, k, chan, idxb;
	double dt, am, ph, scl, steptimephi;
	int dim = wsp->ham_blk->dim;
	int dim2 = dim*dim;
	int Nsh = LEN(OCchanmap);
	mat_complx *kom = NULL, *cdum = NULL, *cdumptr;
	mat_complx *Hleft, *Hmiddle, *Hright;
	mat_complx *grxL, *gryL, *grxM, *gryM, *grxR, *gryR;
	complx cscl;
	// variables for rfmap
	int iz = 0, iphi = 0, iphi_shift = 0;
	double *bx, *by;
	int Nchan = sim->ss->nchan;
	// duration is steptime

	assert(wsp->U != NULL);
	assert(wsp->U->Nblocks == 1); /* no block_diag */
	assert(sim->rfmap != NULL);  // rfmap is defined, no check for proper timings though

	if (wsp->Uisunit != 1) {
		/* there is some pending propagator, store it but make no gradient */
		incr_OCmx_pos(wsp);
		store_OCprop(wsp);
		store_OCdens(sim,wsp); /* store sigma from previous step */
		_evolve_with_prop(sim,wsp); /* get sigma ready for next store */
		_reset_prop(sim,wsp);
	}

	int iter;
	int maxiter = -1; /* exact gradients calculated via exponential */
	double tol;
	double maxtol = 1e-6;
	if (fabs(OCpar.grad_level)<0.5) {
		/* gradient with higher order corrections up to grad_level accuracy */
		maxiter = 10;
		maxtol = fabs(OCpar.grad_level);
	} else if (OCpar.grad_level > 1.5){
		/* gradient with higher order corrections up to grad_level order */
		maxiter = (int)round(OCpar.grad_level);
		maxtol = 1e-10;
	}
	if (maxiter<0) {
		fprintf(stderr,"_pulse_shapedOC_3pointrulerfmap error: option for negative maxiter not implemented\n");
		exit(1);
	}
	//printf("pulse_shaped_OC: maxiter = %d, tol = %g\n",maxiter,maxtol);

	Hleft = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	Hmiddle = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	Hright = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	bx = (double*)malloc(3*Nchan*sizeof(double));
	by = (double*)malloc(3*Nchan*sizeof(double));
	if ( (bx == NULL) || (by == NULL) ) {
		fprintf(stderr,"_pulse_shapedOC_3pointrulerfmap error: bx/by allocation failed\n");
		exit(1);
	}
	if (wsp->Hcplx == NULL) {
		wsp->Hcplx = create_blk_mat_complx_copy2(wsp->ham_blk);
	}
	grxL = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	gryL = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	grxM = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	gryM = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	grxR = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
	gryR = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);

	dt = duration*1.0e-6; // element duration in seconds
	steptimephi = duration*0.5;  // assuming duration = 2*steptimephi !!!
	iz = sim->rfmap->loop[(wsp->rf_idx)*2+0]; // z coil coordinate index
	iphi = sim->rfmap->loop[(wsp->rf_idx)*2+1];  // coil initial "phase" index
	iphi_shift = 0;

	/* do pulsing, gradients, evolving, storing the results */
	for (j=1; j<=Nelem; j++) {  // main loop over all shape elements
		// preparing LEFT Hamiltonian
		for (chan=1; chan<=Nchan; chan++) {
			if (mask[chan] == -1) { // no rf applied on this channel
				_rf(wsp,chan,0.0);
				_ph(wsp,chan,0.0);
			} else {
				idxb = (chan-1)*3;
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,chan-1,bx+idxb,by+idxb); // get distortion coefs
				am = RFshapes[mask[chan]][j].ampl;  // original ampl/phase element
				ph = RFshapes[mask[chan]][j].phase;
				am *= sqrt(bx[idxb]*bx[idxb]+by[idxb]*by[idxb]);  // modify rf parameters at this element
				ph += RAD2DEG*atan2(by[idxb],bx[idxb]);
				_rf(wsp,chan,am);  // set them on
				_ph(wsp,chan,ph);
			}
		}
		ham_hamilton(sim, wsp); // creates interaction Hamiltonian in wsp->ham_blk
		ham_hamrf_complex(Hleft, sim, wsp); // creates final Hleft (interactions and rf complex matrix
			//printf("   ham left\n");
		// preparing MIDDLE Hamiltonian
		wsp->t += steptimephi; // move time forward
		iphi_shift++;         // move rotor forward
		for (chan=1; chan<=Nchan; chan++) {
			if (mask[chan] == -1) { // no rf applied on this channel
				_rf(wsp,chan,0.0);
				_ph(wsp,chan,0.0);
			} else {
				idxb = (chan-1)*3+1;
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,chan-1,bx+idxb,by+idxb); // get distortion coefs
				am = RFshapes[mask[chan]][j].ampl;  // original ampl/phase element
				ph = RFshapes[mask[chan]][j].phase;
				am *= sqrt(bx[idxb]*bx[idxb]+by[idxb]*by[idxb]);  // modify rf parameters at this element
				ph += RAD2DEG*atan2(by[idxb],bx[idxb]);
				_rf(wsp,chan,am);  // set them on
				_ph(wsp,chan,ph);
			}
		}
		ham_hamilton(sim, wsp); // creates interaction Hamiltonian in wsp->ham_blk
		ham_hamrf_complex(Hmiddle, sim, wsp); // creates final Hmiddle (interactions and rf complex matrix
			//printf("   ham middle\n");
		// preparing RIGHT Hamiltonian
		wsp->t += steptimephi; // move time forward
		iphi_shift++;         // move rotor forward
		for (chan=1; chan<=Nchan; chan++) {
			if (mask[chan] == -1) { // no rf applied on this channel
				_rf(wsp,chan,0.0);
				_ph(wsp,chan,0.0);
			} else {
				idxb = (chan-1)*3+2;
				rfmap_get(sim->rfmap,iz,iphi+iphi_shift,chan-1,bx+idxb,by+idxb); // get distortion coefs
				am = RFshapes[mask[chan]][j].ampl;  // original ampl/phase element
				ph = RFshapes[mask[chan]][j].phase;
				am *= sqrt(bx[idxb]*bx[idxb]+by[idxb]*by[idxb]);  // modify rf parameters at this element
				ph += RAD2DEG*atan2(by[idxb],bx[idxb]);
				_rf(wsp,chan,am);  // set them on
				_ph(wsp,chan,ph);
			}
		}
		ham_hamilton(sim, wsp); // creates interaction Hamiltonian in wsp->ham_blk
		ham_hamrf_complex(Hright, sim, wsp); // creates final Hright (interactions and rf complex matrix
			//printf("   ham right\n");
		// combine all Hamiltonians into the complex exponent stored in wsp->Hcplx
		blk_cm_zero(wsp->Hcplx);
		cm_multoc(wsp->Hcplx->m, Hleft, Complx(1.0/6.0,0));
		cm_multoc(wsp->Hcplx->m, Hmiddle, Complx(4.0/6.0,0));
		cm_multoc(wsp->Hcplx->m, Hright, Complx(1.0/6.0,0));
		cdum = cm_commutator(Hleft,Hright); // NOTE: cdum is allocated here!!!
		cm_multoc(wsp->Hcplx->m, cdum, Complx(0.0,dt/12.0));
		// calculate propagator of rf element into wsp->dU
		local_prop_complx(wsp->U->m,wsp->Hcplx->m,dt); // this is STEP propagator
		wsp->Uisunit = 0;
		// prepare pointer where to store OC_deriv
		incr_OCmx_pos(wsp);
			//printf("   propagator\n");
		for (k=1; k<=Nsh; k++) { /* loop over all gradshapes */
			//printf("OCchanmap[%d] = %d\n",k,OCchanmap[k]);
			if (OCchanmap[k] < 0) { /* this shape is not active, erase matrices if exist */
				if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] != NULL) {
					free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)]);
					wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = NULL;
					//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
				}
				if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] != NULL) {
					free_complx_matrix(wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1]);
					wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = NULL;
					//printf("free and NULL (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
				}
				continue;
			}
			// allocate matrices for active rf shapes
			if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] == NULL) {
				wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
				//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1));
			}
			if (wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] == NULL) {
				wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1] = complx_matrix(dim,dim,MAT_DENSE,0,wsp->ham_blk->basis);
				//printf("alloc (%d,%d)\n",wsp->OC_mxpos,2*(k-1)+1);
			}
			/* omega left x channel */
			  // prepare operator connected with omega left x rf parameter
			dm_copy2cm(wsp->chan_Ix[OCchanmap[k]], cdum); // cdum holds Ix
			kom = cm_commutator(cdum, Hright); // NOTE: kom is allocated here!!
			cm_multoc(cdum, kom, Complx(0,dt/2.0));
			cm_muld(cdum, 1.0/6.0); // operator is complete now
			cm_zero(grxL); // holds the gradient
			cm_multoc(grxL, cdum, Complx(0,-1));  // first order gradient
			cm_copy(kom, cdum); // kom holds the operator
			iter=1;
			tol = 1000;
			scl = dt/2.0;
			while (iter < maxiter && tol > maxtol) {
				OC_commutator(wsp->Hcplx->m,kom,iter,cdum, &tol); // kom is updated and tol is calculated
				switch (iter % 4) {
					case 0: cscl.re = 0; cscl.im = -scl;
						break;
					case 1: cscl.re = scl; cscl.im = 0;
						break;
					case 2: cscl.re = 0; cscl.im = scl;
						break;
					case 3: cscl.re = -scl; cscl.im = 0;
						break;
				}
				cblas_zaxpy(dim2,&cscl, kom->data,1,grxL->data,1);
				tol *= scl;
				scl *= (dt/((double)(iter)+2.0));
				iter++;
			}
			cm_multo_rev(grxL,wsp->U->m); // gr is done here (= prop*"commutator series")
			free_complx_matrix(kom);
				//printf("   --> omega left x, chan %d\n",OCchanmap[k]);
			/* omega left y channel */
			  // prepare operator connected with omega left x rf parameter
			dm_copy2cm_imag(wsp->chan_Iy[OCchanmap[k]], cdum); // cdum holds Iy
			kom = cm_commutator(cdum, Hright); // NOTE: kom is allocated here!!
			cm_multoc(cdum, kom, Complx(0,dt/2.0));
			cm_muld(cdum, 1.0/6.0); // operator is complete now
			cm_zero(gryL);  // holds the gradient
			cm_multoc(gryL, cdum, Complx(0,-1));  // first order gradient
			cm_copy(kom, cdum); // kom holds the operator
			iter=1;
			tol = 1000;
			scl = dt/2.0;
			while (iter < maxiter && tol > maxtol) {
				OC_commutator(wsp->Hcplx->m,kom,iter,cdum, &tol); // kom is updated and tol is calculated
				switch (iter % 4) {
					case 0: cscl.re = 0; cscl.im = -scl;
						break;
					case 1: cscl.re = scl; cscl.im = 0;
						break;
					case 2: cscl.re = 0; cscl.im = scl;
						break;
					case 3: cscl.re = -scl; cscl.im = 0;
						break;
				}
				cblas_zaxpy(dim2,&cscl, kom->data,1,gryL->data,1);
				tol *= scl;
				scl *= (dt/((double)(iter)+2.0));
				iter++;
			}
			cm_multo_rev(gryL,wsp->U->m); // gr is done here (= prop*"commutator series")
			free_complx_matrix(kom);
				//printf("   --> omega left y\n");
			/* omega MIDDLE x channel */
			  // prepare operator connected with omega left x rf parameter
			dm_copy2cm(wsp->chan_Ix[OCchanmap[k]], cdum); // cdum holds Ix
			cm_muld(cdum, 4.0/6.0); // operator is complete now
			cm_zero(grxM); // holds the gradient
			cm_multoc(grxM, cdum, Complx(0,-1));  // first order gradient
			kom = cm_dup(cdum); // kom allocated, holds the operator
			iter=1;
			tol = 1000;
			scl = dt/2.0;
			while (iter < maxiter && tol > maxtol) {
				OC_commutator(wsp->Hcplx->m,kom,iter,cdum, &tol); // kom is updated and tol is calculated
				switch (iter % 4) {
					case 0: cscl.re = 0; cscl.im = -scl;
						break;
					case 1: cscl.re = scl; cscl.im = 0;
						break;
					case 2: cscl.re = 0; cscl.im = scl;
						break;
					case 3: cscl.re = -scl; cscl.im = 0;
						break;
				}
				cblas_zaxpy(dim2,&cscl, kom->data,1,grxM->data,1);
				tol *= scl;
				scl *= (dt/((double)(iter)+2.0));
				iter++;
			}
			cm_multo_rev(grxM,wsp->U->m); // gr is done here (= prop*"commutator series")
			free_complx_matrix(kom);
				//printf("   --> omega middle x\n");
			/* omega MIDDLE y channel */
			  // prepare operator connected with omega left x rf parameter
			dm_copy2cm_imag(wsp->chan_Iy[OCchanmap[k]], cdum); // cdum holds Iy
			cm_muld(cdum, 4.0/6.0); // operator is complete now
			cm_zero(gryM);  // holds the gradient
			cm_multoc(gryM, cdum, Complx(0,-1));  // first order gradient
			kom = cm_dup(cdum); // kom allocated, holds the operator
			iter=1;
			tol = 1000;
			scl = dt/2.0;
			while (iter < maxiter && tol > maxtol) {
				OC_commutator(wsp->Hcplx->m,kom,iter,cdum, &tol); // kom is updated and tol is calculated
				switch (iter % 4) {
					case 0: cscl.re = 0; cscl.im = -scl;
						break;
					case 1: cscl.re = scl; cscl.im = 0;
						break;
					case 2: cscl.re = 0; cscl.im = scl;
						break;
					case 3: cscl.re = -scl; cscl.im = 0;
						break;
				}
				cblas_zaxpy(dim2,&cscl, kom->data,1,gryM->data,1);
				tol *= scl;
				scl *= (dt/((double)(iter)+2.0));
				iter++;
			}
			cm_multo_rev(gryM,wsp->U->m); // gr is done here (= prop*"commutator series")
			free_complx_matrix(kom);
				//printf("   --> omega middle y\n");
			/* omega RIGHT x channel */
			  // prepare operator connected with omega left x rf parameter
			dm_copy2cm(wsp->chan_Ix[OCchanmap[k]], cdum); // cdum holds Ix
			kom = cm_commutator(Hleft, cdum); // NOTE: kom is allocated here!!
			cm_multoc(cdum, kom, Complx(0,dt/2.0));
			cm_muld(cdum, 1.0/6.0); // operator is complete now
			cm_zero(grxR); // holds the gradient
			cm_multoc(grxR, cdum, Complx(0,-1));  // first order gradient
			cm_copy(kom, cdum); // kom holds the operator
			iter=1;
			tol = 1000;
			scl = dt/2.0;
			while (iter < maxiter && tol > maxtol) {
				OC_commutator(wsp->Hcplx->m,kom,iter,cdum, &tol); // kom is updated and tol is calculated
				switch (iter % 4) {
					case 0: cscl.re = 0; cscl.im = -scl;
						break;
					case 1: cscl.re = scl; cscl.im = 0;
						break;
					case 2: cscl.re = 0; cscl.im = scl;
						break;
					case 3: cscl.re = -scl; cscl.im = 0;
						break;
				}
				cblas_zaxpy(dim2,&cscl, kom->data,1,grxR->data,1);
				tol *= scl;
				scl *= (dt/((double)(iter)+2.0));
				iter++;
			}
			cm_multo_rev(grxR,wsp->U->m); // gr is done here (= prop*"commutator series")
			free_complx_matrix(kom);
				//printf("   --> omega right x\n");
			/* omega RIGHT y channel */
			  // prepare operator connected with omega left x rf parameter
			dm_copy2cm_imag(wsp->chan_Iy[OCchanmap[k]], cdum); // cdum holds Iy
			kom = cm_commutator(Hleft, cdum); // NOTE: kom is allocated here!!
			cm_multoc(cdum, kom, Complx(0,dt/2.0));
			cm_muld(cdum, 1.0/6.0); // operator is complete now
			cm_zero(gryR);  // holds the gradient
			cm_multoc(gryR, cdum, Complx(0,-1));  // first order gradient
			cm_copy(kom, cdum); // kom holds the operator
			iter=1;
			tol = 1000;
			scl = dt/2.0;
			while (iter < maxiter && tol > maxtol) {
				OC_commutator(wsp->Hcplx->m,kom,iter,cdum, &tol); // kom is updated and tol is calculated
				switch (iter % 4) {
					case 0: cscl.re = 0; cscl.im = -scl;
						break;
					case 1: cscl.re = scl; cscl.im = 0;
						break;
					case 2: cscl.re = 0; cscl.im = scl;
						break;
					case 3: cscl.re = -scl; cscl.im = 0;
						break;
				}
				cblas_zaxpy(dim2,&cscl, kom->data,1,gryR->data,1);
				tol *= scl;
				scl *= (dt/((double)(iter)+2.0));
				iter++;
			}
			cm_multo_rev(gryR,wsp->U->m); // gr is done here (= prop*"commutator series")
			free_complx_matrix(kom);
				//printf("   --> omega right y\n");
			// assemble TOTAL x gradient
			  // will use cdumptr as a pointer to OC_deriv
			cdumptr = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)];
			//printf("  = total grad x start\n");
			cm_zero(cdumptr);
			cm_multod(cdumptr, grxL, bx[(OCchanmap[k]-1)*3]);
			cm_multod(cdumptr, gryL, by[(OCchanmap[k]-1)*3]);
			cm_multod(cdumptr, grxM, bx[(OCchanmap[k]-1)*3+1]);
			cm_multod(cdumptr, gryM, by[(OCchanmap[k]-1)*3+1]);
			cm_multod(cdumptr, grxR, bx[(OCchanmap[k]-1)*3+2]);
			cm_multod(cdumptr, gryR, by[(OCchanmap[k]-1)*3+2]);
			// assemble TOTAL y gradient
			cdumptr = wsp->OC_deriv[wsp->OC_mxpos][2*(k-1)+1];
			//printf("  = total grad y start\n");
			cm_zero(cdumptr);
			cm_multod(cdumptr, grxL, -by[(OCchanmap[k]-1)*3]);
			cm_multod(cdumptr, gryL, bx[(OCchanmap[k]-1)*3]);
			cm_multod(cdumptr, grxM, -by[(OCchanmap[k]-1)*3+1]);
			cm_multod(cdumptr, gryM, bx[(OCchanmap[k]-1)*3+1]);
			cm_multod(cdumptr, grxR, -by[(OCchanmap[k]-1)*3+2]);
			cm_multod(cdumptr, gryR, bx[(OCchanmap[k]-1)*3+2]);
			//printf("  = total grad done\n");
		}/* loop over all gradshapes */
		free_complx_matrix(cdum); // cdum gets allocated for each shape element
		store_OCprop(wsp);
		store_OCdens(sim,wsp); /* sigma of previous step */
		_evolve_with_prop(sim,wsp); /* get sigma ready for the next step */
		_reset_prop(sim,wsp);
	} /* end for j over Nelem */

	// cleaning memory allocations!!!
	free_complx_matrix(Hleft);
	free_complx_matrix(Hright);
	free_complx_matrix(Hmiddle);
	free_complx_matrix(grxL);
	free_complx_matrix(gryL);
	free_complx_matrix(grxM);
	free_complx_matrix(gryM);
	free_complx_matrix(grxR);
	free_complx_matrix(gryR);
	free((char*)bx);
	free((char*)by);
}


/****
 * function calculating gradients for Hermitian state to state transfer
 ****/
 void gradOC_hermit(Sim_info *sim, Sim_wsp *wsp)
{
  int i,j,k,Nelem,Nch,idx,chan,dum, nn, NN;
  mat_complx **lam, *cmdum=NULL;
  complx cc;
  int *gridx, *chnl, *slt;
  char code, *dumstr, mxcd[128], *chptr;

  //printf("in gradOC_hermit\n");
  NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);

  /* check if there is pending propagator and save and evolve with it */
  if ( wsp->Uisunit != 1) {
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      _evolve_with_prop(sim,wsp);
      store_OCdens(sim,wsp);
      _reset_prop(sim,wsp);
      set_OCmx_code(wsp,"P");
   }

  /* initialize info where to store gradients */
  Nelem = OCpar.grad_shapes[0];
  gridx = int_vector(Nelem);
  dum=0;
  for (j=1;j<=Nelem;j++) {
     dum += RFshapes_len(OCpar.grad_shapes[j]);
     gridx[j]=dum;
  }
  /* check the size of fid */
  if ( dum > sim->ntot ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points (%d > %d)\n",dum,sim->ntot);
    exit(1);
  }
  /* optimistically set this and hope all in this function works... */
  wsp->curr_nsig = dum;

  // calculate also values of Phi (oc_acq_hermit)
  for ( i=0; i<NN; i++) {
	  complx c = cm_trace(wsp->sigma[i % sim->Nfstart],wsp->fdetect[i % sim->Nfdetect]);
	  wsp->OC_phivals[i+1].re += c.re;
  }


  /* final lambda is */
  lam = (mat_complx**)malloc(sim->Nfdetect*sizeof(mat_complx*));
  for (nn=0; nn<sim->Nfdetect; nn++) {
	  lam[nn] = cm_dup(wsp->fdetect[nn]);
  }

  i = wsp->OC_mxpos;
  while (i>0) {
	  strcpy(mxcd,wsp->OC_mxcode[i]);
	  code = mxcd[0];
	  /* printf("%d/ code is '%c'\n",i,code); */
	  if (code == 'G') {
		  /* calculate gradients */
		  /* printf("%d is G: %s\n",i,OCpar.mx_code[i]); */
		  /* disassemble the code */
		  dumstr = my_strtok_r(mxcd," ",&chptr);
		  if (sscanf(dumstr,"G%dE%d", &Nelem, &Nch) != 2) {
			  fprintf(stderr,"gradOC_hermit can't decompose first part of code\n");
			  exit(1);
		  }
		  /* printf("   -> %d elements",Nelem); */
		  chnl = int_vector(Nch);
		  slt = int_vector(Nch);
		  for (j=1; j<=Nch; j++) {
			  dumstr = my_strtok_r(chptr, " ",&chptr);
			  if ( sscanf(dumstr,"I%dC%d",&idx, &chan) != 2 ) {
				  fprintf(stderr,"gradOC_hermit - can't read grad code number %d\n",j);
				  exit(1);
			  }
			  chnl[j] = chan;
			  slt[j] = idx;
			  /* printf(", index %d on channel %d",idx, chan); */
		  }
		  /* printf("\n"); */

		  for (j=Nelem; j>0; j--) {
			  /* printf("     doing elem %d\n",i); */
			  if (j == Nelem) {
				  if (i != wsp->OC_mxpos) {
					  if (strncmp(wsp->OC_mxcode[i+1],"F",1) == 0) {
						  /* filter */
						  /* printf("     execute filter with %d\n",i+1); */
						  for (nn=0; nn<sim->Nfdetect; nn++)
							  lambda_filterOC(wsp,i+1,lam[nn]);
					  } else {
						  /* backevolve */
						  /* printf("     backevolve with %d\n",i+1); */
						  for (nn=0; nn<sim->Nfdetect; nn++)
							  blk_simtrans_adj(lam[nn],wsp->OC_props[i+1],sim);
					  }
				  }
			  } else {
				  /* printf("     backevolve with %d\n",i+1); */
				  for (nn=0; nn<sim->Nfdetect; nn++)
					  blk_simtrans_adj(lam[nn],wsp->OC_props[i+1],sim);
			  }
			  for (k=1; k<=Nch; k++) {
				  /* printf("     gradient calc. for shape idx %d on channel %d",slt[k],chnl[k]); */
				  for (nn=0; nn<NN; nn++) {
					  int is = nn % sim->Nfstart;
					  int id = nn % sim->Nfdetect;
					  if (cmdum == NULL) {
						  cmdum = dm_complx(wsp->chan_Ix[chnl[k]]);
					  } else {
						  dm_copy2cm(wsp->chan_Ix[chnl[k]], cmdum);
					  }
					  cm_multo_rev(cmdum, lam[id]);
					  cc = cm_trace(cmdum,wsp->OC_dens[i+is*MAXOCPROPS]);
					  wsp->fid[gridx[slt[k]]+nn*sim->ntot].re += 2.0*cc.im;
					  dm_copy2cm_imag(wsp->chan_Iy[chnl[k]], cmdum);
					  cm_multo_rev(cmdum,lam[id]);
					  cc = cm_trace(cmdum,wsp->OC_dens[i+is*MAXOCPROPS]);
					  wsp->fid[gridx[slt[k]]+nn*sim->ntot].im += 2.0*cc.im;
					  /* printf(", stored in fid[%d]\n",gridx[slt[k]]); */
				  }
				  (gridx[slt[k]])--;
			  }
			  i--;
		  }
		  free_int_vector(chnl);
		  free_int_vector(slt);

	  } else {
		  if (i != wsp->OC_mxpos) {
			  /* code is either P or F, just evolve or filter */
			  /* printf("%d is not G: %s\n",i,OCpar.mx_code[i]); */
			  if (strncmp(wsp->OC_mxcode[i+1],"F",1) == 0) {
				  /* filter */
				  /* printf("     execute filter with %d\n",i+1); */
				  for (nn=0; nn<sim->Nfdetect; nn++)
					  lambda_filterOC(wsp,i+1,lam[nn]);
			  } else {
				  /* backevolve */
				  /* printf("     backevolve with %d\n",i+1); */
				  for (nn=0; nn<sim->Nfdetect; nn++)
					  blk_simtrans_adj(lam[nn],wsp->OC_props[i+1],sim);
			  }
		  }
		  i--;
	  }

  } /*end while */

  if (cmdum != NULL) free_complx_matrix(cmdum);
  for (nn=0; nn<sim->Nfdetect; nn++) free_complx_matrix(lam[nn]);
  free(lam);
  free_int_vector(gridx);
}
/*****
 * function calculating gradients for Hermitian state to state transfer
 *   --> with advanced gradients
 *   NOTE: does not allow usage of filter command !!!
 *****/
void gradOC_hermit_2(Sim_info *sim, Sim_wsp *wsp)
{
  int i,j,Nsh,dum,dim, nn, NN;
  mat_complx **lam, *cm1, *cm2;
  complx cc;
  int *gridx;
  //printf("in gradOC_hermit_2\n");

  NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);
  /* check if there is pending propagator and save and evolve with it */
  if ( wsp->Uisunit != 1) {
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      store_OCdens(sim,wsp); /* sigma of previous step */
      _evolve_with_prop(sim,wsp); /* next step probably will not happen but used to calculate Phi */
      _reset_prop(sim,wsp);
   }

  /* initialize info where to store gradients */
  Nsh = OCpar.grad_shapes[0];
  gridx = int_vector(Nsh);
  dum=0;
  for (j=1;j<=Nsh;j++) {
     dum += RFshapes_len(OCpar.grad_shapes[j]);
     gridx[j]=dum;
     //printf("\t gridx[%d] = %d\n",j, dum);
  }
  /* check the size of fid */
  if ( dum > sim->ntot ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points (%d > %d)\n",dum,sim->ntot);
    exit(1);
  }
  /* optimistically set this and hope all in this function works... */
  wsp->curr_nsig = dum;

  // calculate also values of Phi (oc_acq_hermit)
  for ( i=0; i<NN; i++) {
	  complx c = cm_trace(wsp->sigma[i % sim->Nfstart],wsp->fdetect[i % sim->Nfdetect]);
	  wsp->OC_phivals[i+1].re += c.re;
	  //cm_print(wsp->sigma[i % sim->Nfstart],"sigma");
	  //printf("\t--> oc_acq_hermit form grad: %.4f\n",c.re);
  }

  /* final lambda is */
  lam = (mat_complx**)malloc(sim->Nfdetect*sizeof(mat_complx*));
  for (nn=0; nn<sim->Nfdetect; nn++)
	  lam[nn] = cm_dup(wsp->fdetect[nn]);
  dim = lam[0]->row;
  cm1 = complx_matrix(dim,dim,MAT_DENSE,0,lam[0]->basis);
  cm2 = complx_matrix(dim,dim,MAT_DENSE,0,lam[0]->basis);

  i = wsp->OC_mxpos;
  while (i>0) {
	  //cm_print(lam,"lam");
	  for (j=1; j<=Nsh; j++) {
		  if (wsp->OC_deriv[i][2*(j-1)] == NULL) continue;
		  assert(wsp->OC_deriv[i][2*(j-1)] != NULL);
		  assert(wsp->OC_deriv[i][2*(j-1)+1] != NULL);
		  for (nn=0; nn<NN; nn++) {
			  int is = nn % sim->Nfstart;
			  int id = nn % sim->Nfdetect;
			  if (wsp->OC_dens[i+is*MAXOCPROPS]->type != MAT_DENSE) cm_dense_full(wsp->OC_dens[i+is*MAXOCPROPS]);
			  	  //if (i==1) printf("\t %d: %d; gridx[%d] = %d",i,j,j,gridx[j]);
			  	  //if (i==1) cm_print(wsp->OC_dens[i],"dens");
			  	  //if (i==1) cm_print(wsp->OC_props[i]->m,"prop");
			  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,&Cunit,wsp->OC_dens[i+is*MAXOCPROPS]->data,dim,wsp->OC_props[i]->m->data,dim,&Cnull,cm1->data,dim);
			  	  //if (i==1) printf(" A");
			  	  //cm_print(wsp->OC_dens[i],"rho");
			  	  //cm_print(wsp->OC_props[i]->m,"prop Uk");
			  	  //cm_print(cm1,"cm1");
			  /* x channel */
			  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,wsp->OC_deriv[i][2*(j-1)]->data,dim,cm1->data,dim,&Cnull,cm2->data,dim);
			  	  //if (i==1) printf(" B");
			  	  //cm_print(wsp->OC_deriv[i][2*(j-1)],"deriv Uk");
			  	  //cm_print(cm2,"cm2");
			  cc = cm_trace(lam[id],cm2);
			  	  //printf("X %d shape %d: cc = %g, %g\n",i,j,cc.re, cc.im);
			  	  //if (i==1) printf(" uf");
			  wsp->fid[gridx[j]+nn*sim->ntot].re += 2.0*cc.re;
			  	  //if (i==1) printf(" gr");
			  /* y channel */
			  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,wsp->OC_deriv[i][2*(j-1)+1]->data,dim,cm1->data,dim,&Cnull,cm2->data,dim);
			  cc = cm_trace(lam[id],cm2);
			  	  //printf("Y %d shape %d: cc = %g, %g\n",i,j,cc.re, cc.im);
			  	  //if (i==1) printf(" bu");
			  wsp->fid[gridx[j]+nn*sim->ntot].im += 2.0*cc.re;
		  }
		  //if (i==1) printf(" ha\n");
		  (gridx[j])--;
	  }
	  for (nn=0; nn<sim->Nfdetect; nn++)
		  blk_simtrans_adj(lam[nn],wsp->OC_props[i],sim);
	  i--;
  } /*end while */


  free_complx_matrix(cm1);
  free_complx_matrix(cm2);
  for (nn=0; nn<sim->Nfdetect; nn++)
	  free_complx_matrix(lam[nn]);
  free(lam);
  free_int_vector(gridx);
  //printf("_gradOC_hermit_2 done\n");
}


/****
 * function calculating gradients for non-Hermitian state to state transfer
 ****/
 void gradOC_nonhermit(Sim_info *sim, Sim_wsp *wsp)
{
  int i,j,k,Nelem,Nch,idx,chan,dum, nn, NN;
  mat_complx **lam, *tmpM=NULL, *tmpMM=NULL;
  complx cc, *cc1;
  int *gridx, *chnl, *slt;
  char code, *dumstr, mxcd[128], *chptr;
  
  /* printf("in gradOC_nonhermit\n"); */
  NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);

  /* check if there is pending propagator and save and evolve with it */
  if ( wsp->Uisunit != 1) {
      incr_OCmx_pos(wsp);
      store_OCprop(wsp);
      _evolve_with_prop(sim,wsp);
      store_OCdens(sim,wsp);
      _reset_prop(sim,wsp);
      set_OCmx_code(wsp,"P");
   }
   
  /* initialize info where to store gradients */
  Nelem = OCpar.grad_shapes[0];
  gridx = int_vector(Nelem);
  dum=0;
  for (j=1;j<=Nelem;j++) {
     dum += RFshapes_len(OCpar.grad_shapes[j]);
     gridx[j] = dum;
  }
  /* check the size of fid */
  if ( dum > sim->ntot ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points\n");
    exit(1);
  }
  /* optimistically set this and hope all in this function works... */
  wsp->curr_nsig = dum;

  // calculate also Phi values (oc_acq_nonhermit)
  for (i=0; i<NN; i++) {
	  complx c = cm_trace_adjoint(wsp->fdetect[i % sim->Nfdetect],wsp->sigma[i % sim->Nfstart]);
	  wsp->OC_phivals[i+1].re += c.re*c.re+c.im*c.im;
  }

  /* here, lam is holding adjoint of lambda throughout the calculations */
  lam = (mat_complx**)malloc(sim->Nfdetect*sizeof(mat_complx*));
  for (nn=0; nn<sim->Nfdetect; nn++)
	  lam[nn] = cm_adjoint(wsp->fdetect[nn]);
  /* prepare constant complex term Tr[rho_N+ lam_N] */
  cc1 = (complx*)malloc(NN*sizeof(complx));
  for (nn=0; nn<NN; nn++)
	  cc1[nn] = cm_trace_adjoint(wsp->OC_dens[wsp->OC_mxpos+MAXOCPROPS*(nn % sim->Nfstart)],wsp->fdetect[nn % sim->Nfdetect]);

  i = wsp->OC_mxpos;
  while (i>0) {
	  strcpy(mxcd,wsp->OC_mxcode[i]);
	  code = mxcd[0];
	  /* printf("%d/ code is '%c'\n",i,code); */
	  if (code == 'G') {
		  /* calculate gradients */
		  /* printf("%d is G: %s\n",i,OCpar.mx_code[i]); */
		  /* disassemble the code */
		  dumstr = my_strtok_r(mxcd," ",&chptr);
		  if (sscanf(dumstr,"G%dE%d", &Nelem, &Nch) != 2) {
			  fprintf(stderr,"gradOC_nonhermit can't decompose first part of code\n");
			  exit(1);
		  }
		  /* printf("   -> %d elements",Nelem); */
		  chnl = int_vector(Nch);
		  slt = int_vector(Nch);
		  for (j=1; j<=Nch; j++) {
			  dumstr = my_strtok_r(chptr, " ", &chptr);
			  if ( sscanf(dumstr,"I%dC%d",&idx, &chan) != 2 ) {
				  fprintf(stderr,"gradOC_hermit - can't read grad code number %d\n",j);
				  exit(1);
			  }
			  chnl[j] = chan;
			  slt[j] = idx;
			  /* printf(", index %d on channel %d",idx, chan); */
		  }
		  /* printf("\n"); */
	

		  /* start loop over elements in shapes */
		  for (j=Nelem; j>0; j--) {
			  /* printf("     doing elem %d\n",i); */
			  if (j == Nelem) {
				  if (i != wsp->OC_mxpos) {
					  if (strncmp(wsp->OC_mxcode[i+1],"F",1) == 0) {
						  /* filter */
						  /* printf("     execute filter with %d\n",i+1); */
						  for (nn=0; nn<sim->Nfdetect; nn++)
							  lambda_filterOC(wsp,i+1,lam[nn]);
					  } else {
						  /* backevolve */
						  /* printf("     backevolve with %d\n",i+1); */
						  for (nn=0; nn<sim->Nfdetect; nn++)
							  blk_simtrans_adj(lam[nn],wsp->OC_props[i+1],sim);
					  }
				  }
			  } else {
				  /* printf("     backevolve with %d\n",i+1); */
				  for (nn=0; nn<sim->Nfdetect; nn++)
					  blk_simtrans_adj(lam[nn],wsp->OC_props[i+1],sim);
			  }
			  /* loop over channels */
			  for (k=1; k<=Nch; k++) {
				  /* printf("     gradient calc. for shape idx %d on channel %d",slt[k],chnl[k]); */
				  /* cc = Cmul( m_trace(puls->chan_Ix[chnl[k]], puls->tmp), cc1); */
				  for (nn=0; nn<NN; nn++) {
					int is = nn % sim->Nfstart;
					int id = nn % sim->Nfdetect;
					/* prepare constant matrix [rho_j, lam_j+], note that lam is already adjoint! */
					tmpM = cm_commutator(wsp->OC_dens[i+is*MAXOCPROPS], lam[id]);
					if (tmpMM == NULL) {
						tmpMM = dm_complx(wsp->chan_Ix[chnl[k]]);
					} else {
						dm_copy2cm(wsp->chan_Ix[chnl[k]],tmpMM);
					}
					cc = Cmul( cm_trace_adjoint(tmpMM,tmpM), cc1[nn]);
					/* I use this trace since it is faster, Ix is hermitian so adjoint does not change it */
					wsp->fid[gridx[slt[k]]+nn*sim->ntot].re += 2.0*cc.im;
					/* cc = Cmul( m_trace(puls->chan_Iy[chnl[k]], puls->tmp), cc1); */
					dm_copy2cm_imag(wsp->chan_Iy[chnl[k]],tmpMM);
					cc = Cmul( cm_trace_adjoint(tmpMM,tmpM), cc1[nn]);
					/* I use this trace since it is faster, Iy is hermitian so adjoint does not change it */
					//gr.im = 2.0*cc.im;
					wsp->fid[gridx[slt[k]]+nn*sim->ntot].im += 2.0*cc.im;
					/* printf(", stored in fid[%d]\n",gridx[slt[k]]); */
					//wsp->fid[gridx[slt[k]]] = gr;
					free_complx_matrix(tmpM);
				  }
				  (gridx[slt[k]])--;
			  }
			  i--;
		  }
		  free_int_vector(chnl);
		  free_int_vector(slt);
	  } else {
		  if (i != wsp->OC_mxpos) {
			  /* code is either P or F, just evolve or filter */
			  /* printf("%d is not G: %s\n",i,OCpar.mx_code[i]); */
			  if (strncmp(wsp->OC_mxcode[i+1],"F",1) == 0) {
				  /* filter */
				  /* printf("     execute filter with %d\n",i+1); */
				  for (nn=0; nn<sim->Nfdetect; nn++)
					  lambda_filterOC(wsp,i+1,lam[nn]);
			  } else {
				  /* backevolve */
				  /* printf("     backevolve with %d\n",i+1); */
				  for (nn=0; nn<sim->Nfdetect; nn++)
					  blk_simtrans_adj(lam[nn],wsp->OC_props[i+1],sim);
			  }
		  }
		  i--;
	  }
  } /*end while */
  
  for (nn=0; nn<sim->Nfdetect; nn++)
	  free_complx_matrix(lam[nn]);
  free(lam);
  free(cc1);
  if (tmpMM != NULL) free_complx_matrix(tmpMM);
  free_int_vector(gridx);
}

 /*****
  * function calculating gradients for non-Hermitian state to state transfer
  *   --> with advanced gradients
  *   NOTE: does not allow usage of filter command !!!
  *****/
 void gradOC_nonhermit_2(Sim_info *sim, Sim_wsp *wsp)
 {
   int i,j,Nsh,dum,dim, nn, NN;
   mat_complx **lam, *cm1, *cm2;
   complx cc, *cc2, grx, gry;
   int *gridx;

   NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);

   /* check if there is pending propagator and save and evolve with it */
   if ( wsp->Uisunit != 1) {
       incr_OCmx_pos(wsp);
       store_OCprop(wsp);
       store_OCdens(sim,wsp); /* sigma of previous step */
       _evolve_with_prop(sim,wsp); /* next step probably will not happen */
       _reset_prop(sim,wsp);
    }

   /* initialize info where to store gradients */
   Nsh = OCpar.grad_shapes[0];
   gridx = int_vector(Nsh);
   dum=0;
   for (j=1;j<=Nsh;j++) {
      dum += RFshapes_len(OCpar.grad_shapes[j]);
      gridx[j]=dum;
      //printf("\t gridx[%d] = %d\n",j, dum);
   }
   /* check the size of fid */
   if ( dum > sim->ntot ) {
     fprintf(stderr,"error: gradient function detected overflow in fid points\n");
     exit(1);
   }
   /* optimistically set this and hope all in this function works... */
   wsp->curr_nsig = dum;

   // calculate also Phi values (oc_acq_nonhermit)
   for (i=0; i<NN; i++) {
 	  complx c = cm_trace_adjoint(wsp->fdetect[i % sim->Nfdetect],wsp->sigma[i % sim->Nfstart]);
 	  wsp->OC_phivals[i+1].re += c.re*c.re+c.im*c.im;
   }

   /* final lambda is adjoint as it is not Hermitian */
   lam = (mat_complx**)malloc(sim->Nfdetect*sizeof(mat_complx*));
   for (nn=0; nn<sim->Nfdetect; nn++)
	   lam[nn] = cm_adjoint(wsp->fdetect[nn]);
   dim = lam[0]->row;
   cm1 = complx_matrix(dim,dim,MAT_DENSE,0,lam[0]->basis);
   cm2 = complx_matrix(dim,dim,MAT_DENSE,0,lam[0]->basis);
   /* constant in front of the gradient */
   cc2 = (complx*)malloc(NN*sizeof(complx));
   for (nn=0; nn<NN; nn++)
	   cc2[nn] = cm_trace_adjoint(wsp->sigma[nn % sim->Nfstart],wsp->fdetect[nn % sim->Nfdetect]);

   i = wsp->OC_mxpos;
   while (i>0) {
 	  //cm_print(lam,"lam");
 	  for (j=1; j<=Nsh; j++) {
 		  if (wsp->OC_deriv[i][2*(j-1)] == NULL) continue;
 		  assert(wsp->OC_deriv[i][2*(j-1)] != NULL);
 		  assert(wsp->OC_deriv[i][2*(j-1)+1] != NULL);
 		 for (nn=0; nn<NN; nn++) {
 			int is = nn % sim->Nfstart;
 			int id = nn % sim->Nfdetect;
 			if (wsp->OC_dens[i+is*MAXOCPROPS]->type != MAT_DENSE) cm_dense_full(wsp->OC_dens[i+is*MAXOCPROPS]);
 				//if (i==1) printf("\t %d: %d; gridx[%d] = %d",i,j,j,gridx[j]);
 				//if (i==1) cm_print(wsp->OC_dens[i],"dens");
 				//if (i==1) cm_print(wsp->OC_props[i]->m,"prop");
 			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,&Cunit,wsp->OC_dens[i+is*MAXOCPROPS]->data,dim,wsp->OC_props[i]->m->data,dim,&Cnull,cm1->data,dim);
 				//if (i==1) printf(" A");
 				//cm_print(wsp->OC_dens[i],"rho");
 				//cm_print(wsp->OC_props[i]->m,"prop Uk");
 				//cm_print(cm1,"cm1");
 			/* x channel - first part */
 			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,wsp->OC_deriv[i][2*(j-1)]->data,dim,cm1->data,dim,&Cnull,cm2->data,dim);
 			grx = cm_trace(lam[id],cm2);
 			/* y channel - first part  */
 			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,wsp->OC_deriv[i][2*(j-1)+1]->data,dim,cm1->data,dim,&Cnull,cm2->data,dim);
 			gry = cm_trace(lam[id],cm2);
 			/* preparing second part */
 			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,&Cunit,wsp->OC_props[i]->m->data,dim,wsp->OC_dens[i+is*MAXOCPROPS]->data,dim,&Cnull,cm1->data,dim);
 			/* x channel - second part */
 			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,&Cunit,cm1->data,dim,wsp->OC_deriv[i][2*(j-1)]->data,dim,&Cnull,cm2->data,dim);
 			cc = cm_trace(lam[id],cm2);
 			grx.re += cc.re; grx.im += cc.im;
 			wsp->fid[gridx[j]+nn*sim->ntot].re += 2.0*(cc2[nn].re*grx.re-cc2[nn].im*grx.im);
 			/* y channel - second part  */
 			cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,&Cunit,cm1->data,dim,wsp->OC_deriv[i][2*(j-1)+1]->data,dim,&Cnull,cm2->data,dim);
 			cc = cm_trace(lam[id],cm2);
 			gry.re += cc.re; gry.im += cc.im;
 			wsp->fid[gridx[j]+nn*sim->ntot].im += 2.0*(cc2[nn].re*gry.re-cc2[nn].im*gry.im);
 		 }
 		 (gridx[j])--;
 	  }
 	  for (nn=0; nn<sim->Nfdetect; nn++)
 		  blk_simtrans_adj(lam[nn],wsp->OC_props[i],sim);
 	  i--;
   } /*end while */


   free_complx_matrix(cm1);
   free_complx_matrix(cm2);
   for (nn=0; nn<sim->Nfdetect; nn++)
	   free_complx_matrix(lam[nn]);
   free(lam);
   free(cc2);
   free_int_vector(gridx);
   //printf("_gradOC_hermit_2 done\n");
 }

/****
 * function for calculation gradients for propagator optimization
 *  original grads ala Khaneja
 ****/
void gradOC_prop(Sim_info *sim, Sim_wsp *wsp, int ud)
{
  int i,j,k,Nelem,Nch,idx,chan,dum,N;
  mat_complx *pr, *cmx=NULL;
  complx cc, cc1;
  int *gridx, *chnl, *slt;
  char code, *dumstr, mxcd[128], *chptr;
  
  /* printf("in gradOC_prop\n"); */
  N = sim->matdim;

  if ( (wsp->Uisunit != 1) && (wsp->OC_propstatus == 1) ) {
     /* there is some pending propagator, store it but set code for not take gradient */
     incr_OCmx_pos(wsp);
     store_OCprop(wsp);
     set_OCmx_code(wsp,"P");
  }

  /* initialize info where to store gradients */
  Nelem = OCpar.grad_shapes[0];
  gridx = int_vector(Nelem);
  dum=1;
  for (j=1;j<=Nelem;j++) {
     gridx[j]=dum;
     dum += RFshapes_len(OCpar.grad_shapes[j]);
  }
  /* check the size of fid */
  if ( dum-1 > LEN(wsp->fid) ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points\n");
    exit(1);
  }
  /* Optimistically set this and hope all in this function works... */
  wsp->curr_nsig = dum-1;

  /* check bases of all propagators and matrices*/
  /* props are cumulative, so put everything into basis of the final one */
  // potrebuju: chan_Ix[], chan_Iy[], STO[ud] a OC_props[]

  //blk_cm_print(wsp->STO[ud],"gradOC_prop Ud");
  //blk_cm_print(wsp->OC_props[wsp->OC_mxpos],"gradOC_prop Shape total prop");
  /* pr is holding U_d+ U_final */
  // NOT EFFICIENT CODE !!!!
  blk_mat_complx *blkdum = blk_cm_adjoint(wsp->STO[ud]);
  blk_mat_complx *blkpr = blk_cm_mul(blkdum,wsp->OC_props[wsp->OC_mxpos]);
  free_blk_mat_complx(blkdum);
  pr = complx_matrix(N,N,(sim->sparse) ? MAT_SPARSE : MAT_DENSE, 1, blkpr->basis);
  blk_cm_copy_2(pr,blkpr);
  free_blk_mat_complx(blkpr);
  /* cc1 is constant Tr { U_final+ U_d } */
  //cc1 = cm_trace_adjoint( wsp->OC_props[wsp->OC_mxpos], wsp->STO[ud]);
  cc1 = cm_true_trace(pr);
  cc1.im = -cc1.im;
  
  i=1;
  while (i <= wsp->OC_mxpos) {
	  strcpy(mxcd,wsp->OC_mxcode[i]);
	  code = mxcd[0];
	  DEBUGPRINT("%d/ code is '%c'\n",i,code);
	  //printf("full code  = '%s'\n",mxcd);
	  if (code == 'G') {
		/* calculate gradients */
		  /* printf("%d is G: %s\n",i,OCpar.mx_code[i]); */
		  /* disassemble the code */
		  dumstr = my_strtok_r(mxcd," ",&chptr);
		  if (sscanf(dumstr,"G%dE%d", &Nelem, &Nch) != 2) {
			  fprintf(stderr,"gradOC_prop can't decompose first part of code\n");
			  exit(1);
		  }
		  /* printf("   -> %d elements",Nelem); */
		  chnl = int_vector(Nch);
		  slt = int_vector(Nch);
		  for (j=1; j<=Nch; j++) {
			  dumstr = my_strtok_r(chptr, " ", &chptr);
			  if ( sscanf(dumstr,"I%dC%d",&idx, &chan) != 2 ) {
				  fprintf(stderr,"gradOC_prop - can't read grad code number %d\n",j);
				  exit(1);
			  }
			  chnl[j] = chan;
			  slt[j] = idx;
			  /* printf(", index %d on channel %d",idx, chan); */
		  }
		  /* printf("\n"); */

		  cmx = complx_matrix(sim->matdim,sim->matdim,(sim->sparse) ? MAT_SPARSE : MAT_DENSE, sim->matdim,0);
		  /* start loop over elements in shapes */
		  for (j=1; j<=Nelem; j++) {
			  /* printf("     doing elem %d\n",i); */
			  /* loop over channels */
			  for (k=1; k<=Nch; k++) {
				  /* printf("     gradient calc. for shape idx %d on channel %d",slt[k],chnl[k]); */
				  dm_copy2cm(wsp->chan_Ix[chnl[k]],cmx); // copy real to real part of complex matrix
				  blk_simtrans_adj(cmx,wsp->OC_props[i],sim);
				  cc = Cmul( cm_trace(pr, cmx), cc1);
				  //gr.re = 2.0*cc.im;
				  wsp->fid[gridx[slt[k]]].re += 2.0*cc.im;
				  dm_copy2cm_imag(wsp->chan_Iy[chnl[k]],cmx); // copy real to imag part of complex matrix
				  blk_simtrans_adj(cmx,wsp->OC_props[i],sim);
				  cc = Cmul( cm_trace(pr, cmx), cc1);
				  //gr.im = 2.0*cc.im;
				  wsp->fid[gridx[slt[k]]].im += 2.0*cc.im;
				  /* printf(", stored in fid[%d]\n",gridx[slt[k]]); */
				  //wsp->fid[gridx[slt[k]]] = gr;
				  (gridx[slt[k]])++;
			  }
			  i++;
		  }
		  free_int_vector(chnl);
		  free_int_vector(slt);
	  } else {
		  i++;
	  }
  } /*end while */
  
  free_complx_matrix(pr);
  free_int_vector(gridx);
  if (cmx) free_complx_matrix(cmx);
}

/****
 * function for calculating gradients for propagator optimization
 *  improved grads ala Kuprov
 ****/
void gradOC_prop_2(Sim_info *sim, Sim_wsp *wsp, int ud)
{
  int i,j,Nelem,dum,N;
  mat_complx *pr, *cmx;
  complx cc, cc1;
  int *gridx;
  
  /* printf("in gradOC_prop\n"); */
  
  if ( (wsp->Uisunit != 1) && (wsp->OC_propstatus == 1) ) {
     /* there is some pending propagator, store it but set code for not take gradient */
     incr_OCmx_pos(wsp);
     store_OCprop(wsp);
  }

  /* initialize info where to store gradients */
  N = sim->matdim;
  Nelem = OCpar.grad_shapes[0];
  gridx = int_vector(Nelem);
  dum=1;
  for (j=1;j<=Nelem;j++) {
     gridx[j]=dum;
     dum += RFshapes_len(OCpar.grad_shapes[j]);
  }
  /* check the size of fid */
  if ( dum-1 > LEN(wsp->fid) ) {
    fprintf(stderr,"error: gradient function detected overflow in fid points\n");
    exit(1);
  }
  /* Optimistically set this and hope all in this function works... */
  wsp->curr_nsig = dum-1;

  /* assumes blockdiag option NOT used and dense matrices */
  /* pr is holding U_d+ U_final */
  pr = complx_matrix(N,N,MAT_DENSE,1,wsp->U->basis);
  cm_dense_full(wsp->STO[ud]->m);
  cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,N,N,N,&Cunit,wsp->STO[ud]->m->data,N,wsp->OC_props[wsp->OC_mxpos]->m->data,N,&Cnull,pr->data,N);
  /* cc1 is constant Tr { U_final+ U_d } */
  cc1 = cm_true_trace(pr);
  cc1.im = -cc1.im;
  cmx = complx_matrix(N,N,MAT_DENSE,1,wsp->U->basis);

  i=1;
  while (i <= wsp->OC_mxpos) {
	  for (j=1; j<=Nelem; j++) {
		  if (wsp->OC_deriv[i][2*(j-1)] == NULL) continue;
		  assert(wsp->OC_deriv[i][2*(j-1)] != NULL);
		  assert(wsp->OC_deriv[i][2*(j-1)+1] != NULL);
		  /* x channel */
		  cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,N,N,N,&Cunit,wsp->OC_props[i]->m->data,N,wsp->OC_deriv[i][2*(j-1)]->data,N,&Cnull,cmx->data,N);
		  cc = cm_trace(pr,cmx);
		  wsp->fid[gridx[j]].re += 2.0*(cc.re*cc1.re-cc.im*cc1.im);
		  /* y channel */
		  cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,N,N,N,&Cunit,wsp->OC_props[i]->m->data,N,wsp->OC_deriv[i][2*(j-1)+1]->data,N,&Cnull,cmx->data,N);
		  cc = cm_trace(pr,cmx);
		  wsp->fid[gridx[j]].im += 2.0*(cc.re*cc1.re-cc.im*cc1.im);
		  (gridx[j])++;
	  }
	  i++;
  } /*end while */
  
  free_complx_matrix(pr);
  free_int_vector(gridx);
  free_complx_matrix(cmx);
}



/**************************************************************
 * ZT: functions for Optimal Control by Conjugated Gradients  *
 **************************************************************/

/****
 * ZT: function to evaluate target function with current RFshapes
 *     this value should be maximized
 ****/
 double evaluate_target_function(Tcl_Interp* interp)
{
  double tfval, pen;
  int i,ii;
  Tcl_Obj* objptr;  

  /* reset grad mode to calculate target function */
  OCpar.gradmode = 0;
  
  if (Tcl_EvalEx( (Tcl_Interp*)interp, "target_function", -1, TCL_EVAL_GLOBAL) != TCL_OK) {
    fprintf(stderr,"error: Unable to execute command 'target_function'\n%s\n",Tcl_GetStringResult(interp));
    exit(1);
  }
  objptr = Tcl_GetObjResult( (Tcl_Interp*)interp);
  if ( Tcl_GetDoubleFromObj( (Tcl_Interp*)interp, objptr, &tfval) != TCL_OK) {
    fprintf(stderr,"error in evaluate_target_function: Unable to convert result into double");
    exit(1);
  }

  // adding penalty terms
  pen = 0.0;
  if (OCpar.ispenalty) {
	  for (i=1; i<=OCpar.var_shapes[0]; i++) {
		  if (OCpar.var_shapes_penalty_order[i] != 0) {
			  double peni = 0.0;
			  double lim = OCpar.var_shapes_max[i] > 1e99 ? 0.0 : OCpar.var_shapes_max[i];
			  for (ii=1;ii<=RFshapes_len(OCpar.var_shapes[i]); ii++) {
				  double am = RFshapes[OCpar.var_shapes[i]][ii].ampl;
				  if (am < OCpar.var_shapes_min[i]) {
					  peni += pow(OCpar.var_shapes_min[i]-am,OCpar.var_shapes_penalty_order[i]);
				  }
				  if (am > lim) {
					  peni += pow(am-lim,OCpar.var_shapes_penalty_order[i]);
				  }
			  }
			  pen += peni*OCpar.var_shapes_penalty_factor[i];
		  }
	  }
  }

  OCpar.fvals[0] = tfval;
  OCpar.fvals[1] = pen;

  return tfval-pen;

}

/****
 * ZT: evaluate gradient procedure with current RFshapes
 ****/
 int evaluate_gradient(Tcl_Interp* interp)
{
  int gr_slot, i, ii;
  Tcl_Obj* objptr;
  extern FD** fd;

  /* set gradient mode for calculating gradients and reset matrix counter */
  OCpar.gradmode = 1;
  if (Tcl_EvalEx((Tcl_Interp*)interp,"gradient",-1,TCL_EVAL_GLOBAL) != TCL_OK) {
    fprintf(stderr,"error: Unable to execute command 'gradient':\n%s\n ",Tcl_GetStringResult(interp));
    exit(1);
  }
  objptr = Tcl_GetObjResult(interp);
  if ( Tcl_GetIntFromObj(interp, objptr, &gr_slot) != TCL_OK) {
    TclError(interp, "error in evaluate_gradient: Unable to convert result into integer");
    exit(1);
  }
  OCpar.gradmode = 0;

  // adding penalty gradients
  if (OCpar.ispenalty) {
	  complx *gr = fd[gr_slot]->data;
	  int idx = 0;
	  /* IMPORTANT: this assumes that gradients in fid are {gr_x, gr_y} !!! */
	  for (i=1; i<=OCpar.var_shapes[0]; i++) {
		  if (OCpar.var_shapes_penalty_order[i] != 0) {
			  int slot = OCpar.var_shapes[i];
			  double lim = OCpar.var_shapes_max[i] > 1e99 ? 0.0 : OCpar.var_shapes_max[i];
			  for (ii=1; ii<=RFshapes_len(slot); ii++) {
				  idx++;
				  double am = RFshapes[slot][ii].ampl;
				  double xx = am*cos( (RFshapes[slot][ii].phase)*DEG2RAD );
				  double yy = am*sin( (RFshapes[slot][ii].phase)*DEG2RAD );
				  if (am < OCpar.var_shapes_min[i]) {
					  gr[idx].re	+= pow(OCpar.var_shapes_min[i]-am,OCpar.var_shapes_penalty_order[i]-1)*OCpar.var_shapes_penalty_order[i]*(xx/am)*OCpar.var_shapes_penalty_factor[i];
					  gr[idx].im	+= pow(OCpar.var_shapes_min[i]-am,OCpar.var_shapes_penalty_order[i]-1)*OCpar.var_shapes_penalty_order[i]*(yy/am)*OCpar.var_shapes_penalty_factor[i];
				  }
				  if (am > lim) {
					  gr[idx].re	-= pow(am-lim,OCpar.var_shapes_penalty_order[i]-1)*OCpar.var_shapes_penalty_order[i]*(xx/am)*OCpar.var_shapes_penalty_factor[i];
					  gr[idx].im	-= pow(am-lim,OCpar.var_shapes_penalty_order[i]-1)*OCpar.var_shapes_penalty_order[i]*(yy/am)*OCpar.var_shapes_penalty_factor[i];
				  }
			  }
		  } else {
			  idx += RFshapes_len(OCpar.var_shapes[i]);
		  }
	  }
  }

  return gr_slot;
}

/****
 * ZT: function to update RFshapes with data moved by a step along dir
 *     IMPORTANT: dir is in x,y format!!!
 ****/
 void update_RFshapes_by_step(double step, double *dir, RFelem *CURRshapes[])
{
  double x, y, am, ph, rms;
  int i,j,ii,Nii;
  int Nshapes = OCpar.var_shapes[0];

  if (step == 0) {
     /* just copy CURRshapes to RFshapes */
     for (i=1; i<=Nshapes; i++) {
        for (ii=1; ii<=RFshapes_len(OCpar.var_shapes[i]); ii++) {
           RFshapes[OCpar.var_shapes[i]][ii].ampl = CURRshapes[i-1][ii].ampl;
	   RFshapes[OCpar.var_shapes[i]][ii].phase = CURRshapes[i-1][ii].phase;
	}
     }
  } else {
     j=0;
     for (i=1; i<=Nshapes; i++) {
    	 rms = 0.0;
    	 Nii = RFshapes_len(OCpar.var_shapes[i]);
    	 for (ii=1; ii<=Nii; ii++) {
    		 /* get original x,y components for this shape */
    		 x = CURRshapes[i-1][ii].ampl*cos(CURRshapes[i-1][ii].phase*DEG2RAD);
    		 y = CURRshapes[i-1][ii].ampl*sin(CURRshapes[i-1][ii].phase*DEG2RAD);
    		 /* move them to new value */
    		 j++;
    		 x += step*dir[j];
    		 j++;
    		 y += step*dir[j];
    		 /* transfer them to ampl, phase */
    		 am = sqrt(x*x+y*y);
    		 ph = RAD2DEG*atan2(y,x);
    		 rms += am*am;
    		 /* store them in relevant RFshape */
    		 RFshapes[OCpar.var_shapes[i]][ii].ampl = am;
    		 if (OCpar.var_shapes_penalty_order[i] == 0) {
    			 if (am > OCpar.var_shapes_max[i]) {
    				 RFshapes[OCpar.var_shapes[i]][ii].ampl = OCpar.var_shapes_max[i];
    			 }
    			 if (am < OCpar.var_shapes_min[i]) {
    				 RFshapes[OCpar.var_shapes[i]][ii].ampl = OCpar.var_shapes_min[i];
    			 }
    		 }
    		 RFshapes[OCpar.var_shapes[i]][ii].phase = ph;
    	 }
    	 rms = sqrt(rms/((double)Nii));
    	 if (rms > OCpar.var_shapes_rmsmax[i]) {
    		 double kkk;
    		 kkk = OCpar.var_shapes_rmsmax[i]/rms;
    		 for (ii=1; ii<=Nii; ii++) {
    			 RFshapes[OCpar.var_shapes[i]][ii].ampl *= kkk;
    		 }
    	 }
     }
  } /* end of else */

}

 void update_GROUPshapes_by_step(double step, double *dir, RFelem *CURRshapes[])
{
  double x, y;
  int j, i, ii, Nii;
  int Nshapes = OCpar.var_shapes[0];

  j=0;
  for (i=1; i<=Nshapes; i++) {
 	 Nii = RFshapes_len(OCpar.var_shapes[i]);
 	 for (ii=1; ii<=Nii; ii++) {
 		 /* get original x,y components for this shape */
 		 x = CURRshapes[i-1][ii].ampl;
 		 y = CURRshapes[i-1][ii].phase;
 		 /* move them to new value */
 		 j++;
 		 x += step*dir[j];
 		 j++;
 		 y += step*dir[j];
 		 /* store them in relevant RFshape */
 		 RFshapes[OCpar.var_shapes[i]][ii].ampl = x;
 		 RFshapes[OCpar.var_shapes[i]][ii].phase = y;
 	 }
  }
}

/****
 * ZT: function to minimize during linesearch
 ****/
 double func_to_min(double step, double *dir, RFelem *CURRshapes[], Tcl_Interp *interp)
{
  double fval;
  
  /* update RFshapes. Assumes RFshapes in ampl/phase and gradients in x/y */
  if (OCpar.group == 1)
	  update_GROUPshapes_by_step(step,dir,CURRshapes);
  else
	  update_RFshapes_by_step(step,dir,CURRshapes);
  /* call target function */
  fval = -evaluate_target_function(interp);
  
  return fval;
}

/****
 * ZT: bracketing minimum with 3 values
 ****/
 int bracketminCG(double *a, double *b, double *c, double *fa, double *fb, double *fc, double *dir, RFelem *CURRshapes[], Tcl_Interp *interp)
{
   double golden = 0.381966;
   double EPS = 1.0e-10 ; /*insignificant difference between two numbers*/
   double f_l, f_r, f_m, x_l, x_r, x_m, p, q, u, umax, v;
   int Neval = 0;
   int status=0;

   x_l = *a;
   x_r = *b;
   f_l = *fa; /* passed from caller, save one evaluation! */
   f_r = func_to_min(x_r, dir, CURRshapes, interp);
   Neval++;
   
   if (f_r >= f_l) {
      /* go into interval x_l,x_r */
      x_m = (x_r -x_l)*golden + x_l;
      f_m = func_to_min(x_m, dir, CURRshapes, interp);
      Neval++;
   } else {
      /* go further to the right */
      x_m = x_r;
      f_m = f_r;
      x_r = (x_m - x_l)/golden + x_l;
      f_r = func_to_min(x_r, dir, CURRshapes, interp);
      Neval++;
   }
   
   while ( (Neval<OCpar.max_brack_eval) && (fabs(x_l-x_r)>EPS) ) {
      if (f_m < f_l) {
         if (f_m < f_r) {
	    *a = x_l;
	    *b = x_m;
	    *c = x_r;
	    *fa = f_l;
	    *fb = f_m;
	    *fc = f_r;
	    if ( OCpar.verb ) {
               printf("bracketing reached after %d iterations\n",Neval);
            }
	    return status;
	 } else {
	    /* golden expansion limit */
	    v = (x_r-x_m)/golden+x_m;
	    /* use parabolic function */
	    p = (x_m-x_l)*(f_m-f_r);
	    q = (x_m-x_r)*(f_m-f_l);
            u = ( (x_m+x_r)*q - (x_m+x_l)*p )/(q-p)/2;
	    if ( OCpar.verb ) { printf("-1-"); }
	    if ( u > v ) {
	       /* accept parabolic fit only when it goes further than golden expansion */
	    /* if ( (u-x_r) > EPS ) {  */
	       /* u is further on the right = GOOD, but don't go too far */
	       umax = x_l+100.0*(x_r-x_l);
	       x_l = x_m;
               f_l = f_m;
               x_m = x_r;
               f_m = f_r;
               x_r = (u < umax) ? u : umax;
	       f_r = func_to_min(x_r, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-2-"); }
	    } else if ( ((u-x_m)>EPS) && ((x_r-u)>EPS) ) {
	       /* u is between x_m, x_r = GOOD */
	       x_l = x_m;
	       f_l = f_m;
	       x_m = u;
	       f_m = func_to_min(x_m, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-3-"); }
	    } else {
	       /* parabolic function was useless, take golden expansion */
	       x_l = x_m;
	       f_l = f_m;
	       x_m = x_r;
	       f_m = f_r;
	       x_r = (x_m-x_l)/golden+x_l; 
	       f_r = func_to_min(x_r, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-4-"); }
	    }
	 }
      } else {
         if ( f_m < f_r ) {
	    /* try parabola fist */
	    p = (x_m-x_l)*(f_m-f_r);
	    q = (x_m-x_r)*(f_m-f_l);
            u = ( (x_m+x_r)*q - (x_m+x_l)*p )/(q-p)/2;
	    if ( OCpar.verb ) { printf("-5-"); }
	    if ( ((u-x_l)>EPS) && ((x_m-u)>EPS) ) {
	       /* u is between x_l, x_m = GOOD */
	       x_r = x_m;
	       f_r = f_m;
	       x_m = u;
	       f_m = func_to_min(x_m, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-6-"); }
	    } else {
	       /* parabola was vaste of time, but min must be somewhere between x_l, x_m */
	       x_r = x_m;
	       f_r = f_m;
	       x_m =(x_r-x_l)*golden+x_l;
	       f_m = func_to_min(x_m, dir, CURRshapes, interp);
	       Neval++;
	       if ( OCpar.verb ) { printf("-7-"); }
	    }
	 } else {
	    /* also here we insist that min is right to x_l and step to right was too big */
            x_r = x_m;
            f_r = f_m;
            x_m =(x_r-x_l)*golden+x_l;
            f_m = func_to_min(x_m, dir, CURRshapes, interp);
            Neval++;
	    if ( OCpar.verb ) { printf("-8-"); }
	 }
      }
      
      if ( OCpar.verb ) { 
         printf("bracketing status:x=(%f, %f, %f), F=(%.15g, %.15g, %.15g)\n",x_l,x_m,x_r,f_l,f_m,f_r);   
      }

   } /* end of while */
   
   *a = x_l;
   *b = x_m;
   *c = x_r;
   *fa = f_l;
   *fb = f_m;
   *fc = f_r;
   printf("Bracketmin: Failed to bracket, maximum function evaluation reached or limits closer than machine precision\n");
   printf("ERROR!!! --> Aborting optimization...\n");
   /* exit(1); */
   status = 1;
   return status;
   
}

/****
 * ZT: find minimum with Brent method
 ****/
 void brentCG(double *a, double *b, double *c, double *fa, double *fb, double *fc, double *dir, RFelem *CURRshapes[], Tcl_Interp *interp)
{
   double golden = 0.3819660;
   double EPS = 1.4901161193847656e-08; /* sqrt(double precision) */
   double f_l, f_r, f_m, x_l, x_r, x_m, v, w, d, e, f_v, f_w, z, f_z;
   double w_lo, w_up, tolerance, p, q, r, midpoint, u, f_u;
   double stoptol;
   int Neval = 0;
   
   x_l = *a;
   f_l = *fa;
   x_m = *b;
   f_m = *fb;
   x_r = *c;
   f_r = *fc;
   
   /* initialization */
   stoptol = ( x_l*x_r > 0.0 ) ? (OCpar.tol*x_m) : OCpar.stepmin; 
   v = x_l + golden*(x_r - x_l);
   w = v;
   d = 0.0;
   e = 0.0;
   f_v = func_to_min(v, dir, CURRshapes, interp);
   Neval++;
   f_w = f_v;
   
   while ( (Neval<OCpar.max_brent_eval) && (fabs(x_r-x_l)>stoptol) ) {
    /* printf("Brent status1: x=<%f, %f>, f=<%f, %f>\n",x_l,x_r,f_l,f_r); */
      z = x_m;
      f_z = f_m;
      d = e;
      e = d;
      w_lo = z -x_l;
      w_up = x_r -z;
      tolerance = EPS*fabs(z);
      p = 0.0;
      q = 0.0;
      r = 0.0;
      midpoint = 0.5*(x_l+x_r);
      
      if ( fabs(e)>tolerance) {
         /* fit parabola */
	 r = (z-w)*(f_z-f_v);
	 q = (z-v)*(f_z-f_w);
	 p = (z-v)*q - (z-w)*r;
	 q = 2*(q-r);
	 if (q>0) {
	    p = -p;
	 } else {
	    q = -q;
	 }
	 r = e;
	 e = d;
      }
      
      if  ( (fabs(p)<fabs(0.5*q*r)) && (p<q*w_lo) && (p<q*w_up) ) {
         double t2 = 2.0*tolerance;
	 d = p/q;
	 u = z + d;
	 if ( (u-x_l<t2) || (x_r-u<t2) )
	    d = (z<midpoint) ? tolerance : -tolerance;
      } else {
         e = ( (z<midpoint) ? x_r : x_l ) - z;
	 d = golden*e;
      }
      
      if ( fabs(d) >= tolerance ) {
         u = z + d;
      } else {
         u = ( (d>0.0) ? tolerance : -tolerance) + z;
      }
      
      f_u = func_to_min(u, dir, CURRshapes, interp);
      Neval++;
      
      if ( f_u<f_z ) {
         if ( u<z) {
	    x_r = z;
	    f_r = f_z;
	 } else {
	    x_l = z;
	    f_l = f_z;
	 }
	 v = w;
	 f_v = f_w;
	 w = z;
	 f_w = f_z;
	 x_m = u;
	 f_m = f_u;
      } else {
         if ( u<z) {
	    x_l = u;
	    f_l = f_u;
	 } else {
	    x_r = u;
	    f_r = f_u;
	 }
	 if ( (f_u <= f_w) || (w == z) ) {
	    v = w;
	    f_v = f_w;
	    w = u;
	    f_w = f_u;
	 } else if ( (f_u<=f_v) || (v==z) || (v==w) ) {
	    v = u;
	    f_v = f_u;
	 }
      }
      if ( OCpar.verb ) { 
         printf("Brent status2: x=<%f, %f>, f=<%f, %f>\n",x_l,x_r,f_l,f_r);
      }
      stoptol = ( x_l*x_r > 0.0 ) ? (OCpar.tol*x_m) : OCpar.stepmin; 
   } /* end of 'while' iteration loop */
   
   if ( OCpar.verb ) {
      printf("Brent used %d iterations\n",Neval);
   }
   *b = x_m;
   *fb = f_m;
   /* other brackets are not changed on exit (since they are not relevant) */
   /* but I added print-out so I rather change them, it looks nicer */
   *a = x_l;
   *fa = f_l;
   *c = x_r;
   *fc = f_r;

}

/****
 * ZT: linesearch envelope
 ****/
 double linesearchCG(double* fcurr, double *dir, RFelem *CURRshapes[], Tcl_Interp* interp)
{
  double stepsize, a, b, c, fa, fb, fc;
  int bkstat;
  
  /* find triplet enclosing minimum */
  a = 0.0;                /* we start here with bracketing */
  b = (double)OCpar.mnbkstep; /* this defines initial step size for bracketing */
  c = 2.0*b;                /* this is not realy used on call, just on exit */
  fa = *fcurr;              /* passed from caller, save one evaluation! */
  bkstat = bracketminCG(&a,&b,&c,&fa,&fb,&fc,dir,CURRshapes,interp);
  if ( OCpar.verb ) {
     printf("BRACKETING minimum results:\n %10.5f %10.5f\n %10.5f %10.5f\n %10.5f %10.5f\n",a,fa,b,fb,c,fc);
  }
  if ( bkstat == 0) {
     /* find local minimum using Brent method */
     brentCG(&a,&b,&c,&fa,&fb,&fc,dir,CURRshapes,interp);
     if ( OCpar.verb ) {
        printf("BRENT minimization results:\n %10.5f %10.5f\n %10.5f %10.5f\n %10.5f %10.5f\n",a,fa,b,fb,c,fc);
     }
     stepsize = b;
     *fcurr = fb;
  } else {
     stepsize = 0;
     /* *fcurr is not changed when bracketing failed... */
  }

  
  /******** D E B U G G I N G ****** R E M O V E    T H I S ******/
  /*    double g1=-500.0, g2=1000.0, gs=2.0;
      double dum;
      printf("----PROFILE START-----\n");
      while (g1<=g2) {
          dum = func_to_min(g1,dir,CURRshapes,interp);
          printf("   %.15g   %.15g\n",g1, dum);
	  fflush(stdout);
          g1 += gs; 
      }
      printf("----PROFILE END------\n");
      fflush(stdout);                                            */
  /******** D E B U G G I N G ****** R E M O V E    T H I S ******/
  
  
  /* make sure that current RFshapes hold this new position */
  if (OCpar.group == 1)
	  update_GROUPshapes_by_step(stepsize,dir,CURRshapes);
  else
	  update_RFshapes_by_step(stepsize,dir,CURRshapes);

  return stepsize;
}

/****
 * ZT: function to store RFshapes, option enabling monitoring of changes to shapes
 *     during optimization
 ****/
 void store_OCshapes(Tcl_Interp *interp)
{
  /* not used when 'none' was specified */
  if ( strcmp(OCpar.VarSaveProc,"none") ) {
     if (Tcl_EvalEx( (Tcl_Interp*)interp, OCpar.VarSaveProc, -1, TCL_EVAL_GLOBAL) != TCL_OK) {
        TclError(interp, "error: Unable to execute command '%s' for storing RF shapes during oc_optimize\n",OCpar.VarSaveProc);
        exit(1);
     }
  }
  
}


/****
 * helper function for freeing fid from memory (when gradients are no longer needed)
 ****/
void free_grad_fid(int gr_slot, Tcl_Interp *interp)
{ 
  char buf1[64];
  extern int funload(int fidN);
  
  if (funload(gr_slot) == -1) {
    fprintf(stderr,"Error in free_grad_fid\n");
    exit(1);
  }
  sprintf(buf1,"%d",gr_slot);
  Tcl_UnsetVar2(interp,"FD_Internal::f",buf1,0);
}

/****
 * ZT: conjugated gradients
 ****/
 double OptimizeCG(Tcl_Interp* interp)
{
 double fold, fcurr, dum, g, beta;
 double step, tfval;
 double *r, *rN, *dir;
 int gr_slot, count, iter, i, j, N;
 int Nshapes = OCpar.var_shapes[0];
 int Nvars=0;
 //RFelem *CURRshapes[Nshapes];
 RFelem **CURRshapes = (RFelem**)malloc(Nshapes*sizeof(RFelem*));
 extern FD** fd;
 
 
 /* initialize function value and gradient */
 tfval = evaluate_target_function(interp); /* target function should be maximzed, here we minimize */
 fold = -tfval;
 if (OCpar.ispenalty) {
	 printf("Initial target function: %.10f (phi = %g, pen = %g)\n", tfval, OCpar.fvals[0], OCpar.fvals[1]);
 } else {
	 printf("Initial target function: %.10f \n", tfval);
 }
 gr_slot = evaluate_gradient(interp);

 /* store current RF shapes at safe place */
 for (i=0; i<Nshapes; i++) {
    N = RFshapes_len(OCpar.var_shapes[i+1]);
    Nvars += N;
    CURRshapes[i]=RFshapes_alloc(N);
    for (j=1; j<=N; j++) {
       CURRshapes[i][j].ampl = RFshapes[OCpar.var_shapes[i+1]][j].ampl;
       CURRshapes[i][j].phase = RFshapes[OCpar.var_shapes[i+1]][j].phase;
    }
 }
 
 /* check length of gradients with total number of variables */
 if (Nvars != (fd[gr_slot]->np)) {
    fprintf(stderr,"oc_optimize error: number of gradients does not match total number of variables in RF shapes (%d versus %d)\n",Nvars,(fd[gr_slot]->np) );
    exit(1);
  }
 /* Nvars so far counts just elements in RF shapes. But vector dimension is double (aml, ph) */
 Nvars *=2;
 
 /* initialize variables for conjugated gradients */
 dir = double_vector(Nvars); 
 r = double_vector(Nvars);
 rN = double_vector(Nvars);
 i = 0;
 for (j=1; j<=fd[gr_slot]->np; j++) {
     i++;
     dum = (fd[gr_slot]->data)[j].re;
     r[i] = dum; /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = r[i];
     i++;
     dum = (fd[gr_slot]->data)[j].im;
     r[i] = dum; /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = r[i];
 }
 /* here all work with gadients in FID is done, free it */
 free_grad_fid(gr_slot,interp);
 
 count = 0;
 fcurr = fold;
 
 /* start iteration loop */
 for (iter=1; iter<=OCpar.nIterations; iter++) {
    /* do linesearch */
    step = linesearchCG(&fcurr, dir, CURRshapes, interp);
    tfval = -fcurr;
    /* store new point RFshapes to safe place (update CURRshapes) */
     for (i=1; i<=Nshapes; i++) {
        for (j=1; j<=RFshapes_len(OCpar.var_shapes[i]); j++) {
           CURRshapes[i-1][j].ampl = RFshapes[OCpar.var_shapes[i]][j].ampl;
	   CURRshapes[i-1][j].phase = RFshapes[OCpar.var_shapes[i]][j].phase;
	}
     }
     count++;
    
    /* print out results */
    if (OCpar.ispenalty) {
    	printf("  Iter %d: tf=%.10f (phi = %g, pen = %g, step=%g)\n",iter, tfval, OCpar.fvals[0], OCpar.fvals[1], step);
    } else {
    	printf("  Iter %d: tf=%.10f (step=%g)\n",iter, tfval, step);
    }
    /* timing - TOC every iteration */
    //QueryPerformanceCounter(&_tv2_);
    //printf("\tTiming: %.9f\n",((float)(_tv2_.QuadPart-_tv1_.QuadPart))/(float)_tickpsec_.QuadPart);

    /* store sequences every nreport steps */
    if ( (iter % OCpar.nreport) < 1 ) {
       printf("   Storing data for iteration %d\n",iter);
       store_OCshapes(interp);
    }
    
    /* check for done */
    if (fabs(fcurr-fold)<=OCpar.eps) {
       printf("Done by reaching target function tolerance limit.\n\n");
       break;
    }
    if( (tfval<OCpar.cut) && (iter>OCpar.ncut) ) {
       printf("Done %d iterations, but target function below cut-off limit.\n\n",iter);
       break;
    }
    if ( fabs(step)<OCpar.stepmin ) {
       printf("Done by reaching minimal step in line-search.\n\n");
       break;
    }              
    
    /* calculate new gradient */
    gr_slot = evaluate_gradient(interp);
    
    /* prepare parameters for conjugation */
    g = 0.0;
    beta = 0.0;
    i = 0;
    for (j=1; j<=fd[gr_slot]->np; j++) {
       i++;
       dum = (fd[gr_slot]->data)[j].re;
       rN[i] = dum; /* gradient calculated for maximizing, here we minimize. rN[i]=-dum */
       g += r[i]*r[i];
       beta += (rN[i]-r[i])*rN[i]; /* this is for Polak-Ribiere */
       i++;
       dum = (fd[gr_slot]->data)[j].im;
       rN[i] = dum; /* gradient calculated for maximizing, here we minimize. rN[i]=-dum */
       g += r[i]*r[i];
       beta += (rN[i]-r[i])*rN[i]; /* this is for Polak-Ribiere */
    }
    /* here all work with gadients in FID is done, free it */
    free_grad_fid(gr_slot,interp);
    
    /* check for exact reaching of extremum (nul-gradient) */
    if (g==0.0) {
       printf("Done by reaching extremum (nul gradient).\n\n");
       break;
    }
    beta = beta/g;
    /* reset search direction if beta<0 or all independent directions were already taken */
    if ( (beta < 0.0) || (count == Nvars) ) {
       printf("(reseting search directions)\n");
       beta = 0.0;
       count = 0;
    }
    
    /* generate new direction by conjugation and prepare for next iteration */
    for (i=1; i<=Nvars; i++) {
       dir[i] = rN[i]+beta*dir[i]; 
       r[i]=rN[i];
    }
    fold = fcurr;
        
 } /* end of iteration loop */
 
 if (iter >= OCpar.nIterations) 
    printf("Done by reaching maximal number of iterations.\n\n");
 
 
 free_double_vector(dir);
 free_double_vector(r);
 free_double_vector(rN);
 for (i=0; i<Nshapes; i++) {
    free((char *)CURRshapes[i]);
 }
 free(CURRshapes);
 /* make sure that correct final target function value is returned */
 tfval = -fcurr;
 
 return tfval;
 
}

/*****
 * 'evaluate' function for L-BFGS, returns target_function and gradient
 *****/
static lbfgsfloatval_t lbfgs_evaluate(void *interp,
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step)
{
	lbfgsfloatval_t fx;
	int i, j, k, N;
	extern FD** fd;

	/* copy variables into shapes */
	if (OCpar.group == 0) { // standard RF shapes with ampl/phase
      k = 0;
      for (i=1; i<=OCpar.var_shapes[0]; i++) {
         N = RFshapes_len(OCpar.var_shapes[i]);
         for (j=1; j<=N; j++) {
    	     double re = x[k++];
    	     double im = x[k++];
    	     RFshapes[OCpar.var_shapes[i]][j].ampl = sqrt(re*re+im*im);
    	     RFshapes[OCpar.var_shapes[i]][j].phase = RAD2DEG*atan2(im,re);
         }
      }
	} else if (OCpar.group == 1) { // GROUP RF shapes with x/y
      k = 0;
      for (i=1; i<=OCpar.var_shapes[0]; i++) {
         N = RFshapes_len(OCpar.var_shapes[i]);
         for (j=1; j<=N; j++) {
    	     RFshapes[OCpar.var_shapes[i]][j].ampl = x[k++];
    	     RFshapes[OCpar.var_shapes[i]][j].phase = x[k++];
         }
      }
	} else {
    	fprintf(stderr,"ERROR: unexpected value of OCpar.group in lbfgs\n");
    	exit(1);
	}
    /* evaluate target function */
    fx = -evaluate_target_function((Tcl_Interp*)interp);
    //printf("ev: tf = %g\n",fx);
	/* get gradients */
    N = evaluate_gradient((Tcl_Interp*)interp);
    //printf("ev: grads done\n");
    /* fill gradients */
    i = 0;
    for (j=1; j<=fd[N]->np; j++) {
    	g[i++] = -(fd[N]->data)[j].re;
    	g[i++] = -(fd[N]->data)[j].im;
    	/* gradient calculated for maximizing, here we minimize */
    }
    free_grad_fid(N,(Tcl_Interp*)interp);
    //printf("ev: grads filled\n\n");

	return fx;
}

/*****
 * 'evaluate' function for L-BFGS and just phase optimization, returns target_function and gradient
 *****/
static lbfgsfloatval_t lbfgs_phase_evaluate(void *interp,
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step)
{
	lbfgsfloatval_t fx;
	int i, j, k, N, slot;
	extern FD** fd;

	/* copy variables (i.e. phase) into shapes */
    k = 0;
    for (i=1; i<=OCpar.var_shapes[0]; i++) {
       N = RFshapes_len(OCpar.var_shapes[i]);
       for (j=1; j<=N; j++) {
    	   RFshapes[OCpar.var_shapes[i]][j].phase = RAD2DEG*x[k++];
       }
    }
    /* evaluate target function */
    fx = -evaluate_target_function((Tcl_Interp*)interp);
    //printf("ev: tf = %g\n",fx);
	/* get gradients */
    slot = evaluate_gradient((Tcl_Interp*)interp);
    //printf("ev: grads done\n");
    /* fill gradients */
    k = 0;
    for (i=1; i<=OCpar.var_shapes[0]; i++) {
       N = RFshapes_len(OCpar.var_shapes[i]);
       for (j=1; j<=N; j++) {
    	   double am = RFshapes[OCpar.var_shapes[i]][j].ampl;
     	   double ph = RFshapes[OCpar.var_shapes[i]][j].phase*DEG2RAD;
     	   g[k++] = -am*(-(fd[slot]->data)[j].re*sin(ph) +(fd[slot]->data)[j].im*cos(ph));
     	  /* gradient calculated for maximizing, here we minimize */
       }
    }
    free_grad_fid(slot,(Tcl_Interp*)interp);
    //printf("ev: grads filled\n\n");

	return fx;
}

/*****
 * 'progress' function for L-BFGS: output, temporary store, early termination, ...
 *****/
static int lbfgs_progress(void *interp,
	const lbfgsfloatval_t *x,
	const lbfgsfloatval_t *g,
	const lbfgsfloatval_t fx,
	const lbfgsfloatval_t xnorm,
	const lbfgsfloatval_t gnorm,
	const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    printf("  Iter. %d: ", k);
    if (OCpar.ispenalty) {
    	printf("  tf = %.10f  (phi = %g, pen = %g, step = %f)\n", -fx, OCpar.fvals[0], OCpar.fvals[1], step);
    } else {
    	printf("  tf = %.10f  (step = %f)\n", -fx, step);
    }

    if (OCpar.verb) {
    	printf("\t\tNumber of evaluations   %d\n",ls);
    	printf("\t\tNorm of variable change %.10f\n",xnorm);
    	printf("\t\tNorm of gradients       %.10f\n",gnorm);
    }
    /* store sequences every nreport steps */
    if ( (k % OCpar.nreport) < 1 ) {
       printf("   Storing data for iteration %d\n",k);
       store_OCshapes((Tcl_Interp*)interp);
    }
    /* stop if search is not promissing */
    if( (-fx < OCpar.cut) && (k > OCpar.ncut) ) {
       printf("Done %d iterations, but target function below cut-off limit.\n\n",k);
       return 1;
    }

    /* timing - TOC every iteration */
    //QueryPerformanceCounter(&_tv2_);
    //printf("\tTiming: %.9f\n",((float)(_tv2_.QuadPart-_tv1_.QuadPart))/(float)_tickpsec_.QuadPart);

    return 0;
}

/******
 * optimization using L-BFGS
 ******/
double OptimizeLBFGS(Tcl_Interp* interp)
{
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x;
    lbfgs_parameter_t param;
    int i, j, k, N, ret;
    int Nshapes = OCpar.var_shapes[0];
    int Nvars = 0;

    /* count variables */
    for (i=1; i<=Nshapes; i++) {
       N = RFshapes_len(OCpar.var_shapes[i]);
       Nvars += N;
    }
    Nvars *= 2;
    /* initialize vector of variables x - rf Ix, rf Iy */
    x = lbfgs_malloc(Nvars);
    if (OCpar.group == 0) { // standard RF shapes with ampl/phase
      k = 0;
      for (i=1; i<=Nshapes; i++) {
         N = RFshapes_len(OCpar.var_shapes[i]);
         for (j=1; j<=N; j++) {
    	     double am = RFshapes[OCpar.var_shapes[i]][j].ampl;
    	     double ph = RFshapes[OCpar.var_shapes[i]][j].phase*DEG2RAD;
    	     x[k++] = am*cos(ph);
    	     x[k++] = am*sin(ph);
         }
      }
    } else if (OCpar.group == 1) { // GROUP shapes with x/y
        k = 0;
        for (i=1; i<=Nshapes; i++) {
           N = RFshapes_len(OCpar.var_shapes[i]);
           for (j=1; j<=N; j++) {
      	     x[k++] = RFshapes[OCpar.var_shapes[i]][j].ampl;
      	     x[k++] = RFshapes[OCpar.var_shapes[i]][j].phase;
           }
        }

    } else {
    	fprintf(stderr,"ERROR: unexpected value of OCpar.group in lbfgs\n");
    	exit(1);
    }
    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.epsilon = OCpar.eps; /* convergence criterion */
    param.max_iterations = OCpar.nIterations;
    /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
    param.ftol = OCpar.tol;
    param.max_linesearch = OCpar.max_brack_eval;
    param.m = OCpar.lbfgs_m;

    fx = evaluate_target_function(interp);
    if (OCpar.ispenalty) {
    	printf("Initial target function: %.10f (phi = %g, pen = %g)\n", fx, OCpar.fvals[0],OCpar.fvals[1]);
    } else {
    	printf("Initial target function: %.10f \n", fx);
    }

    /* do optimization */
    ret = lbfgs(Nvars, x, &fx, lbfgs_evaluate, lbfgs_progress, (void*)interp, &param);

    /* Report the result. */
    printf("\nL-BFGS optimization terminated with status message:\n\t'%s'\n", lbfgs_retcode(ret));
    printf("Terminal target function value: %.10f\n\n", -fx);

    lbfgs_free(x);

    return (double)fx;
}

/******
 * optimization of just phase and using L-BFGS
 ******/
double OptimizePhaseLBFGS(Tcl_Interp* interp)
{
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x;
    lbfgs_parameter_t param;
    int i, j, k, N, ret;
    int Nshapes = OCpar.var_shapes[0];
    int Nvars = 0;

    /* count variables */
    for (i=1; i<=Nshapes; i++) {
       N = RFshapes_len(OCpar.var_shapes[i]);
       Nvars += N;
    }
    /* initialize vector of variables x - rf phase */
    x = lbfgs_malloc(Nvars);
    k = 0;
    for (i=1; i<=Nshapes; i++) {
       N = RFshapes_len(OCpar.var_shapes[i]);
       for (j=1; j<=N; j++) {
    	   x[k++] = RFshapes[OCpar.var_shapes[i]][j].phase*DEG2RAD;
       }
    }
    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.epsilon = OCpar.eps; /* convergence criterion */
    param.max_iterations = OCpar.nIterations;
    /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
    param.ftol = OCpar.tol;
    param.max_linesearch = OCpar.max_brack_eval;
    param.m = OCpar.lbfgs_m;

    fx = evaluate_target_function(interp);
    if (OCpar.ispenalty) {
    	printf("Initial target function: %.10f (phi = %g, pen = %g)\n", fx, OCpar.fvals[0],OCpar.fvals[1]);
    } else {
    	printf("Initial target function: %.10f \n", fx);
    }

    /* do optimization */
    ret = lbfgs(Nvars, x, &fx, lbfgs_phase_evaluate, lbfgs_progress, (void*)interp, &param);

    /* Report the result. */
    printf("\nL-BFGS optimization terminated with status message:\n\t'%s'\n", lbfgs_retcode(ret));
    printf("Terminal target function value: %.10f\n\n", -fx);

    lbfgs_free(x);

    return (double)fx;
}

/****
 * ZT: implementation of oc_optimize
 ****/
 int tclOCoptimize2(ClientData data,Tcl_Interp* interp, int argc, char *argv[])
{
  double tfval, dumfloat;
  Tcl_Obj* objptr;
  int i, slot, ndim, nsh, ish, dumint;
  char strbuf[32];

  /* timings - initial TIC */
  //QueryPerformanceFrequency(&_tickpsec_);
  //QueryPerformanceCounter(&_tv1_);

  /* disable when relaxation is ON */
  TclGetString(interp,strbuf,"par","relax",0,"off");
  DEBUGPRINT("OC detects: relax is %s\n",strbuf);
  if (!strncmp(strbuf,"on",2)) {
     fprintf(stderr,"oc_optimize error: relaxation not supported during optimization.\n");
     exit(1);
  }

  /* read input parameters */
  if ( argc == 1)
    return TclError(interp,"Usage: <double> oc_optimize <RFshape> ?-min XXX -max XXX -rmsmax XXX? ?<RFshape> ...?");
   
  /* go through argv and count shapes */
  nsh=0;
  for (i=1; i<argc; i++) {
     if (strncmp(argv[i],"-",1)) {
        /* this seems to be a shape */
        nsh++;
     } else {
        /* it is a switch, skip over its value */
        i++;
     }
  }
  if (nsh == 0)
     return TclError(interp,"oc_optimize has not detected any variable shape among arguments!");
  if (nsh > MAXOCSHAPES)
     return TclError(interp,"oc_optimize: too many shapes, maximum is %d.",MAXOCSHAPES);
     
  ndim = 0;
  ish = 0;
  OCpar_initialize();
  OCpar.var_shapes = int_vector(nsh); 
  OCpar.var_shapes_min = double_vector(nsh); 
  OCpar.var_shapes_max = double_vector(nsh); 
  OCpar.var_shapes_rmsmax = double_vector(nsh);
  OCpar.grad_shapes = int_vector(nsh);
  OCpar.var_shapes_penalty_order = int_vector(nsh);
  OCpar.var_shapes_penalty_factor = double_vector(nsh);
  for (i=1; i<argc; i++) {
     if (strncmp(argv[i],"-",1)) {
        /* this seems to be a shape */
        if (Tcl_GetInt(interp,argv[i],&slot) != TCL_OK) {
           OCpar_destroy();
           return TclError(interp,"oc_optimize: argument %d must be integer <RFshape>",i);
           }
        if (!RFshapes[slot]) {
           OCpar_destroy();
           return TclError(interp,"oc_optimize: argument %d -> RFshape does not exist", i);
        }
        ish++;
        OCpar.var_shapes[ish] = slot;
        OCpar.var_shapes_min[ish] = 0.0;
        OCpar.var_shapes_max[ish] = 2.0e99;
        OCpar.var_shapes_rmsmax[ish] = 2.0e99;
        OCpar.grad_shapes[ish] = slot;
        OCpar.var_shapes_penalty_order[ish] = 0;
        OCpar.var_shapes_penalty_factor[ish] = 0.0;
        ndim += RFshapes_len(slot);
     } else {
        /*  it seems to be limits */
        if (ish == 0) {
           OCpar_destroy();
           return TclError(interp,"oc_optimize error reading arguments: limits preceed shape!");
        }
        if ( !strcmp(argv[i],"-min") ) {
           if ( Tcl_GetDouble( interp,argv[i+1],&dumfloat) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -min: can not get double from %s",ish,argv[i+1]);
           } 
           OCpar.var_shapes_min[ish] = dumfloat;
           i++;
        } else if ( !strcmp(argv[i],"-max") ) {
           if ( Tcl_GetDouble( interp,argv[i+1],&dumfloat) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -max: can not get double from %s",ish,argv[i+1]);
           } 
           OCpar.var_shapes_max[ish] = dumfloat;
           i++;
        } else if ( !strcmp(argv[i],"-rmsmax") ) {
           if ( Tcl_GetDouble( interp,argv[i+1],&dumfloat) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -rmsmax: can not get double from %s",ish,argv[i+1]);
           } 
           OCpar.var_shapes_rmsmax[ish] = dumfloat;
           i++;
        } else if ( !strcmp(argv[i],"-order") ) {
           if ( Tcl_GetInt( interp,argv[i+1],&dumint) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -order: can not get integer from %s",ish,argv[i+1]);
           }
           OCpar.var_shapes_penalty_order[ish] = dumint;
           OCpar.ispenalty = 1;
           i++;
        } else if ( !strcmp(argv[i],"-scl") ) {
           if ( Tcl_GetDouble( interp,argv[i+1],&dumfloat) != TCL_OK ) {
              OCpar_destroy();
              return TclError(interp,"oc_optimize error: %d. shape, switch -scl: can not get double from %s",ish,argv[i+1]);
           }
           OCpar.var_shapes_penalty_factor[ish] = dumfloat;
           i++;
        } else {
           OCpar_destroy();
           return TclError(interp,"oc_optimize error reading %d.argument '%s' ",i,argv[i]);
        }
     }
  }
  /* debug output */
  //for (i=1; i<=nsh; i++) {
  // DEBUGPRINT("Argument %d: shape slot %d, min=%g, max=%g, RMSmax=%g\n",i,OCpar.var_shapes[i],OCpar.var_shapes_min[i], OCpar.var_shapes_max[i], OCpar.var_shapes_rmsmax[i]);
  //}
  /* reading of input parameters is complete */
 
  test_pulseq_for_acqOC_prop(interp);
  OCpar.ndim = ndim;

  /* read in all parameters from par(OC_?) or set their default values */
  TclGetString(interp,strbuf,"par","oc_method",0,"CG");
  if (!strncmp(strbuf,"CG",2)) {
	  OCpar.method = 1;
	  if (!strncmp(strbuf,"CG-GROUP",8))
		  OCpar.group = 1;
	  else
		  OCpar.group = 0;
  } else if (!strncmp(strbuf,"L-BFGS",6)) {
	  OCpar.method = 0;
	  if (!strncmp(strbuf,"L-BFGS-GROUP",12))
		  OCpar.group = 1;
	  else
		  OCpar.group = 0;
  } else {
     fprintf(stderr,"oc_optimize error: invalid oc_method, can be 'CG' or 'L-BFGS'\n");
     exit(1);
  }
  objptr = Tcl_GetVar2Ex( interp, "par", "oc_verbose", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.verb = 0;
  } else {
     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.verb) != TCL_OK ) {
       OCpar.verb = 0;
       printf("Warning! Could'n read par(oc_verbose), using default value %d\n",OCpar.verb);
     } 
  }
  objptr = Tcl_GetVar2Ex( interp, "par", "oc_max_iter", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.nIterations = 1000;
  } else {
     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.nIterations) != TCL_OK ) {
       OCpar.nIterations = 1000;
       printf("Warning! Could'n read par(oc_max_iter), using default value %d\n",OCpar.nIterations);
     } 
  }
  objptr = Tcl_GetVar2Ex( interp, "par", "oc_cutoff", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.cut = 0.0;
  } else {
     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.cut) != TCL_OK ) {
       OCpar.cut = 0.0;
       printf("Warning! Could'n read par(oc_cutoff), using default value %f\n",OCpar.cut);
     } 
  }
  objptr = Tcl_GetVar2Ex( interp, "par", "oc_cutoff_iter", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.ncut = 1000;
  } else {
     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.ncut) != TCL_OK ) {
       OCpar.ncut = 1000;
       printf("Warning! Could'n read par(oc_cutoff_iter), using default value %d\n",OCpar.ncut);
     } 
  }
  objptr = Tcl_GetVar2Ex( interp, "par", "oc_var_save_iter", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.nreport = 1001;
  } else {
     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.nreport) != TCL_OK ) {
       OCpar.nreport = 1001;
       printf("Warning! Could'n read par(oc_var_save_iter), using default value %d\n",OCpar.nreport);
     } 
  } 
  objptr = Tcl_GetVar2Ex( interp, "par", "oc_var_save_proc", TCL_GLOBAL_ONLY);
  if (!objptr) {
     strcpy(OCpar.VarSaveProc,"none");
     objptr = Tcl_NewStringObj(OCpar.VarSaveProc,strlen(OCpar.VarSaveProc));
     objptr = Tcl_SetVar2Ex( interp,"par", "oc_var_save_proc", objptr, TCL_GLOBAL_ONLY);
  } else {
     strcpy(OCpar.VarSaveProc, Tcl_GetString(objptr));
     if ( OCpar.VarSaveProc  == NULL ) {
       strcpy(OCpar.VarSaveProc,"none");
       objptr = Tcl_NewStringObj(OCpar.VarSaveProc,strlen(OCpar.VarSaveProc));
       objptr = Tcl_SetVar2Ex( interp,"par", "oc_var_save_proc", objptr, TCL_GLOBAL_ONLY);
       printf("Warning! Could'n read par(oc_var_save_proc), using default value %s\n",OCpar.VarSaveProc);
     } 
  }
  objptr = Tcl_GetVar2Ex( interp, "par", "oc_grad_level", TCL_GLOBAL_ONLY);
  if (!objptr) {
     OCpar.grad_level = 1;
  } else {
     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.grad_level) != TCL_OK ) {
       OCpar.grad_level = 1;
       printf("Warning! Could'n read par(oc_grad_level), using default value %g\n",OCpar.grad_level);
     }
  }

  if (OCpar.method == 1) { /* conjugated gradients */
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_tol_cg", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.eps = 1.0e-6;
	  } else {
	     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.eps) != TCL_OK ) {
	       OCpar.eps = 1.0e-6;
	       printf("Warning! Could'n read par(oc_tol_cg), using default value %f\n",OCpar.eps);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_tol_ls", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.tol = 1.0e-3;
	  } else {
	     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.tol) != TCL_OK ) {
	       OCpar.tol = 1.0e-3;
	       printf("Warning! Could'n read par(oc_tol_ls), using default value %f\n",OCpar.tol);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_mnbrak_step", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.mnbkstep = 10.0;
	  } else {
	     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.mnbkstep) != TCL_OK ) {
	       OCpar.mnbkstep = 10.0;
	       printf("Warning! Could'n read par(oc_mnbrak_step), using default value %f\n",OCpar.mnbkstep);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_cg_min_step", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.stepmin = 1.0e-3;
	  } else {
	     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.stepmin) != TCL_OK ) {
	       OCpar.stepmin = 1.0e-3;
	       printf("Warning! Could'n read par(oc_cg_min_step), using default value %f\n",OCpar.stepmin);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_max_brack_eval", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.max_brack_eval = 100;
	  } else {
	     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.max_brack_eval) != TCL_OK ) {
	       OCpar.max_brack_eval = 100;
	       printf("Warning! Could'n read par(oc_max_brack_eval), using default value %d\n",OCpar.max_brack_eval);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_max_brent_eval", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.max_brent_eval = 100;
	  } else {
	     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.max_brent_eval) != TCL_OK ) {
	       OCpar.max_brent_eval = 100;
	       printf("Warning! Could'n read par(oc_max_brent_eval), using default value %d\n",OCpar.max_brent_eval);
	     }
	  }

	  printf("Number of variables is %d\n",OCpar.ndim);
	  printf("Global tolerance on target function is %g\n",OCpar.eps);
	  printf("Maximal number of iterations is %d\n",OCpar.nIterations);
	  printf("Optimization is terminated after %d iterations if target function is not higher than %g\n",OCpar.ncut,OCpar.cut);
	  printf("Every %d iterations procedure '%s' (for storing variables) is executed\n",OCpar.nreport, OCpar.VarSaveProc);
	  printf("Minimal step size along conjugated gradient is %g\n",OCpar.stepmin);
	  printf("Tolerance for line-search (Brent) is %g\n",OCpar.tol);
	  printf("Line-search will terminate after reaching %d evaluations\n",OCpar.max_brent_eval);
	  printf("Initial step for bracketing minimum is %g\n",OCpar.mnbkstep);
	  printf("Bracketing will fail after conducting %d evaluations\n",OCpar.max_brack_eval);
	  printf("Gradient level parameter is %g\n",OCpar.grad_level);
	  printf("Optimization method is '%s'\n\n",OCpar.method ? "CG" : "L-BFGS");

	  tfval = OptimizeCG(interp);
  } else { /* L-BFGS */
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_lbfgs_eps", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.eps = 1.0e-8;
	  } else {
	     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.eps) != TCL_OK ) {
	       OCpar.eps = 1.0e-8;
	       printf("Warning! Could'n read par(oc_lbfgs_eps), using default value %f\n",OCpar.eps);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_lbfgs_tol_ls", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.tol = 1.0e-6;
	  } else {
	     if ( Tcl_GetDoubleFromObj( interp,objptr,&OCpar.tol) != TCL_OK ) {
	       OCpar.tol = 1.0e-6;
	       printf("Warning! Could'n read par(oc_lbfgs_tol_ls), using default value %f\n",OCpar.tol);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_lbfgs_max_ls_eval", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.max_brack_eval = 20;
	  } else {
	     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.max_brack_eval) != TCL_OK ) {
	       OCpar.max_brack_eval = 20;
	       printf("Warning! Could'n read par(oc_lbfgs_max_ls_eval), using default value %d\n",OCpar.max_brack_eval);
	     }
	  }
	  objptr = Tcl_GetVar2Ex( interp, "par", "oc_lbfgs_m", TCL_GLOBAL_ONLY);
	  if (!objptr) {
	     OCpar.lbfgs_m = 6;
	  } else {
	     if ( Tcl_GetIntFromObj( interp,objptr,&OCpar.lbfgs_m) != TCL_OK ) {
	       OCpar.lbfgs_m = 6;
	       printf("Warning! Could'n read par(oc_lbfgs_m), using default value %d\n",OCpar.lbfgs_m);
	     }
	  }

	  printf("Number of variables is %d\n",OCpar.ndim);
	  printf("Epsilon for convergence test is %g\n",OCpar.eps);
	  printf("Maximal number of iterations is %d\n",OCpar.nIterations);
	  printf("Optimization is terminated after %d iterations if target function is not higher than %g\n",OCpar.ncut,OCpar.cut);
	  printf("Every %d iterations procedure '%s' (for storing variables) is executed\n",OCpar.nreport, OCpar.VarSaveProc);
	  printf("The number of corrections to approximate the inverse hessian matrix is %d\n",OCpar.lbfgs_m);
	  printf("Tolerance for line-search (Armijo) is %g\n",OCpar.tol);
	  printf("Line-search will terminate after reaching %d evaluations\n",OCpar.max_brack_eval);
	  printf("Gradient level parameter is %g\n",OCpar.grad_level);
	  printf("Optimization method is '%s'\n\n",OCpar.method ? "CG" : "L-BFGS");

	  tfval = OptimizeLBFGS(interp);
  }

  TclSetResult(interp,"%lf",tfval);
 
  OCpar_destroy();
 
  return TCL_OK;
} 





/****
 * grad_shapes is used to define RFshapes for gradient calculation with respect to them
 *   can be called only from inside of 'gradient' procedure (here is just simple
 *    test on existence of OCpar.gradmode)
 ****/
 int tclGradShapes(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int *shapes;
  int i, slot;

  /* test for existence of gradmode */
  if (!OCpar.gradmode) 
     return TclError(interp,"ERROR: gradient mode not activated! oc_grad_shapes can be called only within 'gradient' procedure");
     
  /* read input parameters (only RFshapes) */
  if ( argc == 1)
    return TclError(interp,"Usage: oc_grad_shapes <RFshape> ?<RFshape> ...?");
    
  shapes = int_vector(argc-1);
  for (i=1; i<argc; i++) {
    if (Tcl_GetIntFromObj(interp,argv[i],&slot) != TCL_OK)
      return TclError(interp,"oc_grad_shapes: argument %d must be integer <RFshape>",i);
    if (!RFshapes[slot])
      return TclError(interp,"oc_grad_shapes: argument %d -> RFshape does not exist", i);
    shapes[i] = slot;
  }
  
  i = shapes[0];
  if (OCpar.grad_shapes != NULL) {
     /* printf("oc_grad_shapes is realocating OCpar.grad_shapes\n"); */
     free_int_vector(OCpar.grad_shapes);
  } 
  OCpar.grad_shapes = shapes;
  
  return TCL_OK;
}

/****
 *  this specifically turns gradient mode on or off
 *  can be used only within 'gradient' procedure, the test is insufficient...
 *  (purpose: to enable target function evaluation from 'gradient' procedure)
 ****/
 int tclGradModeSwitch(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  char *buf;

  /* test for existence of OCpar */
  if (OCpar.isinit != 1)
     return TclError(interp,"ERROR: gradient mode not activated! oc_gradmode can be called only within 'gradient' procedure");
     
  /* read input parameters (only RFshapes) */
  if ( argc != 2)
    return TclError(interp,"Usage: oc_gradmode on | off");

  buf = Tcl_GetString(argv[1]);
  if ( !strcmp(buf,"on") ) {
    OCpar.gradmode = 1;
    if (OCpar.verb) printf("GradMode turned  ON\n");
  } else if ( !strcmp(buf,"off") ) {
    OCpar.gradmode = 0;
    if (OCpar.verb) printf("GradMode turned OFF\n");
  }
  
  return TCL_OK;
}

/****
 * this command serves debugging gradient routine in input file
 * Procedure calculating gradients needs to be called either from oc_optimize 
 * or by help of this command. It makes sure that all necessary switches are set
 ****/ 
 int tclOCevaluate_grad(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{

  if (argc != 2)
    return TclError(interp,"usage: oc_evaluate_grad <name of procedure>");

  OCpar_initialize();
  test_pulseq_for_acqOC_prop(interp);
  OCpar.gradmode = 1;

  if (Tcl_EvalObjEx(interp,argv[1],TCL_EVAL_DIRECT) != TCL_OK) {
    fprintf(stderr,"oc_evaluate_grad error in evaluation of '%s':\n %s\n",Tcl_GetString(argv[1]),Tcl_GetStringResult(interp));
    OCpar_destroy();
    exit(1);
  }
  /* on successful evaluation, interp->result is already set by return value
     of 'gradient' procedure */  

  OCpar_destroy();

  return TCL_OK;
}

/****
 * C implementation of oc_acq_hermit
 *   ( result saved in current FID point )
 ****/
void acqOC_hermit(Sim_info *sim, Sim_wsp *wsp)
{
  complx c;
  int i, NN;
  
  /* disable when Hamiltonian is requested in block-diagonal form */
  /*  can be enabled when missing code is added */
  if (sim->block_diag) {
	  fprintf(stderr,"error: oc_acq_hermit - optimal control works only without blocking\n");
	  exit(1);
  }

  if (wsp->fid == NULL) {
    fprintf(stderr,"error: the 'oc_acq_hermit' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (wsp->curr_nsig + 1 > sim->ntot) {
    fprintf(stderr,"error: oc_acq_hermit overflow in fid points\n");
    exit(1);
  }

  _evolve_with_prop(sim,wsp);
  _reset_prop(sim,wsp);
  (wsp->curr_nsig)++;
  NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);
  for ( i=0; i<NN; i++) {
	  c = cm_trace(wsp->sigma[i % sim->Nfstart],wsp->fdetect[i % sim->Nfdetect]);
	  c.im = 0.0;
	  wsp->fid[wsp->curr_nsig+i*sim->ntot].re += c.re;
  }
}

/****
 * C implementation of oc_acq_nonhermit
 *   ( result saved in current FID point )
 *   normalization removed...
 ****/
void acqOC_nonhermit(Sim_info *sim, Sim_wsp *wsp)
{
  complx c;
  double r;
  int i, NN;

  /* disable when Hamiltonian is requested in block-diagonal form */
  /*  can be enabled when missing code is added */
  if (sim->block_diag) {
	  fprintf(stderr,"error: oc_acq_nonhermit - optimal control works only without blocking\n");
	  exit(1);
  }

  if (wsp->fid == NULL) {
    fprintf(stderr,"error: the 'oc_acq_nonhermit' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (wsp->curr_nsig + 1 > sim->ntot) {
    fprintf(stderr,"error: oc_acq_nonhermit overflow in fid points\n");
    exit(1);
  }

  _evolve_with_prop(sim,wsp);
  _reset_prop(sim,wsp);
  NN = (sim->Nfstart > sim->Nfdetect ? sim->Nfstart : sim->Nfdetect);

  /* normalize with sizes of start and detect operators -should we do?- */
  /* m_adjoint(puls->tmp,puls->fstart);
  c = m_trace(puls->tmp,puls->fstart);
  n1 = c.re;
  m_adjoint(puls->tmp,puls->fdetect);
  c = m_trace(puls->tmp, puls->fdetect);
  n2 = c.re;
  */

  (wsp->curr_nsig)++;
  for (i=0; i<NN; i++) {
	  c = cm_trace_adjoint(wsp->fdetect[i % sim->Nfdetect],wsp->sigma[i % sim->Nfstart]);
	  r = c.re*c.re+c.im*c.im;
  /* r /= sqrt(n1);
     r /= sqrt(n2);
  */
	  wsp->fid[wsp->curr_nsig + i*sim->ntot].re += r;
  }
}

/****
 * C implementation of oc_acq_prop
 *   ( result saved in current FID point )
 ****/
void acqOC_prop(Sim_info *sim, Sim_wsp *wsp, int Ud)
{
  complx c;
  double r;

  /* disable when Hamiltonian is requested in block-diagonal form */
  /*  can be enabled when missing code is added */
  if (sim->block_diag) {
	  fprintf(stderr,"error: oc_acq_prop - optimal control works only without blocking\n");
	  exit(1);
  }

  if (wsp->fid == NULL) {
    fprintf(stderr,"error: the 'oc_acq_prop' command can only be used when the computing method is 'direct'\n");
    exit(1);  
  }
  
  if (wsp->curr_nsig + 1 > LEN(wsp->fid)) {
    fprintf(stderr,"error: oc_acq_prop overflow in fid points\n");
    exit(1);
  }

  /* make the adjoint of desired propagator  */
  /* multiply it with optimized propagator and get the trace */
  c = blk_cm_trace_adjoint(wsp->STO[Ud], wsp->U, sim);
  r = c.re*c.re+c.im*c.im;
  //wsp->fid[++(wsp->curr_nsig)] = Complx(r,0.0);
  wsp->fid[++(wsp->curr_nsig)].re += r;
}


/****
 * Tcl implementation of oc_acq_hermit
 ****/
 int tclAcqOCHermit(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
	  Sim_info *sim = NULL;
	  Sim_wsp *wsp = NULL;

	  read_sim_pointers(interp, &sim, &wsp);

  /* acqOC can be called only once in pulseq, check it */
  if (wsp->curr_nsig != 0)
     return TclError(interp,"oc_acq_hermit can be used only once in pulse sequence!!!");
  
  /* this does not accept any arguments */
  if (argc != 1) 
     return TclError(interp,"usage: oc_acq_hermit\n (oc_acq_hermit does not accept any parameters)");
  
  /* test whether to calculate target function or gradient */
  if (OCpar.gradmode) {
     /* printf("do gradients\n"); */
      if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
    	  gradOC_hermit(sim,wsp);
      } else { /* GRAPE with advanced gradients */
    	  gradOC_hermit_2(sim,wsp);
      }
  } else {
     /* printf("do target function\n"); */
     acqOC_hermit(sim,wsp);
  }
    
  return TCL_OK;

}

/****
 * Tcl implementation of oc_acq_nonhermit
 ****/
 int tclAcqOCnonHermit(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
	  Sim_info *sim = NULL;
	  Sim_wsp *wsp = NULL;

	  read_sim_pointers(interp, &sim, &wsp);

	  /* acqOC can be called only once in pulseq, check it */
  if (wsp->curr_nsig != 0)
     return TclError(interp,"oc_acq_nonhermit can be used only once in pulse sequence!!!");
  
  /* this does not accept any arguments */
  if (argc != 1)
     return TclError(interp,"usage: oc_acq_nonhermit\n (oc_acq_nonhermit does not accept any parameters)");
  
  /* test whether to calculate target function or gradient */
  if (OCpar.gradmode) {
     /* printf("do gradients\n"); */
      if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
    	  gradOC_nonhermit(sim,wsp);
      } else { /* GRAPE with advanced gradients */
    	  gradOC_nonhermit_2(sim,wsp);
      }
  } else {
     /* printf("do target function\n"); */
     acqOC_nonhermit(sim,wsp);
  }
  
  return TCL_OK;

}

/****
 * Tcl implementation of oc_acq_prop
 ****/
 int tclAcqOCProp(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int ud;
  Sim_info *sim = NULL;
  Sim_wsp *wsp = NULL;

  read_sim_pointers(interp, &sim, &wsp);

  /* acqOC can be called only once in pulseq, check it */
  if (wsp->curr_nsig != 0)
     return TclError(interp,"oc_acq_prop can be used only once in pulse sequence!!!");
  
  /* check arguments */
  if (argc != 2)
     return TclError(interp,"usage: oc_acq_prop <desired propagator>");
  if (Tcl_GetIntFromObj(interp,argv[1],&ud) != TCL_OK)
       return TclError(interp,"oc_acq_prop: can't read first argument, must be integer");
  if (wsp->STO[ud] == NULL)
     return TclError(interp,"oc_acq_prop: desired propagator seems not to exist");
     
  /* test whether to calculate target function or gradient */
  if (OCpar.gradmode) {
     /* printf("do gradients\n"); */
      if (fabs(OCpar.grad_level-1)<1e-3) { /* original GRAPE ala Khaneja */
    	  gradOC_prop(sim,wsp,ud);
      } else { /* GRAPE with advanced gradients */
    	  gradOC_prop_2(sim,wsp,ud);
      }
  } else {
     /* printf("do target function\n"); */
     acqOC_prop(sim,wsp,ud);
  }

  return TCL_OK;

}

/****
 * Tcl implementation of oc_grad_add_energy_penalty [ units rad.s-1]
 ****/
 int tclGradPenalty(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int fidN, slot, ndim, Nslot, i, ii, idx;
  int *slotvec;
  double *lamvec, lam, xx, yy;
  extern FD** fd;
  extern int nfd;
  complx *gr;
  
  if (argc < 2 || ((argc % 2) != 0))
     return TclError(interp,"Usage: oc_grad_add_energy_penalty <fid> <RFshape> <lambda> ?<RFshape> <lambda>? ...");
     
  if (Tcl_GetIntFromObj(interp,argv[1],&fidN) == TCL_ERROR)
     return TclError(interp,"oc_grad_add_energy_penalty: argument 1 must be integer <data set>");

  if (fidN < 1 || fidN > nfd || fd[fidN] == NULL) 
    return TclError(interp,"oc_grad_add_energy_penalty: data set %d was not previously loaded\n",fidN);

  Nslot = argc/2-1;
  slotvec = int_vector(Nslot);
  lamvec = double_vector(Nslot);
  ii = 1;
  ndim = 0;
  for (i=1; i<=Nslot; i++) {
     ii++;
     if (Tcl_GetIntFromObj(interp,argv[ii],&slot) == TCL_ERROR) {
        free_int_vector(slotvec);
        free_double_vector(lamvec);
        return TclError(interp,"oc_grad_add_energy_penalty: argument %d must be integer <RFshape>",ii);
     }
     if (!RFshapes[slot]) {
        free_int_vector(slotvec);
        free_double_vector(lamvec);
        return TclError(interp,"oc_grad_add_energy_penalty: RFshape %d does not exist",slot);
     }
     ii++;
     if (Tcl_GetDoubleFromObj(interp,argv[ii],&lam) == TCL_ERROR) {
        free_int_vector(slotvec);
        free_double_vector(lamvec);
        return TclError(interp,"oc_grad_add_energy_penalty: unable to get double <lambda> from %s",argv[ii]);
     }
     slotvec[i] = slot;
     lamvec[i] = lam;
     ndim += RFshapes_len(slot);
  }

  if ( ndim != fd[fidN]->np ) {
     free_int_vector(slotvec);
     free_double_vector(lamvec);
     return TclError(interp,"oc_grad_add_energy_penalty: mismatch in lengths of gradients and RF shapes");
  }
  
  gr = fd[fidN]->data;
  idx = 0;
  /* IMPORTANT: this assumes that gradients in fid are {gr_x, gr_y} !!! */
  for (i=1; i<=Nslot; i++) {
     slot = slotvec[i];
     for (ii=1; ii<=RFshapes_len(slot); ii++) {
        xx = (RFshapes[slot][ii].ampl)*cos( (RFshapes[slot][ii].phase)*DEG2RAD );
        yy = (RFshapes[slot][ii].ampl)*sin( (RFshapes[slot][ii].phase)*DEG2RAD );
        idx++;
        gr[idx].re += lamvec[i]*4.0*xx*M_PI;
        gr[idx].im += lamvec[i]*4.0*yy*M_PI;
     }
  }
  
  free_int_vector(slotvec);
  free_double_vector(lamvec);
  
  return TCL_OK;

}

 /****
  * Tcl implementation of oc_grad_add_rflimit_penalty [ arbitrary units!!! ]
  ****/
  int tclGradLimitPenalty(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
 {
   int fidN, slot, ndim, Nslot, i, ii, idx;
   int *slotvec;
   double *lambdavec, *limitvec, lambda, limit, xx, yy;
   extern FD** fd;
   extern int nfd;
   complx *gr;

   if (argc < 5 || (((argc-2) % 3) != 0))
      return TclError(interp,"Usage: oc_grad_add_rflimit_penalty <fid> <RFshape> <limit> <lambda> ?<RFshape> <limit> <lambda>? ...");

   if (Tcl_GetIntFromObj(interp,argv[1],&fidN) == TCL_ERROR)
      return TclError(interp,"oc_grad_add_rflimit_penalty: argument 1 must be integer <data set>");

   if (fidN < 1 || fidN > nfd || fd[fidN] == NULL)
     return TclError(interp,"oc_grad_add_rflimit_penalty: data set %d was not previously loaded\n",fidN);

   Nslot = (argc-2)/3;
   slotvec = int_vector(Nslot);
   limitvec = double_vector(Nslot);
   lambdavec = double_vector(Nslot);
   ii = 1;
   ndim = 0;
   for (i=1; i<=Nslot; i++) {
      ii++;
      if (Tcl_GetIntFromObj(interp,argv[ii],&slot) == TCL_ERROR) {
         free_int_vector(slotvec);
         free_double_vector(lambdavec);
         free_double_vector(limitvec);
         return TclError(interp,"oc_grad_add_rflimit_penalty: argument %d must be integer <RFshape>",ii);
      }
      if (!RFshapes[slot]) {
         free_int_vector(slotvec);
         free_double_vector(lambdavec);
         free_double_vector(limitvec);
         return TclError(interp,"oc_grad_add_rflimit_penalty: RFshape %d does not exist",slot);
      }
      ii++;
      if (Tcl_GetDoubleFromObj(interp,argv[ii],&limit) == TCL_ERROR) {
         free_int_vector(slotvec);
         free_double_vector(lambdavec);
         free_double_vector(limitvec);
         return TclError(interp,"oc_grad_add_rflimit_penalty: unable to get double <limit> from %s",argv[ii]);
      }
      if (limit < TINY)  {
         free_int_vector(slotvec);
         free_double_vector(lambdavec);
         free_double_vector(limitvec);
         return TclError(interp,"oc_grad_add_rflimit_penalty: limit must be larger than zero (%g)",limit);
      }
      ii++;
      if (Tcl_GetDoubleFromObj(interp,argv[ii],&lambda) == TCL_ERROR) {
         free_int_vector(slotvec);
         free_double_vector(lambdavec);
         free_double_vector(limitvec);
         return TclError(interp,"oc_grad_add_rflimit_penalty: unable to get double <lambda> from %s",argv[ii]);
      }
      slotvec[i] = slot;
      lambdavec[i] = lambda;
      limitvec[i] = limit;
      ndim += RFshapes_len(slot);
   }

   if ( ndim != fd[fidN]->np ) {
      free_int_vector(slotvec);
      free_double_vector(lambdavec);
      free_double_vector(limitvec);
      return TclError(interp,"oc_grad_add_rflimit_penalty: mismatch in lengths of gradients and RF shapes");
   }

   gr = fd[fidN]->data;
   idx = 0;
   /* IMPORTANT: this assumes that gradients in fid are {gr_x, gr_y} !!! */
   for (i=1; i<=Nslot; i++) {
      slot = slotvec[i];
      for (ii=1; ii<=RFshapes_len(slot); ii++) {
    	 idx++;
    	 double am = RFshapes[slot][ii].ampl;
    	 if (am > limitvec[i]) {
    		 xx = am*cos( (RFshapes[slot][ii].phase)*DEG2RAD );
    		 yy = am*sin( (RFshapes[slot][ii].phase)*DEG2RAD );
    		 //double dum = am*am - limitvec[i]*limitvec[i];
             //gr[idx].re += lambdavec[i]*4.0*xx*dum*1e-12;
             //gr[idx].im += lambdavec[i]*4.0*yy*dum*1e-12; // remove factor kHz^4
    		 double dum = lambdavec[i]*4.0*(am - limitvec[i])/am*1e-6;  // remove factor kHz^2
             gr[idx].re += dum*xx;
             gr[idx].im += dum*yy;
    	 }
    	 // zero otherwise (for amplitudes below the limit)
      }
   }

   free_int_vector(slotvec);
   free_double_vector(lambdavec);
   free_double_vector(limitvec);

   return TCL_OK;

 }

  /****
   * Tcl implementation of oc_rflimit_penalty [ arbitrary units!!! ]
   ****/
   int tclLimitPenalty(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
  {
    int slot, ndim, Nslot, i, ii, idx;
    int *slotvec;
    double *lambdavec, *limitvec, lambda, limit, result;

    if (argc < 4 || (((argc-1) % 3) != 0))
       return TclError(interp,"Usage: oc_rflimit_penalty <RFshape> <limit> <lambda> ?<RFshape> <limit> <lambda>? ...");

    Nslot = (argc-1)/3;
    slotvec = int_vector(Nslot);
    limitvec = double_vector(Nslot);
    lambdavec = double_vector(Nslot);
    ii = 0;
    ndim = 0;
    for (i=1; i<=Nslot; i++) {
       ii++;
       if (Tcl_GetIntFromObj(interp,argv[ii],&slot) == TCL_ERROR) {
          free_int_vector(slotvec);
          free_double_vector(lambdavec);
          free_double_vector(limitvec);
          return TclError(interp,"oc_rflimit_penalty: argument %d must be integer <RFshape>",ii);
       }
       if (!RFshapes[slot]) {
          free_int_vector(slotvec);
          free_double_vector(lambdavec);
          free_double_vector(limitvec);
          return TclError(interp,"oc_rflimit_penalty: RFshape %d does not exist",slot);
       }
       ii++;
       if (Tcl_GetDoubleFromObj(interp,argv[ii],&limit) == TCL_ERROR) {
          free_int_vector(slotvec);
          free_double_vector(lambdavec);
          free_double_vector(limitvec);
          return TclError(interp,"oc_rflimit_penalty: unable to get double <limit> from %s",argv[ii]);
       }
       if (limit < TINY)  {
          free_int_vector(slotvec);
          free_double_vector(lambdavec);
          free_double_vector(limitvec);
          return TclError(interp,"oc_rflimit_penalty: limit must be larger than zero (%g)",limit);
       }
       ii++;
       if (Tcl_GetDoubleFromObj(interp,argv[ii],&lambda) == TCL_ERROR) {
          free_int_vector(slotvec);
          free_double_vector(lambdavec);
          free_double_vector(limitvec);
          return TclError(interp,"oc_rflimit_penalty: unable to get double <lambda> from %s",argv[ii]);
       }
       slotvec[i] = slot;
       lambdavec[i] = lambda;
       limitvec[i] = limit;
       ndim += RFshapes_len(slot);
    }

    idx = 0;
    result = 0;
    for (i=1; i<=Nslot; i++) {
       slot = slotvec[i];
       for (ii=1; ii<=RFshapes_len(slot); ii++) {
     	 idx++;
     	 double am = RFshapes[slot][ii].ampl;
     	 if (am > limitvec[i]) {
     		 //result += lambdavec[i]*(am*am - limitvec[i]*limitvec[i])*(am*am - limitvec[i]*limitvec[i]);
     		 result += lambdavec[i]*(am - limitvec[i])*(am - limitvec[i]);
     	 }
     	 // zero otherwise (for amplitudes below the limit)
       }
    }
    result *= 1e-6; //remove factor kHz^2

    free_int_vector(slotvec);
    free_double_vector(lambdavec);
    free_double_vector(limitvec);

    return TclSetResult(interp,"%.15g",result);

  }



/****
 * Tcl implementation of TF_line
 *      ( this is meant for debugging mainly!!! )
 ****/
 int tclTFLine(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int gr_slot, i, j, N, Nsh, slot;
  int Nvars=0;
  double g1, g2, gs, dum;
  double *dir;
  extern FD** fd;
  
  /* check arguments */
  if (argc < 5)
     return TclError(interp,"usage: TF_line <from> <to> <step> <RFshape> ?<RFshape> ...?");
  if (Tcl_GetDoubleFromObj(interp,argv[1],&g1) != TCL_OK)
       return TclError(interp,"TF_line: can't read second argument, must be double <from>");
  if (Tcl_GetDoubleFromObj(interp,argv[2],&g2) != TCL_OK)
       return TclError(interp,"TF_line: can't read third argument, must be double <to>");
  if (Tcl_GetDoubleFromObj(interp,argv[3],&gs) != TCL_OK)
       return TclError(interp,"TF_line: can't read fourth argument, must be double <step>");
  if (gs<=0.0)
     return TclError(interp,"TF_line: wrong step size!");

  Nsh = argc-4;
  if (Nsh<=0)
     return TclError(interp,"TF_line: there are no variable shapes!");
  OCpar_initialize();
  OCpar.var_shapes = int_vector(Nsh);
  OCpar.var_shapes_min = double_vector(Nsh);
  OCpar.var_shapes_max = double_vector(Nsh);
  OCpar.var_shapes_rmsmax = double_vector(Nsh);
  OCpar.grad_shapes = int_vector(Nsh);
  j = 0; 
  for (i=4; i<argc; i++) {
     if (Tcl_GetIntFromObj(interp,argv[i],&slot) != TCL_OK) {
          OCpar_destroy();
          return TclError(interp,"TF_line: can't read argument %d, must be integer <RFshape>",i);
     }
     if (!RFshapes[slot]) {
	  OCpar_destroy();
          return TclError(interp,"TF_line: argument %d -> RFshape does not exist", i);
     } 
     j++;
     OCpar.var_shapes[j] = slot;
     OCpar.var_shapes_min[j] = 0.0;
     OCpar.var_shapes_max[j] = 1.0e6;
     OCpar.var_shapes_rmsmax[j] = 1.0e6;
     OCpar.grad_shapes[j] = slot; 
     printf("reading RFshape slot %d\n",j);
  }

  test_pulseq_for_acqOC_prop(interp);
  
  /* evaluate gradient */
  gr_slot = evaluate_gradient(interp);
 printf("gradient done\n"); 
  /* store current RF shapes at safe place */
  RFelem **CURRshapes = (RFelem**)malloc(Nsh*sizeof(RFelem*));
  for (i=0; i<Nsh; i++) {
     N = RFshapes_len(OCpar.var_shapes[i+1]);
     Nvars += N;
     CURRshapes[i]=RFshapes_alloc(N);
     for (j=1; j<=N; j++) {
        CURRshapes[i][j].ampl = RFshapes[OCpar.var_shapes[i+1]][j].ampl;
        CURRshapes[i][j].phase = RFshapes[OCpar.var_shapes[i+1]][j].phase;
     }
  }
 printf("RFshapes stored separately\n");
 /* check length of gradients with total number of variables */
 if (Nvars != (fd[gr_slot]->np)) {
    fprintf(stderr,"TF_line error: number of gradients does not match total number of variables in RF shapes (%d versus %d)\n",Nvars,(fd[gr_slot]->np) );
    exit(1);
  }
 /* Nvars so far counts just elements in RF shapes. But vector dimension is double (aml, ph) */
 Nvars *=2;
printf("FID length and variables checked\n"); 
 /* initialize variable with gradient direction */
 dir = double_vector(Nvars);
 i = 0;
 for (j=1; j<=fd[gr_slot]->np; j++) {
     i++;
     dum = (fd[gr_slot]->data)[j].re;
     /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = dum;
     i++;
     dum = (fd[gr_slot]->data)[j].im;
     /* gradient calculated for maximizing, here we minimize. r[i]=-dum */
     dir[i] = dum;
 }
printf("gradient moved to direction\n");
 Tcl_ResetResult(interp);
 while (g1<=g2) {
    dum = func_to_min(g1,dir,CURRshapes,interp);
    TclAppendResult(interp,"%.15g",dum); /* doesn't work fo unknown reasons... */
    printf("   %.15g   %.15g\n",g1, dum);
    g1 += gs; 
 }
printf("calculations done\n");

  free_grad_fid(gr_slot,interp);
  free_double_vector(dir);
  for (i=0; i<Nsh; i++) {
    free((char *)CURRshapes[i]);
  } 
  free(CURRshapes);
  OCpar_destroy();
  
  return TCL_OK;

}

 // 1.12.2024 ZT: implementing GROUP optimizations, conversion of gradients
/* Tcl implementation of oc_grad_convert_to_group */
int tclGrad2GROUP(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int xygradslot, gradslot, amplslot, basisslot, row, Nrows, Ntime, k, NN;
  extern FD** fd;
  extern int nfd;
  double *xygradptrre,*xygradptrim, *basiselem;
  complx *gradelem;
  int gradidx, gradidxmax;

  /* check arguments */
  if ( (argc < 5) || ((argc-3)%2) )
    return TclError(interp,"usage: oc_grad_convert_to_group <RF x/y grad> <GROUP grad> <GROUP shape> <GROUP basis> ?<GROUP shape 2> <GROUP basis 2> ...?");
  NN = (argc-3)/2;  // this is number of GROUP shapes
  if (Tcl_GetIntFromObj(interp,argv[1],&xygradslot) == TCL_ERROR)
    return TclError(interp,"oc_grad_convert_to_group: first argument must be integer <RF x/y grad>");
  if (xygradslot < 1 || xygradslot > nfd || fd[xygradslot] == NULL)
    return TclError(interp,"oc_grad_convert_to_group: data set <RF x/y grad> does not exist\n");
  if (Tcl_GetIntFromObj(interp,argv[2],&gradslot) == TCL_ERROR)
    return TclError(interp,"oc_grad_convert_to_group: second argument must be integer <GROUP grad>");
  if (gradslot < 1 || gradslot > nfd || fd[gradslot] == NULL)
    return TclError(interp,"oc_grad_convert_to_group: data set <GROUP grad> does not exist\n");

  // real part of xygrad starts here and has step of 2, contains grad w.r.t. x-component
  xygradptrre = (double*)(fd[xygradslot]->data+1);
  xygradptrim = (double*)(fd[xygradslot]->data+1)+1;
  //  GROUP gradient element - ONE-based indexing!!!
  gradelem = fd[gradslot]->data;
  gradidxmax = fd[gradslot]->np;
  gradidx = 0;

  for (k=0; k<NN; k++) { // loop over GROUP shapes & basis pairs given as arguments
    if (Tcl_GetIntFromObj(interp,argv[3+k*2],&amplslot) == TCL_ERROR)
      return TclError(interp,"oc_grad_convert_to_group: (%d). argument must be integer <GROUP shape>",3+k*2);
    if (!RFshapes[amplslot])
      return TclError(interp,"oc_grad_convert_to_group: (%d). argument - shape does not exist",3+k*2);
    if (Tcl_GetIntFromObj(interp,argv[3+k*2+1],&basisslot) == TCL_ERROR)
      return TclError(interp,"oc_grad_convert_to_group: (%d). argument must be integer <GROUP basis>",3+k*2+1);
    if (!GROUPbasis[basisslot])
      return TclError(interp,"oc_grad_convert_to_group: (%d). argument - GROUP basis does not exist",3+k*2+1);
    Nrows = RFshapes_len(amplslot);
    if (Nrows != GROUPbasis[basisslot]->row)
      return TclError(interp,"oc_grad_convert_to_group: (%d)/(%d) arguments - dimension mismatch",3+k*2,3+k*2+1);
    Ntime = GROUPbasis[basisslot]->col;
    for (row=0; row<Nrows; row++) {
    	basiselem = GROUPbasis[basisslot]->data+row; // step is Nrows
    	gradidx++; // increment grad index and check for overflow
    	if (gradidx > gradidxmax)
    	  return TclError(interp,"oc_grad_convert_to_group: GROUP grad overflow");
    	gradelem[gradidx].re = cblas_ddot(Ntime,xygradptrre,2,basiselem,Nrows);
    	gradelem[gradidx].im = cblas_ddot(Ntime,xygradptrim,2,basiselem,Nrows);
    }
    // shift xygrad pointers to another RF shape
    xygradptrre += 2*Ntime;
    xygradptrim += 2*Ntime;
  } // end of loop over GROUP shapes & basis pairs

  return TCL_OK;
}

/* Tcl implementation of oc_rf_penalty */
int tclRFpenalty(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int pwr, i, slot, len;
  double result, scl;
  char *buf;
  const char *usage = "Usage: <double> oc_rf_penalty <RFshape> ?-power <int>? ?-scale <double>?";

  if ( (argc < 2) || (argc > 6) || (argc%2) )
      return TclError(interp,usage);
  if (Tcl_GetIntFromObj(interp,argv[1],&slot) == TCL_ERROR)
    return TclError(interp,"oc_rf_penalty: argument 1 must be integer <RFshape>");
  /* check for RFshape existence */
  if (!RFshapes[slot])
     return TclError(interp,"oc_rf_penalty: trying to access non-existing RFshape");
  len = RFshapes_len(slot);

  // default values
  pwr = 2;
  scl = 1.0;
  result = 0.0;
  for (i=2; i<argc; i++) {
    buf = Tcl_GetString(argv[i]);
	if (!strcmp(buf,"-power")) {
	  i++;
  	  if (Tcl_GetIntFromObj(interp,argv[i],&pwr) == TCL_ERROR)
		return TclError(interp,"oc_rf_penalty: argument following -power must be integer");
	} else if (!strcmp(buf,"-scale")) {
	  i++;
	  if (Tcl_GetDoubleFromObj(interp,argv[i],&scl) == TCL_ERROR)
		return TclError(interp,"oc_rf_penalty: argument following -scale must be double");
	} else
	  return TclError(interp,usage);
  }
  // calculation
  for (i=1; i<=len; i++) {
    result += pow(RFshapes[slot][i].ampl,pwr);
  }
  result *= scl;

  return TclSetResult(interp,"%.15g",result);
}

/* Tcl implementation of oc_rf_penalty_grad */
int tclRFpenaltyGrad(ClientData data,Tcl_Interp* interp, int argc, Tcl_Obj *argv[])
{
  int i, j, len, slot, gradslot, gradidx, gradidxmax, pwr;
  double scl, am, ph, wx, wy;
  extern FD** fd;
  extern int nfd;
  complx *gradelem;
  char *buf;
  const char *usage = "Usage: oc_rf_penalty_grad <gradFID> <RFshape> ?-power <int>? ?-scale <double>? ?<RFshape> ?-power <int>? ?-scale <double>? ...";

  if (argc < 3)
	  return TclError(interp,usage);
  if (Tcl_GetIntFromObj(interp,argv[1],&gradslot) == TCL_ERROR)
    return TclError(interp,"oc_rf_penalty_grad: first argument must be integer <x/y grad>");
  if (gradslot < 1 || gradslot > nfd || fd[gradslot] == NULL)
    return TclError(interp,"oc_rf_penalty_grad: data set <x/y grad> does not exist\n");
  gradelem = fd[gradslot]->data;
  gradidxmax = fd[gradslot]->np;
  gradidx = 0;

  for (i=2; i<argc; i++) {
    // default values
    pwr = 2;
    scl = 1.0;
    if (Tcl_GetIntFromObj(interp,argv[i],&slot) != TCL_OK)
      return TclError(interp,"oc_rf_penalty_grad: error in argument %d, expecting integer <RFshape>",i);
    // check if there is another argument
    while (i+1 < argc) {
      buf = Tcl_GetString(argv[i+1]);
      if (!strcmp(buf,"-power")) {
        i += 2;
        if (i >= argc)
          return TclError(interp,"oc_rf_penalty_grad: missing argument after -power");
        if (Tcl_GetIntFromObj(interp,argv[i],&pwr) == TCL_ERROR)
      	  return TclError(interp,"oc_rf_penalty_grad: argument %d, following -power, must be integer", i);
      } else if (!strcmp(buf,"-scale")) {
        i += 2;
        if (i >= argc)
          return TclError(interp,"oc_rf_penalty_grad: missing argument after -scale");
        if (Tcl_GetDoubleFromObj(interp,argv[i],&scl) == TCL_ERROR)
      	  return TclError(interp,"oc_rf_penalty_grad: argument %d, following -scale, must be double", i);
      } else {
        // next argument seems to be a shape, forget it
    	break;
	  }
    }
    // process the gradient of this shape
    len = RFshapes_len(slot);
    for (j=1; j<=len; j++) {
      am = RFshapes[slot][j].ampl;
      ph = RFshapes[slot][j].phase;
      wx = am*cos(ph*DEG2RAD);
      wy = am*sin(ph*DEG2RAD);
      gradidx++; // increment grad index and check for overflow
      if (gradidx > gradidxmax)
        return TclError(interp,"oc_rf_penalty_grad: grad overflow");
      gradelem[gradidx].re += scl*pwr*pow(am,pwr-2)*wx;
      gradelem[gradidx].im += scl*pwr*pow(am,pwr-2)*wy;
    }

  }

  return TCL_OK;
}

void tclcmd_OCroutines(Tcl_Interp* interp) {

Tcl_CreateCommand(interp,"oc_optimize",(Tcl_CmdProc *)tclOCoptimize2,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_acq_hermit",(Tcl_ObjCmdProc *)tclAcqOCHermit,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_acq_nonhermit",(Tcl_ObjCmdProc *)tclAcqOCnonHermit,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_acq_prop",(Tcl_ObjCmdProc *)tclAcqOCProp,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_evaluate_grad",(Tcl_ObjCmdProc *)tclOCevaluate_grad,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_grad_shapes",(Tcl_ObjCmdProc *)tclGradShapes,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_grad_add_energy_penalty",(Tcl_ObjCmdProc *)tclGradPenalty,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_gradmode",(Tcl_ObjCmdProc *)tclGradModeSwitch,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"TF_line",(Tcl_ObjCmdProc *)tclTFLine,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_rflimit_penalty",(Tcl_ObjCmdProc *)tclLimitPenalty,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_grad_add_rflimit_penalty",(Tcl_ObjCmdProc *)tclGradLimitPenalty,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
// 1.12.2024 ZT: implementing GROUP optimizations, basis creation
Tcl_CreateObjCommand(interp,"oc_grad_convert_to_group",(Tcl_ObjCmdProc *)tclGrad2GROUP,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_rf_penalty",(Tcl_ObjCmdProc *)tclRFpenalty,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
Tcl_CreateObjCommand(interp,"oc_rf_penalty_grad",(Tcl_ObjCmdProc *)tclRFpenaltyGrad,(ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

}
