/*
    Reading and setup of the spinsystem
    Copyright (C) 1999 Mads Bak, Jimmy T. Rasmussen
    2009 ZT modification with new matrix types
    
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
    
    Reads the spin system from the 'spinsys' section in the input
    file and creates the Hamiltonian for the spin system.
    
     Called from sim.c the setup and performs the simulation
*/

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <assert.h>
#include "matrix.h"
#include "spinsys.h"
#include "wigner.h"
#include "sim.h"
#include "defs.h"
#include "cm.h"
#include "blockdiag.h"
#include "ham.h"


/* based on NRC QuickSort algorithm */
void sort_mzmap(int Nmz, int row, int c_ini, int n_c, int *mzmap, int *permvec)
{
	long i, indxt, ir, itemp, j, k, l, *indx;
	int jstack = 0, *istack, *idum;
	int a;
	const int nstack = 64;

//	printf("sort_mzmap: row %d, from column %d, %d elements ... \n",row,c_ini,n_c);

	istack = malloc(nstack*sizeof(int));
	indx = malloc((n_c)*sizeof(unsigned long));
	for (i=0; i<n_c; i++) indx[i] = i;
	l = 0;
	ir = n_c-1;
	for (;;) {
		if (ir-l < 10) {
			/* insertion sort when subarray small */
			for (j=l+1; j<=ir; j++) {
				indxt = indx[j];
				a = mzmap[row+(indxt+c_ini)*Nmz];
				for (i=j-1; i>=l; i--) {
					if (mzmap[row+(indx[i]+c_ini)*Nmz] <= a) break;
					indx[i+1] = indx[i];
				}
				indx[i+1] = indxt;
			}
			if (jstack == 0) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		} else {
			k = (l +ir) >> 1;
			itemp = indx[k]; indx[k] = indx[l+1]; indx[l+1] = itemp;
			if (mzmap[row+(indx[l]+c_ini)*Nmz] > mzmap[row+(indx[ir]+c_ini)*Nmz]) {
				itemp = indx[l]; indx[l] = indx[ir]; indx[ir] = itemp;
			}
			if (mzmap[row+(indx[l+1]+c_ini)*Nmz] > mzmap[row+(indx[ir]+c_ini)*Nmz]) {
				itemp = indx[l+1]; indx[l+1] = indx[ir]; indx[ir] = itemp;
			}
			if (mzmap[row+(indx[l]+c_ini)*Nmz] > mzmap[row+(indx[l+1]+c_ini)*Nmz]) {
				itemp = indx[l]; indx[l] = indx[l+1]; indx[l+1] = itemp;
			}
			i = l+1;
			j = ir;
			indxt = indx[l+1];
			a = mzmap[row+(indxt+c_ini)*Nmz];
			for (;;) {
				do i++; while (mzmap[row+(indx[i]+c_ini)*Nmz] < a);
				do j--; while (mzmap[row+(indx[j]+c_ini)*Nmz] > a);
				if (j < i) break;
				itemp = indx[i]; indx[i] = indx[j]; indx[j] = itemp;
			}
			indx[l+1] = indx[j];
			indx[j] = indxt;
			jstack += 2;
			if (jstack > nstack) {
				fprintf(stderr,"Error in sort_mzmap: nstack too small\n");
				exit(1);
			}
			if (ir-i+1 >= j-l) {
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j-1;
			} else {
				istack[jstack] = j-1;
				istack[jstack-1] = l;
				l = i;
			}
		}
	}
	idum = malloc(n_c*sizeof(int));
	for (i=0; i<n_c; i++) idum[i] = permvec[i+c_ini+1];
	for (i=0; i<n_c; i++) permvec[i+c_ini+1] = idum[indx[i]];
	for (j=0; j< Nmz; j++) {
		for (i=0; i<n_c; i++) idum[i] = mzmap[j +(i+c_ini)*Nmz];
		for (i=0; i<n_c; i++) mzmap[j +(i+c_ini)*Nmz] = idum[indx[i]];
	}
	free(idum);
	free(indx);
	free(istack);

}

void readsys(Tcl_Interp* interp,Sim_info* s)
{
  SpinSys* ss;
  int ver, Nnuc, Nchan, Nelm, NN, i, j, n1, n2, Jiso, Jani;
  int Nshift=0, Nj=0, Ndip=0, Nquad=0, Nmix=0, Ngten=0, Nhf=0, Nhex=0, Nedd=0, Naniso;
  char *buf;
  const char **names, **spins; //changeable pointer to unchangeable char
  Tcl_Obj **vals;
  Shift *csptr;
  Dipole *ddptr;
  Quadrupole *qptr;
  Jcoupling *jptr;
  Mixing *mptr;
  Gtensor *gptr;
  Hyperfine *hfptr;
  Heisenberg_exchange *hexptr;
  Electron_Dipole *eddptr;
  double data[8];
  int *spin_code, *mzmap=NULL, *permvec, Nmz, reps, mzpos, Ival, Icur, c_ini, n_c, Nblk, *blk_dims=NULL;
  int Hiso_isdiag, HQ_isdiag, Nbasis;


  ss=s->ss;
  ss_initialize(ss);

  ver = (verbose & VERBOSE_SPINSYS);
  
  if (ver) printf("Reading spinsystem:\n");

  /* get number of nuclei and their names */
  buf=Tcl_GetVar2(interp,"spinsysres","nuclei",TCL_GLOBAL_ONLY);
  if (!buf) {
     fprintf(stderr,"Error: the spinsystem must contain 'nuclei' field\n");
     exit(1);
  }
  if (Tcl_SplitList(interp, buf, &Nnuc, &spins) != TCL_OK) {
     fprintf(stderr,"Error: readsys cannot decompose list spinsysres(nuclei)\n");
     exit(1);
  }
  if (ver) printf( "Number of spins : %d\n",Nnuc);
  for (i=0; i<Nnuc; i++) {
	  if (ver) {
		  if (!strncmp(spins[i],"e",1)) {
			  printf("  electron %3d : '%s'\n",i+1,spins[i]);
		  } else {
			  printf("   nucleus %3d : '%s'\n",i+1,spins[i]);
		  }
	  }
	  ss_addspin(ss,spins[i]);
  }

  /* read in channels info */
  buf=Tcl_GetVar2(interp,"spinsysres","channels",TCL_GLOBAL_ONLY);
  if (!buf) {
     fprintf(stderr,"Error: the spinsystem must contain 'channels' field\n");
     exit(1);
  }
  if (Tcl_SplitList(interp, buf, &Nchan, &names) != TCL_OK) {
     fprintf(stderr,"Error: readsys cannot decompose list spinsysres(channels)\n");
     exit(1);
  }
  if (Nchan > MAXCHAN) {
	  fprintf(stderr,"Error: maximal number of channels exceeded (MAXCHAN=%d)\n",MAXCHAN);
	  exit(1);
  }
  ss->nchan = Nchan;
  NN = 0;
  spin_code = int_vector(Nnuc);
  iv_zero(spin_code);
  for (i=1; i<=Nchan; i++) {
	  strcpy(ss->channames[i],names[i-1]);
	  ss->nchanelem[i] = 0;
	  for (j=1; j<=Nnuc; j++) {
		  //DEBUGPRINT("readsys: comparing (%d)-'%s' with (%d)-'%s'\n",i,ss->channames[i],j,spins[j-1]);
		  if (!strcmp(spins[j-1],ss->channames[i])) {
			  ss->nchanelem[i]++;
			  ss->chan[i][ss->nchanelem[i]] = j;
			  NN++;
			  spin_code[j] = i;
		  }
	  }
	  if (ss->nchanelem[i] == 0) {
		  fprintf(stderr,"Error: empty channel '%s'\n",names[i-1]);
		  exit(1);
	  }
	  // rf channel frequency taken from the first nucleus on that channel
	  ss->chanfreq[i] = ss_larmor_frequency(s,ss->chan[i][1]);

  }

  /* gather info for Hamiltonian block diagonalization */
  //s->permvec = NULL;
  //s->blk_dims = NULL;
  //s->Nmz = 0;
  if (s->block_diag) {
	  Nmz = Nchan;
	  if (Nnuc != NN) {
		  /* some spins do not have assigned channel but I want to use them */
		  for (i=1; i<= Nnuc; i++) {
			  if (spin_code[i] != 0) continue;
			  Nmz++;
			  spin_code[i] = Nmz;
			  ss->nchanelem[Nmz] = 1;
			  ss->chan[Nmz][ss->nchanelem[Nmz]] = i;
			  for (j=i+1; j<=Nnuc; j++) {
				  if (!strcmp(spins[j-1],spins[i-1])) {
					  spin_code[j] = Nmz;
					  ss->nchanelem[Nmz]++;
					  ss->chan[Nmz][ss->nchanelem[Nmz]] = j;
				  }
			  }
		  }
	  }
	  /* test output BEGIN
	  printf("readsys block_diag info:\n");
	  for (i=1; i<=Nmz; i++) {
		  printf("channel %d:",i);
		  for (j=1; j<=ss->nchanelem[i]; j++) {
			  printf(" %d",ss->chan[i][j]);
		  }
		  printf("\n");
	  }
	   test output END */
	  /* fill map of quantum numbers mz*2 on each channel */
	  mzmap = calloc(Nmz*ss->matdim,sizeof(int));
	  if (mzmap == NULL) {
		  fprintf(stderr,"Error: readsys can not allocate mzmap\n");
		  exit(1);
	  }
	  reps = 1;
	  for (i=Nnuc; i>0; i--) {
		  Icur = (int)round(2*ss->iso[i]->spin);
		  //printf("nucleus %d has spin %d\n",i,Icur);
		  mzpos = 0;
		  while (mzpos < ss->matdim) {
			  for (Ival=Icur; Ival>=-Icur; Ival-=2) {
				  for (j=0; j<reps; j++) {
					  assert(mzpos < ss->matdim);
					  mzmap[spin_code[i]-1 + mzpos*Nmz] += Ival;
					  //printf("mzmap[%d][%d] += %d\n",spin_code[i]-1,mzpos,Ival);
					  mzpos++;
				  }
			  }
		  }
		  reps *= Icur+1;
	  }
	  /* test output BEGIN
	  printf("readsys block_diag info:\n\tmz map\n\t======\n");
	  for (i=0; i<Nmz; i++) {
		  for (j=0; j<ss->matdim; j++) {
			  printf("%3d ",mzmap[i+j*Nmz]);
		  }
		  printf("\n");
	  }
	   test output END */
	  /* sort mz map and get product basis permutation vector */
	  permvec = int_vector(ss->matdim);
	  for (i=1; i<=ss->matdim; i++) permvec[i] = i;
	  sort_mzmap(Nmz,0,0,ss->matdim,mzmap,permvec);
	  for (i=1; i<Nmz; i++) {
		  c_ini = 0;
		  n_c = 0;
		  while (c_ini<ss->matdim) {
			  do n_c++; while ((c_ini+n_c)<ss->matdim && mzmap[i-1+(c_ini+n_c)*Nmz] == mzmap[i-1+(c_ini+n_c-1)*Nmz]);
			  if (n_c != 1 ) sort_mzmap(Nmz,i,c_ini,n_c,mzmap,permvec);
			  c_ini += n_c;
			  n_c = 0;
		  }
	  }
	  /* test output BEGIN
	  printf("readsys block_diag info:\n\tSORTED mz map\n\t======\n");
	  for (i=0; i<Nmz; i++) {
		  for (j=0; j<ss->matdim; j++) {
			  printf("%3d ",mzmap[i+j*Nmz]);
		  }
		  printf("\n");
	  }
	  printf("permutation vector\n");
	  for (i=1;i<=ss->matdim; i++) printf("%d ",permvec[i]);
	  printf("\n\n");
	   test output END */
	  /* prepare information vector about block structure */
	  Nblk = 1;
	  for (i=1; i<=Nmz; i++) Nblk *= ((int)round(2*ss->iso[ss->chan[i][1]]->spin))*ss->nchanelem[i]+1;
	  blk_dims = int_vector(Nblk);
	  Nblk = 0;
	  c_ini = 0;
	  n_c = 0;
	  while (c_ini<ss->matdim) {
		  do n_c++; while ((c_ini+n_c)<ss->matdim && mzmap[Nmz-1+(c_ini+n_c)*Nmz] == mzmap[Nmz-1+(c_ini+n_c-1)*Nmz]);
		  c_ini += n_c;
		  Nblk++;
		  assert(Nblk <= LEN(blk_dims));
		  blk_dims[Nblk] = n_c;
		  n_c = 0;
	  }
	  /* test output BEGIN
	  printf("Number of blocks = %d\ndims =",LEN(blk_dims));
	  for (i=1; i<=LEN(blk_dims);i++) printf(" %d",blk_dims[i]);
	  printf("\n\n");
	   text output END */
	  //s->permvec = permvec;
	  //s->blk_dims = blk_dims;
	  s->Nmz = Nmz;
	  s->Nbasis = Nbasis = 1<<Nmz;

	  /* generation of permutation table */
	  s->perm_table = (int**)malloc(Nbasis*Nbasis*sizeof(int*));
	  s->dims_table = (int**)malloc(Nbasis*sizeof(int*));
	  if (s->perm_table == NULL) {
		  fprintf(stderr,"Error: readsys can not allocate perm_table\n");
		  exit(1);
	  }
	  if (s->dims_table == NULL) {
		  fprintf(stderr,"Error: readsys can not allocate dims_table\n");
		  exit(1);
	  }
	  int *mzmap_tmp = malloc(Nmz*ss->matdim*sizeof(int));
	  int *Pvec, *blkdm;
	  for (i=1; i<Nbasis-1; i++) {
		  //printf("\n\n  basis index %d\n=======================\n",i);
		  memcpy(mzmap_tmp,mzmap,Nmz*ss->matdim*sizeof(int));
		  int code = i;
		  s->perm_table[(Nbasis-1) + i*Nbasis] = Pvec = int_vector(ss->matdim);
		  for (j=1; j<=ss->matdim; j++) Pvec[j] = j;
		  int r1 = 0, r2 = 0;
		  while ( !(code & (1<<0)) ) {
			  r1++;
			  code >>= 1;
		  }
		  //printf("\t\tr1 = %d\n",r1);
		  sort_mzmap(Nmz,r1,0,ss->matdim,mzmap_tmp,Pvec);
		  Nblk = ((int)round(2*ss->iso[ss->chan[r1+1][1]]->spin))*ss->nchanelem[r1+1]+1;
		  r2 = r1 + 1;
		  code >>= 1;
		  while (code != 0) {
			  while ( !(code & (1<<0)) ) {
				  r2++;
				  code >>= 1;
			  }
			  //printf("\t\tr2 = %d",r2);
			  c_ini = 0;
			  n_c = 0;
			  while (c_ini<ss->matdim) {
				  do n_c++; while ((c_ini+n_c)<ss->matdim && mzmap_tmp[r1+(c_ini+n_c)*Nmz] == mzmap_tmp[r1+(c_ini+n_c-1)*Nmz]);
				  if (n_c != 1 ) sort_mzmap(Nmz,r2,c_ini,n_c,mzmap_tmp,Pvec);
				  c_ini += n_c;
				  n_c = 0;
			  }
			  r1 = r2;
			  code >>= 1;
			  r2++;
			  Nblk *= ((int)round(2*ss->iso[ss->chan[r2][1]]->spin))*ss->nchanelem[r2]+1;
		  }
		  /* prepare information vector about block structure */
		  //printf("\tNumber of blocks = %d\n",Nblk);
		  s->dims_table[i] = blkdm = int_vector(Nblk);
		  Nblk = 0;
		  c_ini = 0;
		  n_c = 0;
		  while (c_ini<ss->matdim) {
			  do n_c++; while ((c_ini+n_c)<ss->matdim && mzmap_tmp[r2-1+(c_ini+n_c)*Nmz] == mzmap_tmp[r2-1+(c_ini+n_c-1)*Nmz]);
			  c_ini += n_c;
			  Nblk++;
			  assert(Nblk <= LEN(blkdm));
			  blkdm[Nblk] = n_c;
			  n_c = 0;
		  }
		  //printf("\tpermvec = "); for (j=1;j<=ss->matdim; j++) printf("%d ",Pvec[j]);
		  //printf("\n\tdims = "); for (j=1;j<=Nblk; j++) printf("%d ",blkdm[j]);

		  /* inverse permutation */
		  int *permvec_inv;
		  s->perm_table[i + (Nbasis-1)*Nbasis] = permvec_inv = int_vector(ss->matdim);
		  for (j=1; j<=ss->matdim; j++) permvec_inv[Pvec[j]] = j;

	  }
	  free(mzmap_tmp);
	  free(mzmap);
	  /* all other permutations */
	  int k;
	  int *permv1, *permv2;
	  for (i=1; i<Nbasis-1; i++) {
		  permv1 = s->perm_table[i + (Nbasis-1)*Nbasis];
		  assert(LEN(permv1) == ss->matdim);
		  for (j=i+1; j<Nbasis-1; j++) {
			  permv2 = s->perm_table[(Nbasis-1) + j*Nbasis];
			  assert(LEN(permv2) == ss->matdim);
			  s->perm_table[i + j*Nbasis] = Pvec = int_vector(ss->matdim);
			  for (k=1; k<=ss->matdim; k++) Pvec[k] = permv1[permv2[k]];
			  s->perm_table[j + i*Nbasis] = permv2 = int_vector(ss->matdim);
			  for (k=1; k<=ss->matdim; k++) permv2[Pvec[k]] = k;
			  /* check inverse perm with combined perm
			  printf("\npermutation %d->%d\n",j,i);
			  for (k=1; k<=ss->matdim; k++) printf("%d ",permv2[k]);
			  permvec = s->perm_table[j+(Nbasis-1)*Nbasis];
			  permv2 = s->perm_table[(Nbasis-1)+i*Nbasis];
			  printf("\npermutation %d->%d->%d\n",j,Nbasis-1,i);
			  for (k=1; k<=ss->matdim; k++) printf("%d ",permvec[permv2[k]]);
			  printf("\n");
			  **********************/
		  }
	  }
	  /* permutations from and to basis 0 should also be included */
	  s->perm_table[0 + (Nbasis-1)*Nbasis] = permv1 = permvec;
	  //memcpy(permv1+1,s->permvec+1,LEN(permv1)*sizeof(int));
	  s->perm_table[(Nbasis-1) + 0*Nbasis] = permv2 = int_vector(ss->matdim);
	  for (k=1; k<=ss->matdim; k++) permv2[permv1[k]] = k;
	  //permv1 = s->permvec;
	  for (i=1; i<Nbasis-1; i++) {
		  permv2 = s->perm_table[(Nbasis-1) + i*Nbasis];
		  assert(LEN(permv2) == ss->matdim);
		  s->perm_table[0 + i*Nbasis] = Pvec = int_vector(ss->matdim);
		  for (k=1; k<=ss->matdim; k++) Pvec[k] = permv1[permv2[k]];
		  s->perm_table[i + 0*Nbasis] = permv2 = int_vector(ss->matdim);
		  for (k=1; k<=ss->matdim; k++) permv2[Pvec[k]] = k;
	  }
	  for (i=0; i<Nbasis; i++) s->perm_table[i+i*Nbasis] = NULL;
	  s->dims_table[Nbasis-1] = blk_dims;
	  //memcpy(s->dims_table[0]+1,s->blk_dims+1,LEN(s->blk_dims)*sizeof(int));
	  s->dims_table[0] = blkdm = int_vector(1);
	  blkdm[1] = ss->matdim;
  } else {
	  s->Nmz = 0;
	  s->Nbasis = 0;
	  s->perm_table = NULL;
	  s->dims_table = NULL;
  }





  Tcl_Free((char*)names);
  Tcl_Free((char*)spins);
  free_int_vector(spin_code);

  if (ver) {
    printf("Number of channels : %d\n",ss->nchan);
    for (i=1;i<=ss->nchan;i++) {
      printf( "   channel %3d  : '%s' contain",i,ss->channames[i]);
      for (j=1;j<=ss->nchanelem[i];j++) {
    	 if (ss->iso[ss->chan[i][j]]->number == 0) {
    		 printf(" electron(%d,%s)",ss->chan[i][j],ss->iso[ss->chan[i][j]]->name);
    	 } else {
    		 printf(" nucleus(%d,%d%s)",ss->chan[i][j],ss->iso[ss->chan[i][j]]->number,ss->iso[ss->chan[i][j]]->name);
    	 }
      }
      printf(", carrier frequency %.2f Hz\n",ss->chanfreq[i]/2/M_PI);
    }
  }


  if (Tcl_EvalEx(interp,"array names spinsysres",-1,TCL_EVAL_GLOBAL) != TCL_OK) {
     fprintf(stderr,"Error: readsys can not get spinsysres fields:\n %s\n",Tcl_GetStringResult(interp));
     exit(1);
  }
  buf = Tcl_GetStringResult(interp);
  if (Tcl_SplitList(interp, buf, &Nelm, &names) != TCL_OK) {
     fprintf(stderr,"Error: readsys cannot decompose spinsysres array names from a list\n");
     exit(1);
  }

  for (i=0; i<Nelm; i++) {
    if (!strncmp(names[i],"nuclei",6)) continue;
    if (!strncmp(names[i],"channels",8)) continue;
    if (!strncmp(names[i],"shift",5)) {
    	Nshift++;
    	continue;
    }
    if (!strncmp(names[i],"jcoupling",9)) {
    	Nj++;
    	continue;
    }
    if (!strncmp(names[i],"dipole",6)) {
    	Ndip++;
    	continue;
    }
    if (!strncmp(names[i],"quadrupole",10)) {
    	Nquad++;
    	continue;
    }
    if (!strncmp(names[i],"mixing",6)) {
    	Nmix++;
    	continue;
    }
    if (!strncmp(names[i],"gtensor",7)) {
    	Ngten++;
    	continue;
    }
    if (!strncmp(names[i],"hyperfine",9)) {
    	Nhf++;
    	continue;
    }
    if (!strncmp(names[i],"heisenberg",10)) {
    	Nhex++;
    	continue;
    }
    if (!strncmp(names[i],"edipole",7)) {
    	Nedd++;
    	continue;
    }
  }
  if ( (s->frame == LABFRAME) || (s->frame == DNPFRAME) ) {
	  Nshift = ss->nspins; // all spins need Larmor frequency added to isotropic shift
  }
  s->CS = (Shift**)malloc(Nshift*sizeof(Shift*));
  s->nCS = 0;
  s->J = (Jcoupling**)malloc(Nj*sizeof(Jcoupling*));
  s->nJ = 0;
  s->DD = (Dipole**)malloc(Ndip*sizeof(Dipole*));
  s->nDD = 0;
  s->Q = (Quadrupole**)malloc(Nquad*sizeof(Quadrupole*));
  s->nQ = 0;
  s->MIX = (Mixing**)malloc(Nmix*sizeof(Mixing*));
  s->nMIX = 0;
  s->G = (Gtensor**)malloc(Ngten*sizeof(Gtensor*));
  s->nG = 0;
  s->HF = (Hyperfine**)malloc(Nhf*sizeof(Hyperfine*));
  s->nHF = 0;
  s->HEx = (Heisenberg_exchange**)malloc(Nhex*sizeof(Heisenberg_exchange*));
  s->nHEX = 0;
  s->EDD = (Electron_Dipole**)malloc(Nhex*sizeof(Electron_Dipole*));
  s->nEDD = 0;


  s->matdim = ss->matdim;
  s->Hint_isdiag = 1;
  Hiso_isdiag = 1;
  HQ_isdiag = 1;
  Naniso = 0;

  for (i=0; i<Nelm; i++) {
    if (!strncmp(names[i],"nuclei",6)) continue;
    if (!strncmp(names[i],"channels",8)) continue;
    if (Tcl_ListObjGetElements(interp,Tcl_GetVar2Ex(interp,"spinsysres",names[i],TCL_GLOBAL_ONLY), &NN, &vals) != TCL_OK) {
        fprintf(stderr,"Error: readsys cannot decompose list of parameters for spinsysres(%s)\n",names[i]);
    	exit(1);
    }
    if (!strncmp(names[i],"shift",5)) {
    	if (NN != 7) {
    		fprintf(stderr,"Error: reading shift - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - shift - int conversion failure\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins) {
    		fprintf(stderr,"Error: spinsys: shift %d  out of defined range\n",n1);
    		exit(1);
    	}
    	/* for electron, throw an error */
    	if (ss->iso[n1]->number == 0) {
    		fprintf(stderr,"Error: shift can not be defined for electron spin %d\n",n1);
    		exit(1);
    	}
    	if (shift_exist(s,n1) >= 0) {
    		fprintf(stderr,"Error: spinsys: shift already exists for nucleus %d\n",n1);
    		exit(1);
    	}
    	for (j=1; j<NN; j++) {
    		if (Tcl_GetDoubleFromObj(interp,vals[j],&(data[j-1])) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - shift - double conversion failure\n");
    			exit(1);
    		}
    	}
    	/* when zero amplitudes, do nothing */
    	if ( (fabs(data[0]) < TINY) && (fabs(data[1]) < TINY) ) continue;
    	s->CS[(s->nCS)++] = csptr = (Shift*)malloc(sizeof(Shift));
    	csptr->nuc = n1; csptr->eta = data[2]; csptr->pas[0] = data[3]; csptr->pas[1] = data[4]; csptr->pas[2] = data[5];
    	if (ss_gamma(ss,n1) > 0) {
    		csptr->iso = -data[0]*2.0*M_PI; csptr->delta = -data[1]*2.0*M_PI;
    	} else {
    		csptr->iso = data[0]*2.0*M_PI; csptr->delta = data[1]*2.0*M_PI;
    	}
    	if (fabs(csptr->delta) > TINY) {
    		Naniso++;
    		csptr->Rmol = Dtensor2(SQRT2BY3*csptr->delta, csptr->eta);
    		csptr->Rmol = wig2roti(csptr->Rmol, csptr->pas[0], csptr->pas[1], csptr->pas[2]);
    	} else {
    		csptr->Rmol = NULL;
    	}
    	csptr->T = NULL;
    	csptr->T2q[0] = csptr->T2q[1] = csptr->T2q[2] = csptr->T2q[3] = csptr->T2q[4] = NULL;
    	csptr->Larmor_frequency = 0;
        if ( ver ) {
          printf( "Chemical shift on nucleus %d (%d%s)\n",n1,ss->iso[n1]->number,ss->iso[n1]->name);
          printf( "  isotropic shift        : %g Hz\n",csptr->iso/2.0/M_PI);
          printf( "  anisotropic shift      : %g Hz\n",csptr->delta/2.0/M_PI);
          printf( "  asymmetry parameter    : %g\n",csptr->eta);
          printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",csptr->pas[0],csptr->pas[1],csptr->pas[2]);
        }
    	continue;
    }
    if (!strncmp(names[i],"jcoupling",9)) {
    	if (NN != 8) {
    		fprintf(stderr,"Error: reading jcoupling - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - jcoupling - int conversion failure\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[1],&n2) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - jcoupling - int conversion failure\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins || n2 < 1 || n2 > ss->nspins) {
    		fprintf(stderr,"Error: spinsys: jcoupling %d %d out of defined nuclei range\n",n1,n2);
    		exit(1);
    	}
    	/* for electron, throw an error */
    	if ( (ss->iso[n1]->number == 0) || (ss->iso[n2]->number ==0) ) {
    		fprintf(stderr,"Error: jcoupling %d %d contains electron spin\n",n1,n2);
    		exit(1);
    	}
    	if (jcoupling_exist(s,n1,n2) >= 0) {
    		fprintf(stderr,"Error: spinsys: jcoupling already exists for nuclei (%d, %d)\n",n1, n2);
    		exit(1);
    	}
    	for (j=2; j<NN; j++) {
    		if (Tcl_GetDoubleFromObj(interp,vals[j],&(data[j-2])) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - jcoupling - double conversion failure\n");
    			exit(1);
    		}
    	}
    	if (fabs(data[0]) > TINY) Jiso = 1; else Jiso = 0;
    	if (fabs(data[1]) > TINY) Jani = 1; else Jani = 0;
    	if ( Jani + Jiso == 0) continue;
    	s->J[(s->nJ)++] = jptr = (Jcoupling*)malloc(sizeof(Jcoupling));
    	if (n1 < n2) {
    		jptr->nuc[0] = n1; jptr->nuc[1] = n2;
    	} else {
    		jptr->nuc[0] = n2; jptr->nuc[1] = n1;
    	}
    	jptr->iso = data[0]*2.0*M_PI; jptr->delta = data[1]*2.0*M_PI; jptr->eta = data[2];
    	jptr->pas[0] = data[3]; jptr->pas[1] = data[4]; jptr->pas[2] = data[5];
    	jptr->blk_Tiso = NULL; jptr->blk_T = NULL;
    	jptr->T2q[0] = jptr->T2q[1] = jptr->T2q[2] = jptr->T2q[3] = jptr->T2q[4] = NULL;
    	if (ss_issame(ss,n1,n2)) {
    		if (ver) printf( "Homonuclear ");
    		if (Jiso) Hiso_isdiag = 0;
    		if (Jani) {
    			Naniso++;
    			HQ_isdiag = 0;
    		}
    		s->Hint_isdiag = 0;
    	} else {
    		if (ver) printf( "Heteronuclear ");
    		if (Jani) {
    			Naniso++;
    		}
    	}
    	if (Jani) {
    		jptr->Rmol = Dtensor2(jptr->delta*2.0, jptr->eta);
    		jptr->Rmol = wig2roti(jptr->Rmol,jptr->pas[0],jptr->pas[1],jptr->pas[2]);
    	} else {
    		jptr->Rmol = NULL;
    	}
        if (ver) {
          printf( "J-coupling between nucleus %d (%d%s) and %d (%d%s)\n",n1,ss->iso[n1]->number,ss->iso[n1]->name,n2,ss->iso[n2]->number,ss->iso[n2]->name);
          printf( "  isotropic value        : %g Hz\n",jptr->iso/2.0/M_PI);
          printf( "  anisotropic value      : %g Hz\n",jptr->delta/2.0/M_PI);
          printf( "  asymmetry              : %g \n",jptr->eta);
          printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",jptr->pas[0],jptr->pas[1],jptr->pas[2]);
        }
    	continue;
    }
    if (!strncmp(names[i],"dipole",6)) {
    	if (NN != 6 && NN != 7) {
    		fprintf(stderr,"Error: reading dipole - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - dipole - int conversion failure\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[1],&n2) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - dipole - int conversion failure\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins || n2 < 1 || n2 > ss->nspins) {
    		fprintf(stderr,"Error: spinsys: dipole %d %d out of defined nuclei range\n",n1,n2);
    		exit(1);
    	}
    	/* for electron, throw an error */
    	if ( (ss->iso[n1]->number == 0) || (ss->iso[n2]->number == 0) ) {
    		fprintf(stderr,"Error: dipole %d %d contains electron spin\n",n1,n2);
    		exit(1);
    	}
    	if (dipole_exist(s,n1,n2) >=0 ) {
    		fprintf(stderr,"Error: spinsys: dipole already exists for nuclei (%d, %d)\n",n1, n2);
    		exit(1);
    	}
    	for (j=2; j<NN; j++) {
    		if (Tcl_GetDoubleFromObj(interp,vals[j],&(data[j-2])) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - dipole - double conversion failure\n");
    			exit(1);
    		}
    	}
    	if (fabs(data[0]) < TINY) continue;
    	s->DD[(s->nDD)++] = ddptr = (Dipole*)malloc(sizeof(Dipole));
    	if (n1 < n2) {
    		ddptr->nuc[0] = n1; ddptr->nuc[1] = n2;
    	} else {
    		ddptr->nuc[0] = n2; ddptr->nuc[1] = n1;
    	}
    	ddptr->delta = data[0];
    	if (NN == 6) {
    		ddptr->eta = 0.0; ddptr->pas[0] = data[1]; ddptr->pas[1] = data[2]; ddptr->pas[2] = data[3];
    	} else {
    		assert(NN==7);
    		ddptr->eta = data[1]; ddptr->pas[0] = data[2]; ddptr->pas[1] = data[3]; ddptr->pas[2] = data[4];
    	}
    	if (ss_issame(ss,n1,n2)) {
    		if (ver) printf( "Homonuclear ");
    		HQ_isdiag = 0;
    		s->Hint_isdiag = 0;
    	} else {
    		if (ver) printf( "Heteronuclear ");
    	}
    	ddptr->blk_T = NULL;
    	ddptr->T2q[0] = ddptr->T2q[1] = ddptr->T2q[2] = ddptr->T2q[3] = ddptr->T2q[4] = NULL;
    	Naniso++;
        if (ver) {
          printf( "dipolar coupling between nucleus %d (%d%s) and %d (%d%s)\n",n1,ss->iso[n1]->number,ss->iso[n1]->name,n2,ss->iso[n2]->number,ss->iso[n2]->name);
          printf( "  dipolar coupling       : %g Hz\n",ddptr->delta);
          if (NN==7) printf("  effective asymmetry    : %g Hz\n",ddptr->eta);
          printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",ddptr->pas[0],ddptr->pas[1],ddptr->pas[2]);
        }
        if (s->dipole_check) {
          int bsign = ( -ss_gamma(ss,n1)*ss_gamma(ss,n2) > 0 ? 1 : -1);
          if (ddptr->delta > 0.0 && bsign < 0) {
            fprintf(stderr,"error: the dipolar coupling between nucleus %d and %d must be negative to comply\n"
                           "       to the conventions used in this program. Set 'dipole_check' to 'false' to\n"
                           "       override this sign check.\n",n1,n2);
            exit(1);
          } else if (ddptr->delta < 0.0 && bsign > 0) {
            fprintf(stderr,"error: the dipolar coupling between nucleus %d and %d must be positive to comply\n"
                           "       to the conventions used in this program. Set 'dipole_check' to 'false' to\n"
                           "       override this sign check.\n",n1,n2);
            exit(1);
          }
        }
    	ddptr->Rmol = Dtensor2(4.0*M_PI*ddptr->delta, ddptr->eta);
    	ddptr->Rmol = wig2roti(ddptr->Rmol, ddptr->pas[0], ddptr->pas[1], ddptr->pas[2]);
    	continue;
    }
    if (!strncmp(names[i],"quadrupole",10)) {
    	if (NN != 7) {
    		fprintf(stderr,"Error: reading quadrupole - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - quadrupole - int conversion failure\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins) {
    		fprintf(stderr,"Error: spinsys: quadrupole %d  out of defined nuclei range\n",n1);
    		exit(1);
    	}
    	/* for electron, throw an error */
    	if (ss->iso[n1]->number == 0) {
    		fprintf(stderr,"Error: quadrupole can not be defined for electron spin %d\n",n1);
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[1],&n2) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - quadrupole - int conversion failure\n");
    		exit(1);
    	}
    	if (fabs(ss_qn(ss,n1)) < 0.99) {
    		fprintf(stderr,"Error: spinsys: quadrupole defined for nuclear spin %d/2\n",(int)(ss_qn(ss,n1)*2));
    		exit(1);
		}

    	if (quadrupole_exist(s,n1) >= 0) {
    		fprintf(stderr,"Error: spinsys: quadrupole already exists for nucleus %d\n",n1);
    		exit(1);
    	}
    	for (j=2; j<NN; j++) {
    		if (Tcl_GetDoubleFromObj(interp,vals[j],&(data[j-2])) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - quadrupole - double conversion failure\n");
    			exit(1);
    		}
    	}
    	if (fabs(data[0]) < TINY) continue;
    	s->Q[(s->nQ)++] = qptr = (Quadrupole*)malloc(sizeof(Quadrupole));
    	qptr->nuc = n1; qptr->order = n2; qptr->pas[0] = data[2]; qptr->pas[1] = data[3]; qptr->pas[2] = data[4];
    	qptr->eta = data[1];
    	qptr->wq = data[0]*2*M_PI/(4.0*ss_qn(ss,n1)*(2.0*ss_qn(ss,n1)-1));
    	//DEBUGPRINT("%f; %f\n",data[0],ss_qn(ss,n1));
    	qptr->w0 = ss_gamma(ss,n1)*s->specfreq/ss_gamma1H()*2*M_PI;
    	//DEBUGPRINT("%f * %f / %f * 2pi/n",ss_gamma(ss,n1),s->specfreq,ss_gamma1H());
    	qptr->Rmol = Dtensor2(2.0*qptr->wq,qptr->eta);
    	qptr->Rmol = wig2roti(qptr->Rmol,qptr->pas[0],qptr->pas[1],qptr->pas[2]);
    	Naniso++;
    	qptr->T = qptr->Ta = qptr->Tb = qptr->T3a = qptr->T3b = qptr->T3c = NULL;
    	qptr->T2q[0] = qptr->T2q[1] = qptr->T2q[2] = qptr->T2q[3] = qptr->T2q[4] = NULL;
    	if (n2 > 3) {
    		fprintf(stderr,"Error: order '%d' of quadrupolar interaction not supported, just up to the third.\n",n2);
    		exit(1);
    	}
        if (ver) {
          printf( "Quadrupolar coupling on nucleus %d (%d%s) up to order %d\n",n1,ss->iso[n1]->number,ss->iso[n1]->name,n2);
          printf( "  quadrupolar constant         :  %g Hz\n",data[0]);
          printf( "  quadrupolar assymetry        :  %g\n",qptr->eta);
          printf( "  euler angles of tensor       :  (%g,%g,%g) degrees\n",qptr->pas[0],qptr->pas[1],qptr->pas[2]);
          printf("  perturbation treatment validity report: %3.2f << 1 is%s fulfilled\n",fabs(qptr->wq/qptr->w0), fabs(10*qptr->wq/qptr->w0) < 1 ? "" : " not");
        }
    	continue;
    }
    if (!strncmp(names[i],"mixing",6)) {
        if (s->frame != ROTFRAME) {
        	fprintf(stderr,"Error: mixing cannot be used with LABframe or DNPframe\n");
        	exit(1);
        }

    	if (NN == 1) {
    		/* syntax "mixing N" */
//    		if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
//    		    fprintf(stderr,"Error: readsys - mixing - Q/CSA int conversion failure\n");
//    		   exit(1);
//    		}
//    		s->MIX[(s->nMIX)++] = mptr = (Mixing*)malloc(sizeof(Mixing));
//    		mptr->type = 1;
//    		mptr->couple[0] = n1;
//    		mptr->couple[1] = n1;
//    		mptr->qidx = mptr->idx = -1;
//    		mptr->T = NULL;
//    		mptr->Ta = mptr->Tb = NULL;
    		fprintf(stderr,"Error: reading mixing - syntax mixing N not supported\n");
    		exit(1);
    	} else if (NN == 2) {
    		/* syntax "mixing Nq Nd" */
    		if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		    fprintf(stderr,"Error: readsys - mixing - int conversion failure\n");
    		    exit(1);
    		}
    		if (Tcl_GetIntFromObj(interp,vals[1],&n2) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - mixing - int conversion failure\n");
    			exit(1);
    		}
    		s->MIX[(s->nMIX)++] = mptr = (Mixing*)malloc(sizeof(Mixing));
    		if (n1 == n2) {
    			mptr->type = 1;  // Q-CSA
    		} else {
    			mptr->type = 2;  // Q-DD
    		}
    		mptr->couple[0] = n1;
    		mptr->couple[1] = n2;
    		mptr->qidx = mptr->idx = -1;
    		mptr->T = NULL;
    		mptr->Ta = mptr->Tb = NULL;
            if (ver) {
            	if (mptr->type == 1) {
            		printf( "Second order mixing quadrupole-CSA\n" );
            		printf( "   quadrupole and shift nucleus     :  %d \n",n1);
            	} else {
            		printf( "Second order mixing quadrupole-dipole\n" );
            		printf( "   quadrupole nucleus     :  %d \n",n1);
            		printf( "   dipole nuclei          :  %d and %d\n",n1,n2);
            	}
            }
    	} else {
    		fprintf(stderr,"Error: reading mixing - parameter count mismatch\n");
    		exit(1);
    	}

        continue;
    }
    if (!strncmp(names[i],"gtensor",7)) {
    	if (NN != 7) {
    		fprintf(stderr,"Error: reading gtensor - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - gtensor - int conversion failure\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins) {
    		fprintf(stderr,"Error: spinsys: gtensor %d  out of defined range\n",n1);
    		exit(1);
    	}
    	if (gtensor_exist(s,n1) >= 0) {
    		fprintf(stderr,"Error: spinsys: gtensor already exists for electron %d\n",n1);
    		exit(1);
    	}
    	// test if it is electron
    	if (ss->iso[n1]->number != 0) {
    		fprintf(stderr,"Error: spinsys: gtensor : spin %d is not electron\n",n1);
    		exit(1);
    	}
    	for (j=1; j<NN; j++) {
    		if (Tcl_GetDoubleFromObj(interp,vals[j],&(data[j-1])) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - gtensor - double conversion failure\n");
    			exit(1);
    		}
    	}
    	/* when zero amplitudes, do nothing */
    	if ( (fabs(data[0]) < TINY) && (fabs(data[1]) < TINY) ) continue;
    	s->G[(s->nG)++] = gptr = (Gtensor*)malloc(sizeof(Gtensor));
    	gptr->nuc = n1; gptr->eta = data[2]; gptr->pas[0] = data[3]; gptr->pas[1] = data[4]; gptr->pas[2] = data[5];
   		gptr->iso = -data[0]*2.0*M_PI; gptr->delta = -data[1]*2.0*M_PI;
    	if (fabs(gptr->delta) > TINY) {
    		Naniso++;
    		gptr->Rmol = Dtensor2(SQRT2BY3*gptr->delta, gptr->eta);
    		gptr->Rmol = wig2roti(gptr->Rmol, gptr->pas[0], gptr->pas[1], gptr->pas[2]);
    	} else {
    		gptr->Rmol = NULL;
    	}
    	gptr->T = NULL;
    	gptr->T2q[0] = gptr->T2q[1] = gptr->T2q[2] = gptr->T2q[3] = gptr->T2q[4] = NULL;
        if ( ver ) {
          printf( "G-tensor for electron %d\n",n1);
          printf( "  isotropic value        : %g Hz\n",gptr->iso/2.0/M_PI);
          printf( "  anisotropic value      : %g Hz\n",gptr->delta/2.0/M_PI);
          printf( "  asymmetry parameter    : %g\n",gptr->eta);
          printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",gptr->pas[0],gptr->pas[1],gptr->pas[2]);
        }
    	continue;
    }
    if (!strncmp(names[i],"hyperfine",9)) {
    	if (s->frame == ROTFRAME) {
    		fprintf(stderr,"Error: hyperfine cannot be treated using method ROTframe\n");
    		exit(1);
    	}
    	if (NN != 7 && NN != 8) { // this allows eta which comes in only with _averaged tensor like for dipole
    		fprintf(stderr,"Error: reading hyperfine - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - hyperfine - int conversion failure\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[1],&n2) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - hyperfine - int conversion failure\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins || n2 < 1 || n2 > ss->nspins) {
    		fprintf(stderr,"Error: spinsys: hyperfine %d %d out of defined elements range\n",n1,n2);
    		exit(1);
    	}
    	if (hyperfine_exist(s,n1,n2) >= 0) {
    		fprintf(stderr,"Error: readsys: hyperfine already exists for elements (%d, %d)\n",n1, n2);
    		exit(1);
    	}
    	if ( (ss->iso[n1]->number != 0) && (ss->iso[n2]->number != 0) ) {
    		fprintf(stderr,"Error: readsys: hyperfine %d %d does not involve any electron spin\n", n1, n2);
    		exit(1);
    	}
    	if ( (ss->iso[n1]->number == 0) && (ss->iso[n2]->number == 0) ) {
    		fprintf(stderr,"Error: readsys: hyperfine %d %d - both spins are electrons\n", n1, n2);
    		exit(1);
    	}
    	for (j=2; j<NN; j++) {
    		if (Tcl_GetDoubleFromObj(interp,vals[j],&(data[j-2])) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - hyperfine - double conversion failure\n");
    			exit(1);
    		}
    	}
    	if (fabs(data[0]) > TINY) Jiso = 1; else Jiso = 0;
    	if (fabs(data[1]) > TINY) {
    		Jani = 1;
    		Naniso++;
    		HQ_isdiag = 0;
    		s->Hint_isdiag = 0;
    	}	else {
    		Jani = 0;
    	}
    	if ( Jani + Jiso == 0) continue;
    	s->HF[(s->nHF)++] = hfptr = (Hyperfine*)malloc(sizeof(Hyperfine));
    	if (ss->iso[n1]->number == 0) {  // n1 is electron
    		hfptr->nuc[0] = n1; hfptr->nuc[1] = n2;
    	} else {                         // n2 is electron
    		hfptr->nuc[0] = n2; hfptr->nuc[1] = n1;
    	}
    	hfptr->iso = data[0]*2.0*M_PI; hfptr->delta = data[1]*2.0*M_PI;
    	if (NN == 7) {
    		hfptr->eta = 0;
    	} else {
    		assert(NN==8);
    		hfptr->eta = data[2];
    	}
    	hfptr->pas[0] = data[NN-5]; hfptr->pas[1] = data[NN-4]; hfptr->pas[2] = data[NN-3];
    	hfptr->blk_Tiso = NULL; hfptr->blk_T = NULL; hfptr->blk_Ta = NULL; hfptr->blk_Tb = NULL;
    	hfptr->T2q[0] = hfptr->T2q[1] = hfptr->T2q[2] = hfptr->T2q[3] = hfptr->T2q[4] = NULL;
    	if (Jani) {
    		hfptr->Rmol = Dtensor2(hfptr->delta*2.0, hfptr->eta);
    		hfptr->Rmol = wig2roti(hfptr->Rmol,hfptr->pas[0],hfptr->pas[1],hfptr->pas[2]);
    	} else {
    		hfptr->Rmol = NULL;
    	}
        if (ver) {
          printf( "Hyperfine coupling between electron %d and nucleus %d (%d%s)\n",hfptr->nuc[0],hfptr->nuc[1],ss->iso[hfptr->nuc[1]]->number,ss->iso[hfptr->nuc[1]]->name);
          printf( "  isotropic value (Fermi): %g Hz\n",hfptr->iso/2.0/M_PI);
          printf( "  anisotropic value      : %g Hz\n",hfptr->delta/2.0/M_PI);
          printf( "  asymmetry              : %g \n",hfptr->eta);
          printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",hfptr->pas[0],hfptr->pas[1],hfptr->pas[2]);
        }
    	continue;
    }
    if (!strncmp(names[i],"heisenberg",10)) {
    	if (NN != 3) {
    		fprintf(stderr,"Error: reading heisenberg - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - heisenberg - int conversion failure (1)\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[1],&n2) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - heisenberg - int conversion failure (2)\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins) {
    		fprintf(stderr,"Error: readsys: heisenberg %d %d - first index out of defined range of spins\n",n1,n2);
    		exit(1);
    	}
    	if (n2 < 1 || n2 > ss->nspins) {
    		fprintf(stderr,"Error: readsys: heisenberg %d %d - second index out of defined range of spins\n",n1,n2);
    		exit(1);
    	}
    	/* if not electrons, throw an error */
    	if ( (ss->iso[n1]->number != 0) || (ss->iso[n2]->number !=0) ) {
    		fprintf(stderr,"Error: heisenberg %d %d can be defined for electron spins only\n",n1,n2);
    		exit(1);
    	}
    	if (heisenberg_exist(s,n1,n2) >= 0) {
    		fprintf(stderr,"Error: readsys: heisenberg %d %d  already exists\n",n1,n2);
    		exit(1);
    	}
   		if (Tcl_GetDoubleFromObj(interp,vals[2],&(data[0])) != TCL_OK) {
   			fprintf(stderr,"Error: readsys - heisenberg %d %d - double conversion failure\n",n1,n2);
   			exit(1);
    	}
    	if (fabs(data[0]) < TINY) continue;
    	s->HEx[(s->nHEX)++] = hexptr = (Heisenberg_exchange*)malloc(sizeof(Heisenberg_exchange));
    	if (n1 < n2) {
    		hexptr->electron[0] = n1; hexptr->electron[1] = n2;
    	} else {
    		hexptr->electron[0] = n2; hexptr->electron[1] = n1;
    	}
    	hexptr->iso = data[0]*2.0*M_PI;
    	hexptr->blk_Tiso = NULL;
    	// we are sure that both n1 and n2 are electrons and this interaction is isotropic
    	Hiso_isdiag = 0;
   		s->Hint_isdiag = 0;
        if (ver) {
          printf( "Heisenberg exchange interaction between electrons %d and %d\n",n1,n2);
          printf( "  isotropic value        : %g Hz\n",hexptr->iso/2.0/M_PI);
        }
    	continue;
    }
    if (!strncmp(names[i],"edipole",7)) {
    	if (NN != 6 && NN != 7) {
    		fprintf(stderr,"Error: reading edipole - parameter count mismatch\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[0],&n1) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - edipole - int conversion failure 1\n");
    		exit(1);
    	}
    	if (Tcl_GetIntFromObj(interp,vals[1],&n2) != TCL_OK) {
    		fprintf(stderr,"Error: readsys - edipole - int conversion failure 2\n");
    		exit(1);
    	}
    	if (n1 < 1 || n1 > ss->nspins || n2 < 1 || n2 > ss->nspins) {
    		fprintf(stderr,"Error: readsys: edipole %d %d out of defined range of spins\n",n1,n2);
    		exit(1);
    	}
    	/* if not electrons, throw an error */
    	if ( (ss->iso[n1]->number != 0) || (ss->iso[n2]->number != 0) ) {
    		fprintf(stderr,"Error: edipole %d %d can be defined for electron spins only\n",n1,n2);
    		exit(1);
    	}
    	if (edipole_exist(s,n1,n2) >=0 ) {
    		fprintf(stderr,"Error: readsys: edipole %d %d already exists \n",n1, n2);
    		exit(1);
    	}
    	for (j=2; j<NN; j++) {
    		if (Tcl_GetDoubleFromObj(interp,vals[j],&(data[j-2])) != TCL_OK) {
    			fprintf(stderr,"Error: readsys - edipole - double conversion failure %d\n",j);
    			exit(1);
    		}
    	}
    	if (fabs(data[0]) < TINY) continue;
    	s->EDD[(s->nEDD)++] = eddptr = (Electron_Dipole*)malloc(sizeof(Electron_Dipole));
    	if (n1 < n2) {
    		eddptr->electron[0] = n1; eddptr->electron[1] = n2;
    	} else {
    		eddptr->electron[0] = n2; eddptr->electron[1] = n1;
    	}
    	eddptr->delta = data[0];
    	if (NN == 6) {
    		eddptr->eta = 0.0; eddptr->pas[0] = data[1]; eddptr->pas[1] = data[2]; eddptr->pas[2] = data[3];
    	} else {
    		assert(NN==7);
    		eddptr->eta = data[1]; eddptr->pas[0] = data[2]; eddptr->pas[1] = data[3]; eddptr->pas[2] = data[4];
    	}
    	// we are sure both n1 and n2 are electrons
   		HQ_isdiag = 0;
   		s->Hint_isdiag = 0;
    	eddptr->blk_T = NULL;
    	eddptr->T2q[0] = eddptr->T2q[1] = eddptr->T2q[2] = eddptr->T2q[3] = eddptr->T2q[4] = NULL;
    	Naniso++;
        if (ver) {
          printf( "Dipolar coupling between electrons %d and %d \n",n1,n2);
          printf( "  dipolar coupling       : %g Hz\n",eddptr->delta);
          if (NN==7) printf("  effective asymmetry    : %g Hz\n",eddptr->eta);
          printf( "  euler angles of tensor :  (%g,%g,%g) degrees\n",eddptr->pas[0],eddptr->pas[1],eddptr->pas[2]);
        }
        if (s->dipole_check) {
          if (eddptr->delta < 0.0 ) {
            fprintf(stderr,"error: the dipolar coupling between electrons %d and %d must be positive to comply\n"
                           "       to the conventions used in this program. Set 'dipole_check' to 'false' to\n"
                           "       override this sign check.\n",n1,n2);
            exit(1);
          }
        }
    	eddptr->Rmol = Dtensor2(4.0*M_PI*eddptr->delta, eddptr->eta);
    	eddptr->Rmol = wig2roti(eddptr->Rmol, eddptr->pas[0], eddptr->pas[1], eddptr->pas[2]);
    	continue;
    }
    fprintf(stderr,"Error: unknown identifier '%s' in spinsys, must be one of\n",names[i]);
    fprintf(stderr,"       channels, nuclei, shift, dipole, jcoupling, quadrupole,\n");
    fprintf(stderr,"       gtensor, hyperfine, heisenberg, edipole or mixing.\n");
    exit(1);
  }
  if ( (s->frame == LABFRAME) || (s->frame == DNPFRAME) ) { // add Zeeman terms to nuclear shifts
	  for (i=1; i<=ss->nspins; i++) {
		  if (ss->iso[i]->number == 0) continue; // no Zeeman term for electrons, shift is defined for nuclei only
		  j = shift_exist(s,i);
		  if (j >= 0) { // add to existing shift structure
			  csptr = s->CS[j];
		  } else { // create additional shift structure
			  csptr = s->CS[(s->nCS)++] = (Shift*)malloc(sizeof(Shift));
			  assert(s->nCS <= Nshift);
		      csptr->nuc = i; csptr->iso = 0; csptr->delta = 0; csptr->eta = 0; csptr->pas[0] = 0; csptr->pas[1] = 0; csptr->pas[2] = 0;
		      csptr->Rmol = NULL; csptr->T = NULL; csptr->T2q[0] = csptr->T2q[1] = csptr->T2q[2] = csptr->T2q[3] = csptr->T2q[4] = NULL;
		  }
		  csptr->Larmor_frequency = ss_larmor_frequency(s, csptr->nuc);
		  if (ver) {
			  printf("Nucleus %d (%d%s): Larmor frequency %.2f Hz\n",csptr->nuc,ss->iso[csptr->nuc]->number,ss->iso[csptr->nuc]->name,(csptr->Larmor_frequency)/2/M_PI);
		  }
	  }
  }

  // in LABframe and DNPframe we do not use H assembly and block_diag
  if ( (s->frame == LABFRAME) || (s->frame == DNPFRAME) ) {
	  s->Hassembly = 0;
	  s->block_diag = 0;
	  s->Hint_isdiag = 0;
	  Naniso = -Naniso;
  }
  if (s->block_diag && !(s->Hint_isdiag)) {
	  Nblk = LEN(blk_dims);
	  s->basis = (1 << s->Nmz) - 1;
  } else {
	  Nblk = 1;
	  s->basis = 0;
  }
  if ( (Naniso > 4) || ( (Naniso != 0) && (s->block_diag) ) ) {
	  s->Hassembly = 1;
  } else {
	  s->Hassembly = 0;
  }

  /****************************
  printf("Hamiltonian block structure will");
  if (s->block_diag != 0) {
	  printf(" be based on %d nuclei types\n",Nmz);
	  printf("\t number of blocks: %d\n",LEN(blk_dims));
	  printf("\t dims:");
	  for (i=1; i<=LEN(blk_dims); i++) printf(" %d",blk_dims[i]);
	  printf("\n");
  } else {
	  printf(" not be used\n");
  }
  exit(1);
  ****************************/


  /* create common Hamiltonians  */
	  int dum_type;
	  if (s->sparse > 0) dum_type = MAT_SPARSE; else dum_type = MAT_DENSE;
	  s->Hiso = create_blk_mat_double(s->matdim,Nblk,blk_dims,Hiso_isdiag ? MAT_DENSE_DIAG : dum_type, s->basis);
	  blk_dm_zero(s->Hiso);
	  // initialize matrices for H assembly in rotating frame
	  if ( s->Hassembly == 1 ) {
		  for (i=0; i<5; i++) {
			  s->HQ[i] = create_blk_mat_double(s->matdim,Nblk,blk_dims,HQ_isdiag ? MAT_DENSE_DIAG : dum_type, s->basis);
			  blk_dm_zero(s->HQ[i]);
		  }
	  } else {
		  for (i=0; i<5; i++) s->HQ[i] = NULL;
	  }
//	  if (s->frame == LABFRAME) { // labframe simulation: add Zeeman terms for all nuclei
//		  for (i=1; i<=s->ss->nspins; i++) {
//			  mat_double *dumIz = Iz_ham(s,i);
//			  blk_dm_multod_diag(s->Hiso,dumIz,ss_gamma(s->ss,i)*s->specfreq/ss_gamma1H()*2*M_PI);
//			  free_double_matrix(dumIz);
//			  //printf("Nucleus %d : Larmor freq. %15.10f\n",i,ss_gamma(s->ss,i)*s->specfreq/ss_gamma1H()*2*M_PI);
//		  }
//	  }
	  /* chemical shift */
	  for (i=0; i<s->nCS; i++) {
		  //printf("creating shift %d\n",i);
		  csptr = s->CS[i];
		  csptr->T = Iz_ham(s,csptr->nuc);
		  if (fabs(csptr->iso) > TINY ) {
			  blk_dm_multod_diag(s->Hiso,csptr->T,csptr->iso);
		  }
		  if (fabs(csptr->Larmor_frequency) > TINY ) {
			  blk_dm_multod_diag(s->Hiso,csptr->T,csptr->Larmor_frequency);
		  }
		  if (fabs(csptr->delta) > TINY) {
			  if (s->Hassembly) {
				  blk_dm_multod_diag(s->HQ[0],csptr->T,csptr->Rmol[3].re);
				  blk_dm_multod_diag(s->HQ[1],csptr->T,csptr->Rmol[4].re);
				  blk_dm_multod_diag(s->HQ[2],csptr->T,csptr->Rmol[4].im);
				  blk_dm_multod_diag(s->HQ[3],csptr->T,csptr->Rmol[5].re);
				  blk_dm_multod_diag(s->HQ[4],csptr->T,csptr->Rmol[5].im);
				  free_double_matrix(csptr->T);
				  csptr->T = NULL;
			  }
			  //printf("CSA of nucleus %d - frame is %d\n",csptr->nuc, s->frame);
			  if ( (s->frame == LABFRAME) || (s->frame == DNPFRAME) ) {
				  // shift is for nuclei only
				  //printf("  adding T2q tensors\n");
				  csptr->T2q[0] = csptr->T2q[4] = NULL; // T2-2 and T22 are zero for CSA
				  csptr->T2q[1] = Im_real(s,csptr->nuc); dm_muld(csptr->T2q[1],0.5);
				  csptr->T2q[2] = csptr->T; csptr->T = NULL; dm_muld(csptr->T2q[2],sqrt(2.0/3.0));
				  csptr->T2q[3] = Ip_real(s,csptr->nuc); dm_muld(csptr->T2q[3],-0.5);
				  // readjust Rmol to agree with T2q scalings
				  cv_muld(csptr->Rmol,sqrt(3.0/2.0));
			  } else {
				  csptr->T2q[0] = csptr->T2q[1] = csptr->T2q[2] = csptr->T2q[3] = csptr->T2q[4] = NULL;
			  }
		  } else {
			  free_double_matrix(csptr->T);
			  csptr->T = NULL;
			  csptr->T2q[0] = csptr->T2q[1] = csptr->T2q[2] = csptr->T2q[3] = csptr->T2q[4] = NULL;
		  }
	  }
	  /* dipole - dipole interactions */
	  for (i=0; i<s->nDD; i++) {
		  ddptr = s->DD[i];
		  n1 = ddptr->nuc[0]; n2 = ddptr->nuc[1];
		  //printf("creating dipole %d %d \n",n1,n2);
		  if (ss_issame(ss,n1,n2)) {
			  ddptr->blk_T = T20(s,n1,n2);
		  } else {
			  ddptr->blk_T = IzIz_sqrt2by3(s,n1,n2);
		  }
		  if (s->Hassembly) {
			  blk_dm_multod(s->HQ[0],ddptr->blk_T,ddptr->Rmol[3].re);
			  blk_dm_multod(s->HQ[1],ddptr->blk_T,ddptr->Rmol[4].re);
			  blk_dm_multod(s->HQ[2],ddptr->blk_T,ddptr->Rmol[4].im);
			  blk_dm_multod(s->HQ[3],ddptr->blk_T,ddptr->Rmol[5].re);
			  blk_dm_multod(s->HQ[4],ddptr->blk_T,ddptr->Rmol[5].im);
			  free_blk_mat_double(ddptr->blk_T);
			  ddptr->blk_T = NULL;
		  }
		  if ( (s->frame == LABFRAME) || (s->frame == DNPFRAME) ) {
			  // dipole is defined for nuclei only
			  mat_double *md1, *md2, *md3, *md4;
			  md1 = Im_real(s,n1); md2 = Im_real(s,n2); md3 = Iz_ham(s,n2);
			  ddptr->T2q[0] = dm_mul(md1,md2); dm_muld(ddptr->T2q[0],0.5);
			  ddptr->T2q[1] = Iz_ham(s,n1); dm_multo(ddptr->T2q[1],md2); dm_mm(0.5,md1,md3,0.5,ddptr->T2q[1]);
			  ddptr->T2q[2] = Iz_ham(s,n1); dm_multo(ddptr->T2q[2],md3);
			  md4 = Ip_real(s,n2); dm_mm(-0.5,md1,md4,2.0,ddptr->T2q[2]);
			  free_double_matrix(md1); md1 = Ip_real(s,n1); dm_mm(-0.5/sqrt(6.0),md1,md2,1.0/sqrt(6.0),ddptr->T2q[2]);
			  ddptr->T2q[3] = Iz_ham(s,n1); dm_multo(ddptr->T2q[3],md4); dm_mm(-0.5,md1,md3,-0.5,ddptr->T2q[3]);
			  ddptr->T2q[4] = dm_mul(md1,md4); dm_muld(ddptr->T2q[4],0.5);
			  free_double_matrix(md1); free_double_matrix(md2); free_double_matrix(md3); free_double_matrix(md4);
//			  dm_print(ddptr->T2q[0],"DD T2-2");
//			  dm_print(ddptr->T2q[1],"DD T2-1");
//			  dm_print(ddptr->T2q[2],"DD T20");
//			  dm_print(ddptr->T2q[3],"DD T21");
//			  dm_print(ddptr->T2q[4],"DD T22");
		  } else {
			  ddptr->T2q[0] = ddptr->T2q[1] = ddptr->T2q[2] = ddptr->T2q[3] = ddptr->T2q[4] = NULL;
		  }
	  }
	  /* J-couplings */
	  double corr_factor = 1.0;
	  for (i=0; i<s->nJ; i++) {
		  jptr = s->J[i];
		  n1 = jptr->nuc[0]; n2 = jptr->nuc[1];
		  if (fabs(jptr->iso)> TINY) {
			  if (ss_issame(ss,n1,n2)) {
				  jptr->blk_Tiso = II(s,n1,n2);
				  corr_factor = 1.0;
			  } else {
				  jptr->blk_Tiso = IzIz_sqrt2by3(s,n1,n2);
				  corr_factor = 1.0/SQRT2BY3;
			  }
			  blk_dm_multod(s->Hiso,(jptr->blk_Tiso),jptr->iso*corr_factor);
			  free_blk_mat_double(jptr->blk_Tiso);
			  jptr->blk_Tiso = NULL;
		  }
		  if (fabs(jptr->delta) > TINY) {
			  if (ss_issame(ss,n1,n2)) {
				  jptr->blk_T = T20(s,n1,n2);
			  } else {
				  jptr->blk_T = IzIz_sqrt2by3(s,n1,n2);
			  }
			  if (s->Hassembly) {
				  blk_dm_multod(s->HQ[0],jptr->blk_T,jptr->Rmol[3].re);
				  blk_dm_multod(s->HQ[1],jptr->blk_T,jptr->Rmol[4].re);
				  blk_dm_multod(s->HQ[2],jptr->blk_T,jptr->Rmol[4].im);
				  blk_dm_multod(s->HQ[3],jptr->blk_T,jptr->Rmol[5].re);
				  blk_dm_multod(s->HQ[4],jptr->blk_T,jptr->Rmol[5].im);
				  free_blk_mat_double(jptr->blk_T);
				  jptr->blk_T = NULL;
			  }
			  if ( (s->frame == LABFRAME) || (s->frame == DNPFRAME) ) {
				  fprintf(stderr,"Error: LABframe/DNPframe simulations do not consider J anisotropy explicitly\n");
				  exit(1);
			  }
		  }
		  jptr->T2q[0] = jptr->T2q[1] = jptr->T2q[2] = jptr->T2q[3] = jptr->T2q[4] = NULL;
	  }
	  /* quadrupolar interactions */
	  for (i=0; i<s->nQ; i++) {
		  qptr = s->Q[i];
		  qptr->T = T20II(s,qptr->nuc);
		  if (s->Hassembly) {
			  blk_dm_multod_diag(s->HQ[0],qptr->T,qptr->Rmol[3].re);
			  blk_dm_multod_diag(s->HQ[1],qptr->T,qptr->Rmol[4].re);
			  blk_dm_multod_diag(s->HQ[2],qptr->T,qptr->Rmol[4].im);
			  blk_dm_multod_diag(s->HQ[3],qptr->T,qptr->Rmol[5].re);
			  blk_dm_multod_diag(s->HQ[4],qptr->T,qptr->Rmol[5].im);
			  free_double_matrix(qptr->T);
			  qptr->T = NULL;
		  }
		  if (qptr->order == 2) {
			  fill_Tquad_2(s,qptr->nuc,qptr);
		  }
		  if (qptr->order == 3) {
			  fill_Tquad_2(s,qptr->nuc,qptr);
			  fill_Tquad_3(s,qptr);
//			  dm_print(qptr->T3a,"T3a");
//			  dm_print(qptr->T3b,"T3b");
//			  dm_print(qptr->T3c,"T3c");
//			  exit(1);
		  }
		  if ( (s->frame == LABFRAME) || (s->frame == DNPFRAME) ) {
			  mat_double *md1, *md2;
			  md1 = Im_real(s,qptr->nuc); md2 = Iz_ham(s,qptr->nuc);
			  qptr->T2q[0] = dm_mul(md1,md1); dm_muld(qptr->T2q[0],0.5);
			  qptr->T2q[1] = dm_mul(md1,md2); dm_mm(0.5,md2,md1,0.5,qptr->T2q[1]);
			  qptr->T2q[2] = qptr->T; qptr->T = NULL;
			  free_double_matrix(md1); md1 = Ip_real(s,qptr->nuc);
			  qptr->T2q[3] = dm_mul(md1,md2); dm_mm(-0.5,md2,md1,-0.5,qptr->T2q[3]);
			  qptr->T2q[4] = dm_mul(md1,md1); dm_muld(qptr->T2q[4],0.5);
			  free_double_matrix(md1); free_double_matrix(md2);
//			  dm_print(qptr->T2q[0],"Q T2-2");
//			  dm_print(qptr->T2q[1],"Q T2-1");
//			  dm_print(qptr->T2q[2],"Q T20");
//			  dm_print(qptr->T2q[3],"Q T21");
//			  dm_print(qptr->T2q[4],"Q T22");
		  } else {
			  qptr->T2q[0] = qptr->T2q[1] = qptr->T2q[2] = qptr->T2q[3] = qptr->T2q[4] = NULL;
		  }
	  }
	  /* mixing terms */
	  for (i=0; i<s->nMIX; i++) {
		  mptr = s->MIX[i];
		  if (mptr->type == 1) {
			  /*  Q-CSA */
			  int qidx = quadrupole_exist(s,mptr->couple[0]);
			  if ( qidx < 0) {
				  fprintf(stderr,"Error: spinsys: mixing - nucleus %d has no defined quadrupole\n",mptr->couple[0]);
				  exit(1);
			  }
			  int idx = shift_exist(s,mptr->couple[1]);
			  if ( idx < 0) {
				  fprintf(stderr,"Error: spinsys: mixing - no shift interaction defined for nucleus %d\n",mptr->couple[1]);
				  exit(1);
			  }
			  mptr->qidx = qidx;
			  mptr->idx = idx;
			  mptr->T = T20II(s,mptr->couple[0]);
			  dm_muld(mptr->T,sqrt(6.0));
			  //fprintf(stderr,"Error: readsys - creation of Hams - mixing type 1 not implemented\n");
			  //exit(1);
		  } else {
			  /*  Q-DD  */
			  int qidx = quadrupole_exist(s,mptr->couple[0]);
			  if ( qidx < 0) {
				  fprintf(stderr,"Error: spinsys: mixing - nucleus %d has no defined quadrupole\n",mptr->couple[0]);
				  exit(1);
			  }
			  int idx = dipole_exist(s,mptr->couple[0],mptr->couple[1]);
			  if ( idx < 0) {
				  fprintf(stderr,"Error: spinsys: mixing - no dipole-dipole interaction defined between nuclei %d and %d\n",mptr->couple[0],mptr->couple[1]);
				  exit(1);
			  }
			  mptr->qidx = qidx;
			  mptr->idx = idx;
			  assert(mptr->type == 2);
			  fill_Tmix_dipole(s,mptr);
			  if (ss_issame(ss,mptr->couple[0],mptr->couple[1])) {
				  fill_Tabmix_dipole(s,mptr);
			  }
		  }
	  }
	  // G-tensor for electrons
	  for (i=0; i<s->nG; i++) {
		  //printf("creating G-tensor %d\n",i);
		  gptr = s->G[i];
		  gptr->T = Iz_ham(s,gptr->nuc);
		  //dm_print(gptr->T,"operator T");
		  if (fabs(gptr->iso) > TINY ) {
			  blk_dm_multod_diag(s->Hiso,gptr->T,gptr->iso);
		  }
		  if (fabs(gptr->delta) > TINY) {
			  if (s->Hassembly) {
				  blk_dm_multod_diag(s->HQ[0],gptr->T,gptr->Rmol[3].re);
				  blk_dm_multod_diag(s->HQ[1],gptr->T,gptr->Rmol[4].re);
				  blk_dm_multod_diag(s->HQ[2],gptr->T,gptr->Rmol[4].im);
				  blk_dm_multod_diag(s->HQ[3],gptr->T,gptr->Rmol[5].re);
				  blk_dm_multod_diag(s->HQ[4],gptr->T,gptr->Rmol[5].im);
				  free_double_matrix(gptr->T);
				  gptr->T = NULL;
			  }
			  if (s->frame == LABFRAME) {
				  fprintf(stderr,"Error: gtensor cannot be combined with method LABframe\n");
				  exit(1);
				  /*****
				  gptr->T2q[0] = gptr->T2q[4] = NULL; // T2-2 and T22 are zero for CSA
				  gptr->T2q[1] = Im_real(s,gptr->nuc); dm_muld(gptr->T2q[1],0.5);
				  gptr->T2q[2] = gptr->T; gptr->T = NULL; dm_muld(gptr->T2q[2],sqrt(2.0/3.0));
				  gptr->T2q[3] = Ip_real(s,gptr->nuc); dm_muld(gptr->T2q[3],-0.5);
				  // readjust Rmol to agree with T2q scalings
				  cv_muld(gptr->Rmol,sqrt(3.0/2.0));
				  *****/
			  }
			  // note T2q[] are NULLed above in creating Rmol
		  } else {
			  free_double_matrix(gptr->T);
			  gptr->T = NULL;
			  //gptr->T2q[0] = gptr->T2q[1] = gptr->T2q[2] = gptr->T2q[3] = gptr->T2q[4] = NULL;
		  }
	  }
	  // Hyperfine couplings for electrons
	  for (i=0; i<s->nHF; i++) {
		  //printf("creating Hyperfine coupling %d\n",i);
		  hfptr = s->HF[i];
		  n1 = hfptr->nuc[0]; // this is electron
		  n2 = hfptr->nuc[1];
		  if (fabs(hfptr->iso)> TINY) {
			  hfptr->blk_Tiso = IzIz_sqrt2by3(s,n1,n2);
			  blk_dm_print(hfptr->blk_Tiso,"oparator Tiso");
			  blk_dm_multod(s->Hiso,hfptr->blk_Tiso,hfptr->iso/SQRT2BY3);
			  free_blk_mat_double(hfptr->blk_Tiso);
			  hfptr->blk_Tiso = NULL;
		  }
		  if (fabs(hfptr->delta) > TINY) {
			  hfptr->blk_T = IzIz_sqrt2by3(s,n1,n2);
			  //blk_dm_print(hfptr->blk_T,"oparator T (=T20)");
			  if ( (s->Hassembly) && (s->frame == ROTFRAME) ) {
				  // NOTE 30.5.2021 - Hassembly is disabled in DNPframe where hyperfine is used
				  //     but we can use this to include the secular T20 and ignore pseudosecular T2+-1 in ROTFRAME
				  blk_dm_multod(s->HQ[0],hfptr->blk_T,hfptr->Rmol[3].re);
				  blk_dm_multod(s->HQ[1],hfptr->blk_T,hfptr->Rmol[4].re);
				  blk_dm_multod(s->HQ[2],hfptr->blk_T,hfptr->Rmol[4].im);
				  blk_dm_multod(s->HQ[3],hfptr->blk_T,hfptr->Rmol[5].re);
				  blk_dm_multod(s->HQ[4],hfptr->blk_T,hfptr->Rmol[5].im);
				  free_blk_mat_double(hfptr->blk_T);
				  hfptr->blk_T = NULL;
			  }
			  if (s->frame == DNPFRAME) fill_hyperfine_Tab(s,hfptr);
			  //blk_dm_print(hfptr->blk_Ta,"oparator Ta (=T2-1)");
			  //blk_dm_print(hfptr->blk_Tb,"oparator Tb (=T2+1)");
			  if (s->frame == LABFRAME) {
				  fprintf(stderr,"Error: LABframe simulations do not support hyperfine interaction\n");
				  exit(1);
			  }
		  }
		  hfptr->T2q[0] = hfptr->T2q[1] = hfptr->T2q[2] = hfptr->T2q[3] = hfptr->T2q[4] = NULL;
	  }
	  // Heisenberg exchange for electrons
	  for (i=0; i<s->nHEX; i++) {
		  hexptr = s->HEx[i];
		  n1 = hexptr->electron[0]; n2 = hexptr->electron[1];
		  if (fabs(hexptr->iso)> TINY) {
			  // we know it is homo-spin
			  hexptr->blk_Tiso = II(s,n1,n2);
			  blk_dm_multod(s->Hiso,hexptr->blk_Tiso,hexptr->iso);
			  free_blk_mat_double(hexptr->blk_Tiso);
			  hexptr->blk_Tiso = NULL;
			  // jak je to s turnoff???
		  }
	  }
	  // Dipole-dipole for electrons
	  for (i=0; i<s->nEDD; i++) {
		  eddptr = s->EDD[i];
		  n1 = eddptr->electron[0]; n2 = eddptr->electron[1];
		  //printf("creating dipole %d %d \n",n1,n2);
		  eddptr->blk_T = T20(s,n1,n2);
		  if (s->Hassembly) {
			  blk_dm_multod(s->HQ[0],eddptr->blk_T,eddptr->Rmol[3].re);
			  blk_dm_multod(s->HQ[1],eddptr->blk_T,eddptr->Rmol[4].re);
			  blk_dm_multod(s->HQ[2],eddptr->blk_T,eddptr->Rmol[4].im);
			  blk_dm_multod(s->HQ[3],eddptr->blk_T,eddptr->Rmol[5].re);
			  blk_dm_multod(s->HQ[4],eddptr->blk_T,eddptr->Rmol[5].im);
			  free_blk_mat_double(eddptr->blk_T);
			  eddptr->blk_T = NULL;
		  }
		  if ( s->frame == LABFRAME ) {
			  // edipole is defined for electrons only
			  mat_double *md1, *md2, *md3, *md4;
			  md1 = Im_real(s,n1); md2 = Im_real(s,n2); md3 = Iz_ham(s,n2);
			  eddptr->T2q[0] = dm_mul(md1,md2); dm_muld(eddptr->T2q[0],0.5);
			  eddptr->T2q[1] = Iz_ham(s,n1); dm_multo(eddptr->T2q[1],md2); dm_mm(0.5,md1,md3,0.5,eddptr->T2q[1]);
			  eddptr->T2q[2] = Iz_ham(s,n1); dm_multo(eddptr->T2q[2],md3);
			  md4 = Ip_real(s,n2); dm_mm(-0.5,md1,md4,2.0,eddptr->T2q[2]);
			  free_double_matrix(md1); md1 = Ip_real(s,n1); dm_mm(-0.5/sqrt(6.0),md1,md2,1.0/sqrt(6.0),eddptr->T2q[2]);
			  eddptr->T2q[3] = Iz_ham(s,n1); dm_multo(eddptr->T2q[3],md4); dm_mm(-0.5,md1,md3,-0.5,eddptr->T2q[3]);
			  eddptr->T2q[4] = dm_mul(md1,md4); dm_muld(eddptr->T2q[4],0.5);
			  free_double_matrix(md1); free_double_matrix(md2); free_double_matrix(md3); free_double_matrix(md4);
		  } else {
			  eddptr->T2q[0] = eddptr->T2q[1] = eddptr->T2q[2] = eddptr->T2q[3] = eddptr->T2q[4] = NULL;
		  }
	  }

	  // end of common Hamiltonian creation


  if (ver) {
	  printf("Done reading spinsystem file.\n");
	  printf("Number of anisotropic interactions: %d\n",abs(Naniso));
	  printf("Spin system matrix dimension: %d\n",ss->matdim);
	  printf("Interaction Hamiltonian is%s diagonal\n", (s->Hint_isdiag ?  " " : " not") );
	  printf("Hamiltonian block structure will");
	  if (s->block_diag != 0) {
		  printf(" be based on %d nuclei types\n",Nmz);
		  printf("\t number of blocks: %d\n",LEN(blk_dims));
		  printf("\t dims:");
		  for (i=1; i<=LEN(blk_dims); i++) printf(" %d",blk_dims[i]);
		  printf("\n");
	  } else {
		  printf(" not be used\n");
	  }
	  if (s->frame == LABFRAME) {
		  printf("Simulation will be carried on in laboratory frame.\n");
	  }
	  if (s->frame == DNPFRAME) {
		  printf("Simulation will be carried on in DNP frame (electron rotating frame and laboratory frame for nuclear interactions).\n");
	  }
  }

  Tcl_Free((char*)names);

  /* test print out BEGIN */
  //blk_dm_print(s->Hiso,"readsys Hiso");
  /* test print out END */


 // exit(1);

}

