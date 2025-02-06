/*
    Reads rf inhomogeneity files
    Copyright (C) 2007 Zdenek Tosner

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
    
    The rf inhomogeneity file must have the following data format:

       <N> <M>
       <scalefactor1 chan1> <scalefactor1 chan2> ... <weight1>
       <scalefactor2 chan1> <scalefactor2 chan2> ... <weight2>
       <scalefactor3 chan1> <scalefactor3 chan2> ... <weight3>
       ...
       <scalefactorN chan1> <scalefactorN chan2> ... <weightN>

    where <N> is the number of rf isochromats, <M> is the number of rf channels (is
    M is set to 1 the same profile is assumed on all channels), scalefactor is the 
    rf scale factor typically ranging from 0 to 2 (nominal field = 1). 
    Weight is relative occurance of the isochromat; sum of all weights doesn't need 
    to sum up to 1, program takes care of that...
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "matrix.h"
#include "defs.h"
#include "rfprof.h"

/****
 * alocate rfprof structure
 ****/
double ** rfprof_alloc(int N, int chan)
{
	int i;
	double ** obj=NULL;

	obj = (double**)malloc((N+1)*sizeof(double*));
	obj[0] = (double*)malloc(sizeof(double));
	*(int*)(obj[0]) = N;
	//obj[0][0] = (double)N;
	//*((int*)*obj) = N;
	for (i=1;i<=N; i++) {
		obj[i] = double_vector(chan+1);
	}
	return obj;
}
/*
 * rfprof structure length
 */
int rfprof_len(double ** obj){
	//return *((int*)*obj);
	//return (int)(obj[0][0]);
	return *(int*)(obj[0]);
}
/*
 * free rfprof structure
 */
void rfprof_free(double ** obj)
{
	int i, N;

	N = rfprof_len(obj);
	for (i=1; i<=N; i++) {
		free_double_vector(obj[i]);
	}
	free(obj[0]);
	free(obj);
}


/****
 * ZT: reads in rf profile from a file
 ****/
double ** read_rfproffile(const char* name,int nchan)
{
  FILE* fp;
  char fname[256];
  int N,i,j,M;
  int ver;
  double dum, **data;
  
  ver= verbose & VERBOSE_RFPROF;
  strcpy(fname,name);
  if (!strcmp(fname,"none")) {
	  data = rfprof_alloc(1, nchan);
      for (i=1; i<=nchan+1; i++) {
        data[1][i]=1.0;
      }
      return data;
  }
  
#ifdef UNIX
  if (name[0] == '~') {
    char* p=getenv("HOME");
    if (p != NULL) {
      strcpy(fname,p);
      strcat(fname,&name[1]);
    }
  }
#endif
  if (!(fp=fopen(fname,"r"))) {
    strcat(fname,".rf");
    fp=fopen(fname,"r");
  }
  
  if (!fp) {
    fprintf(stderr,"error: unable to open file '%s'\n",fname);
    fprintf(stderr,"\n");
    exit(1);
  } 
  
  if (ver) printf("loading external rf profile file '%s'\n",fname);

  if (fscanf(fp,"%d%d",&N,&M) != 2) {
    fprintf(stderr,"unable to read rf profile structure data from file '%s'\n",fname);
    exit(1);
  }
  data = rfprof_alloc(N,nchan);

  for (i=1;i<=N;i++) {
    for (j=1;j<=M;j++) {
       if (fscanf(fp,"%lg",&dum) != 1) {
          fprintf(stderr,"error: unable to read line %d, column %d in file %s\n",i,j,fname);
          exit(1);
       }
       data[i][j] = dum;
    }
    for (j=M+1; j<=nchan; j++) {
    	data[i][j] = 1.0;
    }
    if (fscanf(fp,"%lg",&dum) != 1) {
          fprintf(stderr,"error: unable to read weight on line %d in file %s\n",i,fname);
          exit(1);
    }
    if (dum == 0.0) {
         fprintf(stderr,"error: rf isochromat number %d in file '%s' has zero weight\n",i,name);
         exit(1);
    }
    data[i][nchan+1] = dum;
  }

  fclose(fp);

  if (ver) {
	  for (i=1; i<=N; i++) {
		  printf("%5d ",i);
		  for (j=1; j<=nchan+1; j++)
			  printf("%15g ",data[i][j]);
		  printf("\n");
	  }
  }
  
  return data;
}

double rfprof_sumweight(double **rfdata)
{
	double w=0.0;
	int i, N, NN;

	N = rfprof_len(rfdata);
	DEBUGPRINT("rfprof_sumweight: rfdata length = %d\n",N);
	NN = LEN(rfdata[1]);
	for (i=1; i<=N; i++)
		w += rfdata[i][NN];

	return w;
}

rfmap_struct* read_rfmap(const char* name,int nchan, int avestep)
{
  FILE* fp;
  char fname[256];
  int Nz, Nch, Ntot, iz, iloop, Nphi, iphi, ich, ib1;
  int ver;
  double dum, bx, by;
  rfmap_struct *rfmap;
  //int avestep = 2; // this will be input parameter
  //printf("rfmap reading: avestep = %d\n",avestep);

  ver= verbose & VERBOSE_RFPROF;
  strcpy(fname,name);
  if (!strcmp(fname,"none")) {
	  rfmap = NULL;
      return rfmap;
  }

#ifdef UNIX
  if (name[0] == '~') {
    char* p=getenv("HOME");
    if (p != NULL) {
      strcpy(fname,p);
      strcat(fname,&name[1]);
    }
  }
#endif
  if (!(fp=fopen(fname,"r"))) {
    strcat(fname,".rf");
    fp=fopen(fname,"r");
  }
  if (!fp) {
    fprintf(stderr,"error: rfmap unable to open file '%s' nor '%s'\n",name,fname);
    fprintf(stderr,"\n");
    exit(1);
  }

  if (ver) printf("loading external rfmap file '%s'\n",fname);

  if (fscanf(fp,"%d%d%d",&Nz,&Nch,&Ntot) != 3) {
    fprintf(stderr,"unable to read rfmap structure data from first line of file '%s'\nExpecting 3 integers",fname);
    exit(1);
  }
  if (Nch != nchan) {
	fprintf(stderr,"error: rfmap with wrong number of channels (rfmap: %d, spinsys: %d)\n",Nch,nchan);
	exit(1);
  }

  rfmap = (rfmap_struct*)malloc(sizeof(rfmap_struct));
  rfmap->Nch = Nch;
  rfmap->Nz = Nz;
  rfmap->z = (int*)malloc((Nch+1)*Nz*sizeof(int));
  if (rfmap->z == NULL) {
	  fprintf(stderr,"error: rfmap can not allocate z component\n");
	  exit(1);
  }
  rfmap->b1 = (double*)malloc(Ntot*2*Nch*sizeof(double));
  if (rfmap->b1 == NULL) {
	  fprintf(stderr,"error: rfmap can not allocate b1 component\n");
	  exit(1);
  }
  rfmap->loop_size = Ntot;
  rfmap->loop = (int*)malloc(2*Ntot*sizeof(int));
  if (rfmap->loop == NULL) {
	  fprintf(stderr,"error: rfmap can not allocate loop component\n");
	  exit(1);
  }
  rfmap->weight = (double*)malloc(Ntot*sizeof(double));
  if (rfmap->weight == NULL) {
	  fprintf(stderr,"error: rfmap can not allocate weight component\n");
	  exit(1);
  }

  iloop = 0;
  ib1 = 0;
  for (iz=0; iz<Nz; iz++) {
	  if (fscanf(fp,"%lg",&dum) != 1) {
		  fprintf(stderr,"error: rfmap crashed when reading weight from line %d\n",iz+1);
		  exit(1);
	  }
	  if (fscanf(fp,"%d",&Nphi) != 1) {
		  fprintf(stderr,"error: rfmap crashed when reading Nphi from line %d\n",iz+1);
		  exit(1);
	  }
	  //for (iphi=0; iphi<Nphi; iphi++) {
	  if ( (Nphi > 1) && (Nphi % avestep != 0) ) {
		  fprintf(stderr,"error reading rfmap: number of phi elements (%d) is not multiple of avesteps (%d)\n",Nphi, avestep);
		  exit(1);
	  }
	  for (iphi=0; iphi<Nphi; iphi += avestep) {
		  //rfmap->weight[iloop] = dum/Nphi;
		  rfmap->weight[iloop] = dum/Nphi*(Nphi>1 ? avestep:1);
		  rfmap->loop[iloop*2+0] = iz;
		  rfmap->loop[iloop*2+1] = iphi;
		  iloop++;
	  }
	  rfmap->z[iz*(Nch+1)+Nch] = Nphi;
	  for (ich=0; ich<Nch; ich++) {
		  rfmap->z[iz*(Nch+1)+ich] = ib1;
		  for (iphi=0; iphi<Nphi; iphi++) {
			  if (fscanf(fp,"%lg%lg",&bx,&by) != 2) {
				  fprintf(stderr,"error: rfmap crashed when reading line %d at position chan %d, couple %d\n",iz+1,ich+1,iphi+1);
				  exit(1);
			  }
			  rfmap->b1[ib1*2+0] = bx;
			  rfmap->b1[ib1*2+1] = by;
			  ib1++;
		  }
 	  }
  }
  fclose(fp);
  rfmap->loop_size = iloop; // correct the size of loop (relevant entries)

  /*** control printout
  printf("rfmap: iloop=%d, ib1=%d\n",iloop,ib1);
  printf("rfmap->z\n");
  for (iz=0; iz<Nz; iz++) {
	  for (ich=0; ich<=Nch; ich++) {
		  printf(" %3d",rfmap->z[iz*(Nch+1)+ich]);
	  }
	  printf("\n");
  }
  printf("rfmap->loop, weight\n");
  for (iloop=0; iloop<Ntot; iloop++) {
	  printf("%3d %3d %g\n",rfmap->loop[iloop*2+0],rfmap->loop[iloop*2+1],rfmap->weight[iloop]);
  }
  printf("Averaging points in rfmap->b1 \n");
  for (iloop=0; iloop<rfmap->loop_size; iloop++) {
	  iz = rfmap->loop[iloop*2+0];
	  iphi = rfmap->loop[iloop*2+1];
	  printf("%d) %d %d ",iloop,iz,iphi);
	  for (ich=0; ich<Nch; ich++) {
		  rfmap_get(rfmap,iz,iphi,ich,&bx,&by);
		  printf("ch %d (%g %g) ",ich, bx, by);
	  }
	  printf("\n");
  }
  printf("rfmap->b1\n");
  for (ib1=0; ib1<Ntot*Nch; ib1++) {
	  printf("%d) %g %g\n",ib1*2,rfmap->b1[ib1*2+0],rfmap->b1[ib1*2+1]);
  }
  printf("rfmap->b1 again and again\n");
  for (iz=0; iz<Nz; iz++) {
	  Nphi = rfmap->z[iz*(Nch+1)+Nch];
	  printf("Row (%d)",iz+1);
	  for (ich=0; ich<Nch; ich++) {
		  printf("\n chan %d (%d): ",ich+1,Nphi);
		  for (iphi=0; iphi<Nphi; iphi++) {
			  rfmap_get(rfmap,iz,iphi,ich,&bx,&by);
			  printf("(%g %g / %d %d) ",bx,by,rfmap->z[iz*(Nch+1)+ich],(rfmap->z[iz*(Nch+1)+ich]+iphi)*2);
		  }
	  }
	  printf("\n");
  }
  ***/

  return rfmap;
}

void rfmap_get(rfmap_struct *rfmap, int iz, int iphi, int ich, double *bx, double *by) {
	int Nch = rfmap->Nch;
	int Nphi = rfmap->z[iz*(Nch+1)+Nch];
	iphi %= Nphi;
	int ib1 = (rfmap->z[iz*(Nch+1)+ich]+iphi)*2;
	*bx = rfmap->b1[ib1];
	*by = rfmap->b1[ib1+1];
}

void rfmap_free(rfmap_struct *rfmap) {

	if (rfmap == NULL) return;

	free(rfmap->b1);
	free(rfmap->loop);
	free(rfmap->weight);
	free(rfmap->z);
	free(rfmap);

}
