/*
 * rfprof.h
 *
 *  Created on: Jun 1, 2010
 *      Author: zdenek
 */

#ifndef RFPROF_H_
#define RFPROF_H_

typedef struct _rfmap_struct {
	int loop_size, Nz, Nch;
	int *loop, *z;
	double *weight, *b1;
} rfmap_struct;

double ** rfprof_alloc(int N, int chan);
void rfprof_free(double ** obj);
double ** read_rfproffile(const char * fname, int chan);
int rfprof_len(double ** obj);
double rfprof_sumweight(double **rfdata);
//rfmap_struct* read_rfmap(const char *fname, int chan);
rfmap_struct* read_rfmap(const char *fname, int chan, int avestep); //updated version allows to skip some iphi values
void rfmap_free(rfmap_struct *rfmap);
void rfmap_get(rfmap_struct *rfmap, int iz, int iphi, int ich, double *bx, double *by);

#endif /* RFPROF_H_ */
