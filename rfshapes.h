#ifndef __RFSHAPES_H
#define __RFSHAPES_H


/* maximal number of rf shapes */
#define MAXRFSHAPES 128

 /* Define rf shape element structure */
 typedef struct _RFelem {
     double ampl;
     double phase;
 } RFelem;


extern RFelem* RFshapes[];

void RFshapes_init(void);
void RFshapes_reset(void);
int RFshapes_len(int slot);
RFelem* RFshapes_alloc(int len);
int RFshapes_slot(void);
void free_RFshapes(int a);

// 1.12.2024 ZT: implementing GROUP optimizations, basis definition
#include "matrix.h"
extern mat_double* GROUPbasis[];
void GROUPbasis_init(void);
int GROUPbasis_slot(void);
mat_double* GROUPbasis_alloc(int Nrows, int Ntimes);
void free_GROUPbasis(int slot);
void GROUPbasis_reset(void);







#endif
