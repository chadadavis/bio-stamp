#include <stdio.h>
#include <math.h>
#include "stamp.h"
#include "rbr_vector.h"
#include "rbr_newlib.h"
#include "rbr_aadat.h"

#define MAX_SEQ_LEN 100000
#define PRECIS 1000

float RBR_ipowdist(int *x1, int *x2, int PRECISION, int power);
float RBR_r_of_gyration_pep(int *Ro, int **coords, int ncoords, char *seq, int w, int PRECISION, int DIST_POWER);
float RBR_r_max_pep(int *Ro, int **coords, int ncoords, char *seq, int PRECISION, int DIST_POWER);
float RBR_r_of_gyration(int *Ro, int **coords, int ncoords, int *masses, int w,   int PRECISION, int DIST_POWER);

int *RBR_c_of_m_pep(int **coords, int ncoords, char *seq, int w, float *total_mass);
int *RBR_c_of_m(int **coords, int ncoords, int *masses, int w);

float RBR_ran3(int *idum);
