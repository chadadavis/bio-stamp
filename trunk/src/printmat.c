#include <stdio.h>

int printmat(float **R, float *V, int n, FILE *OUTPUT) {

	int i,j;
	for(i=0; i<n; ++i) {
   	   for(j=0; j<n; ++j) fprintf(OUTPUT,"%10.5f ",R[i][j]);
	   fprintf(OUTPUT,"         %10.5f \n",V[i]);
	}
	return 0;
}
