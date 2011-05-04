#include <stdio.h>
#include <math.h>

void fmatmult(float **r, float *v, float **coord, int n) {

	int i,j,k;
	float temp;
	float t[3];


        for(k=0; k<n; ++k) {
	 for(i=0; i<3; ++i) 
	   t[i]=(float)coord[k][i];

	  for(i=0; i<3; ++i) {
	   temp=0.0;
	   for(j=0; j<3; ++j) 
	      temp+=r[i][j]*t[j]; 
	   coord[k][i]=(float)temp+v[i];
	  }  /* End of for(i=... */

        }  /* End of for(k.... */


        return;

} /* End of function */
	       
