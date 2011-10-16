#include <stdio.h>
#include <math.h>

/* special integer atoms version for STAMP */
int matmult(float **r, float *v, int **coord, int n, int PRECISION) {

      int i,j,k;
      float temp;
      float t[3];
	
/*	printmat(r,v,3,stdout); */

      for(k=0; k<n; ++k) {
	for(i=0; i<3; ++i) 
	   t[i]=(float)coord[k][i]/(float)PRECISION;
/*	for(i=0; i<3; ++i)
	   printf("%8d ",coord[k][i]);
	printf("--> "); */

	for(i=0; i<3; ++i) {
	   temp=0.0;
	   for(j=0; j<3; ++j) 
	      temp+=r[i][j]*t[j]; 
	   coord[k][i]=(int)( (float)PRECISION*((float)temp+v[i]));
	}  /* End of for(i=... */
/*	for(i=0; i<3; ++i)
           printf("%8d ",coord[k][i]);
	printf("\n"); */

      }  /* End of for(k.... */

      return 0;

} /* End of function */
	       
