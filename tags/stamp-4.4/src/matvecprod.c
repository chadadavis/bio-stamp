#include <stdio.h>
#include <math.h>

/* Multiplies a vector by a matrix.  The
 *  sense is: C = A x  B */
int matvecprod(float **A, float *C, float *B, FILE *OUTPUT) {

      	int i,j,k;
      	float temp;
      	float t[3];

/*	printf("matvecprod:\n");
	printf("A x V:\n");
	printmat(A,B,3,OUTPUT); */

	for(j=0; j<3; ++j) {
	   t[j]=0.0;
	   for(k=0; k<3; ++k) {
	      t[j]+=A[j][k]*B[k];
	   }
	}
	for(i=0; i<3; ++i) C[i]=t[i];
	
/*	printf("product:\n");
	for(i=0; i<3; ++i) 
	 printf("%8.5f\n",C[i]); */

   	return 0;

} /* End of function */
	       
