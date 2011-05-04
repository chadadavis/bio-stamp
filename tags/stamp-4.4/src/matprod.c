#include <stdio.h>
#include <stdlib.h>

/* multiplies two 3x3 matrices 
 * the sense is P = A x B */

int matprod(float **P, float **A, float **B, FILE *OUTPUT) {

	int i,j,k,l;

	float **T;
	float *V;
	T=(float**)malloc(3*sizeof(float*));
	V=(float*)malloc(3*sizeof(float));
	for(i=0; i<3; ++i) {
	  V[i]=0;
	  T[i]=(float*)malloc(3*sizeof(float));
	  for(j=0; j<3; ++j) T[i][j]=0.0;
	}
	for(i=0; i<3; ++i) {
	  for(j=0; j<3; ++j) {
	    for(k=0; k<3; ++k) {
	        T[i][j]+=A[i][k]*B[k][j];
	      }
	  }
	}
/*	printf("matprod:\n");
	printf("A:\n");
	printmat(A,V,3,OUTPUT);
	printf("x B:\n");
	printmat(B,V,3,OUTPUT);
	printf("is C:\n");
	printmat(T,V,3,OUTPUT);
*/
	for(i=0; i<3; ++i) {
	   for(j=0; j<3; ++j) {
	      P[i][j]=T[i][j];
	   }
	}
	free(V);
	for(i=0; i<3; ++i) free(T[i]);
	free(T);
	return 0;
}
