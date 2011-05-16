#include <stdio.h>

/* This function takes in a change matrix and vector (dR, dV) and 
 *   applies it to an 'old' matrix and vector (R,V).  The result
 *   is copied into the 'old' R and V.  */

void update(float **dR, float **R, float *dV, float *V) {

	int i,j,k;
	float rtemp[3][3],vtemp[3];
/*
	printf("dR,dV:\n");
	printmat(dR,dV,3,stdout);
	printf("R,V:\n");
	printmat(R,V,3,stdout);
*/
	/* Matrix multiplication and vector addition */
	for(i=0; i<3; ++i) {
	   vtemp[i]=0;
	   for(j=0; j<3; ++j) {
	      rtemp[i][j]=0;
	      vtemp[i]+=dR[i][j]*V[j];
	      for(k=0; k<3; ++k)
		 rtemp[i][j]+=dR[i][k]*R[k][j];
	    } /* End of for(j=... */
	} /* End of for(i=.... */

	/* Transfering rtemp into R */
	for(i=0; i<3; ++i) {
	   V[i]=vtemp[i]+dV[i];
	   for(j=0; j<3; ++j) R[i][j]=rtemp[i][j];
	   } /* End of for(i=... */
/*	printf("After combining:\nR,V:\n");
	printmat(R,V,3,stdout); */
}

