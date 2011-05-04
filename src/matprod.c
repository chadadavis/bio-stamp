/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current contact addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Bioinformatics                      EMBL-European Bioinformatics Institute
 SmithKline Beecham Pharmaceuticals  Wellcome Trust Genome Campus
 New Frontiers Science Park (North)  Hinxton, Cambridge, CB10 1SD U.K.
 Harlow, Essex, CM19 5AW, U.K.       
 Tel: +44 1279 622 884               Tel: +44 1223 494 414
 FAX: +44 1279 622 200               FAX: +44 1223 494 468
 e-mail: russelr1@mh.uk.sbphrd.com   e-mail geoff@ebi.ac.uk
                                     WWW: http://barton.ebi.ac.uk/

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
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
