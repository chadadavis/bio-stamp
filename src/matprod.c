/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 Appendix A of the file LICENSE.  Conditions and restrictions for use of
 this package are also in this file.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Biomolecular Modelling Laboratory   Laboratory of Molecular Biophysics
 Imperial Cancer Research Fund       The Rex Richards Building
 Lincoln's Inn Fields, P.O. Box 123  South Parks Road
 London, WC2A 3PX, U.K.              Oxford, OX1 3PG, U.K.
 Tel: +44 171 269 3583               Tel: +44 865 275368
 FAX: +44 171 269 3417               FAX: 44 865 510454
 e-mail: russell@icrf.icnet.uk       e-mail gjb@bioch.ox.ac.uk
 WWW: http://bonsai.lif.icnet.uk/    WWW: http://geoff.biop.ox.ac.uk/

 The WORK is Copyright (1992,1993,1995,1996) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>

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
