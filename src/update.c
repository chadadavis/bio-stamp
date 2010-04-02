/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)	            Prof. Geoffrey J. Barton (GJB)
 EMBL Heidelberg                            School of Life Sciences
 Meyerhofstrasse 1                          University of Dundee
 D-69117 Heidelberg                         Dow Street
 Germany                                    Dundee, DD1 5EH
                                          
 Tel: +49 6221 387 473                      Tel: +44 1382 345860
 FAX: +44 6221 387 517                      FAX: +44 1382 345764
 E-mail: russell@embl-heidelberg.de         E-mail geoff@compbio.dundee.ac.uk
 WWW: http://www.russell.emb-heidelberg.de  WWW: http://www.compbio.dundee.ac.uk

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
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

