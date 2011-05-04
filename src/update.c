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

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>

/* This function takes in a change matrix and vector (dR, dV) and 
 *   applies it to an 'old' matrix and vector (R,V).  The result
 *   is copied into the 'old' R and V.  */

void update(dR,dV,R,V)
float **dR,**R;
float *dV,*V;
{
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

