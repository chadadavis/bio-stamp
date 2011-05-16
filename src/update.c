/*
Copyright (1997,1998,1999,2010) Robert B. Russell & Geoffrey J. Barton

This file is part of STAMP.

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. A copy of the license
can be found in the LICENSE file in the STAMP installation directory.

STAMP was developed by Robert B. Russell and Geoffrey J. Barton of
current addresses:

 Prof. Robert B. Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
 Cell Networks, University of Heidelberg            College of Life Sciences
 Room 564, Bioquant                                 University of Dundee
 Im Neuenheimer Feld 267                            Dow Street
 69120 Heidelberg                                   Dundee DD1 5EH
 Germany                                            UK
                                                
 Tel: +49 6221 54 513 62                            Tel: +44 1382 385860
 Fax: +49 6221 54 514 86                            FAX: +44 1382 385764
 Email: robert.russell@bioquant.uni-heidelberg.de   E-mail g.j.barton@dundee.ac.uk
 WWW: http://www.russell.embl-heidelberg.de         WWW: http://www.compbio.dundee.ac.uk

 All use of STAMP must cite: 

 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*/
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

