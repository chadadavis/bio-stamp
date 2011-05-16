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
#include <math.h>
#include "stamp.h"

/* REVerse of MATMULT
 *
 * Since with matmult we assume the following is true:
 *
 *     X2 = R*X1 + V
 *
 *  to get X1 back we must do:
 *
 *     X1 = R^-1*(X2 - V)
 *
 *  where R^-1 is the inverse of R 
 *
 * That is what this routine does.  It makes use of a routine called MATINV from
 *  Numerical recipes.
 *
 * "R" and "V" above are the "r" and "v" supplied below.
 * X2 and X1 are both *integer* coordinate sets.
 *
 * RBR November 1992
 */
int revmatmult(float **r, float *v, int **coord, int n, int PRECISION) {

	int i,j,k;
	int *indx;

	float temp;
	float t[3];
	float **R,**RI;
	float pos;


	pos=1.0;
	/* getting the inverse matrix */
	indx=(int*)malloc(4*sizeof(int));
	RI=(float**)malloc(4*sizeof(float*));
	R=(float**)malloc(4*sizeof(float*));
	for(i=0; i<4; ++i) {
	   RI[i]=(float*)malloc(4*sizeof(float));
	   R[i]=(float*)malloc(4*sizeof(float));
	}
	for(i=0; i<3; ++i)  {
	  for(j=0; j<3; ++j) {
	    R[i+1][j+1]=r[i][j];
	  }
	}

	matinv(R,RI,pos,indx);


      for(k=0; k<n; ++k) {
	for(i=0; i<3; ++i) 
	   t[i]=((float)coord[k][i]/(float)PRECISION) - v[i];

	for(i=0; i<3; ++i) {
	   temp=0.0;
	   for(j=0; j<3; ++j) 
	      temp+=RI[i+1][j+1]*t[j]; 
	   coord[k][i]=(int)( (float)PRECISION*((float)temp));
	}  /* End of for(i=... */

      }  /* End of for(k.... */

      for(i=0; i<4; ++i) {
	free(RI[i]);
	free(R[i]);
      }
      free(R); free(RI); free(indx);
      return 0;

} /* End of function */
	       
