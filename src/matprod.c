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
