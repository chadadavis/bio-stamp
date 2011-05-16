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
#include <math.h>
#include "igetcb.h"
#define PI 3.141592653589793

int RBR_print_vector(float *V) {
	int i;
	printf("(");
	for(i=0; i<3; ++i) {
	   printf("%8.5f",V[i]);
	   if(i<2) printf(", ");
	}
	printf(")");
	printf(" length: %8.5f\n",
	   sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]));
        return 1;
}

int RBR_vector_unify(float *V) {

	int i;
	float magnitude;

	magnitude=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
	for(i=0; i<3; ++i) 
	   V[i]=((1/magnitude)*(float)V[i]);
	
	return 0;
}

float *RBR_vector_diff(float *V1, float *V2) {
	int i;
	float *V3;

	V3=(float*)malloc(3*sizeof(float));

	for(i=0; i<3; ++i) {
	   V3[i]=V1[i]-V2[i];
	}

	return V3;
}

float *RBR_vector_ave(float *V1, float *V2) {
	int i;
	float *V3;

	V3=(float*)malloc(3*sizeof(float));

	for(i=0; i<3; ++i) {
	   V3[i]=((float)(V1[i]+V2[i])/2);
	}

	return V3;
}

float *RBR_vector_cross(float *V1, float *V2) {

	float *V3;

/* SMJS Removed V3=V3 */
	V3=(float*)malloc(3*sizeof(float));

	V3[0]=V1[1]*V2[2]-V1[2]*V2[1];
	V3[1]=V1[2]*V2[0]-V1[0]*V2[2];
	V3[2]=V1[0]*V2[1]-V1[1]*V2[0];

	return V3;
}

float *RBR_vector_set_dist(float *V1, float R) {

	/* given a vector and a distance, update the vector to have a new length */
	int i;
	float *V2;
	float r;

	V2=(float*)malloc(3*sizeof(float));
	r=sqrt(V1[0]*V1[0]+V1[1]*V1[1]+V1[2]*V1[2]);
	for(i=0; i<3; ++i) 
	   V2[i]=(R/r)*V1[i];

	return V2;
}
