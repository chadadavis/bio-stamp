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
#include <math.h>
#include <igetcb.h>
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

	int i;
	float *V3;

	V3=V3=(float*)malloc(3*sizeof(float));

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
