#include <stdio.h>
#include <malloc.h>
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
	return 0;
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
