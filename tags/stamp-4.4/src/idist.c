#include <stdio.h>
#include <math.h>

/* Finds the square distance between to ints *atm1, *atm2 */

float idist(int *atm1, int *atm2, int PRECISION) {

	int i;
	float D;
	
	D=0.0;
	for(i=0; i<3; ++i) 
	   D+=((float)(atm1[i]-atm2[i])/(float)PRECISION)*((float)(atm1[i]-atm2[i])/(float)PRECISION);
	return sqrt(D);
}

