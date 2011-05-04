#include <stdio.h>
#include <math.h>
#include <stamp.h>

/* This function calculates the Rossmann and Argos probability for a 
 *  single pair of (integer) coordinates.  Pij is returned, as
 *  are Dij and Cij (the distance and conformational components of
 *  Pij respectively */

float rossmann(int **atoms1, int **atoms2, int start, int end,
	float const1, float const2, float *Dij, float *Cij, int PRECISION) {

	   int i,j;
   	   float dijsq,Sijsq,t1,t2,t3,t4,t5,t6;
   	   float dx,dy,dz;
	   float Pij;
	
           if(start) t1=t2=t3=0;
           else {
               t1=(float)(abs(atoms1[0][0] - atoms2[0][0])- abs(atoms1[-1][0] - atoms2[-1][0]));
               t2=(float)(abs(atoms1[0][1] - atoms2[0][1])- abs(atoms1[-1][1] - atoms2[-1][1]));
               t3=(float)(abs(atoms1[0][2] - atoms2[0][2])- abs(atoms1[-1][2] - atoms2[-1][2]));
           } /*End of if(start)... */
               t1*=t1; t2*=t2; t3*=t3;
 
           if(end) t4=t5=t6=0;
           else {
               t4=(float)(abs(atoms1[0][0] - atoms2[0][0])- abs(atoms1[1][0] - atoms2[1][0]));
               t5=(float)(abs(atoms1[0][1] - atoms2[0][1])- abs(atoms1[1][1] - atoms2[1][1]));
               t6=(float)(abs(atoms1[0][2] - atoms2[0][2])- abs(atoms1[1][2] - atoms2[1][2]));
               t4*=t4; t5*=t5; t6*=t6;
            } /* End of if(end).... */
 
            dx=(float)(atoms1[0][0] - atoms2[0][0]);
            dy=(float)(atoms1[0][1] - atoms2[0][1]);
            dz=(float)(atoms1[0][2] - atoms2[0][2]);

            dijsq=dx*dx+dy*dy+dz*dz;
            Sijsq=t1+t2+t3+t4+t5+t6;

	    /* conversion to floating point occurs at this stage */
 	    (*Dij)=(dijsq/(float)(PRECISION*PRECISION))/const1;
            (*Cij)=(Sijsq/(float)(PRECISION*PRECISION))/const2;
	    Pij=(float)exp((double)((*Dij)+(*Cij)));
/*	    printf("Pij=%f\n",Pij); */
	    return Pij;
}
