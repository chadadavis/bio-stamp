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
