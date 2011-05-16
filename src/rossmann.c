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

/* This function calculates the Rossmann and Argos probability for a 
 *  single pair of (integer) coordinates.  Pij is returned, as
 *  are Dij and Cij (the distance and conformational components of
 *  Pij respectively */

/* SMJS Changed abs to ABS */
#define ABS(x) (x>=0 ? x : -(x))
#define fexp expf

/* SMJS Changed to use inverse const1 and const2 and precalculated */
/*      inverse squared PRECISION */
/* SMJS Removed precision argument */
float rossmann(int **atoms1, int **atoms2, int start, int end,
	float const1, float const2, float *Dij, float *Cij) {

   	   float dijsq,Sijsq,t1,t2,t3,t4,t5,t6;
   	   float dx,dy,dz;
	   float Pij;
/* SMJS Added vars to hopefully speed up routine */
           int *down1,*down2;
           int *up1,*up2;
           int *this1=atoms1[0],*this2=atoms2[0];
           int thisdiff[3];
	
           thisdiff[0]=ABS(this1[0] - this2[0]);
           thisdiff[1]=ABS(this1[1] - this2[1]);
           thisdiff[2]=ABS(this1[2] - this2[2]);

           if(start) t1=t2=t3=0;
           else {
               down1=atoms1[-1];down2=atoms2[-1];
               t1=(float)(thisdiff[0] - ABS(down1[0] - down2[0]));
               t2=(float)(thisdiff[1] - ABS(down1[1] - down2[1]));
               t3=(float)(thisdiff[2] - ABS(down1[2] - down2[2]));
               t1*=t1; t2*=t2; t3*=t3;
           } /*End of if(start)... */
 
           if(end) t4=t5=t6=0;
           else {
               up1=atoms1[1];up2=atoms2[1];
               t4=(float)(thisdiff[0] - ABS(up1[0] - up2[0]));
               t5=(float)(thisdiff[1] - ABS(up1[1] - up2[1]));
               t6=(float)(thisdiff[2] - ABS(up1[2] - up2[2]));
               t4*=t4; t5*=t5; t6*=t6;
            } /* End of if(end).... */

#ifdef DBGPREC
            printf("t1 = %f t2 = %f t3 = %f t4 = %f t5 = %f t6 = %f\n",t1,t2,t3,t4,t5,t6);
#endif
 
/* SMJS Use the precalculated values even though they are absolute 
            dx=(float)(atoms1[0][0] - atoms2[0][0]);
            dy=(float)(atoms1[0][1] - atoms2[0][1]);
            dz=(float)(atoms1[0][2] - atoms2[0][2]);
*/
            dx=(float)(thisdiff[0]);
            dy=(float)(thisdiff[1]);
            dz=(float)(thisdiff[2]);

            dijsq=dx*dx+dy*dy+dz*dz;
            Sijsq=t1+t2+t3+t4+t5+t6;

	    /* conversion to floating point occurs at this stage */
/* SMJS Changed to use inverse for PRECISION and const1 and const2 */
 	    (*Dij)=/*(*/dijsq/**PREC2I)*/*const1;
            (*Cij)=/*(*/Sijsq/**PREC2I)*/*const2;
/* SMJS	    Pij=(float)exp((double)((*Dij)+(*Cij))); */
	    Pij=fexp((*Dij)+(*Cij));
/*	    printf("Pij=%f\n",Pij); */
	    return Pij;
}


#define DBGPREC
float rossmanni(int **atoms1, int **atoms2, 
	const float const1, const float const2, float * const Dij, float * const Cij) {

   	   float dijsq,Sijsq,t1,t2,t3,t4,t5,t6;
   	   float dx,dy,dz;
	   float Pij;
/* SMJS Added vars to hopefully speed up routine */
           int *down1,*down2;
           int *up1,*up2;
           int *this1=atoms1[0],*this2=atoms2[0];
           int thisdiff[3];
	
           thisdiff[0]=ABS(this1[0] - this2[0]);
           thisdiff[1]=ABS(this1[1] - this2[1]);
           thisdiff[2]=ABS(this1[2] - this2[2]);

           down1=atoms1[-1];down2=atoms2[-1];
           t1=(float)(thisdiff[0] - ABS(down1[0] - down2[0]));
           t2=(float)(thisdiff[1] - ABS(down1[1] - down2[1]));
           t3=(float)(thisdiff[2] - ABS(down1[2] - down2[2]));
           t1*=t1; t2*=t2; t3*=t3;

           up1=atoms1[1];up2=atoms2[1];
           t4=(float)(thisdiff[0] - ABS(up1[0] - up2[0]));
           t5=(float)(thisdiff[1] - ABS(up1[1] - up2[1]));
           t6=(float)(thisdiff[2] - ABS(up1[2] - up2[2]));
           t4*=t4; t5*=t5; t6*=t6;
 
           dx=(float)(thisdiff[0]);
           dy=(float)(thisdiff[1]);
           dz=(float)(thisdiff[2]);

           dijsq=dx*dx+dy*dy+dz*dz;
           Sijsq=t1+t2+t3+t4+t5+t6;

	    /* conversion to floating point occurs at this stage */
 	   (*Dij)=dijsq*const1;
           (*Cij)=Sijsq*const2;
	   Pij=fexp((*Dij)+(*Cij));
	   return Pij; 
}
