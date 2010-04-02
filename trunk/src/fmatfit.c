/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)	            Prof. Geoffrey J. Barton (GJB)
 EMBL Heidelberg                            School of Life Sciences
 Meyerhofstrasse 1                          University of Dundee
 D-69117 Heidelberg                         Dow Street
 Germany                                    Dundee, DD1 5EH
                                          
 Tel: +49 6221 387 473                      Tel: +44 1382 345860
 FAX: +44 6221 387 517                      FAX: +44 1382 345764
 E-mail: russell@embl-heidelberg.de         E-mail geoff@compbio.dundee.ac.uk
 WWW: http://www.russell.emb-heidelberg.de  WWW: http://www.compbio.dundee.ac.uk

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <math.h>

/*  This function mimics the first part of the FORTRAN routine
 *    'matfit.f'.  It makes use of the FORTRAN routine 'qkfit.f'.
 *     This keeps complications involving passing great heaps
 *     of memory to FORTRAN from arising.  The resulting transformation
 *     is in the sense of:
 *	   atoms1 = R*atoms2 + V     */

/* SMJS Added */
#include "stamp.h"

float fmatfit(float **atoms1, float **atoms2, float **R, float *V, int nats, int entry) {

	/* Note that everything that is passed to the routine
	 *  qkfit must be double precision. */

/* SMJS Commented
	void qkfit();
*/
	double cma[3], cmb[3];
	double umat[3][3];
	double xasq, xbsq, xni,t,rtsum,rmse;
	double r[3][3];
	int i,j,k,xn;

/* SMJS Added */
        long long_entry;

	xn=nats;
	xasq=0.0;
	xbsq=0.0;
	xni=1.0/((double)xn);

	/* Accumulate uncorrected (for c.m.) sums and squares */

	for(i=0; i<3; ++i) {
	   cma[i]=cmb[i]=0.0;
	   for(j=0; j<3; ++j) {
		umat[j][i]=0.0;
	   }

           for(j=0; j<nats; ++j) {
	      for(k=0; k<3; ++k) {
	         umat[k][i]+=(double)atoms1[j][i]*(double)atoms2[j][k];
	      }

/*
	      printf("Matfit...\n");
	      for(k=0; k<3; ++k) printf("%10.5f ",atoms1[j][k]);
	      printf("and ");
	      for(k=0; k<3; ++k) printf("%10.5f ",atoms2[j][k]);
	      printf("\n");
*/

	      t=(double)atoms1[j][i];
	      cma[i]+=t;
	      xasq+=t*t;
	      t=(double)atoms2[j][i];
	      cmb[i]+=t;
	      xbsq+=t*t;

	   } 
	}  

	/* Subtract cm offsets */

	for(i=0; i<3; ++i) {
	   xasq-=cma[i]*cma[i]*xni;
	   xbsq-=cmb[i]*cmb[i]*xni;
	   for(j=0; j<3; ++j)  {
	       umat[j][i]=(umat[j][i]-cma[i]*cmb[j]*xni)*xni;
	       r[i][j] = 0.0;
	   }
	} 

	/* Fit it */

/* SMJS Changed entry to long_entry */
        long_entry = entry;
/*	printf("Calling qkfit...\n"); */
	qkfit(&umat[0][0],&rtsum,&r[0][0],&long_entry);
/* SMJS Added */
        entry = (int)long_entry;

	rmse=(xasq+xbsq)*xni - 2.0*rtsum;
	if(rmse<0.0) rmse=0.0;
	rmse=sqrt(rmse);

	/* The matrix obtained from the FORTRAN routine (r) must be transfered to
	 *  the 'malloc'ed C matrix (R).  */
        for(i=0; i<3; ++i) {
           for(j=0; j<3; ++j)  {
/*               printf("%10f => %10f \n",R[i][j],r[j][i]);  */
               R[i][j]=(float)r[j][i];
           }
/*	   printf("\n");  */
        }


	/* Calculate offset if entry=1 */

	if(!entry) return rmse;
	for(i=0; i<3; ++i) {
	   t=0.0;
	   for(j=0; j<3; ++j)  {
		t+=((double)R[i][j])*cmb[j];
	   }
	   V[i]=(float)(cma[i]-t)*xni;
	} 

	return rmse;

} 
