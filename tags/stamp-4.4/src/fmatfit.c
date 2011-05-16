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
