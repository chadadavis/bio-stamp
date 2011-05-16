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

#include "stamp.h"

int ccprobcalc(int **atoms1, int **atoms2, int **prob, int lena, int lenb, 
	struct parameters *parms)

/* Calculates the probability matrix, prob, after Rossmann and Argos, given
 * the two sets of coordinates, (float)atoms1 and (float)atoms2. 
 *
 * Modified version of probcalc to corner cut.  This avoids calculating
 *  part of the matrix to save time */
{

        int i,j,k,ii,jj;
	int start,end,ll;
	int upper,uval,lstrt,lstart;

	float sum,sumsq;
        float Dij,Cij,mean,sd,const1,const2;
	float fcut,fm,fn,jt,fk1;
/* SMJS Added inverse precision squared constant */
        float prec2i;

/* SMJS Use inverse constants in rossmann */
        prec2i=1.0/(float)(parms[0].PRECISION*parms[0].PRECISION);
	const1=(1.0/parms[0].const1)*prec2i;
	const2=(1.0/parms[0].const2)*prec2i;

	mean=((float)parms[0].PRECISION*parms[0].NMEAN); 
	sd=((float)parms[0].PRECISION*parms[0].NSD); 
        sum=sumsq=0.0;
	

	if(parms[0].CCADD) {
          /* auto corner cutting defined so derive fk1 from lengths */
	  fk1 = parms[0].CCFACTOR + fabs((float)(lena-lenb-4));
	} else {
	  fk1 = parms[0].CCFACTOR;
	}
	fn = (float) lena+1;  /* len + 1 */
	fm = (float) lenb+1;
	fcut = (fk1/sqrt(2.0))*(((fm/fn)+(fn/fm))/sqrt(fm*fm+fn*fn));
	upper = 0;
/*	fprintf(parms[0].LOG,"CC factor=%f, |lena-lenb-4|=%f\n",parms[0].CCFACTOR,fabs((float)(lena-lenb-4)));
	fprintf(parms[0].LOG,"CC variables: m=%f, n=%f, fk1=%f, distance=%f\n",fm,fn,fk1,fcut);
*/

	
	lstart=0;
	for(j=0; j<lena; ++j) {
	   /* corner cutting */
	   jt=j/fn; 
	   jj=j+1;
           for(i=lstart; i<lenb; ++i)  {
	     ii=i+1;
	     if(fabs((i/fm)-jt) > fcut) {
	       if(upper) {
		 goto l1;
	       } else {
		 goto l2;
	       }
	    } else {
	      if(!upper) {
		lstrt = i;
		upper = 1;
	      }
	    }
	    if(!parms[0].BOOLEAN) {
/* Use prec2i instead of parms[0].PRECISION */
/* Now incorporate prec2i into const1 and const2 */
               prob[jj][ii]=(int)((float)parms[0].PRECISION*rossmann(&atoms1[j],&atoms2[i],
			   (i==0 || j==0),(j==lena-1 || i==lenb-1),
			   const1,const2,&Dij,&Cij)); 
               prob[jj][ii]=(int)( (float)parms[0].PRECISION*((float)((float)prob[jj][ii]-mean)/(float)(sd)));
	       sum+=(float)prob[jj][ii]; sumsq+=(float)(prob[jj][ii]*prob[jj][ii]); 
/*		prob[jj][ii]=10*parms[0].PRECISION; */
	    } else {
/* Use prec2i instead of parms[0].PRECISION */
/* Now incorporate prec2i into const1 and const2 */
	       prob[jj][ii]=(rossmann(&atoms1[j],&atoms2[i], (i==0 || j==0),(j==lena-1 || i==lenb-1),
			       const1,const2,&Dij,&Cij)>=parms[0].BOOLCUT);
	    }
	    l2:; /* jump to here if corner cutting */
         }  /* for (i=... */
	 l1: /* jump to here if corner cutting */

	 upper = 0;
	 lstart=lstrt;
     	}  /* for (j=0... */

	return 0;

} /* End of function ccprobcalc */
