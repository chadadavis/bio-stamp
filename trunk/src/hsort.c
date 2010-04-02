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

#include "array.h"
#include "pathdef.h"

/* modified heap sort from nr routines */

void hsort(n,sortarr)
int n;
/*float ra[];*/
struct path *sortarr;
{
	int l,j,ir,i;
/*	float rra;*/
	struct path rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
/*			rra=ra[--l];*/
			rra=sortarr[--l];
		else {
/*			rra=ra[ir];*/
			rra=sortarr[ir];
/*			ra[ir]=ra[1];*/
			sortarr[ir]=sortarr[1];
			if (--ir == 1) {
/*				ra[1]=rra;*/
				sortarr[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
/*			if (j < ir && ra[j] < ra[j+1]) ++j;*/
			if (j < ir && sortarr[j].score > sortarr[j+1].score) ++j;
/*			if (rra < ra[j]) {*/
			if (rra.score > sortarr[j].score){
/*				ra[i]=ra[j];*/
				sortarr[i]=sortarr[j];
				j += (i=j);
			}
			else j=ir+1;
		}
/*		ra[i]=rra;*/
		sortarr[i]=rra;
	}
}
