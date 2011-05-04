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

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include "dstamp.h"

/* Given a structure of STAMP values, and specifications about
 *  which data to use, this returns an array of  booleans where 
 *  1 indicates that the corresponding blocfile position is
 *  part of a reliably aligned region based on the supplied criteria.
 *
 * 'stamp' is an array of type structdat
 * 'nval' is the number of STAMP parameters in this array
 * 'npos' is the number of positions in the alignment
 * 'type' is a character, which must be one of the 
 *	data types ('G' == > Pij') read in from the file
 * 'cuttoff' is a float telling what values to consider
 * 'window' is the minimum stretch of values above 'cutoff' which
 *	can constitute a reliable region 
 *
 * The routine returns an integer array, the memory for which is
 *  allocated within this routine */

int *getstamprel(stamp,nval,npos,type,cutoff,window)
struct stampdat *stamp;
int nval,npos;
char type;
float cutoff;
int window;
{
	int i,j,which,neighbors;
	int *rel;

	/* find which of nval to use according to the supplied type */
	which=-1;
	for(i=0; i<nval; ++i) 
	   if(stamp[i].what==type) which=i;
	
	if(which==-1) {
	   fprintf(stderr,"error: specified STAMP type not found\n");
	   return NULL;
	}

	rel=(int*)malloc((npos+1)*sizeof(int));

	/* first run through and set values to 1 or 0 according to
	 *  whether they are greater than or equal to the cutoff */
	
	for(i=0; i<npos; ++i) {
	   if( (stamp[which].n[i]>=cutoff && type!='A') ||
	       (stamp[which].n[i]<=cutoff && type=='A' && stamp[which].n[i]>0.0000)) rel[i]=1;
	   else rel[i]=0;
	  /* 'A' corresponds to a distance in Angtroms, so we must consider all
	   *  positions less than this */
	}

	/* now smooth the array out according to the window */
	for(i=0; i<npos; ++i) {
	 if(rel[i]==1) {
	   neighbors=0;
	   for(j=1; j<=window; ++j) {
	      if(i+j>(npos-1) || rel[i+j]==0) break;
	      else neighbors++;
	   }
	   for(j=1; j<=window; ++j) {
	      if((i-j)<0 || rel[i-j]==0) break;
	      else neighbors++;
	   }
	   if(neighbors>=(window-1)) rel[i]=1;
	   else rel[i]=0;
	 }
	}
/*	printf("rel: ");
	for(i=0; i<npos; ++i) printf("%1d",rel[i]);
	printf("\n"); */
	return rel;
}
