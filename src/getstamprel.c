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

int *getstamprel(struct stampdat *stamp, int nval, int npos, char type, float cutoff, int window) {

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

/* SMJS Changed malloc to calloc */
	rel=(int*)calloc((npos+1),sizeof(int));

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
