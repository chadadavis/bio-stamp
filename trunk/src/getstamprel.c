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
