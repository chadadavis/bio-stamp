#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* returns the sequence identity for a pair of sequences 
 *
 * nid 		is the number of identical positions
 * align_len	is the corrected alignment length (ie. after 
 *                removing N- and C- tails 
 *
 * the routine will return -1.0 if an error is encountered */

float seq_identity(char *seq1, char *seq2, int *nid, int *align_len, FILE *OUTPUT) {

	
	int i,len;

	len=strlen(seq1);

	(*nid)=0;
	(*align_len)=len;

	if(len!=strlen(seq2)) {
	  fprintf(stderr,"error: sequences are of different lengths (%d and %d)\n",
	     len,strlen(seq2));
	  return -1.0;
	}

	/* get number of identities */
	for(i=0; i<len; ++i) {
	   if(seq1[i]!=' ' && seq2[i]!=' ' && seq1[i]==seq2[i]) (*nid)++;
	}

	/* correct the alignment length */

	/* N-terminal overhang */
	for(i=0; i<len; ++i) {
	   if(seq1[i]!=' ' && seq2[i]!=' ') break;
	    (*align_len)--;
	}

	/* C-terminal overhang */
	for(i=len-1; i>=0; --i) {
	   if(seq1[i]!=' ' && seq2[i]!=' ') break;
           (*align_len)--;
	}

	return (float)(100*(float)(*nid)/(float)(*align_len));

}
