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
