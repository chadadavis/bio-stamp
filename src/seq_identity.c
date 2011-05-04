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

 The WORK is Copyright (1992,1993,1995,1996) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>

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
