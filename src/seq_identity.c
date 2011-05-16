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
