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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_STAMP_NUM 10
#define MAX_N_SEQ 1000
struct parameters{
   char filename[200];
   char TYPE;
   float CUTOFF;
   int WINDOW;
   int SECDISP;
   int SMALLSEQ,SMALLSEC;
   int SECSUM;
   int CASESEQ,CASESEC;
   int VERBOSE;
};
struct stampdat{
   char what;
   char *title;
   float *n;
};
struct seqdat { /* all lengths include char terminator and [0] */
    int ilen;   /* length of identifier*/
    char *id;   /* identifier */
    int tlen;   /* length of title */
    char *title;        /* title */
    int slen;   /* length of sequence*/
    char *seq;  /* sequence */
};
int getstampdat(struct stampdat *stamp, FILE *IN, int *nstamp, int *nseq, int *npos, int maxpos);
int *getstamprel(struct stampdat *stamp, int nval, int npos, char type, float cutoff, int window);
int Agetbloc(FILE *bfile, struct seqdat *bloc, int *nbloc);
char utol(char c);
char ltou(char c);
void exit_error();
