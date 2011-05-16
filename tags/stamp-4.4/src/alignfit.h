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

#define MAXslen 10000
#define MAXnbloc 1000
/* SMJS Increased MAXtlen from 200 */
#define MAXtlen 1000

struct brookn { /* structure to represent brookhaven residue numbers */
    int n;      /* numerical part of number */
    char cid;   /* chain identifier */
    char in;    /* insertion code */
   };

struct domain_loc {		/* This structure allows rather complex domains to be described */
   char filename[100];
   char id[100];
   int nobj;			/* The number of objects considered within the named file */
   int *type;			/* The type that each object is:
					0 ==> an error
					1 ==> All of the residues in the file
					2 ==> A particular chain
					3 ==> A particular named region (eg. A 25 _ to B 10 _ ) */
   struct brookn *start;	/* There will be a start and end for each 'object' */
   struct brookn *end;
   int *reverse;
   float **coords;
   char *aa;
   char *align;
   char *oldalign;
   char *sec;
   struct brookn *numb;
   int ncoords;
   int *use;
   float **R;			/* Initial transformation */
   float *V;
   float **r;			/* current transformation, when STAMP is done, we must update(r,v,R,V) to get the final transformation */
   float *v;
   float value;
   };

struct parameters {
   char bloc_file[200];
   char dom_file[200];
   char parm_file[200];
   char prefix[200];
   int MAX_SEQ_LEN;
   int PAIRWISE;
   int TREEWISE;
   int OLDFORMAT;
   char MATFILE[100];
   char TREEFILE[100];
   char ORDFILE[100];
   char TRANSFILE[100];
   char STAMPDIR[200];
   };
/* Standard structure for storing protein sequence data */

struct seqdat {	/* all lengths include char terminator and [0] */
    int ilen;	/* length of identifier*/
    char *id;	/* identifier */
    int tlen;	/* length of title */
    char *title;/* title */
    int slen;	/* length of sequence*/
    char *seq;	/* sequence */
};
/* tree structures */

#if !defined(CLUST_STRUCT) 
struct indclust {
	int number;
	int *member;
	};
struct cluster {
	struct indclust a;
	struct indclust b;
	};
#define CLUST_STRUCT
#endif
#include "gjutil.h"

struct cluster *readtree(char *tordfile, char *treefile, int *number, int method, FILE *OUT);
float fmatfit(float **atoms1, float **atoms2, float **R, float *V, int nats, int entry);
struct cluster *get_clust(double **matrix, char **ids, int ndomain, char *noc_parms);
int closefile(FILE *handle,char *filename);
FILE *openfile(char *filename, char *type);
char ltou(char c);
char utol(char c);
void fmatmult(float **r, float *v, float **coord, int n);
void update(float **dR, float **R, float *dV, float *V);
int count_domain(FILE *IN);
void rmsp(char *c);
int getca(FILE *IN, float **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct brookn end, int type,  int MAXats, int REVERSE, FILE *OUTPUT);
int newoutput(FILE *TRANS, struct domain_loc *domain, int ndomain, int writetrans);
void exit_error();
int getpars(FILE *fp, struct parameters *var);
int Agetbloc(FILE *bfile, struct seqdat *bloc, int *nbloc);
int getdomain(FILE *IN, struct domain_loc *domains, int *ndomain, int maxdomain, 
        int *gottrans, char *env, int DSSP, FILE *OUTPUT);

int getdomain_error(char *buff); 
