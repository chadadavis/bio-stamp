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

 The WORK is Copyright (1997,1998) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "stamp.h"
#include "stamprel.h"

/* Clean up a block file by attaching lone isolated
 *  residues to the nearest continuous segment */

main(int argc, char *argv[]) {

	int i,j;
	int nbloc,minlen,bloclen;
	int nstamp,nstamppos,nstampseq;
	int ndomain,gottrans;

	char c;
	char *env;

	FILE *BLOC;

	struct seqdat *bloc;
	struct stampdat *stamp;
	struct domain_loc *domain;

	if(argc!=3) {
	  printf("format: stamp_clean (block file) (minimum continuous segment length) > (output file)\n");
	  exit(-1);
	}

	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
        }


	sscanf(argv[2],"%d",&minlen);


	if((BLOC=fopen(argv[1],"r"))==NULL) {
		fprintf(stderr,"error opening file %s\n",argv[1]);
		exit(-1);
	}
	printf("%% ALIGN_CLEAN, R.B. Russell, 1995\n%% Searching for domain descriptors...\n");
	ndomain=count_domain(BLOC);
	rewind(BLOC);
	if(ndomain!=0) {
	  domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	  if(getdomain(BLOC,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
/*	  int getdomain(FILE *IN, struct domain_loc *domains, int *ndomain, int maxdomain, int *gottrans, char *env, int DSSP, FILE *OUTPUT); */

	  printf("%%   %d domain descriptions read in\n",ndomain);
	  for(i=0; i<ndomain; ++i) printdomain(stdout,domain[i],gottrans);
	} else {
	  printf("%%   no domain information found\n");
	}
	printf("%% Reading block file...\n");
	nbloc=0;
	while((c=getc(BLOC))!=(char)EOF) nbloc+=(c=='>');
	rewind(BLOC);
	bloc=(struct seqdat*)malloc((nbloc+1)*sizeof(struct seqdat));
	if(Agetbloc(BLOC,bloc,&nbloc)==-1) exit(-1);
	rewind(BLOC);
	bloclen=strlen(&bloc[1].seq[1]);
	printf("%% Searching for STAMP data...\n");
	nstamp=0;
	while((c=getc(BLOC))!=(char)EOF) nstamp+=(c=='#');
	rewind(BLOC);
	if(nstamp>0) {
	   stamp=(struct stampdat*)malloc(nstamp*sizeof(struct stampdat));
	   if(getstampdat(stamp,BLOC,&nstamp,&nstampseq,&nstamppos,bloclen)==-1) exit(-1);
	   if(nstamppos!=bloclen) {
	      fprintf(stderr,"error: STAMP and sequence data disagree\n");
	      exit(-1);
	   }
	   printf("%%   %d STAMP fields found: ",nstamp);
	   for(i=0; i<nstamp; ++i) printf("%c ",stamp[i].what);
	   printf("\n");
	} else {
	  printf("%%   no STAMP data found\n");
	}
	fclose(BLOC);
	printf("%% Block file contains %d sequences; the alignment length is %d\n",
	       nbloc,strlen(&bloc[1].seq[1]));
	printf("%% Cleaning up allowing continuous segments of %d or greater...\n",minlen);
	if(nstamp==0) clean_block(bloc,nbloc,minlen);
	else stamp_clean_block(bloc,nbloc,minlen,stamp,nstamp);
	printf("%%  Cleaning done.\n");
	bloclen=strlen(&bloc[1].seq[1]);
	printf("%% The final alignment length is %d\n",bloclen);
	printf("%% The alignment:\n");
	for(i=0; i<nbloc; ++i)  {
	   printf(">%s %s\n",bloc[i+1].id,bloc[i+1].title);
	}
	for(i=0; i<nstamp; ++i) {
/* SMJS Removed \n  (its included in the title) */
	   printf("#%c %s",stamp[i].what,stamp[i].title);
	}
	printf("*\n");
	for(i=0; i<bloclen; ++i) {
	   for(j=0; j<nbloc; ++j) 
	      printf("%c",bloc[j+1].seq[i+1]);
	   if(nstamp>0 && stamp[0].n[i]>-0.001) {
	     printf(" ");
	     for(j=0; j<nstamp; ++j) {
	      if(stamp[j].what=='T') printf("%1.0f ",stamp[j].n[i]);
	      else printf("%10.5f ",stamp[j].n[i]);
             }
	   }
	   printf("\n");
	}
	printf("*\n");
	exit(0);
}
