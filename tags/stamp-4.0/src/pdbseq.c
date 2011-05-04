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
#include <stdlib.h>
#include "include.h"
#define MAX_SEQ_LEN 3000

/* Reads a list of domains from the standard input (or a file),
 *  and  outputs a list of sequences in PIR format */

/* Modification 17 October 1995, now outputs FASTA format as well
 *  with the option -format fasta */

main(argc,argv)
int argc;
char *argv[];
{

	char c;
	char *env;

	char infile[200];
	char format[50];
	int ftype;
	int i,j,k,nbloc,ndomain,total,add;
	int min_len,max_len;
	int *skip;
	float **coords1;
	FILE *IN,*OUT,*PDB;
	struct domain_loc *domain;

	min_len = 20;
	max_len = 5000;
	strcpy(&format[0],"NBRF(PIR)");
	ftype=0;

	if(argc<3) exit_error();
	for(i=1; i<argc; ++i) {
           if(argv[i][0]!='-') exit_error();
           if(strcmp(&argv[i][1],"f")==0)  {
              /* file name */
              if((IN=fopen(argv[i+1],"r"))==NULL) {
                fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
                exit(-1);
              }
              strcpy(&infile[0],argv[i+1]);
              i++;
           } else if(strcmp(&argv[i][1],"min")==0) {
	      if((i+1)>=argc) exit_error();
              sscanf(argv[i+1],"%d",&min_len);
	      i++;
	   } else if(strcmp(&argv[i][1],"max")==0) {
              if((i+1)>=argc) exit_error();
              sscanf(argv[i+1],"%d",&max_len);
              i++;
	   } else if(strcmp(&argv[i][1],"format")==0) {
	      if((i+1)>=argc) exit_error();
	      strcpy(&format[0],argv[i+1]);
	      for(j=0; j<strlen(format); ++j) format[j]=ltou(format[j]);
	      if(strcmp(format,"NBRF(PIR)")==0 || strcmp(format,"NBRF")==0 || strcmp(format,"PIR")==0) {
		  ftype=0;
	      } else if(strcmp(format,"FASTA")==0) {
		   ftype=1;
	      } else {
		fprintf(stderr,"error: format %s not recognised\n",format);
	      }
	      i++;
           } else {
	      exit_error();
	   }
        }


	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
        }

	if((IN=fopen(infile,"r"))==NULL) {
	     fprintf(stderr,"error: file %s not found\n",infile);
	     exit(-1);
	}
	if(ftype==0) {
  	  printf("\nPDBSEQ, R.B. Russell 1995\n Extracts amino acid sequence from PDB files\n\n");
  	  printf("Min sequence length %4d, Maximum %4d\n",min_len, max_len);
	  printf("Sequence format will be %s\n",format);
	}
	/* read in list of domains */
	nbloc=0;
	nbloc=count_domain(IN);
/*	while((c=getc(IN))!=(char)EOF) nbloc+=(c=='{'); */
	rewind(IN);
	if(ftype==0) printf("Reading in domain descriptions...\n");
	domain=(struct domain_loc*)malloc(nbloc*sizeof(struct domain_loc));
	if(getdomain(IN,domain,&ndomain,nbloc,&i,env,0,stdout)==-1) exit(-1);
	if(ndomain!=nbloc) {
	   fprintf(stderr,"error: something wrong with input file %s\n",infile);
	   exit(-1);
	}
	if(ftype==0) printf("Reading sequence...\n");
	/* get the sequences from the brookhaven files */
	skip=(int*)malloc(ndomain*sizeof(int));
	for(i=0; i<ndomain; ++i) {
	   skip[i]=0;
	   if(ftype==0) printf("Domain %3d %s %s\n   ",i+1,domain[i].filename,domain[i].id);
	   if((PDB=fopen(domain[i].filename,"r"))==NULL) {
	      if(ftype==0) printf("\nError: file %s does not exist\n",domain[i].filename);
	      if(ftype==0) printf("\nSkipping this domain...\n");
	      skip[i]=1;
	   }
	   if(skip[i]==0) {
	     domain[i].ncoords=0;
	     domain[i].aa=(char*)malloc((MAX_SEQ_LEN+1)*sizeof(char)); 
	     total=0;
	     if(ftype==0) printf("   ");
	     for(j=0; j<domain[i].nobj; ++j) {
	      if(getca(PDB,&domain[i].coords,&domain[i].aa[total],&domain[i].numb[total],
	         &add,domain[i].start[j],domain[i].end[j],
		 domain[i].type[j],(MAX_SEQ_LEN-total),domain[i].reverse[j],1)==-1) {
		if(ftype==0) printf("\nSkipping this domain...\n");
		skip[i]=1;
		break;
	      }
	      if(ftype==0) switch(domain[i].type[j]) {
	  	  case 1: printf(" all residues"); break;
		  case 2: printf(" chain %c",domain[i].start[j].cid); break;
		  case 3: printf(" from %c %4d %c to %c %4d %c",
		   domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
		   domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
	       }
	       if(ftype==0) printf("%4d CAs ",add);
	         total+=add;
	         rewind(PDB);
	       }
	     domain[i].ncoords=total;
	     if(ftype==0) printf("= %4d CAs in total\n",domain[i].ncoords);
	     if(domain[i].ncoords<min_len || domain[i].ncoords>max_len) {
		if(ftype==0) printf("Sequence is too short or too long, will ignore\n");
		skip[i]=1;
	     }
	     fclose(PDB);
	   }
	}
	if(ftype==0) printf("  ...done\n");
	if(ftype==0) printf("\n\n");
	/* Modification, now outputs a few different formats */
	

	for(i=0; i<nbloc; ++i) if(skip[i]==0) {
	  printf(">%s ",domain[i].id);
	  if(ftype==0) printf("\n");
	  for(j=0; j<domain[i].nobj; ++j) { 
	     if(j>0) printf("and ");
	     switch(domain[i].type[j]) {
		case 1: printf("All residues "); break;
		case 2: printf("Chain %c ",domain[i].start[j].cid); break;
		case 3: printf("From %c %d %c to %c %d %c ",
			  domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
			  domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
	     }
	  }
	  printf("taken from PDB file %s\n",domain[i].filename);
	  for(j=0; j<strlen(domain[i].aa); ++j) {
	     printf("%c",domain[i].aa[j]);
	     if(j>0 && j%80==0) printf("\n");
	  }
	  if(ftype==0) printf("*");
	  printf("\n");
	}
	exit(0);
}
int exit_error() {

	   fprintf(stderr,"format: pdbseq -f <domain descriptor file> [-min <val> -max <val>] > <output file>\n");
	   fprintf(stderr,"               -format <pir, fasta>\n");
	   exit(-1);
}
