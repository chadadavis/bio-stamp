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
#include "stamp.h"

#define chainstring "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

/* Given a domain descriptor file, this
 *   program outputs a series of PDB format files
 *   called <ident>.pdb (where ident is from the
 *   domain_loc structure */

void exit_error();

main(int argc, char *argv[]) {

	int i,j,k;
	char c;
	char *env;
	int ndomain;
	int gottrans;
	int count;
	int hetero,dssp,graphics,waters,nucleic;
	int verbose;
	int tmp_hetero,tmp_waters,tmp_nucleic;
	char chainlabel;
	char filename[200];
	char infile[200];
	char outfile[200];
	struct domain_loc *domain;
	FILE *IN,*OUT;
	
	if(argc<3) exit_error();

	dssp=0;
	hetero=0;
	nucleic=0;
	waters=0;
	graphics=0;
	verbose=0;
	strcpy(&outfile[0],"all.pdb");

	
	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(argv[i][1]=='f') {
	      /* file name */
	      if((IN=fopen(argv[i+1],"r"))==NULL) {
	        fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
	        exit(-1);
	      }
	      strcpy(&infile[0],argv[i+1]);
	      i++;
	   } else if(argv[i][1]=='d') {
	      dssp=1;
	   } else if(strcmp(&argv[i][1],"het")==0) {
	      hetero=1;
	   } else if(strcmp(&argv[i][1],"nuc")==0) {
	      nucleic = 1;
	   } else if(argv[i][1]=='g') {
	      graphics=1;
	   } else if(strcmp(&argv[i][1],"hoh")==0 || strcmp(&argv[i][1],"HOH")==0) {
	      waters=1;
	   } else if(argv[i][1]=='v' || argv[i][1]=='V') {
	      verbose=1;
	   } else if(argv[i][1]=='o') {
             if((i+1)>=argc) exit_error();
	     strcpy(&outfile[0],argv[i+1]);
	     i++;
	   } else exit_error();
	}
	printf("TRANSFORM R.B. Russell, 1995\n");

	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
      	}

	/* count the number of domains */
	ndomain=0; 
	ndomain=count_domain(IN);
	rewind(IN);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	
	/* read in the domains */
	if(getdomain(IN,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	fclose(IN);

	if(dssp) printf(" Using DSSP files\n");
	else	 printf(" Using PDB files\n"); 

	if(hetero) printf(" Files will include heteroatoms\n");
	else	printf(" Files will not include heteroatoms\n");

	if(nucleic) printf(" Files will include DNA/RNA \n");
        else    printf(" Files will not include DNA/RNA \n");

	if(waters) printf(" Files will include waters\n");
	else printf(" Files will not include waters\n");
	
	if(graphics) {
	    if((OUT=fopen(outfile,"w"))==NULL) {
 		fprintf(stderr,"error opening file %s \n",outfile);
	        exit(-1);
	    }	
 	    printf(" All coordinates will be in file %s\n",outfile);
	    chainlabel='A';
	    /* write some initial comments in the file outfile */
	    fprintf(OUT,"REMARK Output from transform\n");
	    fprintf(OUT,"REMARK  STAMP Package (Russell and Barton Proteins, 14, 309-323, 1992)\n");
	    fprintf(OUT,"REMARK Domains were read from the file %s\n",infile);
	    fprintf(OUT,"REMARK Chains are labelled sequentially starting with 'A' and\n");
	    fprintf(OUT,"REMARK  after the order given in the file %s\n",infile);
	    fprintf(OUT,"REMARK The domains in this file are:\n");
	    count = 0;
	    for(i=0; i<ndomain; ++i) {
		fprintf(OUT,"REMARK       %s  chain %c \n",domain[i].id,chainstring[count]);
		count++;
                if(count>=36) { count = 0; }
	    }
	    count = 0;
	    if(hetero) fprintf(OUT,"REMARK Includes heteroatoms\n");
	    else fprintf(OUT,"REMARK Does not include heteroatoms\n");
	    if(nucleic) fprintf(OUT,"REMARK  Includes DNA/RNA \n");
            else    fprintf(OUT,"REMARK  Does not include DNA/RNA \n");
	    if(waters) fprintf(OUT,"REMARK Includes waters\n");
            else fprintf(OUT,"REMARK Does not include waters\n");
	} else {
 	    chainlabel='\0';
	}

	for(i=0; i<ndomain; ++i) {
	   if(!dssp) sprintf(&filename[0],"%s.pdb",domain[i].id);
	   else      sprintf(&filename[0],"%s.dssp",domain[i].id);
	   if(graphics==0) {
	      if((OUT=fopen(filename,"r"))!=NULL) {
	         fprintf(stderr,"error file %s already exists \n",filename);
		 printf(" you should either delete the existing file, \n");
		 printf(" or change the identifier in your input file\n");
	         exit(-1);
	      }
	      if((OUT=fopen(filename,"w"))==NULL) {
		 fprintf(stderr,"error opening file %s\n",filename);
	      }
	      printf(" Domain %3d, %6s => to %s\n",i+1,domain[i].id,filename);
	      fprintf(OUT,"REMARK Output from transform\n");
              fprintf(OUT,"REMARK  STAMP Package (Russell and Barton Proteins, 14, 309-323, 1992)\n");
              fprintf(OUT,"REMARK Domain/transformation was read from the file %s\n",infile);
              fprintf(OUT,"REMARK Domain name is %s\n",domain[i].id);
              if(hetero) fprintf(OUT,"REMARK Includes heteroatoms\n");
              else fprintf(OUT,"REMARK Does not include heteroatoms\n");
              if(waters) fprintf(OUT,"REMARK Includes waters\n");
              else fprintf(OUT,"REMARK Does not include waters\n");


	   } else { 	
              printf(" Domain %3d, %6s => to %s (chain %c)\n",
	       i+1,domain[i].id,outfile,chainlabel);
	   }
	   
	   tmp_hetero=hetero; tmp_waters=waters; tmp_nucleic=nucleic;
	   if((IN=openfile(domain[i].filename,"r"))==NULL) {
	      fprintf(stderr,"error: PDB file %s does not exist.  Skipping this domain.\n",domain[i].filename);
	   } else {
	      for(j=0; j<domain[i].nobj; ++j) {
		 if((j==0 && chainlabel=='\0') || (chainlabel=='A')) k=1;
		 else k=0;
		 if(graphics) k=0;
	         if(!dssp) extract_pdb(IN,domain[i].start[j],domain[i].end[j],domain[i].type[j],
			     domain[i].R,domain[i].V,k,tmp_hetero,tmp_nucleic,tmp_waters,chainlabel,verbose,filename,OUT);
		 else extract_dssp(IN,domain[i].start[j],domain[i].end[j],domain[i].type[j],
			     domain[i].R,domain[i].V,k,chainlabel,OUT);
	   	 closefile(IN,domain[i].filename);
	   	 IN=openfile(domain[i].filename,"r");
		 tmp_hetero=0;
		 tmp_nucleic=0;
		 tmp_waters=0; /* only output hetero-atoms, waters, nucleic acid */
	      }
	      closefile(IN,domain[i].filename);
	   }
	   if(graphics==0) {
              fclose(OUT);
	   } else {
	      count++;
	      if(count>=36) {
			count = 0;
	         	printf("Warning: Chains starting from 'A' again\n");
	      }
	      chainlabel=chainstring[count];
	   }
	}

	exit(0);
}
void exit_error()
{
	  fprintf(stderr,"format: transform -f <domain descriptor file> [ -het -hoh -nuc -d -g -v ]\n");
	  fprintf(stderr,"        -o <combined output file> (-g only) \n");
	  fprintf(stderr,"        -het ==> include all heteroatoms\n");
	  fprintf(stderr,"        -hoh ==> include all waters\n");
	  fprintf(stderr,"        -d ==> DSSP file (PDB is default)\n");
	  fprintf(stderr,"        -g ==> label chains sequentially for graphical output\n");
	  fprintf(stderr,"	  -v ==> include all non ATOM/HETATM records in the PDB files\n");
	  exit(-1);
}
