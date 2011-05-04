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

/* Given a domain descriptor file, this
 *   program outputs a series of PDB format files
 *   called <ident>.pdb (where ident is from the
 *   domain_loc structure */

main(argc,argv)
int argc;
char *argv[];
{
	int i,j,k;
	char c;
	char *env;
	int ndomain;
	int gottrans;
	int hetero,dssp,graphics,waters,verbose;
	char chainlabel;
	char filename[200];
	char infile[200];
	char outfile[200];
	struct domain_loc *domain;
	FILE *IN,*OUT;
	
	if(argc<3) exit_error();

	dssp=0;
	hetero=0;
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
	if(dssp) printf(" Using DSSP files\n");
	else	 printf(" Using PDB files\n"); 
	if(hetero) printf(" Files will include heteroatoms\n");
	else	printf(" Files will not include heteroatoms\n");
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
	} else {
 	    chainlabel='\0';
	}
	/* count the number of domains */
	ndomain=0; 
	ndomain=count_domain(IN);
	rewind(IN);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	
	/* read in the domains */
	if(getdomain(IN,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	fclose(IN);

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

	   } else { 	
              printf(" Domain %3d, %6s => to %s (chain %c)\n",
	       i+1,domain[i].id,outfile,chainlabel);
	   }
	   
	   if((IN=fopen(domain[i].filename,"r"))==NULL) {
	      fprintf(stderr,"error: PDB file %s does not exist.  Skipping this domain.\n",domain[i].filename);
	   } else {
	      for(j=0; j<domain[i].nobj; ++j) {
		 if((j==0 && chainlabel=='\0') || (chainlabel=='A')) k=1;
		 else k=0;
		 if(graphics) k=0;
	         if(!dssp) extract_pdb(IN,domain[i].start[j],domain[i].end[j],domain[i].type[j],
			     domain[i].R,domain[i].V,k,hetero,waters,chainlabel,verbose,filename,OUT);
		 else extract_dssp(IN,domain[i].start[j],domain[i].end[j],domain[i].type[j],
			     domain[i].R,domain[i].V,k,chainlabel,OUT);
	   	 rewind(IN);
	      }
	      fclose(IN);
	   }
	   if(graphics==0) {
              fclose(OUT);
	   } else {
	      chainlabel+=1;
	   }
	}

	exit(0);
}
int exit_error()
{
	  fprintf(stderr,"format: transform -f <domain descriptor file> -het -d -g\n");
	  fprintf(stderr,"        -o <combined output file> (-g only) \n");
	  fprintf(stderr,"        -het ==> include all heteroatoms\n");
	  fprintf(stderr,"        -hoh ==> include all waters\n");
	  fprintf(stderr,"        -d ==> DSSP file (PDB is default)\n");
	  fprintf(stderr,"        -g ==> label chains sequentially for graphical output\n");
	  fprintf(stderr,"	  -v ==> include all non ATOM/HETATM records in the PDB files\n");
	  exit(-1);
}
