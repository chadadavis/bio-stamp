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
#include <stdlib.h>
#include <stamp.h>

/* EXTRANS
 * Given a file of transformations and a set of domain IDs, this program
 *  simply extracts the various chosen transformations out from the file
 *  (probably one could do this in perl, but hell)
 *
 * RBR 20 June 1996 */

main(int argc, char *argv[]) {
	
	
	int i,j,k;
	int ndomain,gotdomain;
	int got_file;
	int got_id;
	int n_ids,len;
	int gottrans;

	int *indx;
	int *which;

	float sign;

	float *negvec;

	float **invmat;
	float **R,**RI;

	char c;

	char *id;
	char **id_list;
	char *buff;
	char *infile;
	char *env;
	
	FILE *IN;

	struct domain_loc *domain;

	id=(char*)malloc(100*sizeof(char));
	buff=(char*)malloc(1000*sizeof(char));
	infile=(char*)malloc(1000*sizeof(char));
	indx=(int*)malloc(100*sizeof(int));
	invmat=(float**)malloc(3*sizeof(float));
	R=(float**)malloc(4*sizeof(float));
	RI=(float**)malloc(4*sizeof(float));
	negvec=(float*)malloc(3*sizeof(float));
	for(i=0; i<4; ++i) { 
	   R[i]=(float*)malloc(4*sizeof(float));
	   RI[i]=(float*)malloc(4*sizeof(float));
	}
	for(i=0; i<3; ++i) 
	  invmat[i]=(float*)malloc(3*sizeof(float));


	if(argc<3) exit_error();

	got_file=got_id=0;
	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(argv[i][1]=='f') { 
	      if((i+1)>=argc) exit_error();
	      if((IN=fopen(argv[i+1],"r"))==NULL) {
		 fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
		 exit(-1);
	      }
	      got_file=1;
	      strcpy(infile,argv[i+1]);
	      i++;
	   } else if(argv[i][1]=='i') {
	      if((i+1)>=argc) exit_error(); /* Must give at least on ID */
	      j=i+1; n_ids=0; got_id=1;
	      id_list = (char**)malloc(sizeof(char*));
	      while(j<argc && argv[j][0]!='-') { /* Read in the ids sought */
		 len=strlen(argv[j]);
		 id_list[n_ids]=(char*)malloc((len+1)*sizeof(char));
	         strcpy(id_list[n_ids],argv[j]);
/*		 printf("ID %d is %s\n",n_ids+1,id_list[n_ids]); */
		 n_ids++;
		 id_list=(char**)realloc(id_list,(n_ids+1)*sizeof(char*));
		 j++;
	      }
	      i=j-1;
	   } else {
		exit_error();
	   }
	}

	if(!got_file) {
	  fprintf(stderr,"error: must specify file name\n");
	  exit(-1);
	}

	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
        }

	ndomain=count_domain(IN);
	rewind(IN);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	if(getdomain(IN,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	if(!gottrans) {
	  fprintf(stderr,"error: file does not contain transformations\n");
	  exit(-1);
	}
	rewind(IN);

	if(!got_id) {
	  fprintf(stderr,"error: must specify on or more identifiers from domain file\n");
	  fprintf(stderr," The identifiers in the file are:\n");
	  for(i=0; i<ndomain; ++i) fprintf(stderr,"          %s\n",domain[i].id);
	  exit(-1);
	}
	
	for(i=0; i<n_ids; ++i) {
	    for(j=0; j<strlen(id_list[i]); ++j) id_list[i][j]=utol(id_list[i][j]);
	}
	   
	/* find which domain ID to use */
	which=(int*)malloc(ndomain*sizeof(int));
	for(i=0; i<ndomain; ++i) {
	  which[i]=0;
	  for(j=0; j<n_ids; ++j) { 
	     if(strcmp(id_list[j],domain[i].id)==0) {
	      which[i]=1;
	     }
	  }
	}

	/* put some comments in to tell which domain is centered */
	printf("%%\n%%\n%% EXTRANS R.B. Russell, 1995\n");
	printf("%%  The input file was %s\n",infile);
	printf("%%  the following domains were extracted from the file\n");
	for(j=0; j<n_ids; ++j) {
		printf("%% %s \n",id_list[j]);
	}
	printf("%%  The original comments have been removed.\n");
	printf("%%\n%%\n");
	
	for(i=0; i<ndomain; ++i) if(which[i]==1) {
	   printdomain(stdout,domain[i],1);
	}

	exit(0);
}

int exit_error()
{
	fprintf(stderr,"format: extrans -f <transformation file> -i <ids to extract> \n");
	exit(-1);
}
