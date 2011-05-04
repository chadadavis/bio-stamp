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
#define MAX_SEQ_LEN 3000

/* Reads a list of domains from the standard input (or a file),
 *  and  outputs a list of sequences in PIR format */

/* Modification 17 October 1995: now outputs FASTA format as well
 *  with the option -format fasta */
/* Modification  7 February 1996: now can output each sequence to a
 *  separate file (<id>.seq) if the flag '-separate' is used 
 * Modification 7 March 1996 Fixed naff tendancy to put an extra
 *  residue on the first line (call me anal-retentive call me what you will) 
 * Modification 15 January 1997 Now outputs a more useful summary on
 *  the command line rather than "Chain B taken from PDB file /disk3/pdb/pdb1gsf.ent" */

main(int argc, char *argv[]) {

	char c;
	char *env;

	char infile[200];
	char outfilename[200];
	char format[50];
	char buff[1000];
	int ftype;
	int i,j,k,nbloc,ndomain,total,add;
	int min_len,max_len;
	int sepfiles;
	int verbose;
	int title_limit;
	int *skip;
	FILE *IN,*OUT,*PDB;
	struct domain_loc *domain;

	verbose=0;
	min_len = 20;
	max_len = 5000;
	title_limit = -1; /* All of the title will be printed out (MAX 1000 hardwired) */
	strcpy(&format[0],"NBRF(PIR)");
	ftype=0; sepfiles=0;

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
	   } else if(strcmp(&argv[i][1],"title_limit") ==0 || strcmp(&argv[i][1],"tl")==0) {
	      if((i+1)>=argc) exit_error();
              sscanf(argv[i+1],"%d",&title_limit);
	      i++;
	   } else if(strcmp(&argv[i][1],"separate")==0) {
	      sepfiles=1;
	   } else if(strcmp(&argv[i][1],"V")==0 || strcmp(&argv[i][1],"v")==0) {
	      verbose=1;
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
	if(verbose==1 && ftype==0) {
  	  printf("\nPDBSEQ, R.B. Russell 1995\n Extracts amino acid sequence from PDB files\n\n");
  	  printf("Min sequence length %4d, Maximum %4d\n",min_len, max_len);
	  printf("Sequence format will be %s\n",format);
	}
	/* read in list of domains */
	nbloc=0;
	nbloc=count_domain(IN);
/*	while((c=(char)getc(IN))!=(char)EOF) nbloc+=(c=='{'); */
	rewind(IN);
	if(verbose==1 && ftype==0) printf("Reading in domain descriptions...\n");
	domain=(struct domain_loc*)malloc(nbloc*sizeof(struct domain_loc));
	if(getdomain(IN,domain,&ndomain,nbloc,&i,env,0,stdout)==-1) exit(-1);
	if(ndomain!=nbloc) {
	   fprintf(stderr,"error: something wrong with input file %s\n",infile);
	   exit(-1);
	}
	if(ftype==0 && verbose==1) printf("Reading sequences...\n");
	/* get the sequences from the brookhaven files */
	skip=(int*)malloc(ndomain*sizeof(int));
	for(i=0; i<ndomain; ++i) {
	   skip[i]=0;
	   if(ftype==0 && verbose==1) printf("Domain %3d %s %s ",i+1,domain[i].filename,domain[i].id); 
	   if((PDB=fopen(domain[i].filename,"r"))==NULL) {
	      if(ftype==0 && verbose==1) printf("\nError: file %s does not exist\n",domain[i].filename);
	      if(ftype==0 && verbose==1) printf("\nSkipping this domain...\n");
	      skip[i]=1;
	   }
	   if(skip[i]==0) {
	     /* Whiz through the PDB file and try to get a name for the protein */
	     /* For the moment will read HEADER/SOURCE/TITLE lines,
	      *  TITLE will take precedence */
	     domain[i].align=(char*)malloc(1000*sizeof(char));
	     domain[i].align[0]='\0';
	     while(fgets(buff,79,PDB)!=NULL && strncmp(buff,"ATOM  ",6)!=0 && strlen(domain[i].align)<920) {
		 if(((strncmp(buff,"TITLE ",6)==0) && (strncmp(&buff[10],"MOL_ID:",7)!=0)) || 
		    ((strncmp(buff,"COMPND",6)==0) && (strncmp(&buff[10],"MOL_ID:",7)!=0)) || 
		    ((strncmp(buff,"SOURCE",6)==0) && (strncmp(&buff[10],"MOL_ID:",7)!=0))) {
			k=71;
			if(strlen(buff)<=71) { k=strlen(buff)-2; } 
			for(j=k; j>0; --j) {
			   if(buff[j] == ' ') { buff[j] = '\0'; }
			   else { break; } 
		        }
		/* TITLE     APOSTREPTAVIDIN, PH 5.6, TWO MOLECULES OF (SO4)2 BOUND AT     1SLF   3 
		   TITLE    2 THE BIOTIN BINDING SITE                                      1SLF   4
		   COMPND    MOL_ID: 1;                                                    1SLF   5
		   COMPND   2 MOLECULE: STREPTAVIDIN;                                      1SLF   6
		   COMPND   3 CHAIN: B, D                                                  1SLF   7 */
		/* 012345678901234567890123456789012345678901234567890123456789012345678901234567890 */
			sprintf(&domain[i].align[strlen(domain[i].align)],"%s ",&buff[10]);
		 }
	     }
	     rewind(PDB);
	     if(title_limit!=-1) domain[i].align[title_limit]='\0';
	     for(j=0; j<strlen(domain[i].align); ++j) {
		if(domain[i].align[j]=='\n') domain[i].align[j]=' ';
	     }
	     domain[i].ncoords=0;
	     domain[i].aa=(char*)malloc((MAX_SEQ_LEN+1)*sizeof(char)); 
	     total=0;
             domain[i].coords=NULL;
	     for(j=0; j<domain[i].nobj; ++j) {
	      if(igetca(PDB,domain[i].coords,&domain[i].aa[total],&domain[i].numb[total],
	                &add,domain[i].start[j],domain[i].end[j],
		        domain[i].type[j],(MAX_SEQ_LEN-total),	
		        domain[i].reverse[j],1000,stdout)==-1) {
		   if(ftype==0 && verbose==1) fprintf(stderr,"Domain %s skipped, couldn't get sequence\n",domain[i].id);
		   skip[i]=1;
		   break;
	      }
	      if(ftype==0 && verbose==1) switch(domain[i].type[j]) {
	  	  case 1: printf(" all residues"); break;
		  case 2: printf(" chain %c",domain[i].start[j].cid); break;
		  case 3: printf(" from %c %4d %c to %c %4d %c",
		   domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
		   domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
	       }
	       if(ftype==0 && verbose==1) printf("%4d CAs ",add);
	       total+=add;
	       rewind(PDB);
	     }
	     domain[i].ncoords=total;
	     if(ftype==0 && verbose==1) printf("= %4d CAs in total\n",domain[i].ncoords);
	     if(domain[i].ncoords<min_len || domain[i].ncoords>max_len) {
		if(ftype==0 && verbose==1) printf("Sequence is too short or too long, will ignore\n");
		skip[i]=1;
	     }
	     fclose(PDB);
	     if(ftype==0 && verbose==1) { printf("Descriptor is %s\n",domain[i].align); }
	   }
	}
	if(ftype==0 && verbose==1) printf("  ...done\n\n\n");
	/* Modification, now outputs a few different formats */
	

	for(i=0; i<nbloc; ++i) if(skip[i]==0) {
	  if(sepfiles==1) {
		sprintf(outfilename,"%s.seq",domain[i].id);
		if((OUT=fopen(outfilename,"r"))!=NULL) {
		   fprintf(stderr,"Error: file %s already exists - delete first\n",outfilename);
		} 
		fclose(OUT);
		if((OUT=fopen(outfilename,"w"))==NULL) {
		   fprintf(stderr,"Error opening file %s\n",outfilename);
		   exit(-1);
		}
		printf("Writing sequence of domain %s to file %s\n",domain[i].id,outfilename);
	  } else {
		OUT=stdout;
	  }
	  if(ftype==0) fprintf(OUT,">P1;%s ",domain[i].id);
          else         fprintf(OUT,">%s ",domain[i].id);
	  if(ftype==0) fprintf(OUT,"\n");
	  fprintf(OUT,"%s : ",domain[i].align);
	  for(j=0; j<domain[i].nobj; ++j) { 
	     if(j>0) fprintf(OUT,"& ");
	     switch(domain[i].type[j]) {
		case 1: fprintf(OUT,"All "); break;
		case 2: fprintf(OUT,"Chain %c ",domain[i].start[j].cid); break;
		case 3: fprintf(OUT,"%c%d%c-%c%d%c ",
			  domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
			  domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
	     }
	  }
	  fprintf(OUT,"\n");
	  for(j=0; j<strlen(domain[i].aa); ++j) {
	     fprintf(OUT,"%c",domain[i].aa[j]);
	     if(((j+1)%80)==0) fprintf(OUT,"\n");
	  }
	  if(ftype==0) fprintf(OUT,"*");
	  fprintf(OUT,"\n");
	  if(sepfiles==1) fclose(OUT);
	}
	exit(0);
}
int exit_error() {

	   fprintf(stderr,"format: pdbseq -f <domain descriptor file> [-min <val> -max <val>] > <output file>\n");
	   fprintf(stderr,"               -format <pir, fasta> -separate -tl <title string max length>\n");
	   fprintf(stderr,"               -separate => write each sequence to it's own file\n");
	   exit(-1);
}