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
#include <math.h>
#include <malloc.h>
#include <time.h>

#include "general.h"
#include "domain.h"
#include "stamprel.h"

#define MAX_SEQ_LEN 2000
#define PRECISION 1000

/* avestruc: reads in a stamp output file and generates an average structure file */
main(argc, argv)
int argc; 
char *argv[];
{
	int i,j,k,l,test;
	int ndomain,total,add,nbloc;
	int nstamp,nseq,nstamppos;
	int gottrans;
	int bloclen;
	int stampwindow;
	int polyA,nrel;
	int ncoords,nres;
	int all_aligned;

	int *pointer;
	int *counter;
	int *all;
	int *rel;
	int **avecoords;
	int ***n,***c,***o,***cb;
	

	char filename[200],outfile[200];
	char stampchar;

	char *env;
	char *buff;
	char *aa;
	char *N,*C,*O;

	float stampthresh;
	float *average;

	FILE *TRANS,*PDB,*OUT;

	struct domain_loc *domain;
	struct seqdat *bloc;
	struct stampdat *stamp;
	struct brookn *numb;

	if(argc<3) exit_error();
	
	/* defaults */
	polyA=0;
	stampchar='G';
	stampwindow=3;
	stampthresh=6.0;
	strcpy(&outfile[0],"average.pdb");
	N=(char*)malloc(5*sizeof(char));
	C=(char*)malloc(5*sizeof(char));
	O=(char*)malloc(5*sizeof(char));
	strcpy(N," N  ");
	strcpy(C," C  ");
	strcpy(O," O  ");
	
	average=(float*)malloc(3*sizeof(float));
	buff=(char*)malloc(100*sizeof(char));
	all_aligned=0;

	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(strcmp(&argv[i][1],"f")==0) {
	      if((i+1)>=argc) exit_error();
	      strcpy(&filename[0],argv[i+1]);
	      i++;
	   } else if(strcmp(&argv[i][1],"c")==0) {
	     if((i+1)>=argc) exit_error();
             stampchar=argv[i+1][0];
             i++;
 	   } else if(strcmp(&argv[i][1],"t")==0) {
	     if((i+1)>=argc) exit_error();
	     sscanf(argv[i+1],"%f",&stampthresh);
	     i++;
	   } else if(strcmp(&argv[i][1],"w")==0) {
	     if((i+1)>=argc) exit_error();
             sscanf(argv[i+1],"%d",&stampwindow);		
             i++;
	   } else if(strcmp(&argv[i][1],"polyA")==0) {
	     polyA=1;
	   } else if(strcmp(&argv[i][1],"o")==0) {
	     if((i+1)>=argc) exit_error();
             strcpy(&outfile[0],argv[i+1]);
             i++;
	   } else if(strcmp(&argv[i][1],"aligned")==0) {
	     all_aligned=1;
	   } else {
	     exit_error();
	   }
	}

	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
      	}
	/* read in coordinate locations and initial transformations */
	if((TRANS = fopen(filename,"r")) == NULL) {
	   fprintf(stderr,"error: file %s does not exist\n",argv[1]);
	   exit(-1);
	}
	if((OUT=fopen(outfile,"w"))==NULL) {
	  fprintf(stderr,"error opening file %s\n",outfile);
	  exit(-1);
	}
	printf("AVESTRUC R.B. Russell, 1995\n");
	printf(" Using STAMP file %s\n",filename);
	if(polyA) printf(" Will make an averaged poly-alanine structure\n");
	else printf(" Will make an averaged C-alpha-only structure\n");
	printf(" Averaged coordinates will be in %s\n",outfile);
	fprintf(OUT,"REMARK  Averaged coordinates \n");
	fprintf(OUT,"REMARK  Calculated using the transformations and alignments\n");
	fprintf(OUT,"REMARK   in the STAMP file %s\n",filename);
	fprintf(OUT,"REMARK   R. B. Russell, 1994\n");
	fprintf(OUT,"REMARK   \n");

	/* determine the number of domains specified */
	printf(" Reading domain descriptors...\n");
	ndomain=count_domain(TRANS);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	rewind(TRANS);
	if(getdomain(TRANS,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	pointer=(int*)malloc(ndomain*sizeof(int));
	if(polyA) {
	  n=(int***)malloc(ndomain*sizeof(int**));
	  c=(int***)malloc(ndomain*sizeof(int**));
	  o=(int***)malloc(ndomain*sizeof(int**));
	  cb=(int***)malloc(ndomain*sizeof(int**));
	  aa=(char*)malloc(MAX_SEQ_LEN*sizeof(char));
	  numb=(struct brookn*)malloc(MAX_SEQ_LEN*sizeof(struct brookn));
	}
	for(i=0; i<ndomain; ++i) {
	  if(polyA) {
	    n[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	    o[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	    c[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	    cb[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	  }
	  for(j=0; j<domain[i].nobj; ++j) {
	     if(domain[i].start[j].cid=='_') domain[i].start[j].cid=' ';
	     if( domain[i].start[j].in=='_') domain[i].start[j].in=' ';
	     if(  domain[i].end[j].cid=='_') domain[i].end[j].cid=' '; 
	     if(   domain[i].end[j].in=='_') domain[i].end[j].in=' ';
	  }
	  fprintf(OUT,"REMARK Domain %4d Name is %s\n",i+1,domain[i].id);
	}

	/* read in the alignment */
	printf(" Reading alignment...\n");
	rewind(TRANS);
	bloc=(struct seqdat*)malloc((ndomain*2+2)*sizeof(struct seqdat));
	printf(" ");
	Agetbloc(TRANS,bloc,&nbloc);
	bloclen=strlen(&bloc[1].seq[1]);
	if(nbloc!=(ndomain*2+1)) {
	   fprintf(stderr,"error: number of domains and alignment disagree\n");
	   exit(-1);
	}
	counter=(int*)malloc(bloclen*sizeof(int));
	all=(int*)malloc(bloclen*sizeof(int));

	fprintf(OUT,"REMARK   \n");
	/* read in the STAMP output and get reliable regions */
	rewind(TRANS);
	if(all_aligned) {
	   rel=(int*)malloc(MAX_SEQ_LEN*sizeof(int));
	   for(i=0; i<bloclen; ++i) {
	     k=0;
	     for(j=0; j<ndomain; ++j) {
		if(bloc[j+1].seq[i+1]!=' ') k++;
	     }
	     if(k==ndomain) rel[i]=1;
	     else rel[i]=0;
	   }
        } else { 
	   stamp=(struct stampdat*)malloc(100*sizeof(struct stampdat));
	   if(getstampdat(stamp,TRANS,&nstamp,&nseq,&nstamppos,MAX_SEQ_LEN)==-1) exit(-1);
	   if((rel=getstamprel(stamp,nstamp,bloclen,stampchar,stampthresh,stampwindow))==NULL) exit(-1);
	}
	nrel=0;
	   for(i=0; i<bloclen; ++i) {
             if(rel[i]==1) nrel++;
	}
	fprintf(OUT,"REMARK %4d STAMP reliable or core positions found\n",nrel);
	fprintf(OUT,"REMARK The characters to the far right of the coordinates show\n");
	fprintf(OUT,"REMARK  the residues occuring at each position within the alignment\n");
	fclose(TRANS);

	printf(" Reading coordinates...\n");
	for(i=0; i<ndomain; ++i) {
	   /* output the domain descriptors again */
	   printf(" Domain %3d %s %s\n   ",i+1,domain[i].filename,domain[i].id); 
	   if((PDB=fopen(domain[i].filename,"r"))==NULL) {
	      fprintf(stderr,"error: file %s does not exist\n",domain[i].filename);
	      exit(-1);
	   }
	   domain[i].ncoords=0;
	   domain[i].coords=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	   domain[i].aa=(char*)malloc((MAX_SEQ_LEN+1)*sizeof(char)); 
	   domain[i].numb=(struct brookn*)malloc((MAX_SEQ_LEN)*sizeof(struct brookn));
	   total=0;
	   printf("    "); 
	   for(j=0; j<domain[i].nobj; ++j) {
	       if(igetca(PDB,&domain[i].coords[total],&domain[i].aa[total],&domain[i].numb[total],
		  &add,domain[i].start[j],domain[i].end[j],
		  domain[i].type[j],(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
	       if(polyA) {
		  rewind(PDB);
		  if(igetgen(PDB,&n[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],domain[i].type[j],N,(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
		  if(test!=add) {
			fprintf(stderr,"error: different number of N atoms in %s\n",domain[i].id);
		        exit(-1);
		  }
		  rewind(PDB);
                  if(igetgen(PDB,&c[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],
                    domain[i].type[j],C,(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
                  if(test!=add) {
                        fprintf(stderr,"error: missing C atoms in %s\n",domain[i].id);
                        exit(-1);
                  }     
		  rewind(PDB);
                  if(igetgen(PDB,&o[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],
                     domain[i].type[j],O,(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
                  if(test!=add) {
                        fprintf(stderr,"error: missing O atoms in %s\n",domain[i].id);
                        exit(-1);
                  }
		  rewind(PDB);
                  if(igetcb(PDB,&cb[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],
                    domain[i].type[j],(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
                  if(test!=add) {
                        fprintf(stderr,"error: missing CB atoms in %s\n",domain[i].id);
                        exit(-1);
                  }
	       }
	       switch(domain[i].type[j]) {
	 	  case 1: printf(" all residues"); break; 
		  case 2: printf(" chain %c",domain[i].start[j].cid); break;
		  case 3: printf(" from %c %4d %c to %c %4d %c",
			 domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
			 domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
		} 
		if(polyA) printf("   %4d main chain/C-betas\n",add);
		else printf("   %4d C-alphas\n",add); 
	        total+=add;
	        rewind(PDB);
	    }
	    domain[i].ncoords=total;
	    if(polyA) printf("  A total of %4d sets of main chain/C-beta atoms\n",domain[i].ncoords);
	    else printf("  A total of %4d C-alphas\n",domain[i].ncoords);
	    /* disp(domain[i],stdout); */
	    printf(" Applying transformation... \n");
/*	    printmat(domain[i].R,domain[i].V,3,stdout);
	    printf("      ...to these coordinates.\n");  */
/*	    if(polyA) {
               for(j=0; j<5; ++j) {
                  printf("Res %3d %c N : %8d %8d %8d\n",j+1,domain[i].aa[j],n[i][j][0],n[i][j][1],n[i][j][2]);
                  printf("Res %3d %c CA: %8d %8d %8d\n",j+1,domain[i].aa[j],domain[i].coords[j][0],domain[i].coords[j][1],domain[i].coords[j][2]);
                  printf("Res %3d %c C : %8d %8d %8d\n",j+1,domain[i].aa[j],c[i][j][0],c[i][j][1],c[i][j][2]);
                  printf("Res %3d %c O : %8d %8d %8d\n",j+1,domain[i].aa[j],o[i][j][0],o[i][j][1],o[i][j][2]);
                  printf("Res %3d %c CB: %8d %8d %8d\n",j+1,domain[i].aa[j],cb[i][j][0],cb[i][j][1],cb[i][j][2]);
               }
	    } */
	    matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,PRECISION);  
	    if(polyA) {
	       matmult(domain[i].R,domain[i].V,n[i],domain[i].ncoords,PRECISION);  
               matmult(domain[i].R,domain[i].V,c[i],domain[i].ncoords,PRECISION);  
               matmult(domain[i].R,domain[i].V,o[i],domain[i].ncoords,PRECISION);  
	       matmult(domain[i].R,domain[i].V,cb[i],domain[i].ncoords,PRECISION);
/*	       for(j=0; j<5; ++j) {
		  printf("Res %3d %c N : %8d %8d %8d\n",j+1,domain[i].aa[j],n[i][j][0],n[i][j][1],n[i][j][2]);
		  printf("Res %3d %c CA: %8d %8d %8d\n",j+1,domain[i].aa[j],domain[i].coords[j][0],domain[i].coords[j][1],domain[i].coords[j][2]);
 		  printf("Res %3d %c C : %8d %8d %8d\n",j+1,domain[i].aa[j],c[i][j][0],c[i][j][1],c[i][j][2]);
		  printf("Res %3d %c O : %8d %8d %8d\n",j+1,domain[i].aa[j],o[i][j][0],o[i][j][1],o[i][j][2]);
		  printf("Res %3d %c CB: %8d %8d %8d\n",j+1,domain[i].aa[j],cb[i][j][0],cb[i][j][1],cb[i][j][2]); 
	       } */
	    }
/*	    printmat(domain[i].R,domain[i].V,3,stdout); */
	    /* disp(domain[i],stdout); */
	    fclose(PDB);
	}
	
	/* Now proceed through the alignment and calculate an average coordinate using the
	 *  reliable positions */
	for(i=0; i<ndomain; ++i) pointer[i]=0;
	ncoords=0; nres=0;
	for(i=0; i<bloclen; ++i) {
	  if(rel[i]==1) {
	    for(j=0; j<3; ++j) average[j]=0.0;
	    if(polyA) {
		/* N atoms */
		for(j=0; j<3; ++j) average[j]=0.0;
		for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)n[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
		fprintf(OUT,"ATOM  %5d  N   ALA Z%4d    %8.3f%8.3f%8.3f  0.00  0.00   0 \n",
                  ncoords+1,nres+1,average[0],average[1],average[2]);
		ncoords++;
	    }
	    /* CA atoms */	
 	    for(j=0; j<3; ++j) average[j]=0.0;
	    for(j=0; j<ndomain; ++j) {	
		for(k=0; k<3; ++k) average[k]+=((float)domain[j].coords[pointer[j]][k]/(float)ndomain)/(float)PRECISION;
	    }
	    /* Output the average coordiates */ 
	    sprintf(buff,"ATOM  %5d  CA  ALA Z%4d    %8.3f%8.3f%8.3f  0.00  0.00   0 ",
		ncoords+1,nres+1,average[0],average[1],average[2]);
	    fprintf(OUT,"%s",buff);
	    for(j=0; j<ndomain; ++j) 
               if(j<13) fprintf(OUT,"%c",bloc[j+1].seq[i+1]); 
	    fprintf(OUT,"\n");
	    ncoords++;
	    if(polyA) {
		/* C O and CB atoms */
                for(j=0; j<3; ++j) average[j]=0.0;
                for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)c[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
                fprintf(OUT,"ATOM  %5d  C   ALA Z%4d    %8.3f%8.3f%8.3f  0.00  0.00   0 \n",
                  ncoords+1,nres+1,average[0],average[1],average[2]);
                ncoords++;
                for(j=0; j<3; ++j) average[j]=0.0;
                for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)o[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
                fprintf(OUT,"ATOM  %5d  O   ALA Z%4d    %8.3f%8.3f%8.3f  0.00  0.00   0 \n",
                  ncoords+1,nres+1,average[0],average[1],average[2]);
                ncoords++;
                for(j=0; j<3; ++j) average[j]=0.0;
                for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)cb[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
                fprintf(OUT,"ATOM  %5d  CB  ALA Z%4d    %8.3f%8.3f%8.3f  0.00  0.00   0 \n",
                  ncoords+1,nres+1,average[0],average[1],average[2]);
                ncoords++;
	    }
	    nres++; 
	  }
	  for(j=0; j<ndomain; ++j) {
	      if(bloc[j+1].seq[i+1]!=' ') pointer[j]++;
	  }
	}
	free(pointer);
	printf(" ...done.\n");
	exit(0);
}
int exit_error()
{
	 fprintf(stderr,"format:  avestruc -f <file> -c <char> -w <win> -t <thresh> -polyA -aligned\n");
	 fprintf(stderr,"          where -polyA ==> poly alanine structure (default is CAs only)\n");
 	 fprintf(stderr,"                -aligned ==> consider all aligned positions\n");
	 exit(-1);
}
