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
#include <math.h>
#include <time.h>

#include <general.h>
#include <domain.h>
#include <stamprel.h>

#define MAX_SEQ_LEN 2000
#define PRECISION 1000

/* STAMP2CASP - converts STAMP alignment to CASP format */
main(int argc, char *argv[]) {

	int i,j,k,l,test;
	int ndomain,total,add,nbloc;
	int nstamp,nseq,nstamppos;
	int gottrans;
	int bloclen,maxlen;
	int stampwindow;
	int nrel;
	int start1,start2;
	int ncoords,nres1,nres2;
	int all_aligned;
	int gen_average;
	int helix,strand,coil;
	int cons_sec;
	int verbose;

	int *pointer;
	int *counter;
	int *all;
	int *rel;
	int **avecoords;
	int ***n,***c,***o,***cb;
	

	char filename[200];
	char id[200];
	char chain;
	char stampchar;
	char type;
	char *summary;

	char *env;
	char *buff;
	char *aa;
	char *sec;
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
	stampchar='G';
	stampwindow=3;
	stampthresh=6.0;
	verbose=0;
	
	buff=(char*)malloc(100*sizeof(char));
	all_aligned=0; gen_average=0;
	cons_sec=0; 

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
	   } else if(strcmp(&argv[i][1],"aligned")==0) {
	     all_aligned=1;
	   } else if(strcmp(&argv[i][1],"a")==0) {
	     gen_average=1;
	   } else if(strcmp(&argv[i][1],"cons")==0) {
	     cons_sec=1;
	   } else if(strcmp(&argv[i][1],"V")==0) {
	      verbose=1;
	   } else {
	     exit_error();
	   }
	}

	printf("STAMP2CASP R.B. Russell, 1996\n");
        printf(" Converts STAMP alignments to CASP(2) format\n");
	printf(" You must run `transform' or `avestruc' on the file:\n   %s\n",filename);
	printf("  in order to run MOLSCRIPT subsequently\n\n");

	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
      	}
	/* read in coordinate locations and initial transformations */
	if((TRANS = fopen(filename,"r")) == NULL) {
	   fprintf(stderr,"error: file %s does not exist\n",filename);
	   exit(-1);
	}

	/* determine the number of domains specified */
	printf(" Reading domain descriptors...\n");
	ndomain=count_domain(TRANS);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	rewind(TRANS);
	if(getdomain(TRANS,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	for(i=0; i<ndomain; ++i) {
	  for(j=0; j<domain[i].nobj; ++j) {
	     if(domain[i].start[j].cid=='_') domain[i].start[j].cid=' ';
	     if( domain[i].start[j].in=='_') domain[i].start[j].in=' ';
	     if(  domain[i].end[j].cid=='_') domain[i].end[j].cid=' '; 
	     if(   domain[i].end[j].in=='_') domain[i].end[j].in=' ';
	  }
	} 

	bloc=(struct seqdat*)malloc((ndomain*2+2)*sizeof(struct seqdat));
	/* read in the alignment */
	printf(" Reading alignment...\n");
	rewind(TRANS);
	printf(" ");
	if(Agetbloc(TRANS,bloc,&nbloc)==-1) {
	  fprintf(stderr,"error: alignment file %s appears to be incomplete %s\n",filename); 
	  exit(-1); 
        } 
	bloclen=strlen(&bloc[1].seq[1]);
	if(nbloc!=(ndomain*2+1)) {
	     fprintf(stderr,"error: number of domains and alignment disagree\n");
	     exit(-1);
	}
	counter=(int*)malloc(bloclen*sizeof(int));
	all=(int*)malloc(bloclen*sizeof(int));
	for(i=0; i<ndomain; ++i) {
	   /* three state assignment */
	   if(threestate(&bloc[i+ndomain+2].seq[1],"HG","EB","TS-Ic")==-1) exit(-1);
	   if(smoothsec(&bloc[i+ndomain+2].seq[1],4,2)==-1) exit(-1);
	}

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
	fclose(TRANS);

	/* display the secondary structure summary */
	if(verbose) {
	  j=0;
	  printf(" The alignment looks like this...\n");
	  while(j<bloclen) {
	      for(l=0; l<nbloc; ++l) if(strncmp(bloc[l+1].id,"space",5)!=0) {
		printf(" %20s ",bloc[l+1].id);
	        k=0; 
		while((k+j)<bloclen && k<60) { 
		   printf("%c",bloc[l+1].seq[j+k+1]); 
		   k++;
	        } 
	        printf("\n");
	      }
              k=0; 
	      printf("      Reliablity      ");
	      while((k+j)<bloclen && k<60) { 
		printf("%1d",rel[j+k]); 
		k++; 
	      } 
	      printf("\n\n");
	      j+=k;
	  }
	}

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
		  domain[i].type[j],(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) {
	          fprintf(stderr,"Error in domain %s object %d \n",domain[i].id,j+1);
                    exit(-1);
               }

	       switch(domain[i].type[j]) {
	 	  case 1: printf(" all residues"); break; 
		  case 2: printf(" chain %c",domain[i].start[j].cid); break;
		  case 3: printf(" from %c %4d %c to %c %4d %c",
			 domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
			 domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
		} 
		printf("%4d CAs ",add); 
	        total+=add;
	        rewind(PDB);
	    }
	    domain[i].ncoords=total;
	    printf(" A total of %4d CAs in total\n",domain[i].ncoords);
	    /* disp(domain[i],stdout); */
	    printf(" Applying transformation... \n");
/*	    printmat(domain[i].R,domain[i].V,3,stdout);
	    printf("      ...to these coordinates.\n");  */
	    matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,PRECISION);  
	    fclose(PDB);
	}
        summary=(char*)malloc(bloclen*sizeof(char));

	printf(" Processing alignment...\n");
	for(i=0; i<bloclen; ++i) {
	   helix=strand=coil=0;
	   for(j=0; j<ndomain; ++j) {
	      /* change '_' to ' ' in brookhaven numbering */
	      if(domain[j].numb[i].cid=='_') domain[j].numb[i].cid=' ';
	      if(domain[j].numb[i].in=='_') domain[j].numb[i].in=' ';
	      if(bloc[j+ndomain+2].seq[i+1]=='H') helix++;
	      if(bloc[j+ndomain+2].seq[i+1]=='B') strand++;
	      if(bloc[j+ndomain+2].seq[i+1]=='c') coil++;
	   }
	   if(cons_sec==1) {
	     /* if we are after a consensus, then the majority secondary
	      *  structure will dominate */
	     if(helix>=(int)((float)ndomain/2)) summary[i]='H';	
             else if(strand>=(int)((float)ndomain/2)) summary[i]='B';
	     else summary[i]='c';
	   } else {
	     if(helix==ndomain) summary[i]='H';
	     else if(strand==ndomain) summary[i]='B';
	     else summary[i]='c';
	   }
	}
	printf(" Generating TSCORE file(s)...\n");
	for(i=1; i<ndomain; ++i) { /* Assume first is target */
	    sprintf(filename,"%s.talign",domain[i].id);
	    printf(" Domain %s to file %s\n",domain[i].id,filename);
/*	    if((OUT=fopen(filename,"r"))!=NULL) {
		fprintf(stderr,"error: file %s already exists\n",filename);
		exit(-1);
	    }
	    fclose(OUT); 
*/
	    OUT=fopen(filename,"w");
	    /* Preamble */
	    fprintf(OUT,"REMARK CASP(2) format alignment file generated from STAMP output file\n");
	    fprintf(OUT,"REMARK  R.B. Russell 1994\n");
	    fprintf(OUT,"REMARK \nREMARK \n");
	    fprintf(OUT,"REMARK  A total of %d equivalent residues were found\n",nrel);
	    fprintf(OUT,"REMARK \nREMARK \n");
	    j=0;
	    nres1=0; nres2=0;
	    while(j<bloclen) {
		if(rel[j]==1) {
		   start1=nres1; start2=nres2;
		   while(j<bloclen && rel[j]==1) { 
			j++;
			nres1++;
		        nres2++;
		   }
		   if(j>=bloclen) { j--; nres1--; nres2--; } 
		   strncpy(&id[0],domain[i].id,4); id[4]='\0';
		   chain=domain[i].id[4];
		   if(chain>='a' && chain <='z') chain = 'A' + chain - 'a';
		   if(chain<'A' || chain > 'Z') chain = '-';
		   /* format is - |TALIGN T0004 0    49    60 1CSP - 0    40    51 1.0 1| */
		   fprintf(OUT,"TALIGN %s 0 %4d %4d %s %c 0 %4d %4d 1.0 1\n",
			  domain[0].id,domain[0].numb[start1].n,domain[0].numb[nres1].n,
			  id,chain,domain[i].numb[start2].n,domain[i].numb[nres2].n);
		} else {
		   if(bloc[1].seq[j+1]!=' ') nres1++;
		   if(bloc[i+1].seq[j+1]!=' ') nres2++;
		}
	        j++;
	    }
	    fclose(OUT);
	}
	printf(" ...done.\n");
	exit(0);
}
int exit_error()
{
	 fprintf(stderr,"format:  stamp2casp -f <file> -c <char> -w <win> -t <thresh> -a -aligned -cons\n");
 	 fprintf(stderr,"          -aligned ==> consider all aligned positions\n");
	 fprintf(stderr,"          -a ==> generate average structure file\n");
	 fprintf(stderr,"          -cons ==> generate consensus secondary structure summary\n");	
	 fprintf(stderr,"           (default is all or nothing, e.g. all helix all strand or coil)\n");
	 fprintf(stderr,"          -V ==> verbose output\n");
	 exit(-1);
}
