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

#include <poststamp.h>

#define MAX_SEQ_LEN 2000
#define PRECISION 1000
#define E1 3.8
#define E2 3.8

/* poststamp: a program to perform a secondary analysis on stamp output 
 *  see ~/distribute/stamp/doc/stamp.doc for details 
 *
 * Modification 21 March 1996
 *  Now reads more things from the command line (- options etc.)
 * Also will now provide STAMP like equivalences from normal
 *  alignments (you have to give a domain/transformation file using -d)
 */
main(int argc, char *argv[]) {

	int i,j,k,l,test;
	int ndomain,total,add,nbloc;
	int nstamp,nseq,nstamppos;
	int gottrans;
	int bloclen;
	int aligned;
	int domfile;
	int atstart,atend;
	int *pointer;
	int *counter;
	int *all;

	char c,filename[200],filename2[200],outfile[200];
	char *env;

	float min_Pij,Dij,Cij;
	float const1,const2;
	float **Pij;

	FILE *DOM,*ALIGN,*PDB,*OUT;

	struct domain_loc *domain;
	struct seqdat *bloc;
	struct stampdat *stamp;

	printf("POSTSTAMP, R.B. Russell 1995\n");

	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
      	}
	min_Pij=0.5;
	aligned=0;
	filename2[0]='\0';
	filename[0]='\0';

	if(argc<3) exit_error();
	for(i=1; i<argc; ++i) {
           if(argv[i][0]!='-') exit_error();
           if(strcmp(&argv[i][1],"f")==0) {
              if((i+1)>=argc) exit_error();
              strcpy(&filename[0],argv[i+1]);
              i++;
           } else if(strcmp(&argv[i][1],"d")==0) {
             if((i+1)>=argc) exit_error();
             strcpy(&filename2[0],argv[i+1]);
             i++;
           } else if(strcmp(&argv[i][1],"min")==0) {
             if((i+1)>=argc) exit_error();
             sscanf(argv[i+1],"%f",&min_Pij);
             i++;
           } else if(strcmp(&argv[i][1],"aligned")==0) {
	     aligned=1;
           } else {
             exit_error();
           }
        }
	const1=-2*E1*E1;
	const2=-2*E2*E2;
	if(min_Pij<0.0 || min_Pij>1.0) {
	   fprintf(stderr,"error: minimum Pij value must be between 0 and 1\n");
	   exit(-1);
	}

	if(strlen(filename)==0) {
	   fprintf(stderr,"error: you must specify a file name\n");
	   exit(-1);
	}
	sprintf(&outfile[0],"%s.post",filename);
	if((OUT=fopen(outfile,"w"))==NULL) {
	  fprintf(stderr,"error opening file %s \n",outfile);
	  exit(-1);
	}
	printf(" New output will be in file %s\n",filename);
	printf(" E1 = %7.3f, E2 = %7.3f\n",E1,E2); 
	printf(" Minimum Pij set to %5.3f\n",min_Pij);

	/* read in coordinate locations and initial transformations */
	printf(" Reading domain descriptors/transformations from the file ");
	if(filename2[0]=='\0') {
	  if((DOM = fopen(filename,"r")) == NULL) {
	   fprintf(stderr,"error: file %s does not exist\n",filename);
	   exit(-1);
	  }
	  printf("%s\n",filename);
	} else {
	   if((DOM = fopen(filename2,"r")) == NULL) {
           fprintf(stderr,"error: file %s does not exist\n",filename2);
           exit(-1);
          }
	  printf("%s\n",filename2);
	}
	/* determine the number of domains specified */
	ndomain=count_domain(DOM);
	if(ndomain==0) {
	   fprintf(stderr,"Error no domain descriptors found\n");
	   exit(-1);
	}
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	rewind(DOM);
	if(getdomain(DOM,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	pointer=(int*)malloc(ndomain*sizeof(int));
	for(i=0; i<ndomain; ++i) {
	  for(j=0; j<domain[i].nobj; ++j) {
	     if(domain[i].start[j].cid=='_') domain[i].start[j].cid=' ';
	     if( domain[i].start[j].in=='_') domain[i].start[j].in=' ';
	     if(  domain[i].end[j].cid=='_') domain[i].end[j].cid=' '; 
	     if(   domain[i].end[j].in=='_') domain[i].end[j].in=' ';
	  }
	}
	fclose(DOM);

	if((ALIGN = fopen(filename,"r")) == NULL) {
           fprintf(stderr,"error: file %s does not exist\n",filename);
           exit(-1);
        }


	/* read in the alignment */
	printf(" Reading alignment...\n");
	rewind(ALIGN);
	bloc=(struct seqdat*)malloc((ndomain*2+2)*sizeof(struct seqdat));
	printf(" ");
	Agetbloc(ALIGN,bloc,&nbloc);
	bloclen=strlen(&bloc[1].seq[1]);
	if(nbloc!=(ndomain*2+1)) {
	   fprintf(stderr,"error: number of domains and alignment disagree\n");
	   exit(-1);
	}
	Pij=(float**)malloc((ndomain*(ndomain-1)/2)*sizeof(float*));
	for(i=0; i<(ndomain*(ndomain-1)/2); ++i) 
	  Pij[i]=(float*)malloc(bloclen*sizeof(float));
	counter=(int*)malloc(bloclen*sizeof(int));
	all=(int*)malloc(bloclen*sizeof(int));

	/* read in the STAMP output */
	rewind(ALIGN);
	stamp=(struct stampdat*)malloc(100*sizeof(struct stampdat));
	if(getstampdat(stamp,ALIGN,&nstamp,&nseq,&nstamppos,MAX_SEQ_LEN)==-1) exit(-1);
	fclose(ALIGN);

	printf(" Reading coordinates...\n");
	for(i=0; i<ndomain; ++i) {
	   /* output the domain descriptors again */
	   printdomain(OUT,domain[i],1);
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
	    printf("=> %4d CAs in total\n",domain[i].ncoords);
	    printf(" Transforming coordinates...\n");
	    matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,PRECISION);  
	    fclose(PDB);
	}
	
	/* output ">" descriptors */
	fprintf(OUT,"\n\n");
	for(i=0; i<nbloc; ++i) 
	  fprintf(OUT,">%s %s",bloc[i+1].id,bloc[i+1].title);
	/* output "#" descriptors */
	for(i=0; i<nstamp; ++i) 
	  fprintf(OUT,"#%c %s",stamp[i].what,stamp[i].title);
	fprintf(OUT,"#B 1 if all pairiwse Pij greater than %5.3f\n",min_Pij);
	fprintf(OUT,"#R total number of pairwise comparisons having Pij greater than %5.3f (out of %4d)\n",
	   min_Pij,ndomain*(ndomain-1)/2);
	/* Now proceed through the alignment calculating all pairwise Pij values */
	for(i=0; i<ndomain; ++i) pointer[i]=0;
	fprintf(OUT,"*\n");
	for(i=0; i<bloclen; ++i) {
	  counter[i]=0;
	  for(j=0; j<(ndomain*2+1); ++j) {
	     fprintf(OUT,"%c",bloc[j+1].seq[i+1]);
	  }
	  all[i]=1; l=0;
/*	  fprintf(OUT," %3d",i+1);  */
	  for(j=0; j<ndomain; ++j) {
	    for(k=j+1; k<ndomain; ++k) {
	       if(bloc[j+1].seq[i+1]!=' ' && bloc[k+1].seq[i+1]!=' ') {
		  if(pointer[j]==0 || pointer[k]==0) atstart=1;
		  if(pointer[j]>=(domain[j].ncoords-1) || pointer[k]>=(domain[k].ncoords-1)) atend=1;
		  Pij[l][i]=rossmann(&domain[j].coords[pointer[j]],&domain[k].coords[pointer[k]],
		   atstart,atend,const1,const2,&Dij,&Cij,PRECISION);
		} else Pij[l][i]=0.0;
	  	counter[i]+=(Pij[l][i]>=min_Pij);
		all[i]*=(Pij[l][i]>=min_Pij);
/*		fprintf(OUT," %4.2f",Pij[l][i]);     */
	        l++;
	     }
	  }
	  if(aligned || (nstamp>0 && (stamp[0].n[i]>-0.001))) {
	    fprintf(OUT,"  ");
	    for(j=0; j<nstamp; ++j) 
	       if(stamp[j].what=='T') fprintf(OUT,"%1.0f ",stamp[j].n[i]);
	       else fprintf(OUT,"%10.5f ",stamp[j].n[i]);
	    fprintf(OUT," %1d",all[i]); 
	    fprintf(OUT," %3d",counter[i]); 
	  }
	  for(j=0; j<ndomain; ++j) {
	     pointer[j]+=(bloc[j+1].seq[i+1]!=' ');
/*	     fprintf(OUT," %3d",pointer[j]);
	     fprintf(OUT," %c",bloc[j+1].seq[i+1]);  */
	  }
	  fprintf(OUT,"\n");
	}
	fprintf(OUT,"*\n");
	printf(" ...done.\n");
	free(pointer);
	exit(0);
}
int exit_error() {
    fprintf(stderr,"format: poststamp -f <alignment file> -min <minimum Pij> -aligned\n");
    fprintf(stderr,"        -d <domains file> [use only if using a non-STAMP alignment]\n");
    fprintf(stderr,"        -aligned will consider all aligned positions\n");
    exit(-1);
}
