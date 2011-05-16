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

/* MERGETRANS
 * Given two file of transformations and a domain ID, this program
 *  centres all transformations such that all transformations are
 *  expressed as a transformation of all domains onto the first
 *  ID common to both files 
 *
 * RBR 20 June 1996 
 *
 * Modification 9 July 1996 - no ignores duplicates in the second file (i.e. 
 *  those not involved in indexing */

void exit_error();

/* TPW : Added int as main() return type; removed unused variables */
int main(int argc, char *argv[]) {
	
	
	int i,j;
	int ndomain;
	int ndomain2;
	int got_file,got_file2;
	int got_id;
	int which,which2;
	int gottrans,gottrans2;
	int ignore;

	int *indx;

	float sign;

	float *negvec;

	float **invmat;
	float **R,**RI;


	char *id;
	char *buff;
	char infile[200],infile2[200];
	char *env;
	
	FILE *IN,*IN2;

	struct domain_loc *domain,*domain2;

	id=(char*)malloc(100*sizeof(char));
	buff=(char*)malloc(1000*sizeof(char));
	indx=(int*)malloc(100*sizeof(int));
/* Was sizeof(float) */
	invmat=(float**)malloc(3*sizeof(float *));
/* Was sizeof(float) */
	R=(float**)malloc(4*sizeof(float *));
/* Was sizeof(float) */
	RI=(float**)malloc(4*sizeof(float *));
	negvec=(float*)malloc(3*sizeof(float));
	for(i=0; i<4; ++i) { 
	   R[i]=(float*)malloc(4*sizeof(float));
	   RI[i]=(float*)malloc(4*sizeof(float));
	}
	for(i=0; i<3; ++i) 
	  invmat[i]=(float*)malloc(3*sizeof(float));


	if(argc<3) exit_error();

	got_file=got_file2=got_id=0;
	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(strcmp(&argv[i][1],"f1")==0) { 
	      if((i+1)>=argc) exit_error();
	      if((IN=fopen(argv[i+1],"r"))==NULL) {
		 fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
		 exit(-1);
	      }
	      got_file=1;
	      strcpy(infile,argv[i+1]);
	      i++;
	   } else if(strcmp(&argv[i][1],"f2")==0) {
              if((i+1)>=argc) exit_error();
              if((IN2=fopen(argv[i+1],"r"))==NULL) {
                 fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
                 exit(-1);
              }
              got_file2=1;
              strcpy(infile2,argv[i+1]);
              i++;
	   } else if(argv[i][1]=='i') {
	      if((i+1)>=argc) exit_error();
	      strcpy(id,argv[i+1]);
	      got_id=1;
	      i++;
	   } else exit_error();
	}

	if(!got_file || !got_file2) {
	  fprintf(stderr,"One or both files not found\n");
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
/*	if(!gottrans) {
	  fprintf(stderr,"error: file %s does not contain transformations\n",infile);
	  exit(-1);
	} 
*/
	fclose(IN);

	ndomain2=count_domain(IN2);
        rewind(IN2);
        domain2=(struct domain_loc*)malloc(ndomain2*sizeof(struct domain_loc));
        if(getdomain(IN2,domain2,&ndomain2,ndomain2,&gottrans2,env,0,stdout)==-1) exit(-1);
/*        if(!gottrans) {
          fprintf(stderr,"error: file %s does not contain transformations\n",infile2);
          exit(-1);
        }
*/
	fclose(IN2);


	if(got_id) { /* User has specified ID */
  	  /* find which domain ID to use */
	  which=-1;
	  for(i=0; i<ndomain; ++i) {
	     if(strcmp(id,domain[i].id)==0) {
	        which=i;
/*	        printf("which = %4d\n",which); */
	        break;
	     }
	  }
	  if(which==-1) {
	     fprintf(stderr," Error: ID %s not found in file %s\n",id,infile);
	     fprintf(stderr," The identifiers in the file are:\n"); 
	     for(i=0; i<ndomain; ++i) fprintf(stderr,"          %s\n",domain[i].id);
	     exit(-1);
	  }
	  which2=-1;
          for(i=0; i<ndomain2; ++i) {
             if(strcmp(id,domain2[i].id)==0) {
                which2=i;
/*              printf("which = %4d\n",which2); */
                break;
             }
          }
          if(which2==-1) {
             fprintf(stderr," Error: ID %s not found in file\n",infile2);
             fprintf(stderr," The identifiers in the file are:\n"); 
             for(i=0; i<ndomain2; ++i) fprintf(stderr,"          %s\n",domain2[i].id);
             exit(-1);
          }
	} else if(ndomain2>0) {
	  /* Just look for the first ID in common */
	  which = which2 = -1;
	  for(i=0; i<ndomain; ++i) {
		for(j=0; j<ndomain2; ++j) {
		    if(strcmp(domain[i].id,domain2[j].id)==0) {
			which = i; which2 = j;
			printf("%% Using domain %s to centre (%d in %s; %d in %s)\n",
			   domain[i].id,i+1,infile,j+1,infile2);
			break;
		    }
	         }
		 if(which !=-1) break;
	   }
	} else {
	   fprintf(stderr,"Warning no domains found in file %s - ignored in output\n",infile2);
	   which = 0;
	}
	if(which==-1 || which2==-1) {
	   fprintf(stderr,"Warning couldn't index the two files.  File %s will be ignored\n",infile2);
	   ndomain2 = 0;
	   which = 0;
	}
	if(ndomain==0) {
	    fprintf(stderr,"Warning: ignoring file %s in output\n",infile);
	    which2 =0;
	}
	if(ndomain2==0) {
	     fprintf(stderr,"Warning: ignoring file %s in ouptut\n",infile2);
	     which = 0;
	}



	/* put some comments in to tell which domain is centered */
	printf("%%\n%%\n%% MERGETRANS R.B. Russell, 1996\n");
	printf("%%  The input files were %s and %s\n",infile,infile2);
	printf("%%  the domains in this file have been centred on domain %s\n",domain[which].id);
	printf("%%  The original comments are have been removed.\n");
	printf("%%\n%%\n");
	/* N.B. Comments are lost */


	if(ndomain>0) {
	  /* Do each file one at a time, then just ignore the second copy (i.e. which2) */
	  /* First domains first */
	  /* get the inverse */
	  for(i=0; i<3; ++i) {
	    for(j=0; j<3; ++j) {
	     R[i+1][j+1]=domain[which].R[i][j];
	    }
	  }
	  matinv(R,RI,sign,indx);
	  for(i=0; i<3; ++i) { 
	   for(j=0; j<3; ++j) {
	      invmat[i][j]=RI[i+1][j+1];
	   }
	  }
	  /* store the negative of domain[which].V in negvec */
	  for(i=0; i<3; ++i) 
	    negvec[i]=-1*domain[which].V[i];

	  /* See PICKFRAME for an explanation */
	  for(i=0; i<ndomain; ++i) {
	   /* first apply the inverse to the translation */
	   for(j=0; j<3; ++j) {
	     domain[i].V[j]=domain[i].V[j]+negvec[j];
	   }
	   matvecprod(invmat,domain[i].V,domain[i].V,stdout);
	   /* now apply the inverse to the old matrix */
	   matprod(domain[i].R,invmat,domain[i].R,stdout);
	   printdomain(stdout,domain[i],1);
	  }
	} else {
	   printf("%% First domain file ignore as it was empty\n");
	}

	/* Second set of domains */
	if(ndomain2>0) {
	  for(i=0; i<3; ++i) {
            for(j=0; j<3; ++j) {
               R[i+1][j+1]=domain2[which2].R[i][j];
            }
          }
          matinv(R,RI,sign,indx);
          for(i=0; i<3; ++i) {
             for(j=0; j<3; ++j) {
                invmat[i][j]=RI[i+1][j+1];
             }
          }
          /* store the negative of domain[which].V in negvec */
          for(i=0; i<3; ++i)
            negvec[i]=-1*domain2[which2].V[i];

          /* See PICKFRAME for an explanation */
          for(i=0; i<ndomain2; ++i) {
	   ignore=0;
           for(j=0; j<ndomain; ++j) if(strcmp(domain[j].id,domain2[i].id)==0) ignore=1;
	   if(!ignore) { 
              /* first apply the inverse to the translation */
              for(j=0; j<3; ++j) {
                domain2[i].V[j]=domain2[i].V[j]+negvec[j];
              }
              matvecprod(invmat,domain2[i].V,domain2[i].V,stdout);
              /* now apply the inverse to the old matrix */
              matprod(domain2[i].R,invmat,domain2[i].R,stdout);
              printdomain(stdout,domain2[i],1);
           } else {
	      printf("%% Second copy of domain %s ignored\n",domain2[i].id);
	   }
	 }
	} else {
	  printf("%% Second domain file ignored as it was empty\n");
	}
    return 0;
}

void exit_error()
{
	fprintf(stderr,"format: mergetrans -f1 <trans file> -f2 <trans file> \n");
	fprintf(stderr,"                  [ -i <id to centre on> ]\n");
	exit(-1);
}
