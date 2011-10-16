#include <stdio.h>
#include <stdlib.h>
#include <stamp.h>

/* Given a file of transformations and a domain ID, this program
 *  centres all transformations such that all transformations are
 *  expressed as a transformation of all domains onto the specified
 *  ID */

void exit_error();

int main(int argc, char *argv[]) {
	
	
	int i,j,k;
	int ndomain,gotdomain;
	int got_file;
	int got_id;
	int which;
	int gottrans;

	int *indx;

	float sign;

	float *negvec;

	float **invmat;
	float **R,**RI;

	char c;

	char *id;
	char *buff;
	char *infile;
    char *stampdir = AM_STAMPDIR;

	FILE *IN;

	struct domain_loc *domain;

	id=(char*)malloc(100*sizeof(char));
	buff=(char*)malloc(1000*sizeof(char));
	infile=(char*)malloc(1000*sizeof(char));
	indx=(int*)malloc(100*sizeof(int));
/* SMJS Was sizeof(float) */
	invmat=(float**)malloc(3*sizeof(float *));
/* SMJS Was sizeof(float) */
	R=(float**)malloc(4*sizeof(float *));
/* SMJS Was sizeof(float) */
	RI=(float**)malloc(4*sizeof(float *));
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
	      if((i+1)>=argc) exit_error();
	      strcpy(id,argv[i+1]);
	      got_id=1;
	      i++;
	   } else exit_error();
	}

	if(!got_file) {
	  fprintf(stderr,"error: must specify file name\n");
	  exit(-1);
	}

    if(getenv("STAMPDIR")!=NULL) {
      /* Allow environment variable to override config setting */
      stampdir=getenv("STAMPDIR");
    }

	ndomain=count_domain(IN);
	rewind(IN);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	if(getdomain(IN,domain,&ndomain,ndomain,&gottrans,stampdir,0,stdout)==-1) exit(-1);
	if(!gottrans) {
	  fprintf(stderr,"error: file does not contain transformations\n");
	  exit(-1);
	}
	rewind(IN);

	if(!got_id) {
	  fprintf(stderr,"error: must specify identifier from domain file\n");
	  fprintf(stderr," The identifiers in the file are:\n");
	  for(i=0; i<ndomain; ++i) fprintf(stderr,"          %s\n",domain[i].id);
	  exit(-1);
	}
	
	/* find which domain ID to use */
	for(i=0; i<ndomain; ++i) {
	   if(strcmp(id,domain[i].id)==0) {
	      which=i;
/*	      printf("which = %4d\n",which); */
	      break;
	   }
	}

	/* put some comments in to tell which domain is centered */
	printf("%%\n%%\n%% PICKFRAME R.B. Russell, 1995\n");
	printf("%%  The input file was %s\n",infile);
	printf("%%  the domains in this file have been centred on domain %s\n",domain[which].id);
	printf("%%  The original comments are included.\n");
	printf("%%\n%%\n");
	
	/* put some of the comments back into the output */
	while(fgets(buff,599,IN)!=NULL && (buff[0]=='%' || buff[0]=='#')) 
	   printf("%s",buff);

	/* get the inverse */
	for(i=0; i<3; ++i) {
	  for(j=0; j<3; ++j) {
	     R[i+1][j+1]=domain[which].R[i][j];
	  }
	}
/*	printf("The matrix:\n");
	printmat(domain[which].R,domain[which].V,3,stdout);
*/
	matinv(R,RI,&sign,indx);

	for(i=0; i<3; ++i) { 
	   for(j=0; j<3; ++j) {
	      invmat[i][j]=RI[i+1][j+1];
	   }
	}

	/* store the negative of domain[which].V in negvec */
	for(i=0; i<3; ++i) 
	  negvec[i]=-1*domain[which].V[i];

/*	printf("The inverse:\n");
	printmat(invmat,domain[which].V,3,stdout);
*/
	/* Now apply the reverse transformation to each of the matrices 
	 *
	 * Think, if X is the 
	 *    X' = RxX + Vx
	 * then,
	 *    X = Rx^-1(X'-Vx)
	 *
	 * where X and X' are the old and new coordinates.  Now we
	 * want to apply these to all other transformations, Y, where:
	 *
	 *  Y' = RyY + Vy
	 *
	 * Therefor,
	 * 
	 *  Y'' = Rx^-1(Y'-Vx)
	 *      = Rx^-1(RyY + Vy - Vx)
	 *
	 * We must thus first subtract the translation, then apply the 
	 *  inverse matrix, giving the new transformations as:
	 *
	 *  Ry' = Rx^-1 Ry
	 *  Vy' = Rx^-1(Vy - Vx) */

	for(i=0; i<ndomain; ++i) {
	   /* first apply the inverse to the translation */
/*	   printf("subtracting:\n");
	   printf("Vector from %s from vector from %s\n",
	      domain[which].id,domain[i].id); */
	   for(j=0; j<3; ++j) {
/*	     printf("%8.5f +  %8.5f\n",domain[i].V[j],negvec[j]); */
	     domain[i].V[j]=domain[i].V[j]+negvec[j];
	   }
	   matvecprod(invmat,domain[i].V,domain[i].V,stdout);
	   /* now apply the inverse to the old matrix */
	   matprod(domain[i].R,invmat,domain[i].R,stdout);
	   printdomain(stdout,domain[i],1);
	   /* now put some more of the comments back in 
	    *  first read till the closing brace */
	   while((c=getc(IN))!=(char)EOF && c!='}');
	   if(c==(char)EOF) break;
	   /* read till end of line */
	   while((c=getc(IN))!=(char)EOF && c!='\n');
	   if(c==(char)EOF) break;
	   while(fgets(buff,599,IN)!=NULL && (buff[0]=='%' || buff[0]=='#'))
	     printf("%s",buff);

	}

	exit(0);
}

void exit_error()
{
	fprintf(stderr,"format: pickframe -f <transformation file> -i <id to centre on> \n");
	exit(-1);
}
