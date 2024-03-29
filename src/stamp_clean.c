#include <stdio.h>
#include <stdlib.h>

#include <stamp.h>
#include <stamprel.h>

/* Clean up a block file by attaching lone isolated
 *  residues to the nearest continuous segment */

int main(int argc, char *argv[]) {

	int i,j;
	int nbloc,minlen,bloclen;
	int nstamp,nstamppos,nstampseq;
	int ndomain,gottrans;
        float minfrac;

	char c;
    char *stampdir = AM_STAMPDIR;

	FILE *BLOC;

	struct seqdat *bloc;
	struct stampdat *stamp;
	struct domain_loc *domain;

	if(argc!=4) {
	  printf("format: stamp_clean (block file) (minimum continuous segment length) (min frac) > (output file)\n");
	  exit(-1);
	}

    if(getenv("STAMPDIR")!=NULL) {
      /* Allow environment variable to override config setting */
      stampdir=getenv("STAMPDIR");
    }


	sscanf(argv[2],"%d",&minlen);
	sscanf(argv[3],"%f",&minfrac);


	if((BLOC=fopen(argv[1],"r"))==NULL) {
		fprintf(stderr,"error opening file %s\n",argv[1]);
		exit(-1);
	}
	printf("%% ALIGN_CLEAN, R.B. Russell, 1995\n%% Searching for domain descriptors...\n");
	ndomain=count_domain(BLOC);
	rewind(BLOC);
	if(ndomain!=0) {
	  domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	  if(getdomain(BLOC,domain,&ndomain,ndomain,&gottrans,stampdir,0,stdout)==-1) exit(-1);
/*	  int getdomain(FILE *IN, struct domain_loc *domains, int *ndomain, int maxdomain, int *gottrans, char *env, int DSSP, FILE *OUTPUT); */

	  printf("%%   %d domain descriptions read in\n",ndomain);
	  for(i=0; i<ndomain; ++i) printdomain(stdout,domain[i],gottrans);
	} else {
	  printf("%%   no domain information found\n");
	}
	printf("%% Reading block file...\n");
	nbloc=0;
	while((c=getc(BLOC))!=(char)EOF) nbloc+=(c=='>');
	rewind(BLOC);
	bloc=(struct seqdat*)malloc((nbloc+1)*sizeof(struct seqdat));
	if(Agetbloc(BLOC,bloc,&nbloc)==-1) exit(-1);
	rewind(BLOC);
	bloclen=strlen(&bloc[1].seq[1]);
	printf("%% Searching for STAMP data...\n");
	nstamp=0;
	while((c=getc(BLOC))!=(char)EOF) nstamp+=(c=='#');
	rewind(BLOC);
	if(nstamp>0) {
	   stamp=(struct stampdat*)malloc(nstamp*sizeof(struct stampdat));
	   if(getstampdat(stamp,BLOC,&nstamp,&nstampseq,&nstamppos,bloclen)==-1) exit(-1);
	   if(nstamppos!=bloclen) {
	      fprintf(stderr,"error: STAMP and sequence data disagree\n");
	      exit(-1);
	   }
	   printf("%%   %d STAMP fields found: ",nstamp);
	   for(i=0; i<nstamp; ++i) printf("%c ",stamp[i].what);
	   printf("\n");
	} else {
	  printf("%%   no STAMP data found\n");
	}
	fclose(BLOC);
	printf("%% Block file contains %d sequences; the alignment length is %d\n",
	       nbloc,strlen(&bloc[1].seq[1]));
	printf("%% Cleaning up allowing continuous segments of %d or greater...\n",minlen);
	if(nstamp==0) clean_block(bloc,nbloc,minlen);
	else stamp_clean_block(bloc,nbloc,minlen,minfrac,stamp,nstamp);
	printf("%%  Cleaning done.\n");
	bloclen=strlen(&bloc[1].seq[1]);
	printf("%% The final alignment length is %d\n",bloclen);
	printf("%% The alignment:\n");
	for(i=0; i<nbloc; ++i)  {
           /* Fix the stupid title problem with spaces */
           j=strlen(bloc[i+1].title)-1;
           while(((bloc[i+1].title[j]==' ') || (bloc[i+1].title[j]=='\n')) && (j>=0)) { 
              bloc[i+1].title[j]='\0'; 
              j--; 
           }
	   printf(">%s %s\n",bloc[i+1].id,bloc[i+1].title);
	}
	for(i=0; i<nstamp; ++i) {
	   printf("#%c %s\n",stamp[i].what,stamp[i].title);
	}
	printf("*\n");
	for(i=0; i<bloclen; ++i) {
	   for(j=0; j<nbloc; ++j)  {
	      printf("%c",bloc[j+1].seq[i+1]);
/*              if(bloc[j+1].seq[i+1] == '\n') {
                 fprintf(stderr,"Error: newline character found in output - exiting\n");
                 exit(-1); 
              }
*/

           }
	   if(nstamp>0 && stamp[0].n[i]>-0.001) {
	     printf(" ");
	     for(j=0; j<nstamp; ++j) {
	      if(stamp[j].what=='T') printf("%1.0f ",stamp[j].n[i]);
	      else printf("%10.5f ",stamp[j].n[i]);
             }
	   }
	   printf("\n");
	}
	printf("*\n");
	exit(0);
}
