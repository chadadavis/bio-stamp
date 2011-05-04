#include <stdio.h>
#include <dstamp.h>

/* This routine works in conjunction with Agetbloc.
 *
 *  Given a blockfile containing a list of '>' which specifiy 
 *   character sequences, and a list of "#" which specify
 *   numeric data, this reads them in and returns them as
 *   a string of floats called 'stamp' */

int getstampdat(struct stampdat *stamp, FILE *IN, int *nstamp, int *nseq, int *npos, int maxpos) {
	
	int i,j,k,n;
	float value;

	char c;
	char *tmp;

	tmp=(char*)malloc(1000*sizeof(char));

	(*nstamp)=0;
	(*nseq)=0;

	/* read in the descriptors */
	while((c=getc(IN))!=(char)EOF) {
	   (*nseq)+=(c=='>');
	   if(c=='#') {
	      fgets(tmp,100,IN);
	      stamp[(*nstamp)].what=tmp[0];
	      stamp[(*nstamp)].title=(char*)malloc(100*sizeof(char));
	      stamp[(*nstamp)].n=(float*)malloc(maxpos*sizeof(float));
	      strncpy(stamp[(*nstamp)].title,&tmp[1],99);
	      stamp[(*nstamp)].title[99]='\0';
	      (*nstamp)++;
	   }
	}
	rewind(IN);

	/* now find the first '*' and which column it is in */
	n=0;
	while((c=getc(IN))!=(char)EOF && c!='*') n=(n+(c!='\n'))*(c!='\n');
	if(c==(char)EOF) return -1;
	while((c=getc(IN))!=(char)EOF && c!='\n'); /* read to the end of the line */
	if(c==(char)EOF) return -1;

	/* now read in the file line by line */
	(*npos)=0;
	while(tmp[n]!='*') {
	  fgets(tmp,900,IN);
	  for(i=0; i<strlen(tmp); ++i) if(tmp[i]=='\n') tmp[i]='\0';
	  if(tmp[n]=='*') break;
	  /* Lets allow for missing values, and set them to zero.
	   * In other words, when the length of the string read in is
	   * equal to the number of '>' characters + n, ignore the line */
	  if(strlen(tmp)<=((*nseq)+n+1)) {
	     for(i=0; i<(*nstamp); ++i) 
		stamp[i].n[(*npos)]=-1.0;
	  } else { 
	     j=((*nseq)+n);
	     for(i=0; i<(*nstamp); ++i) {
		while(tmp[j]==' ') ++j;  /* move to the next space in the string */
		sscanf(&tmp[j],"%f",&value);
		stamp[i].n[(*npos)]=value; /* read in the next float */
		while(tmp[j]!=' ' && tmp[j]!='\0') ++j; /* move to the next space */
	     }
	  }
	  (*npos)++;
	}
	free(tmp);
	return 0;
}
