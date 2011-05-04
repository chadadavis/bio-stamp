#include <stdio.h>
#include "stamp.h"
#define MAXL 200    /* max length for input line */

/* get_dssp_sum:  return the secondary structure summary
	 from a dssp file */

int get_dssp_sum(FILE *DSSP, struct brookn begin, struct brookn end, int type,
	char *sec, int REVERSE,	int maxlen, int *count, FILE *OUT) {

    char start,datstart;
    char cid,in;
    char css;

    char number[6];
    char *line;

    int n,i;

    datstart=start=0;
    *count=0;

    line = (char *) malloc(sizeof(char)*MAXL);

    while(fgets(line,MAXL,DSSP) != NULL) {
	strncpy(&number[0],&line[5],5); number[5]='\0';
	sscanf(&number[0],"%d",&n);
	cid=line[11];
	in=line[10];
	if(datstart && (line[13]!='!') &&
	   ( (datstart && type==1) ||
	     (cid==begin.cid && type==2) ||
	     (cid==begin.cid && in==begin.in && n==begin.n && type==3) 
	    )) {
	      start=1;
	}
	if(datstart && start && type==2 && (line[13]!='!') && begin.cid!=cid) break;
	if(start && (line[13]!='!')) { /* skips ill placed chain breaks */
	    if((*count)==maxlen-1) {
	       fprintf(stderr,"error: maximum sequence length surpassed when attempting to get DSSP summary\n");
	       return -1;
	    }
	    sec[(*count)] = line[16];
	    if(sec[(*count)]==' ') sec[(*count)]='-';
	    (*count)++;
	}
	if(start && line[13]!='!' && cid==end.cid && in==end.in && n==end.n && type==3) break;
	if(line[2] == '#') datstart=1;
    }
    sec[(*count)]='\0';
    if(REVERSE) {
       for(i=0; i<(int)(strlen(sec)/2); ++i) {
	  css=sec[i];
	  sec[i]=sec[strlen(sec)-i-1];
	  sec[strlen(sec)-i-1]=css;
	}
    }
    free(line);

    return 0;
}
