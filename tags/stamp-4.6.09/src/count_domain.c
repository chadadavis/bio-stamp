#include <stdio.h>
#include <stdlib.h>


/* Counts the number of domain descriptors in a file ignoring comments */
int count_domain(FILE *IN) {

	int N,col,comment;
	char c;

	N=0;
	comment=0;
	col=0;
	while((c=getc(IN))!=(char)EOF) {
	  if((col==0) && (c=='#' || c=='%')) {
	     comment=1;
	  }
	  if(col!=0 && !comment && c=='{') 
	     N++;
	  col++;
	  if(c=='\n') { 
	     col=0; 
	     comment=0; 
	  }
	}
	return N;
}
