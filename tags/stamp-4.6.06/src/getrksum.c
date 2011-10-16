#include <stdio.h>

char *getrksum(FILE *IN) { 
/* returns the summary from a DEFINE file.  The DEFINE output is from a modified
 *  version of the Richards and Kundrot program.  Refer to RBR for details */

	char *buff;
	char *sec;
	int start,len;

	buff=(char*)malloc(1000*sizeof(char));
	sec=(char*)malloc(sizeof(char));
	start=len=0;

	while(fgets(buff,100,IN)!=NULL) {
	   if(start) {
	     if(buff[0]=='*') 
		break;
	     if(buff[20]!='-') 
		sec[len]='H';
	     else if(buff[22]!='-') 
		sec[len]='3';
	     else if(buff[26]!='-') 
		sec[len]='B';
	     else if(buff[24]!='-') 
		sec[len]='E';
	     else  
		sec[len]='-';
	     len++;
	     sec=(char*)realloc(sec,(len+1)*sizeof(char));
	    }
	    if(!start && buff[0]=='*') 
	       start=1;
	}
	sec[len]='\0';
	free(buff);
	return sec;
}
