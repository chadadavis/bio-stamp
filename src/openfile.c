
/* This checks for compression, etc and runs open/popen as appropriate */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FILE *openfile(char *filename, char *type) {

	FILE *handle;
	char command[1000];

	if(strcmp(type,"r")==0) { /* A read file */
		if(strcmp(&filename[strlen(filename)-2],".Z")==0) { /* UNIX compression */
/*		     sprintf(command,"zcat %s 2> /dev/null",filename);  */
		     sprintf(command,"zcat %s",filename);  
		     handle=popen(command,"r");
		} else if(strcmp(&filename[strlen(filename)-3],".gz")==0) { /* gzipped  */
/*		     sprintf(command,"gunzip -c %s 2> /dev/null",filename);  */
		     sprintf(command,"gunzip -c %s",filename); 
		     handle=popen(command,"r");
		} else { /* presume no compression */
		    handle=fopen(filename,"r");
		}
	} else {
	  	    handle=fopen(filename,"w");
	}
	return handle;
}

	
