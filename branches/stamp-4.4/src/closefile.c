/* Close for openfile() */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void closefile(FILE *handle, char *filename) {


	char b[1000];
	if((strcmp(&filename[strlen(filename)-2],".Z")==0) || 
           (strcmp(&filename[strlen(filename)-3],".gz")==0)) { 

	    /* This fgets prevents "Broken pipe" messages from appearing when running under Solaris,
	     *  and silly newlines to stderr when running under OSF - something about a wrong implementation
	     *  of pclose()? */
	    while(fgets(b,999,handle)!=NULL);

	    pclose(handle);
	} else { 
	    fclose(handle);
	}
}
