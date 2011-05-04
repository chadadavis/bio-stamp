#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This routine removes spaces from a string and returns the result. */

void rmsp(char *c) {

	int i,j;
	char *temp;

	temp=(char*)malloc(strlen(c)*sizeof(char));
	i=0; j=0;
	while(c[i] !='\0') {
	   if(c[i]!=' ') temp[j++]=c[i];
	   i++;
	} /* End of while. */
	temp[j]='\0';
	strcpy(c,temp);
	free(temp);
	return;
} /* End of routine. */

