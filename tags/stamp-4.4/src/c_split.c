#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char **RBR_c_split(char *str, int *n,  char delimiter) {
	/* breaks a string up into substrings using a delimiter */

	int i,j,k,len;
	char **sstr;

	sstr=(char**)malloc(sizeof(char*));
	
	(*n)=0;
	i=0;
	len=strlen(str);

 	while(i<len && str[i]!='\0' && str[i]!='\n') {
	   while(i<len && str[i]==delimiter) i++;
	   if(i>=len || str[i]=='\0' || str[i]=='\n') break;
	   /* we are at a new sub-string */
	   sstr[(*n)]=(char*)malloc(sizeof(char));
	   j=0;
	   while(i<len && str[i]!='\0' && str[i]!='\n' && str[i]!=delimiter) {
		sstr[(*n)][j]=str[i];
		j++; i++;
	        sstr[(*n)]=(char*)realloc(sstr[(*n)],(j+1)*sizeof(char));
	   }
	   sstr[(*n)][j]='\0';
	   (*n)++;
	   sstr=(char**)realloc(sstr,((*n)+1)*sizeof(char*));
	   i++;
	}
	return sstr;
}
