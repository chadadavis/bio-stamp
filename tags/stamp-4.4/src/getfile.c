#include <stdlib.h>
#include <math.h>

#include <stamp.h>

#define MAXLINE 200 

/* getfile: 
 *  opens a file containing a list of appropriate directories, suffixes and
 *  prefixes.  It then makes combinates of these with the supplied code
 *  to see if a file exists 
 *
 * Allocates memory for a string corresponding to the filename read in,
 *  then returns this string or '\0' if none is found 
 *
 * Modification March 1999 - now tries to expand into a "distr" style hierarchy */

char *getfile(char *code_safe, char *dirfile, int code_length, FILE *OUTPUT) {

	char *temp;
	char *dir;
	char *buff;
	char prefix[10],suffix[10];
	char *code;
	char *CODE;
	char *two_letter; 
	char a;

	int i,j;
	int found;

	FILE *DIR;

	dir=(char*)malloc(1000*sizeof(char)); 
	buff=(char*)malloc(1000*sizeof(char)); 
	temp=(char*)malloc(1000*sizeof(char)); 
	CODE=(char*)malloc((strlen(code_safe)+1)*sizeof(char));
	code=(char*)malloc((strlen(code_safe)+1)*sizeof(char));
	two_letter = (char*)malloc(3*sizeof(char));

	/* copy code_safe into code to prevent the changing case from
	 *  carrying on through the program */

	a = code[code_length];
	code[code_length]='\0'; /* truncates the code if necessary */
	for(i=0; i<strlen(code_safe); ++i) code[i]=utol(code_safe[i]);
	for(i=0; i<strlen(code_safe); ++i) CODE[i]=ltou(code_safe[i]); /* upper case code */
	code[code_length]='\0';
	CODE[code_length]='\0';
	two_letter[0] = code[1]; two_letter[1] = code[2]; two_letter[2] = '\0'; /* For PDB distr style */

	if((DIR=fopen(dirfile,"r"))==NULL) {
	   fprintf(stderr,"error: no directories file found\n");
	   temp[0]='\0';
	   return temp;
	}
	found=0;
	while((fgets(buff,MAXLINE,DIR))!=NULL) {

	  sscanf(buff,"%s %s %s",dir,&prefix[0],&suffix[0]);
	  if(prefix[0] == '_') { prefix[0] = '\0'; }
	  if(suffix[0] == '_') { suffix[0] = '\0'; }

	  /* Lower case, all files one directory */
	  sprintf(temp,"%s/%s%s%s",dir,prefix,code,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Lower case, files in 'distr' type directories */
	  sprintf(temp,"%s/%s/%s%s%s",dir,two_letter,prefix,code,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Upper case, all files one directory */
	  sprintf(temp,"%s/%s%s%s",dir,prefix,CODE,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Upper case, files in 'distr' type directories */
	  sprintf(temp,"%s/%s/%s%s%s",dir,two_letter,prefix,CODE,suffix);
	  if(testfile(temp)) { found=1; break; } 

	}

	if(found==0) {
	   temp[0]='\0';
	}
	fclose(DIR);
	free(dir); 
	free(buff);
	free(code);
	free(CODE);
	free(two_letter);

	return temp;
}
