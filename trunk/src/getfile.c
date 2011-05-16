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

	/*
	 * ensure that code_length is always <= strlen(code_safe).
	 * - sometimes code_length is fixed at four, when it is assumed
	 *   that the first four characters give the PDB idcode, but
	 *   the actual code may have fewer than four characters
	 */
	code_length = (code_length > strlen(code_safe)) ? strlen(code_safe) : code_length;

	temp=(char*)malloc(1000*sizeof(char)); 
	if(temp == NULL) {
	    fprintf(stderr, "error: getfile: malloc failed on temp\n");
	    return temp;
 	}
	memset(temp, '\0', 1000);

	dir=(char*)malloc(1000*sizeof(char)); 
	if(dir == NULL) {
	    fprintf(stderr, "error: getfile: malloc failed on dir\n");
	    return temp;
 	}

	buff=(char*)malloc(1000*sizeof(char)); 
	if(buff == NULL) {
	    fprintf(stderr, "error: getfile: malloc failed on buff\n");
	    return temp;
 	}

	CODE=(char*)malloc((code_length + 1) * sizeof(char));
	if(CODE == NULL) {
	    fprintf(stderr, "error: getfile: malloc failed on CODE\n");
	    return temp;
 	}

	code=(char*)malloc((code_length + 1) * sizeof(char));
	if(code == NULL) {
	    fprintf(stderr, "error: getfile: malloc failed on code\n");
	    return temp;
 	}

	two_letter = (char*)malloc(3 * sizeof(char));
	if(two_letter == NULL) {
	    fprintf(stderr, "error: getfile: malloc failed on two_letter\n");
	    return temp;
 	}

	/* copy code_safe into code to prevent the changing case from
	 *  carrying on through the program */

	/*
	  a = code[code_length];
	  code[code_length] = '\0';
	  CODE[code_length] = '\0';
	*/
	memset(code, '\0', code_length + 1);
	for(i = 0; i < code_length; ++i) {
	    code[i] = utol(code_safe[i]);
	}

	memset(CODE, '\0', code_length + 1);
	for(i = 0; i < code_length; ++i) {
	    CODE[i] = ltou(code_safe[i]); /* upper case code */
	}

	/* For PDB distr style */
	/*
	  two_letter[0] = code[1];
	  two_letter[1] = code[2];
	*/
	memset(two_letter, '\0', 3);
	 /* code_length can sometimes be < 3, though only for ids that don't include the PDB idcode where the following isn't relevant */
	for(i = 1; (i < code_length) && (i < 3); ++i) {
	    two_letter[i - 1] = code[i];
	}

	if((DIR=fopen(dirfile,"r"))==NULL) {
	   fprintf(stderr,"error: no directories file found\n");
	   return temp;
	}
	found=0;

	while((fgets(buff,MAXLINE,DIR))!=NULL) {

	  sscanf(buff,"%s %s %s",dir,&prefix[0],&suffix[0]);
	  if(prefix[0] == '_') { prefix[0] = '\0'; }
	  if(suffix[0] == '_') { suffix[0] = '\0'; }
	  
	  /* Lower case, all files one directory */
	  memset(temp, '\0', 1000);
	  sprintf(temp,"%s/%s%s%s",dir,prefix,code,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Lower case, files in 'distr' type directories */
	  memset(temp, '\0', 1000);
	  sprintf(temp,"%s/%s/%s%s%s",dir,two_letter,prefix,code,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Upper case, all files one directory */
	  memset(temp, '\0', 1000);
	  sprintf(temp,"%s/%s%s%s",dir,prefix,CODE,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Upper case, files in 'distr' type directories */
	  memset(temp, '\0', 1000);
	  sprintf(temp,"%s/%s/%s%s%s",dir,two_letter,prefix,CODE,suffix);
	  if(testfile(temp)) { found=1; break; } 

	}

	if(found == 0) {
	    memset(temp, '\0', 1000);
	}
	fclose(DIR);
	free(dir); 
	free(buff);
	free(code);
	free(CODE);
	free(two_letter);

	return temp;
}
