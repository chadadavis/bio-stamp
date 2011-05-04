#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define acids3 "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASX GLX UNK CYH"
#define acids1 "ARNDCQEGHILKMFPSTWYVBZXc"

/* Converts three letter amino acid code to one letter
 *  amino acid code 
 * Note: CYS refers to cystine, CYH refers to cysteine */

char a3to1(char *a3) {
	int i;
	char new;

	new='X';
	for(i=0;i<24; i++) {
	   if (strncmp(&acids3[i*4],a3,3) == 0) 
	      new=(char)acids1[i];
	} 

	return(new);
} 
