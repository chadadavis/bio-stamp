#include <stdio.h>


/* This function changes a character to uppercase. */
char ltou(char c) {

	if (c>='a' && c<='z')
	   return c - 'a' + 'A';
	else return c;
}

