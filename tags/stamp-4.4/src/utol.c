#include <stdio.h>

/* This function changes a character to lowercase. */
char utol(char c) {

	if (c>='A' && c<='Z')
	   return c - 'A' + 'a';
	else return c;
}
