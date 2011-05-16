#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* takes a three state SS assignment (H,B or c) and removes those H or B runs that are less than
 *  minhelixlen or minstrandlen (replacing them by c) */
int smoothsec(char *sec, int minhelixlen, int minstrandlen) {

	int i,j,k,l;
	int len,add;

	len=strlen(sec);
	for(j=0; j<len; ++j) {
          if(sec[j]=='H') {
	     add=0;
	     /* a helix */
	     for(k=j; k<len; ++k) {
	        if(sec[k]=='H') add++;
	        else break;
   	     }
	     if(add<minhelixlen) { /* too short -- remove */
	        for(l=j; l<k; ++l) sec[l]='c';
	     }
	     j=k-1;
	   } else  if(sec[j]=='B') {
	      add=0;
	      /* a strand */
	      for(k=j; k<len; ++k) {
                  if(sec[k]=='B') add++;
                  else break;
              }    
              if(add<minstrandlen) { /* too short -- remove */ 
                 for(l=j; l<k; ++l) sec[l]='c';
              }
	      j=k-1;
           } else if(sec[j]!='c' && sec[j]!=' ') {
		fprintf(stderr,"error: unrecognised secondary structure character found: %c\n",sec[j]);
		return -1;
	   }
        }
	return 0;
}
