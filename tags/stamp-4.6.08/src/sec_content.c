#include <stdio.h>

/* return the percent helix sheet and coil for a given secondary
 *  structure string */

int sec_content(char *sec, int npos, int type, int *pera, int *perb, int *perc) {

	int i;
	int a,b,c;
	
	a=b=c=0;
	for(i=0; i<npos; ++i) {
   	   if((type==1 && (sec[i]=='H' || sec[i]=='G')) ||
	      (type==2 && (sec[i]=='H' || sec[i]=='3')) ||
	      (type==3 && (sec[i]=='H')) ) a++; 
	   else if ((type==1 && sec[i]=='E') ||
	           ((type==2 || type==3) && sec[i]=='B') ) b++; 
	   else c++;
	}
	(*pera)=(int)(100*(float)a/(float)npos);  /* Alpha */
        (*perb)=(int)(100*(float)b/(float)npos);  /* Beta */
	(*perc)=(int)(100*(float)c/(float)npos);  /* Coil */
	return 0;
}
