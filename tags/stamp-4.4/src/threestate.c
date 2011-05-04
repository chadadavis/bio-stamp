#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* threestate -- converts secondary structures into three state 
 *  you must give the routine four strings:
 *
 * sec -- the secondary structure assignment 
 * helix -- a string of helical definitions (eg. DSSP "HGI")
 * extended -- a string of extended definitions (eg. DSSP "EB")
 * coil     -- a string of coil definitions (eg. DSSP "-TS") 
 *
 * The string supplied is altered to contain only H, B and c
 *  according to the supplied definitions 
 *
 * It will leave un-recognised characters (ie. spaces) unaltered 
 *
 */

int threestate(char *sec,char *helix,char *extended,char *coil) {

	int i,j;
	int slen,hlen,elen,clen;

	slen=strlen(sec);
	hlen=strlen(helix);
	elen=strlen(extended);
	clen=strlen(coil);

	for(i=0; i<strlen(sec); ++i) {
	     /* coil is done first, since this will keep changed 'H' or 'B' from being coil */
	     for(j=0; j<clen; ++j) {
		if(sec[i]==coil[j]) {
		   sec[i]='c';
		   break;
	   	}
	     }
	     for(j=0; j<hlen; ++j) {
		if(sec[i]==helix[j]) {
		   sec[i]='H'; 
		   break;
		}
	     }
	     for(j=0; j<elen; ++j) {
		if(sec[i]==extended[j]) {
		   sec[i]='B';
		   break;
		}
	     }
	}

	return 0;
}
