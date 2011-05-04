#include <stdio.h>
#include <stdlib.h>

#include "stamp.h"

/* This routine looks for a DSSP file, if found it reads in the DSSP summary
 *  and stores it in domain[i].sec for future use */

int getks(struct domain_loc *domain, int ndomain, struct parameters *parms) {

	FILE *DSSP;
	int i,j,k;
	int include,count;
	int total,add;
	int retval;
	char c,chain;
	char *filename;
	char *label;
	char *empty = NULL; 
	
	count=0;
	retval=0;


	for(i=0; i<ndomain; ++i) {
	   fprintf(parms[0].LOG,"%s -- ",domain[i].id);
	   /* first check to see if there is a DSSP file that uses the whole ID name */
           filename=getfile(domain[i].id,parms[0].dsspfile,strlen(domain[i].id),parms[0].LOG);
	   /* if there is, use it, otherwise try to get one using the four letter code */
	   if(filename[0]=='\0') filename=getfile(domain[i].id,parms[0].dsspfile,4,parms[0].LOG);
           if(filename[0]=='\0') {
               free(filename);
               filename=getfile(domain[i].id,parms[0].dsspfile,4,parms[0].LOG);
           }
	   if(filename[0]=='\0') {
	      fprintf(parms[0].LOG," no DSSP file found for %s\n",domain[i].id);
	      for(j=0; j<domain[i].ncoords; ++j) domain[i].sec[j]='?';
	      domain[i].sec[j]='\0';
	      retval=-1; /* if any of the sequences have missing secondary structures */
	   } else {
	      DSSP=openfile(filename,"r");
	      total=0;
	      fprintf(parms[0].LOG," using file %s\n",filename);
	      for(j=0; j<domain[i].nobj; ++j) {
	         if(get_dssp_sum(DSSP,domain[i].start[j],domain[i].end[j],
		    domain[i].type[j],&domain[i].sec[total],domain[i].reverse[j],
		    (parms[0].MAX_SEQ_LEN-total),&add,parms[0].LOG)==-1) 
		    retval=-1;
	         total+=add;
		 closefile(DSSP,filename);
	         DSSP=openfile(filename,"r");
	      }
	      closefile(DSSP,filename);
	      free(filename);
	      if(total!=domain[i].ncoords) {
		 fprintf(parms[0].LOG,"warning: DSSP summary found was incomplete -- the results may have errors\n");
		 for(j=total; j<domain[i].ncoords; ++j) domain[i].sec[j]='?';
		 domain[i].sec[j]='\0';
	      }
	      display_align(&domain[i].aa,1,&domain[i].sec,1,&domain[i].aa,&domain[i].sec,empty,empty,parms[0].COLUMNS,0,0,parms[0].LOG);
	   }
	}
	return retval;
	      
} 
	  


