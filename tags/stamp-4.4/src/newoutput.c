#include <stdio.h>
#include <stamp.h>

int newoutput(FILE *TRANS, struct domain_loc *domain, int ndomain, int writetrans) {

	int i,j,k;
        for(i=0; i<ndomain; ++i) {
	   if(printdomain(TRANS,domain[i],writetrans)==-1) return -1;
        } 
        return 0;
}

int printdomain(FILE *TRANS, struct domain_loc domain, int writetrans) {

	int i,j,k,l,m;
   	fprintf(TRANS,"%s ",domain.filename);
	fprintf(TRANS,"%s {",domain.id);
	for(j=0; j<domain.nobj; ++j) {
	   if(domain.start[j].cid==' ') domain.start[j].cid='_';
	   if(domain.start[j].in==' ') domain.start[j].in='_';
	   if(domain.end[j].cid==' ') domain.end[j].cid='_';
	   if(domain.end[j].in==' ') domain.end[j].in='_';

	   switch(domain.type[j]) {
	      case 1: fprintf(TRANS," ALL "); break;
	      case 2: fprintf(TRANS," CHAIN %c ",domain.start[j].cid); break;
	      case 3: fprintf(TRANS," %c %d %c to %c %d %c ",
		domain.start[j].cid,domain.start[j].n,
		domain.start[j].in,
		domain.end[j].cid,domain.end[j].n,
		domain.end[j].in); break;
	   } 
	}
	if(writetrans) {
	   fprintf(TRANS,"\n");
	   for(j=0; j<3; ++j) {
	     for(k=0; k<3; ++k) fprintf(TRANS,"%10.5f ",domain.R[j][k]); 
	     fprintf(TRANS,"     %10.5f ",domain.V[j]);
	     if(j!=2) fprintf(TRANS,"\n");
	   }
	} 
	fprintf(TRANS," }\n");

	return 0;
}
