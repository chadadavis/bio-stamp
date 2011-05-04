#include <stamp.h>

int disp(struct domain_loc domain, FILE *OUTPUT) {

	int i,j,k;
	fprintf(OUTPUT,"The first ten coordinates\n");
        for(i=0; i<10; ++i) 
	   fprintf(OUTPUT,"%c %8d %8d %8d\n",domain.aa[i],domain.coords[i][0],domain.coords[i][1],domain.coords[i][2]);
	fprintf(OUTPUT,"R,V:\n");
	printmat(domain.R,domain.V,3,OUTPUT);
	fprintf(OUTPUT,"r,v:\n");
	printmat(domain.r,domain.v,3,OUTPUT);
	return 0;
}
