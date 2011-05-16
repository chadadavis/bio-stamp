#include "array.h"
#include "pathdef.h"

/* modified heap sort from nr routines */

void hsort(n,sortarr)
int n;
/*float ra[];*/
struct path *sortarr;
{
	int l,j,ir,i;
/*	float rra;*/
	struct path rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
/*			rra=ra[--l];*/
			rra=sortarr[--l];
		else {
/*			rra=ra[ir];*/
			rra=sortarr[ir];
/*			ra[ir]=ra[1];*/
			sortarr[ir]=sortarr[1];
			if (--ir == 1) {
/*				ra[1]=rra;*/
				sortarr[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
/*			if (j < ir && ra[j] < ra[j+1]) ++j;*/
			if (j < ir && sortarr[j].score > sortarr[j+1].score) ++j;
/*			if (rra < ra[j]) {*/
			if (rra.score > sortarr[j].score){
/*				ra[i]=ra[j];*/
				sortarr[i]=sortarr[j];
				j += (i=j);
			}
			else j=ir+1;
		}
/*		ra[i]=rra;*/
		sortarr[i]=rra;
	}
}
