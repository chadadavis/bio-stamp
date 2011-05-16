/*
Copyright (1997,1998,1999,2010) Robert B. Russell & Geoffrey J. Barton

This file is part of STAMP.

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. A copy of the license
can be found in the LICENSE file in the STAMP installation directory.

STAMP was developed by Robert B. Russell and Geoffrey J. Barton of
current addresses:

 Prof. Robert B. Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
 Cell Networks, University of Heidelberg            College of Life Sciences
 Room 564, Bioquant                                 University of Dundee
 Im Neuenheimer Feld 267                            Dow Street
 69120 Heidelberg                                   Dundee DD1 5EH
 Germany                                            UK
                                                
 Tel: +49 6221 54 513 62                            Tel: +44 1382 385860
 Fax: +49 6221 54 514 86                            FAX: +44 1382 385764
 Email: robert.russell@bioquant.uni-heidelberg.de   E-mail g.j.barton@dundee.ac.uk
 WWW: http://www.russell.embl-heidelberg.de         WWW: http://www.compbio.dundee.ac.uk

 All use of STAMP must cite: 

 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*/
#include <stdio.h>
#include <stdlib.h>
#include "gjutil.h"
#include "gjnoc.h"

/* Given an upper diagonal matrix, return a tree in stamp format via Geoff's
 *  OC routine */

struct cluster *get_clust(double **matrix, char **ids, int ndomain, char *noc_parms) {

	int i,j,k;
	
	int nclust;

	struct gjnoc *gjclust;
	struct cluster *cl;

/*	printf("Parameters are: %s\n",noc_parms);  */

	nclust=ndomain-1; 
	gjclust=GJnoc(matrix,ids,ndomain,noc_parms);
	cl=(struct cluster*)malloc((nclust)*sizeof(struct cluster));

/*	printf("Order:\n");
	for(i=0; i<ndomain; ++i) {
		printf("%4d\n",gjclust->order[i]);
	}
	for(i=0; i<nclust; ++i) {
		printf("\nCluster %4d\nA:",i+1);
		for(j=0; j<gjclust->clusters[i].a.number; ++j) {
		printf("%4d ",gjclust->clusters[i].a.member[j]);
		}
		printf("\nB:");
		for(j=0; j<gjclust->clusters[i].b.number; ++j) {
                  printf("%4d ",gjclust->clusters[i].b.member[j]);
                }
	}
*/

	for(i=0; i<nclust; ++i) { /* take Geoff's structure and convert it to a tree */
	      cl[i].a.member=(int*)malloc(gjclust->clusters[i].a.number*sizeof(int));
	      cl[i].a.number=gjclust->clusters[i].a.number;
	      cl[i].b.member=(int*)malloc(gjclust->clusters[i].b.number*sizeof(int));
              cl[i].b.number=gjclust->clusters[i].b.number;
	      /* Now sort out the clusters */
	      for(j=0; j<gjclust->clusters[i].a.number; ++j) {
			cl[i].a.member[j]=gjclust->order[gjclust->clusters[i].a.member[j]]; 
	      }
	      for(j=0; j<gjclust->clusters[i].b.number; ++j) {
                        cl[i].b.member[j]=gjclust->order[gjclust->clusters[i].b.member[j]]; 
              }
	}
/*	printf("Returned from GJnoc...\n");
	for(i=0; i<nclust; ++i) {
	   printf("Cluster: %4d (",i+1);
              for(k=0; k<cl[i].a.number; ++k) printf("%d ",cl[i].a.member[k]);
	   printf("   and ");
	   for(k=0; k<cl[i].b.number; ++k) printf("%d ",cl[i].b.member[k]);
	   printf(") \n\n");
	}
*/
	for(i=0; i<nclust; ++i) {
	  free(gjclust->clusters[i].a.member);
	  free(gjclust->clusters[i].b.member);
	}
	free(gjclust->clusters);
	free(gjclust->order);
	free(gjclust);

	return cl;
}
