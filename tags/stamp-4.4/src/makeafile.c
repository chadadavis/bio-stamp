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
#include "stamp.h"

#define VSMALL 0.000001

int makefile(struct domain_loc *domain, int ndomain, struct cluster cl, 
	int nclust, float score, float rms, int length, int nfit, 
	float *Pij, float *Dij, float *dist, float *Pijp,
	int PAIRWISE, struct parameters *parms) {

	char *tmp;

	int i,j,k,l,m;
	int nogap;

	int *asec,*bsec;

	float Cij;

	FILE *OUT;

	


	tmp=(char*)malloc(1000*sizeof(char));


	/* making the filename */
	if(!PAIRWISE) sprintf(tmp,"%s.%d",parms[0].transprefix,nclust+1);
	else sprintf(tmp,"%s.pairs.%d",parms[0].transprefix,nclust+1);
	fprintf(parms[0].LOG,"File is: %s\n",tmp);
	if((OUT=fopen(tmp,"w"))==NULL) {
	   fprintf(stderr,"error: opening file %s\n",tmp);
	   return -1;
	} 

	asec=(int*)malloc(cl.a.number*sizeof(int));
	bsec=(int*)malloc(cl.b.number*sizeof(int));

	/* outputing the domain descriptors */
	for(i=0; i<cl.a.number; ++i) {
	   if(printdomain(OUT,domain[cl.a.member[i]],1)==-1) return -1;
        }
   	for(i=0; i<cl.b.number; ++i) {
           if(printdomain(OUT,domain[cl.b.member[i]],1)==-1) return -1;
        }

	/* now for some comments and junk */
	fprintf(OUT,"\n\n");
	fprintf(OUT,"Alignment score  Sc = %f\n",score);
	fprintf(OUT,"Alignment length Lp = %d\n",length);
	fprintf(OUT,"RMS deviation after fitting on %d atoms =  %f\n",nfit,rms);
	fprintf(OUT,"Secondary structures are ");
	switch(parms[0].SECTYPE) {
	   case 0: fprintf(OUT,"not specified\n"); break;
	   case 1: fprintf(OUT,"from DSSP\n"); break;
	   case 2: fprintf(OUT,"read in from the file %s\n",parms[0].secfile);
	   default: fprintf(OUT,"not specified\n");
	}
	fprintf(OUT,"\n\n\n");

	/* now for the alignment in blocfile format 
	 * '>' identifiers first */
	for(i=0; i<cl.a.number; ++i) 
	   fprintf(OUT,">%s   (cluster A) sequence\n",domain[cl.a.member[i]].id);
	for(i=0; i<cl.b.number; ++i)
	   fprintf(OUT,">%s   (cluster B) sequence\n",domain[cl.b.member[i]].id);
	fprintf(OUT,">space \n");
	for(i=0; i<cl.a.number; ++i) {
	   fprintf(OUT,">%s_",domain[cl.a.member[i]].id);
	   switch(parms[0].SECTYPE) {
	    case 0: fprintf(OUT,"none (cluster A) secondary structure not specified\n"); break;
            case 1: fprintf(OUT,"dssp (cluster A) secondary structure from DSSP\n"); break;
	    case 2: fprintf(OUT,"spec read in from the file %s\n",parms[0].secfile);
	    default: fprintf(OUT,"none  (cluster A) secondary structure not specified\n");
           }
	}

	for(i=0; i<cl.b.number; ++i) {
	   fprintf(OUT,">%s_",domain[cl.b.member[i]].id);
	   switch(parms[0].SECTYPE) {
	      case 0: fprintf(OUT,"none  (cluster B) secondary structure not specified\n"); break;
	      case 1: fprintf(OUT,"dssp  (cluster B) secondary structure from DSSP\n"); break;
	      case 2: fprintf(OUT,"spec  (cluster B) secondary structure read in from the file %s\n",parms[0].secfile);
	      default: fprintf(OUT,"none  (cluster B) secondary structure not specified\n");
	   }
	}
	/* now for the '#' descriptors */
	fprintf(OUT,"#T -- '1' = equivalenced residues \n");

	if(!parms[0].BOOLEAN) {
	  if(!PAIRWISE) {
	     fprintf(OUT,"#P -- averaged Pij\n");
	     fprintf(OUT,"#A -- distance between averaged CA atoms in angtroms\n");
	  }
	  fprintf(OUT,"#G -- Pij' value\n");
	} else {
	  if(!PAIRWISE) fprintf(OUT,"#A -- distance between averaged CA atoms in angtroms\n");
	}
	/* header */
	for(i=0; i<cl.a.number; ++i) fprintf(OUT,"A");
	for(i=0; i<cl.b.number; ++i) fprintf(OUT,"B");
	fprintf(OUT," ");
	for(i=0; i<cl.a.number; ++i) fprintf(OUT,"A");
	for(i=0; i<cl.b.number; ++i) fprintf(OUT,"B");
	/* if(!parms[0].BOOLEAN) fprintf(OUT," use  Pij      Cij      Dij      Distance Pij'\n"); */
	if(!parms[0].BOOLEAN) {
	   fprintf(OUT,"equiv Pij      Distance Pij'\n"); 
	} else {
	   fprintf(OUT," use Distance\n");
	}

	/* now the alignment is output vertically */
	fprintf(OUT,"* iteration 1\n");
	for(i=0; i<cl.a.number; ++i) asec[i]=0;
	for(i=0; i<cl.b.number; ++i) bsec[i]=0;  /* counters for secondary structure strings */
	for(i=0; i<strlen(domain[cl.a.member[0]].align); ++i) { 
	   nogap=1;
	   for(j=0; j<cl.a.number; ++j) {
	      fprintf(OUT,"%c",domain[cl.a.member[j]].align[i]);
	      nogap*=(domain[cl.a.member[j]].align[i]!=' ');
	   }
	   for(j=0; j<cl.b.number; ++j) {
	      fprintf(OUT,"%c",domain[cl.b.member[j]].align[i]);
	      nogap*=(domain[cl.b.member[j]].align[i]!=' ');
	   }
	   fprintf(OUT," ");
	   for(j=0; j<cl.a.number; ++j) 
	      if(domain[cl.a.member[j]].align[i]!=' ') {
		 fprintf(OUT,"%c",domain[cl.a.member[j]].sec[asec[j]]);
		 asec[j]++;
	      } else fprintf(OUT," ");
	   for(j=0; j<cl.b.number; ++j) 
	      if(domain[cl.b.member[j]].align[i]!=' ') { 
		 fprintf(OUT,"%c",domain[cl.b.member[j]].sec[bsec[j]]);
	         bsec[j]++; 
  	      } else fprintf(OUT," ");
	   fprintf(OUT," ");


	if(!parms[0].BOOLEAN) {
	   if(nogap) {
	      fprintf(OUT,"%1d ",(Pijp[i]>=parms[0].CUTOFF));
	      fprintf(OUT,"%8.5f ",Pij[i]*(Pij[i]>VSMALL)); /* Zero values imply rounding error */
	      if(!PAIRWISE) {
	         fprintf(OUT,"%9.5f ",dist[i]*(dist[i]>VSMALL));
	         fprintf(OUT,"%9.5f ",Pijp[i]*(Pijp[i]>VSMALL));
	      }
	   }
	} else {
	   if(nogap) {
	      fprintf(OUT,"%1d ",(int)Pij[i]);
	      if(!PAIRWISE) fprintf(OUT,"%9.5f ",dist[i]*(dist[i]>VSMALL));
	   }
	}
	   
	   fprintf(OUT,"\n");
	}
	fprintf(OUT,"*\n");

	free(asec); free(bsec); free(tmp);

	fclose(OUT);

	return 0;
}
