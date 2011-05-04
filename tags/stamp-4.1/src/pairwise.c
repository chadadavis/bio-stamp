/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 Appendix A of the file LICENSE.  Conditions and restrictions for use of
 this package are also in this file.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Biomolecular Modelling Laboratory   Laboratory of Molecular Biophysics
 Imperial Cancer Research Fund       The Rex Richards Building
 Lincoln's Inn Fields, P.O. Box 123  South Parks Road
 London, WC2A 3PX, U.K.              Oxford, OX1 3PG, U.K.
 Tel: +44 171 269 3583               Tel: +44 865 275368
 FAX: +44 171 269 3417               FAX: 44 865 510454
 e-mail: russell@icrf.icnet.uk       e-mail gjb@bioch.ox.ac.uk
 WWW: http://bonsai.lif.icnet.uk/    WWW: http://geoff.biop.ox.ac.uk/

 The WORK is Copyright (1992,1993,1995,1996) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stamp.h>

int pairwise(struct domain_loc *domain, int ndomain, struct parameters *parms) {

	int i,j,k,l;
	int length,nfit;
	int start1,end1;
	int start2,end2;
	int a1,b1,c1,a2,b2,c2;
 	int nsec,nequiv;

	int **hbcmat;

	float rms;
	float score;
	float seqid,secid;
	float secdist,ratio;
	float **pairmat;

	FILE *MAT;

	fprintf(parms[0].LOG,"\n\nPAIRWISE comparisons\n");

	pairmat=(float**)malloc(ndomain*sizeof(float*));
	for(i=0; i<ndomain; ++i) 
	   pairmat[i]=(float*)malloc(ndomain*sizeof(float));

	hbcmat=(int**)malloc(3*sizeof(int*));
	for(i=0; i<3; ++i) {
	   hbcmat[i]=(int*)malloc(3*sizeof(int));
	}
	k=0;
	if(ndomain<2) {
	   fprintf(stderr,"error you can't run PAIRWISE mode on one domain\n");
	   fprintf(stderr,"  have you forgot -s ?\n");
	   exit(-1);
	}
	if(strcmp(parms[0].logfile,"silent")==0) {
	    printf("\n    Sc = STAMP score, RMS = RMS deviation, Align = alignment length\n");
	    printf("    Len1, Len2 = length of domain, Nfit = residues fitted\n");
	    printf("    Secs = no. equivalent sec. strucs. Eq = no. equivalent residues\n");
	    printf("     %%I = seq. identity, %%S = sec. str. identity\n\n");
	    printf("     No.   Domain1  Domain2 Sc     RMS    Len1 Len2  Align NFit Eq. Secs.   %%I   %%S \n");
/* 	           "Pair   1      4mbn    2hhba 7.71 100.00    153   141        133   7  132"  */  

	}
	for(i=0; i<ndomain; ++i) 
	    for(j=i+1; j<ndomain; ++j) {
/*	      fprintf(parms[0].LOG,"before fitting:\n"); disp(domain[j],parms[0].LOG); */
	      fprintf(parms[0].LOG,"\n\nComparison  %d, %s and %s\n",k+1,domain[i].id,domain[j].id);
	      /* now perform either one or two comparisons as requested *
	       * first fit */
	      if(parms[0].NPASS==2) {
		 if(parms[0].BOOLEAN) 
	      	   fprintf(parms[0].LOG,"First fit: BOOLCUT = %5.3f\n",parms[0].first_BOOLCUT);
		 else 
	           fprintf(parms[0].LOG,"First fit: E1 = %5.2f, E2 = %5.2f, CUT=%5.2f, PEN = %5.2f\n",
		      parms[0].first_E1,parms[0].first_E2,parms[0].first_CUTOFF,parms[0].first_PAIRPEN);
	        parms[0].const1=-2*parms[0].first_E1*parms[0].first_E1;
	        parms[0].const2=-2*parms[0].first_E2*parms[0].first_E2;
		parms[0].PAIRPEN=parms[0].first_PAIRPEN;
	        parms[0].CUTOFF=parms[0].first_CUTOFF;
		parms[0].BOOLCUT=parms[0].first_BOOLCUT;		
		if(pairfit(&domain[i],&domain[j],&score,&rms,&length,&nfit,parms,0,&start1,&end1,&start2,&end2,&seqid,&secid,&nequiv,&nsec,hbcmat,0,-1,0)==-1) return -1;
		fprintf(parms[0].LOG,"Second fit: ");
	      } else fprintf(parms[0].LOG,"Fitting with: ");
	      if((score>=parms[0].first_THRESH) || parms[0].NPASS==1) {
		 if(parms[0].BOOLEAN)
	     	   fprintf(parms[0].LOG,"BOOLCUT = %5.3f\n",parms[0].second_BOOLCUT);
		 else
		   fprintf(parms[0].LOG,"E1 = %5.2f, E2 = %5.2f, CUT=%5.2f, PEN = %5.2f\n",
		      parms[0].second_E1,parms[0].second_E2,parms[0].second_CUTOFF,parms[0].second_PAIRPEN);
	         parms[0].const1=-2*parms[0].second_E1*parms[0].second_E1;
	         parms[0].const2=-2*parms[0].second_E2*parms[0].second_E2;
	         parms[0].PAIRPEN=parms[0].second_PAIRPEN;
	         parms[0].CUTOFF=parms[0].second_CUTOFF;
	         parms[0].BOOLCUT=parms[0].second_BOOLCUT;
	         if(pairfit(&domain[i],&domain[j],&score,&rms,&length,&nfit,parms,1,&start1,&end1,&start2,&end2,&seqid,&secid,&nequiv,&nsec,hbcmat,parms[0].PAIRALIGN,k,1)==-1) return -1;
	      } else {
		 fprintf(parms[0].LOG,"not performed since Sc < %7.3f\n",parms[0].first_THRESH);
		 score/=10.0;
	      }
	      if(parms[0].SEC!=0) { /* secondary structure distance */
		sec_content(domain[i].sec,domain[i].ncoords,parms[0].SEC,&a1,&b1,&c1);
		sec_content(domain[j].sec,domain[j].ncoords,parms[0].SEC,&a2,&b2,&c2);
		fprintf(parms[0].LOG,"secondary structure: %10s,  H: %3d, S: %3d (C: %3d)\n",
			domain[i].id,a1,b1,c1);
	        fprintf(parms[0].LOG,"secondary structure: %10s,  H: %3d, S: %3d (C: %3d)\n",
			domain[j].id,a2,b2,c2);
		secdist=sqrt((float)((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)));
		fprintf(parms[0].LOG,"distance = %6.2f\n",secdist);
	      } else { 
		secdist=0.0;
	      }
	      if(strcmp(parms[0].logfile,"silent")!=0) {
	        fprintf(parms[0].LOG,"Sum: %s & %s, Sc: %7.3f, RMS: %7.3f, Len: %d, maxlen: %d nfit: %d nsec: %d nequiv %d\n",
		      domain[i].id,domain[j].id,score,rms,length,max(domain[i].ncoords,domain[j].ncoords),nfit,nsec,nequiv);
	      } else {
	         printf("Pair %3d  %8s %8s %4.2f %6.2f ",
                      k+1,domain[i].id,domain[j].id,score,rms);
		 printf("  %4d %4d  %4d %4d %3d %4d ",
                      domain[i].ncoords,domain[j].ncoords,length,nfit,nequiv,nsec);
		 printf("%6.2f %6.2f ",seqid,secid);
		 if(score<2.0) {
  		    printf(" LOW SCORE");
		 }
	         printf("\n");
		 if(parms[0].PAIROUTPUT==1) printf("See file %s.pairs.%d for alignment and transformations\n",parms[0].transprefix,k+1);
	         fflush(stdout);
	      }
	      if(domain[i].ncoords>domain[j].ncoords) 
		 ratio=(float)domain[i].ncoords/(float)domain[j].ncoords;
	      else 
		 ratio=(float)domain[j].ncoords/(float)domain[i].ncoords;
	      fprintf(parms[0].LOG,"Sum:     lena: %4d, lenb: %4d, ratio: %7.4f, secdist: %6.2f, seqid: %6.2f, secid: %6.2f\n",
		 domain[i].ncoords,domain[j].ncoords,ratio,secdist,seqid,secid);
	      k++;
	      /* create the similarity matrix */
	      if(parms[0].CLUSTMETHOD==0) pairmat[i][j]=pairmat[j][i]=1/rms;
	      if(parms[0].CLUSTMETHOD==1) pairmat[i][j]=pairmat[j][i]=(float)score;
	      l=clock();
	      parms[0].CPUtime+=(float)l/(60000000);
/*	      fprintf(parms[0].LOG,"after fitting:\n"); disp(domain[j],parms[0].LOG);  */
	    } /* end of for (j=.... */ 
	
	fprintf(parms[0].LOG,"\nPairwise calculations done.\n\n");
	/* Now we write pairmat to the file parms[0].matfile. */
	if((MAT=fopen(parms[0].matfile,"w"))==NULL) {
	   fprintf(stderr,"error opening file %s \n",parms[0].matfile);
	   return -1;
	}
	fprintf(MAT,"%d\n",ndomain);
	for(i=0; i<ndomain; ++i) {
		for(j=i+1; j<ndomain; ++j) 
		    fprintf(MAT,"%f  ",pairmat[i][j]);
		fprintf(MAT,"\n");
		}
	fclose(MAT);

	for(i=0; i<ndomain; ++i) free(pairmat[i]);
	free(pairmat);

	return 0;
} /* end of pairwise */
