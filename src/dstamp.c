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

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include "dstamp.h"

/* Reads in a STAMP file, and an optional parameter file and
 *  produces a file for ALScript (G.J. Barton, Protein Eng. 6, 37-40, 1993) */

main(argc,argv)
int argc;
char *argv[];
{
	int i,j,k;
	int nbloc,nseq;
	int nstamp,nstampseq,nstamppos;
	int startrel,endrel;
	int startnonrel,endnonrel;
	int Helix,Strand,Other,Gap;
	int bloclen;

	int *reliable;

	char tmp[100];
	char c;

	FILE *IN,*OUT,*PARM;

	struct parameters *parms;
	struct seqdat *bloc;
	struct stampdat *stampstuff;


	parms=(struct parameters*)malloc(sizeof(struct parameters));
	/* open files */
	/* set default parameter values
	 *  these are changed if a paramter file is specified */
	parms[0].TYPE='G';  	/*  Use Pij' values */
	parms[0].CUTOFF=6.0;	/*  Pij' >= 6.0 only */
	parms[0].WINDOW=3;		/*  Stretches of three or more only */
	parms[0].POINTSIZE=8.0;	/*  Pointsize */
	parms[0].SEC=1;		/*  Display secondary structure */
	parms[0].BOXSEQ=1;		/*  Box secondary structre based on STAMP reliability */
	parms[0].BOXSEC=1;		/*  Box secondary structre as above */
	parms[0].SMALLSEQ=1;	/*  1 ==> put non-reliable regions in small */
	parms[0].SMALLSEC=1;	/*  as above for secondary structure */
	parms[0].CASESEQ=1;	/*  1 ==> non-reliable regions in lower case */
	parms[0].CASESEC=1;	/*  as above for secondary structure */
	parms[0].BOLDSEQ=1;	/*  1 ==> reliable regions in bold */
	parms[0].BOLDSEC=1;	/*  as above for secondary structure */ 
	parms[0].SECSUM=0;		/*  1 ==> display secondary structure  summary only */
	parms[0].VERBOSE=0;	/*  Run ALScript in silent mode */
	strcpy(&parms[0].prefix[0],"dstamp");

	readcom(argc,argv,parms);

	printf("DSTAMP R.B. Russell, 1995\n");
	printf(" Makes input for ALSCRIPT, to display STAMP\n");
	printf("  alignments in PostScript format\n\n");
	if((IN=fopen(parms[0].filename,"r"))==NULL) {
	   fprintf(stderr,"error opening file %s\n",parms[0].filename);
	   exit(-1);
	}

	/* read the input file */
	printf(" Reading Alignment...\n");
	bloc=(struct seqdat*)malloc(MAX_N_SEQ*sizeof(struct seqdat));
	printf(" ");
	if(Agetbloc(IN,bloc,&nbloc)==-1) exit(-1); rewind(IN);

	bloclen=strlen(&bloc[1].seq[1]);

	printf(" Getting STAMP information...\n");
	stampstuff=(struct stampdat*)malloc(MAX_STAMP_NUM*sizeof(struct stampdat));
	if(getstampdat(stampstuff,IN,&nstamp,&nstampseq,&nstamppos,bloclen)==-1) 
	    exit(-1);
	
	/* determine the reliable regions */
	if((reliable=getstamprel(stampstuff,nstamp,nstamppos,parms[0].TYPE,parms[0].CUTOFF,parms[0].WINDOW))==NULL) 
	   exit(-1);

	if(nstamppos!=strlen(&bloc[1].seq[1]) || nstampseq != nbloc) {
	   fprintf(stderr,"error: something wrong with STAMP file\n");
	   fprintf(stderr,"	  STAMP length is %d, Alignment length is %d\n",nstamppos,strlen(&bloc[1].seq[1]));
	   fprintf(stderr,"       STAMP nseq is %d, Alignment nseq is %d\n",nstampseq,nbloc);
	   exit(-1);
	}
	nseq=(nbloc-1)/2;

	/* First we make another bloc file, according to the 
	 *  parameters */
	sprintf(&tmp[0],"%s_als.bloc",parms[0].prefix);
	if((OUT=fopen(tmp,"w"))==NULL) {
	   fprintf(stderr,"error opening file %s \n",tmp);
	   exit(-1);
	}
	printf(" Writing new block file to %s...\n",tmp);
	fprintf(OUT,"This is a blockfile for use in conjuction with the program\n DSTAMP.  There ought to be a corresponding ALScript file\n in this directory.\n\n");
	/* sequences first */
	for(i=0; i<nseq; ++i) 
	  fprintf(OUT," %4d >%s \n%s\n",i+1,bloc[i+1].id,bloc[i+1].title);
	
	fprintf(OUT,"      > \n");
	/* DSSP summary if necessary */
	if(parms[0].SEC) {
	   if(!parms[0].SECSUM)  {
	     for(i=0; i<nseq; ++i) 
	        fprintf(OUT," %4d >%s \n%s\n",i+1,bloc[nseq+i+2].id,bloc[nseq+i+2].title);
	   } else {
	     fprintf(OUT,"      >DSSP Summary\n");
	   }
	}
	/* now the sequence alignment */
	fprintf(OUT,"*\n");
	for(i=0; i<nstamppos; ++i) {
	   for(j=0; j<nseq; ++j) 
	      if(parms[0].CASESEQ && !reliable[i]) 
		fprintf(OUT,"%c",utol(bloc[j+1].seq[i+1]));
	      else 
		fprintf(OUT,"%c",ltou(bloc[j+1].seq[i+1]));
	   fprintf(OUT," ");
	   if(parms[0].SEC) {
	      if(!parms[0].SECSUM) {
		 for(j=0; j<nseq; ++j) {
		   if(parms[0].CASESEC && !reliable[i]) 
		      fprintf(OUT,"%c",utol(bloc[nseq+j+2].seq[i+1]));
		   else 
		      fprintf(OUT,"%c",ltou(bloc[nseq+j+2].seq[i+1]));
		}
	      } else {
		 Helix=Strand=Other=Gap=0;
		 for(j=0; j<nseq; ++j) {
		    switch(bloc[nseq+j+2].seq[i+1]) {
		       case 'H': case 'G': Helix++; break;
		       case 'E': case 'B': Strand++; break;
		       case ' ': Gap++; break;
		       default: Other++;
		    }
		 }
		 if(Helix>Strand && Helix >Other && Helix>Gap) 
		    c='H';
		 else if(Strand>Helix && Strand>Other && Strand>Gap)
		    c='E';
		 else if(Gap>Helix && Gap>Strand && Gap>Other)
		    c=' ';
		 else 
		    c='-';
		 if(parms[0].CASESEC && !reliable[i]) 
		   c=utol(c);
		 else 
		   c=ltou(c);
		 fprintf(OUT,"%c",c);
	     }
	  }
	  fprintf(OUT,"\n");
	}
	fprintf(OUT,"*\n");
	fclose(OUT);

	
	/* AlScript output */
	sprintf(&tmp[0],"%s_als.com",parms[0].prefix);
	if((OUT=fopen(tmp,"w"))==NULL) {
	   fprintf(stderr,"error opening file %s\n",tmp);
	   exit(-1);
	}
	printf(" Writing ALSCRIPT command file to %s...\n",tmp);
	/*  introductory stuff first */
	if(!parms[0].VERBOSE) fprintf(OUT,"SILENT_MODE\n");
	fprintf(OUT,"#\n# STAMP output for AlScript\n#\n# RBR October 1992\n#\n");
	fprintf(OUT,"# Input and output files \n");
	fprintf(OUT,"BLOCK_FILE	%s_als.bloc\n",parms[0].prefix);
	fprintf(OUT,"OUTPUT_FILE %s_als.ps\n",parms[0].prefix);
	fprintf(OUT,"#\n#\n");
	fprintf(OUT,"LANDSCAPE\n");
	fprintf(OUT,"POINTSIZE %f\n",parms[0].POINTSIZE);
	fprintf(OUT,"IDENT_WIDTH 10\n");
	fprintf(OUT,"VERTICAL_SPACING 5\n");
	fprintf(OUT,"#\n#\n#\n");
	fprintf(OUT,"DEFINE_FONT 0 Helvetica      DEFAULT\n");
	fprintf(OUT,"DEFINE_FONT 1 Helvetica REL  0.75  \n");
	fprintf(OUT,"DEFINE_FONT 7 Helvetica REL 0.5\n");
	fprintf(OUT,"DEFINE_FONT 3 Helvetica-Bold DEFAULT \n");
	fprintf(OUT,"DEFINE_FONT 4 Times-Bold     DEFAULT \n");
	fprintf(OUT,"DEFINE_FONT 5 Helvetica-BoldOblique  DEFAULT\n");
	fprintf(OUT,"SETUP\n");
	fprintf(OUT,"SUB_ID %d \"\" \n",nseq+1); /* removes the word 'space' */


	/* Finding reliable regions and boxing, etc them appropriately */
	i=0;
	while(i<=nstamppos) { 
	   if(reliable[i]) { /* reliable region */
	      startrel=i;
	      while(reliable[i] && i<=nstamppos) i++;
	      endrel=i-1;
	      if(parms[0].BOXSEQ) fprintf(OUT,"BOX_REGION %d %d %d %d\n", startrel+1,1,endrel+1,nseq);
	      if(parms[0].BOXSEC) {
	         if(!parms[0].SECSUM)
		    fprintf(OUT,"BOX_REGION %d %d %d %d\n", startrel+1,nseq+2,endrel+1,nseq*2+1);
		 else 
		    fprintf(OUT,"BOX_REGION %d %d %d %d\n", startrel+1,nseq+2,endrel+1,nseq+3);
	      }
	      if(parms[0].BOLDSEQ) fprintf(OUT,"FONT_REGION %d %d %d %d 3\n", startrel+1,1,endrel+1,nseq);
	         if(parms[0].BOLDSEC) {
	            if(!parms[0].SECSUM) 
	   	       fprintf(OUT,"FONT_REGION %d %d %d %d 3\n", startrel+1,nseq+2,endrel+1,nseq*2+1);
		    else
		       fprintf(OUT,"FONT_REGION %d %d %d %d 3\n", startrel+1,nseq+2,endrel+1,nseq+3);
		 }
	   }
	   if(!reliable[i]) { /* non reliable region */
	     startnonrel=i;
	     while(!reliable[i] && i<=nstamppos) i++;
	     endnonrel=i-1;
	     if(parms[0].SMALLSEQ) fprintf(OUT,"FONT_REGION %d %d %d %d 1\n",startnonrel+1,1,endnonrel+1,nseq);
	     if(parms[0].SMALLSEC) {
	        if(!parms[0].SECSUM)
	   	   fprintf(OUT,"FONT_REGION %d %d %d %d 1\n",startnonrel+1,nseq+2,endnonrel+1,nseq*2+1);
		else
		   fprintf(OUT,"FONT_REGION %d %d %d %d 1\n",startnonrel+1,nseq+2,endnonrel+1,nseq+3);
		}
	   }
	}
	fclose(OUT);
	printf(" ...done.\n You need to run ALSCRIPT now by typing `alscript %s'\n",tmp);
	printf(" A PostScript alignment will then be in %s_als.ps\n",parms[0].prefix);
	exit(0);
}
