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
#include "include.h"

/* Takes a DSSP file, a description of a domain, and a transformation.
 *  It applies the transformation to the coordinates outputs the 
 *  corresponding porition of the DSSP file in DSSP format */
int extract_dssp(IN,start,end,type,R,V,startats,chainlabel,OUT)
FILE *IN;
struct brookn start,end;
int type; /* 1 = all CA atoms in the file, 2 = single chain, 3 = specific start and end */
float **R;
float *V;
int startats;	/* 1 => write the transformation (ie. this is the first object in the file) */
char chainlabel;
FILE *OUT;
{
	int i,j,k;
	int begin,endnext;
	int startdat;
	int number;
	int count;

	char cid,in;
	char c;

	char *buff;
	char tmp[10];

	float *coord;

	struct brookn old,this;

	coord=(float*)malloc(3*sizeof(float));
	buff=(char*)malloc(200*sizeof(char));
	begin=endnext=count=startdat=0;
	
	/* get to start of data */
	while((c=getc(IN))!='#' && c!=(char)EOF) fprintf(OUT,"%c",c);
	fprintf(OUT,"#");
	if(c==(char)EOF) {
	  fprintf(stderr,"error: no starting # found in DSSP file \n");
	  return -1;
	}
	while((c=getc(IN))!='\n') fprintf(OUT,"%c",c); /* read till the end of line */
	fprintf(OUT,"\n");

	while((buff=fgets(buff,199,IN))!=NULL) {
	   buff[strlen(buff)-2]='\0';
	      /* get chain, number and insertion code */
	      this.cid=buff[11];
	      strncpy(tmp,&buff[6],4); tmp[4]='\0';
	      sscanf(tmp,"%d",&this.n);
	      this.in=buff[10];
	      if(endnext && type==3 && (this.cid!=old.cid || this.in!=old.in || this.n!=old.n)) 
		 break;
	      if(!begin && 
		 ((start.cid==this.cid && start.n==this.n && start.in==this.in && type==3) ||
		  (start.cid==this.cid && type==2) ||
		  (type==1) )) {
		  begin=1;
	      }
	      if(type==2 && start.cid!=this.cid && begin) break;
	      if(begin) {
		 for(i=0; i<3; ++i) {
		   strncpy(&tmp[0],&buff[107+i*7],7); tmp[7]='\0'; 
		   sscanf(tmp,"%f",&coord[i]);
		 }
		 buff[107]='\0';
		 /* transform coordinates */
/*		 if(count<10) 
		    printf("%d: %5.1f %5.1f %5.1f\n",count+1,coord[0],coord[1],coord[2]);  */
		 matmult(R,V,&coord,1);
/*		 if(count<10)
		    printf("%d: %8.4f %8.4f %8.4f\n",count+1,coord[0],coord[1],coord[2]);  */
		 count++;	
                 /* if a chain label is given write this in the file */
		 if(chainlabel!='\0') 
		    buff[12]=chainlabel;
		 fprintf(OUT,"%s",buff);
 		 fprintf(OUT,"%7.1f%7.1f%7.1f",coord[0],coord[1],coord[2]); 
/*		 fprintf(OUT,"%s",&buff[54]); 
		 for(i=0; i<(26-strlen(&buff[54])); ++i) fprintf(OUT," ");  */
		 fprintf(OUT,"\n");
	      } /* end of if(begin... */
	      if(begin && type==3 && end.cid==this.cid && end.n==this.n && end.in==this.in) 
	          endnext=1;
	        /* this residing after the last "if" makes the set of type==3 atoms inclusive */
	      old.cid=this.cid; old.n=this.n; old.in=this.in;
	} /* end of while((buff=... */
	free(buff);
	if(!begin) {
	   fprintf(stderr,"error: begin of sequence not found in DSSP file\n");
	   return -1;
	} else
	   return 0;
}
