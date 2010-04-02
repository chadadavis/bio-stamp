/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current contact addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Bioinformatics                      EMBL-European Bioinformatics Institute
 SmithKline Beecham Pharmaceuticals  Wellcome Trust Genome Campus
 New Frontiers Science Park (North)  Hinxton, Cambridge, CB10 1SD U.K.
 Harlow, Essex, CM19 5AW, U.K.       
 Tel: +44 1279 622 884               Tel: +44 1223 494 414
 FAX: +44 1279 622 200               FAX: +44 1223 494 468
 e-mail: russelr1@mh.uk.sbphrd.com   e-mail geoff@ebi.ac.uk
                                     WWW: http://barton.ebi.ac.uk/

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <stamp.h>

/* Takes a DSSP file, a description of a domain, and a transformation.
 *  It applies the transformation to the coordinates outputs the 
 *  corresponding porition of the DSSP file in DSSP format */
int extract_dssp(FILE *IN, struct brookn start, struct brookn end,int type,
	float **R, float *V, int startats, char chainlabel, FILE *OUT) {

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
		 fmatmult(R,V,&coord,1);
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
