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
#include <dstamp.h>

/***********************************************************************
getpars():
Read pairs of Command, Dimension values from the file and assign these 
values to the globally accessible variables.
Example parfile:

MAXSEQ 1000
MAX_SEQ_LEN 5000
etc.
-----------------------------------------------------------------------*/

int getpars(FILE *fp, struct parameters *var) {

    int i,T_FLAG;
    char c;
    char *parm;		/* name of following dimension */
    char *dim;		/* dimension */
	

   /* set default parameter values
    *  these are changed if a paramter file is specified */
    var[0].TYPE='G';  	/*  Use Pij' values */
    var[0].CUTOFF=6.0;	/*  Pij' >= 6.0 only */
    var[0].WINDOW=3;	/*  Stretches of three or more only */
    var[0].POINTSIZE=8.0;/*  Pointsize */
    var[0].SEC=1;	/*  Display secondary structure */
    var[0].BOXSEQ=1;	/*  Box secondary structre based on STAMP reliability */
    var[0].BOXSEC=1;	/*  Box secondary structre as above */
    var[0].SMALLSEQ=1;	/*  1 ==> put non-reliable regions in small */
    var[0].SMALLSEC=1;	/*  as above for secondary structure */
    var[0].BOLDSEQ=1;	/*  1 ==> reliable regions in bold */
    var[0].BOLDSEC=1;	/*  as above for secondary structure */ 
    var[0].SECSUM=0;	/*  1 ==> display a summary of DSSP only */
    var[0].VERBOSE=0;	/*  Run ALScript in silent mode */
    var[0].CASESEC=1;
    var[0].CASESEQ=1;

    parm = (char*)malloc(200*sizeof(char));
    dim  = (char*)malloc(200*sizeof(char));

    while(fscanf(fp,"%s%s",parm,dim) != (int)EOF){
	T_FLAG=(dim[0]=='y' || dim[0]=='Y' || dim[0]=='1' || dim[0]=='t' || dim[0]=='T');
	/* enables one to write '1', 'YES', 'Yes', 'yes', 'T_FLAG', 'True' or 'true' to 
	 *  set any boolean var[0]iable to one */
	if(strcmp(parm,"CUTOFF")==0)
		sscanf(dim,"%f",&var[0].CUTOFF);
	else if(strcmp(parm,"WINDOW")==0)
		sscanf(dim,"%d",&var[0].WINDOW);
	else if(strcmp(parm,"TYPE")==0)
		sscanf(dim,"%s",&var[0].TYPE);
	else if(strcmp(parm,"POINTSIZE")==0)
		sscanf(dim,"%f",&var[0].POINTSIZE);
	else if(strcmp(parm,"BOXSEQ")==0)
		var[0].BOXSEQ=T_FLAG;
	else if(strcmp(parm,"BOXSEC")==0)
		var[0].BOXSEC=T_FLAG;
	else if(strcmp(parm,"SMALLSEQ")==0)
		var[0].SMALLSEQ=T_FLAG;
	else if(strcmp(parm,"SMALLSEC")==0)
		var[0].SMALLSEC=T_FLAG;
	else if(strcmp(parm,"BOLDSEQ")==0)
		var[0].BOLDSEQ=T_FLAG;
	else if(strcmp(parm,"BOLDSEC")==0)
		var[0].BOLDSEC=T_FLAG;
	else if(strcmp(parm,"SEC")==0)
		var[0].SEC=T_FLAG;
	else if(strcmp(parm,"SECSUM")==0)
		var[0].SECSUM=T_FLAG;
	else if(strcmp(parm,"VERBOSE")==0)
		var[0].VERBOSE=T_FLAG;
	else if(strcmp(parm,"CASESEQ")==0)
 		var[0].CASESEQ=T_FLAG;
	else if(strcmp(parm,"CASESEC")==0) 
		var[0].CASESEC=T_FLAG;
	else{
    	    printf("Unrecognised Dimension Command\n");
	    printf("%s %s\n",parm,dim);
	    return -1;
	}
	while((c=getc(fp))!=(char)EOF && c!='\n'); /* read the end of the line, allows for comments */
	if(c==(char)EOF) break;
    }
    free(parm); free(dim); 
    return 0;

}
