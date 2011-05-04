/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)	            Prof. Geoffrey J. Barton (GJB)
 EMBL Heidelberg                            School of Life Sciences
 Meyerhofstrasse 1                          University of Dundee
 D-69117 Heidelberg                         Dow Street
 Germany                                    Dundee, DD1 5EH
                                          
 Tel: +49 6221 387 473                      Tel: +44 1382 345860
 FAX: +44 6221 387 517                      FAX: +44 1382 345764
 E-mail: russell@embl-heidelberg.de         E-mail geoff@compbio.dundee.ac.uk
 WWW: http://www.russell.emb-heidelberg.de  WWW: http://www.compbio.dundee.ac.uk

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "alignfit.h"

/* Reads in paramters from a specified file
 *
 * for example
 *
 *  WINDOW 3   		Window length of three
 *  THRESTATE yes	Use three states only
 *
 * N.B. Comments may only be one line long (no return characters)
 */

int getpars(FILE *fp, struct parameters *var) {

    int i,BOOLEAN;
    char c;
    char *parm;		/* name of following dimension */
    char *dim;		/* dimension */

	

   /* set default parameter values
    *  these are changed if a paramter file is specified */
    var[0].MAX_SEQ_LEN=3000;
    var[0].PAIRWISE=1;
    var[0].TREEWISE=1;
    var[0].OLDFORMAT=0;
    strcpy(&var[0].MATFILE[0],"alignfit.mat");
    strcpy(&var[0].TREEFILE[0],"alignfit.tree");
    strcpy(&var[0].ORDFILE[0],"alignfit.ord");
    strcpy(&var[0].TRANSFILE[0],"align.trans");

    parm = (char*)malloc(200*sizeof(char));
    dim  = (char*)malloc(200*sizeof(char));

    while(fscanf(fp,"%s%s",parm,dim) != (int)EOF){
	for(i=0; i<strlen(parm); ++i) parm[i]=ltou(parm[i]); /* change to upper case */
	BOOLEAN=(dim[0]=='Y' ||dim[0]=='y' || dim[0]=='1' || dim[0]=='T' || dim[0]=='t');
	/* enables one to write '1', 'YES', 'Yes', 'yes', 'BOOLEAN', 'True' or 'true' to 
	 *  set any boolean var[0]iable to one */
	if(strcmp(parm,"MAX_SEQ_LEN")==0)
		sscanf(dim,"%d",&var[0].MAX_SEQ_LEN);
	else if(strcmp(parm,"PAIRWISE")==0)
		var[0].PAIRWISE=BOOLEAN;
	else if(strcmp(parm,"TREEWISE")==0)
		var[0].TREEWISE=BOOLEAN;
	else if(strcmp(parm,"MATFILE")==0)
		strcpy(&var[0].MATFILE[0],dim);
	else if(strcmp(parm,"ORDFILE")==0)
		strcpy(&var[0].ORDFILE[0],dim);
	else if(strcmp(parm,"TREEFILE")==0)
		strcpy(&var[0].TREEFILE[0],dim);
	else if(strcmp(parm,"TRANSFILE")==0)
		strcpy(&var[0].TRANSFILE[0],dim);
	else if(strcmp(parm,"OLDFORMAT")==0)
		var[0].OLDFORMAT=BOOLEAN;
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
