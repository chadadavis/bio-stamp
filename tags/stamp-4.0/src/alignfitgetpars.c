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

int getpars(fp,var)
FILE *fp;
struct parameters *var;
{
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
