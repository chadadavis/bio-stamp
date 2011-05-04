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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <alignfit.h>

/*
11th June 1992.
Agetbloc:  like getbloc, but does not require that every character read into
the seqs structure is an alphabetic character.  Also does not contain the 
option to convert the "sequences" read in into integer format


getbloc:  Read an AMPS style block file into the seqs array
the nbloc aligned sequences are stored in positions 1-nbloc.

This is a straight-ish translation of the fortran routine fbloc.f, hence
the non-C like goto's...

Sequence lengths are actual length +1 + 1.  This allows position 0 to 
be reserved for future use, and preserves the '\0' for output.

Personalised by RBR October 1992
*/

int Agetbloc(FILE *bfile, struct seqdat *bloc, int *nbloc) {

    int i,llen;
    char *buff;

    char *idstart, *idend, *bstart, sident = 0;
    int idlen,totseq = 0,k,j;
	
    buff = malloc(sizeof(char) * MAXtlen);

l1: 
    buff = fgets(buff,MAXtlen,bfile);
    if(buff == NULL){
	printf("Premature end of BLOCK FILE\n");
	return -1;
    }
    if((idstart = strchr(buff,'>')) != NULL){
	if(++totseq == MAXnbloc){
	    printf("Max Number of block file sequences exceeded: %d\n", totseq);
	    printf("Use MAX_NSEQ command to increase value");
	    return -1;
	}
	sident = 1;
	idend = strchr(idstart,' ');
	if(idend == NULL){
	  idend = strchr(idstart,'\0');
	}
	if(idend == NULL){
	  printf("Error reading identifier:%s\n",idstart);
	  printf("Exiting\n");
	}
	idlen = (idend - idstart) + 1;
	bloc[totseq].id = malloc(sizeof(char) * idlen);
	bloc[totseq].id = GJstrblank(bloc[totseq].id,idlen);
	strncpy(bloc[totseq].id,idstart+1,idlen-1);   /* don't copy the ">" symbol */
	bloc[totseq].ilen = idlen-1;
	bloc[totseq].id[idlen-1] = '\0';

	bloc[totseq].tlen = strlen(idend)+1;
	bloc[totseq].title = malloc(sizeof(char) * bloc[totseq].tlen);
	bloc[totseq].title = GJstrblank(bloc[totseq].title,bloc[totseq].tlen);
	strcpy(bloc[totseq].title,idend);

	bloc[totseq].seq = (char *) malloc(sizeof(char) * MAXslen);
        bloc[totseq].seq[0] = ' ';
	goto l1;
    } else if(sident){
	if((idstart = strchr(buff,'*')) != NULL){
	    i = 0;
	    while((buff = fgets(buff,MAXtlen,bfile)) != NULL){
		if(*idstart == '*'){
		    printf("Blocfile read: Length: %d\n",i);
		    ++i;
		    for(k=1;k<totseq+1;++k){
			bloc[k].slen = i;
			bloc[k].seq[i] = '\0';
			bloc[k].seq = realloc(bloc[k].seq,sizeof(char)*(i+1));
		    }
		    *nbloc = totseq;
		    free(buff);
		    return 0;
		}
		bstart = idstart;
		++i;
		if(i==MAXslen) printf("Max Sequence length exceeded - use MAX_SEQ_LEN command to increase");
		for(j=1;j<totseq+1;++j){
		    /*cope with short lines */
		    bloc[j].seq[i] = *bstart++;
		}
	    }
	    printf("No terminating * in blocfile\n");
	    return -1;
	} else goto l1;
    } else {
	goto l1;
    }
}

/*char *GJstrblank(string,len)
char *string;
int len;
{
  --len;
  string[len] = '\0';
  --len;
  while(len > -1){
    string[len] = ' ';
    --len;
  }
  return string;
}
  

*/ 
