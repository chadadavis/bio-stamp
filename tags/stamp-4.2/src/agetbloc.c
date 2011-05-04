/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Bioinformatics                      EMBL-European Bioinformatics Institute
 SmithKline Beecham Pharmaceuticals  Wellcome Trust Genome Campus
 New Frontiers Science Park (North)  Hinxton, Cambridge, CB10 1SD U.K.
 Harlow, Essex, CM19 5AW, U.K.       
 Tel: +44 1279 622 884               Tel: +44 1223 494 414
 FAX: +44 1279 622 200               FAX: +44 1223 494 468
 e-mail: russelr1@mh.uk.sbphrd.com   e-mail geoff@ebi.ac.uk
                                     WWW: http://barton.ebi.ac.uk/

 The WORK is Copyright (1997,1998) Robert B. Russell & Geoffrey J. Barton
	
	
	

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
	  idend = strchr(idstart,'\n');
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
	for(i=0; i<strlen(bloc[totseq].title); ++i) {
		if(bloc[totseq].title[i]=='\n') bloc[totseq].title[i]='\0';
	}

	bloc[totseq].seq = (char *) malloc(sizeof(char) * MAXslen);
        bloc[totseq].seq[0] = ' ';
	goto l1;
    } else if(sident){
	if((idstart = strchr(buff,'*')) != NULL){
	    i = 0;
	    while((buff = fgets(buff,MAXtlen,bfile)) != NULL){
		if(*idstart == '*'){
/*		    printf("Blocfile read: Length: %d\n",i); */
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
		if(i==MAXslen) {
			fprintf(stderr,"GETBLOC: Max Sequence length exceeded - use MAX_SEQ_LEN command to increase");
			exit(-1);
		}
		for(j=1;j<totseq+1;++j){
		    /*cope with short lines */
		    bloc[j].seq[i] = *bstart++;
		}
	    }
	    fprintf(stderr,"No terminating * in blocfile\n");
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
