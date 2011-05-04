#include <stdio.h>
#include <msf2blc.h>

struct seqdat *RBR_getmsfalign(FILE *IN, int *nseq);


/* Rob's quick and dirty MSF to BLOCK conversion */
/* The format is msf2blc < MSF-file > BLOCK-file */

main() {


	int i,j;
	int nseq;
	int alignlen;

	struct seqdat *seq;


	seq=RBR_getmsfalign(stdin,&nseq);

	for(i=0; i<nseq; ++i) {
	   for(j=0; j<seq[i].slen; ++j) {
		if(seq[i].seq[j]=='.') seq[i].seq[j]=' ';
	   }
	}
	for(i=0; i<nseq; ++i) {
	   printf(">%s\n",seq[i].id);
	}
	printf("\n*\n");
	for(i=0; i<seq[0].slen; ++i) {
	   for(j=0; j<nseq; ++j) {
		printf("%c",seq[j].seq[i]);
	   }
	   printf("\n");
	}
	printf("*\n");
}



struct seqdat *RBR_getmsfalign(FILE *IN, int *nseq) {

	int i,j,k;
	int astart;
	int which;
	int nbits;

	char **bits;
	char buff[1000];
	struct seqdat *seq;

        seq=(struct seqdat*)malloc(sizeof(struct seqdat));

	(*nseq)=0;
	astart=0;

	while(fgets(buff,999,stdin)!=NULL) {
	  bits=RBR_c_split(buff, &nbits, ' ');
	  if(astart==0 && strstr(buff,"Name:")!=NULL /*  && strstr(buff,"Len:")!=NULL 
		&& strstr(buff,"Check:")!=NULL && strstr(buff,"Weight:")!=NULL */) { /* Name line */
	     seq[(*nseq)].id=(char*)malloc((strlen(bits[1])+1)*sizeof(char));
	     strcpy(seq[(*nseq)].id,bits[1]);
	     seq[(*nseq)].seq=(char*)malloc(sizeof(char));
	     seq[(*nseq)].slen=0;
	     (*nseq)++;
             seq=(struct seqdat*)realloc(seq,((*nseq)+1)*sizeof(struct seqdat));
	  } else if(astart==0 && strstr(buff,"//")!=NULL) { 
		astart=1; 
	  } else if(astart==1) {
	     which=-1;
	     for(i=0; i<(*nseq); ++i) {
		if(strcmp(seq[i].id,bits[0])==0) {
			which=i;
			break;
		}
	    }
	    if(which!=-1) {
		for(i=1; i<nbits; ++i) {
	           seq[which].seq=(char *)realloc(seq[which].seq,(seq[which].slen+strlen(bits[i])+2)*sizeof(char));
		   sprintf(&seq[which].seq[seq[which].slen],"%s",bits[i]);
		   seq[which].slen+=strlen(bits[i]);
		   seq[which].seq[seq[which].slen]='\0';
		}
	     }
	  }
	  for(i=0; i<nbits; ++i) free(bits[i]); 
	  free(bits);
	}

	return seq;

}
