#include <ver2hor.h>

/* Reads in a STAMP file, and an optional parameter file and
 *  produces a  horizontal alignment.  Pirated from DSTAMP.
 *
 * RBR September 1995 */

int main(int argc, char *argv[]) {

	int i,j,k,n;
	int nbloc,nseq;
	int nstamp,nstampseq,nstamppos;
	int startrel,endrel;
	int startnonrel,endnonrel;
	int Helix,Strand,Other,Gap;
	int bloclen;
	int count;
	int post;
	int ON;
	int user;
	int columns;

	int *reliable;
	int *reliable2;
	int *reliable3;
	int *reliable4;
	char *summary;
	int *gapped;

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
	parms[0].SECDISP=1;		/*  Display secondary structure */
	parms[0].SMALLSEQ=1;	/*  1 ==> put non-reliable regions in small */
	parms[0].SMALLSEC=1;	/*  as above for secondary structure */
	parms[0].CASESEQ=1;	/*  1 ==> non-reliable regions in lower case */
	parms[0].CASESEC=1;	/*  as above for secondary structure */
	parms[0].SECSUM=0;		/*  1 ==> display secondary structure  summary only */
	parms[0].VERBOSE=0;	/*  Run ALScript in silent mode */
	columns=70;


	n=argc;
	user=0;
	if(n<2) exit_error();
	for(i=1; i<n; ++i) {
           if(argv[i][0]!='-') exit_error();
	   if((i+1)<n && (argv[i+1][0]=='n' || argv[i+1][0]=='N' || 
		argv[i+1][0]=='T' || argv[i+1][0]=='t' || 
		argv[i+1][0]=='1' || argv[i+1][0]=='O' || 
		argv[i+1][0]=='o')) ON=1;
	   else ON=0;
	   for(j=1; j<strlen(argv[i]); ++j) argv[i][j]=utol(argv[i][j]);
           if(strcmp(&argv[i][1],"f")==0) {
	     if((i+1)>=n) exit_error();
	     strcpy(&parms[0].filename[0],argv[i+1]);
             i++;
           } else if(strcmp(&argv[i][1],"c")==0) {
             if((i+1)>=n) exit_error();
             parms[0].TYPE=argv[i+1][0];
             i++; user=1;
           } else if(strcmp(&argv[i][1],"t")==0) {
             if((i+1)>=n) exit_error();
             sscanf(argv[i+1],"%f",&parms[0].CUTOFF);
             i++; user=1;
	   } else if((strcmp(&argv[i][1],"columns")==0) || (strcmp(&argv[i][1],"col")==0)) {
	      if((i+1)>=n) exit_error();
             sscanf(argv[i+1],"%d",&columns);
	     i++;
           } else if(strcmp(&argv[i][1],"w")==0) {
             if((i+1)>=n) exit_error();
             sscanf(argv[i+1],"%d",&parms[0].WINDOW);
             i++; user=1;
	   } else if(strcmp(&argv[i][1],"cutoff")==0) {
		if((i+1)>=n) exit_error();
		sscanf(argv[i+1],"%f",&parms[0].CUTOFF);
		i++; user=1;
	   } else if(strcmp(&argv[i][1],"window")==0) {
		if((i+1)>=n) exit_error();
		sscanf(argv[i+1],"%d",&parms[0].WINDOW);
		i++; user=1;
	   } else if(strcmp(&argv[i][1],"type")==0) {
		if((i+1)>=n) exit_error();
		sscanf(argv[i+1],"%s",&parms[0].TYPE);
		i++; user=1;
	   } else if(strcmp(&argv[i][1],"smallseq")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SMALLSEQ=ON;
		i++;
	   } else if(strcmp(&argv[i][1],"smallsec")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SMALLSEC=ON;
		i++;
	   } else if(strcmp(&argv[i][1],"sec")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SECDISP=ON;
		i++;
	   } else if(strcmp(&argv[i][1],"secsum")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SECSUM=ON;
		i++;
	   } else if(strcmp(&argv[i][1],"verbose")==0) {
		if((i+1)>=n) exit_error();
		parms[0].VERBOSE=ON;
		i++;
	   } else if(strcmp(&argv[i][1],"caseseq")==0) {
		if((i+1)>=n) exit_error();
 		parms[0].CASESEQ=ON;
		i++;
	   } else if(strcmp(&argv[i][1],"casesec")==0) { 
		if((i+1)>=n) exit_error();
		parms[0].CASESEC=ON;
		i++;
           } else {
	     printf("unrecognised command: %s \n",&argv[i][1]);
             exit_error();
           }
        }
	
	if(columns<20) {
	  fprintf(stderr,"error: column value of %d is too narrow\n",columns);
	  exit(-1);
	}
	columns-=13;
	printf("VER2HOR R.B. Russell, 1995\n");
	printf(" Prints STAMP alignments in horizontal format\n");
	printf("  for quick viewing\n");
	if((IN=fopen(parms[0].filename,"r"))==NULL) {
	   fprintf(stderr,"error opening file %s\n",parms[0].filename);
	   exit(-1);
	}

	/* read the input file */
	post=0;
	printf(" Reading Alignment...\n");
	bloc=(struct seqdat*)malloc(MAX_N_SEQ*sizeof(struct seqdat));
	printf(" ");
	if(Agetbloc(IN,bloc,&nbloc)==-1) exit(-1); rewind(IN);

	bloclen=strlen(&bloc[1].seq[1]);
	summary=(char*)malloc(bloclen*sizeof(char));
	gapped=(int*)malloc(bloclen*sizeof(int));

	printf(" Getting STAMP information...\n");
	stampstuff=(struct stampdat*)malloc(MAX_STAMP_NUM*sizeof(struct stampdat));
	if(getstampdat(stampstuff,IN,&nstamp,&nstampseq,&nstamppos,bloclen)==-1) 
	    exit(-1);

	post=-1;
	for(i=0; i<nstamp; ++i) if(stampstuff[i].what=='B') post=i;

	/* determine the reliable regions */
        if((reliable4=getstamprel(stampstuff,nstamp,nstamppos,parms[0].TYPE,parms[0].CUTOFF,parms[0].WINDOW))==NULL) exit(-1);

	if((reliable=getstamprel(stampstuff,nstamp,nstamppos,'G',6.0,3))==NULL) exit(-1);
	if((reliable2=getstamprel(stampstuff,nstamp,nstamppos,'G',4.5,3))==NULL) exit(-1);
	if(post!=-1) {
	    if((reliable3=getstamprel(stampstuff,nstamp,nstamppos,'B',1,3))==NULL) exit(-1);
	}

	    
	printf(" %d STAMP fields read in for %d positions \n",nstamp,nstamppos);
	if(nstamppos!=strlen(&bloc[1].seq[1]) || nstampseq != nbloc) {
	   fprintf(stderr,"error: something wrong with STAMP file\n");
	   fprintf(stderr,"	  STAMP length is %d, Alignment length is %d\n",nstamppos,strlen(&bloc[1].seq[1]));
	   fprintf(stderr,"       STAMP nseq is %d, Alignment nseq is %d\n",nstampseq,nbloc);
	   exit(-1);
	}
	nseq=(nbloc-1)/2;

	for(i=0; i<bloclen; ++i) {
	  gapped[i]=0;
	  for(j=0; j<nseq; ++j) if(bloc[j+1].seq[i+1]==' ') gapped[i]=1;
	}
	for(i=0; i<nbloc; ++i) {
	  bloc[i+1].id[10]='\0';
	}
	
	/* Modify the sequence alignment as for DSTAMP */
	printf(" Processing the alignment...\n");
	for(i=0; i<nstamppos; ++i) {
	   for(j=0; j<nseq; ++j) 
	      if(parms[0].CASESEQ && !reliable[i]) 
		bloc[j+1].seq[i+1]=utol(bloc[j+1].seq[i+1]);
	      else 
		bloc[j+1].seq[i+1]=ltou(bloc[j+1].seq[i+1]);
	   if(parms[0].SECDISP) {
	      if(!parms[0].SECSUM) {
		 for(j=0; j<nseq; ++j) {
		   if(parms[0].CASESEC && !reliable[i]) 
		      bloc[nseq+j+2].seq[i+1]=utol(bloc[nseq+j+2].seq[i+1]);
		   else 
		      bloc[nseq+j+2].seq[i+1]=ltou(bloc[nseq+j+2].seq[i+1]);
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
		    summary[i]='H';
		 else if(Strand>Helix && Strand>Other && Strand>Gap)
		    summary[i]='E';
		 else if(Gap>Helix && Gap>Strand && Gap>Other)
		    summary[i]=' ';
		 else 
		   summary[i]='-';
	
		 if(parms[0].CASESEC && !reliable[i]) 
		   summary[i]=utol(summary[i]);
		 else 
		   summary[i]=ltou(summary[i]);
/*		 printf("Helix %d Strand %d Other %d Gap %d => %c \n",
		 	Helix, Strand, Other, Gap, summary[i]); */
	     }
	  }
	}
/*	printf("Summary %s\n",summary); */
	count=0;
	printf(" Output:\n");
	printf(" Very reliable => Pij' >=6 for stretches of >=3\n");
	printf(" Less reliable => Pij' >=4.5 for stretches of >=3\n");
	if(user) printf(" User reliable => STAMP field %c thresh %5.2f window %4d\n",
		parms[0].TYPE,parms[0].CUTOFF,parms[0].WINDOW);
	if(post!=-1) printf(" Post reliable => All Pij' > stamp_post parameter for stretches >=3\n");
	printf("\n");
	while(count<bloclen) {
	  /* Display a numbering scheme */
	  j=0;
	  printf("Number    ");
	  while(j<columns && (count+j)<bloclen) {
	      if(((j+count+1)%10)==0) {
		 printf("%4d",(j+count+1));
		 j+=4;
	      } else {
		printf(" ");
		j++;
	      }
	  }
	  printf("\n");
	
	  for(i=0; i<nseq; ++i) {
		if(strncmp(bloc[i+1].id,"space",5)!=0) printf("%10s   ",bloc[i+1].id);
		else printf("             ");
		j=0;
		while(j<columns && (count+j)<bloclen) {
		   printf("%c",bloc[i+1].seq[count+j+1]);
		   j++;
		}
		printf("\n");
	  }
	  printf("\n");
	  if(parms[0].SECDISP) {
	    if(parms[0].SECSUM) {
		j=0;
		printf("Sec. Summary ");
                while(j<columns && (count+j)<bloclen) {
                   printf("%c",summary[count+j]);
		   j++;
		}
		printf("\n");
	     } else {
	       for(i=0; i<nseq; ++i) {
                  if(strncmp(bloc[nseq+i+2].id,"space",5)!=0) printf("%10s   ",bloc[nseq+i+2].id);
                  else printf("             ");
                  j=0;
                  while(j<columns && (count+j)<bloclen) {
                     printf("%c",bloc[nseq+i+2].seq[count+j+1]);
                     j++;
                  }
                  printf("\n");
               }
	    }
	  }
	  printf("\nVery similar ");
	  j=0; while(j<columns && (count+j)<bloclen) {
                   if(gapped[count+j]==0) printf("%1d",reliable[count+j]);
                   else printf("-");
		   j++;
	  }
	  printf("\n");
	  printf("Less similar ");
	  j=0; while(j<columns && (count+j)<bloclen) {
		   if(gapped[count+j]==0) printf("%1d",reliable2[count+j]);
                   else printf("-");

                   j++;
          }
          printf("\n");
	  if(post!=-1) {
	   printf("Post similar ");
	   j=0; while(j<columns && (count+j)<bloclen) {
                   if(gapped[count+j]==0) printf("%1d",reliable3[count+j]);
		   else printf("-");
                   j++;
           }
           printf("\n");
	  }
	  if(user==1) {
           printf("User similar ");
           j=0; while(j<columns && (count+j)<bloclen) {
		   if(gapped[count+j]==0) printf("%1d",reliable4[count+j]);
                   else printf("-");
                   j++;
           }
           printf("\n");
          }

	  count+=j;
	  printf("\n");
	}
	
	free(summary);
	exit(0);
}
	
void exit_error()
{
	  fprintf(stderr,"format: ver2hor -f <STAMP file> \n");
 	  fprintf(stderr,"       -<parameter> <value>\n");
	  exit(-1);
}
