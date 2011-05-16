#include <dstamp.h>

/* DSTAMP 96 - New version of DSTAMP that makes prettier alignments
 * Does a prettification of an alignment file with or without secondary structures
 *  residues to the nearest continuous segment 
 *
 * (c) R.B. Russell 1996
 */

main(int argc, char *argv[]) {

	int i,j;
	int nbloc,minlen,bloclen;
	int nstamp,nstamppos,nstampseq;
	int ndomain,gottrans;
	int nseq,nsec,nacc;
	int n_hydr, n_polar, n_cys, n_small, n_interesting;
	int n_pos, n_neg;
	int space;
	int n_buried, n_exposed, n_half;
	int total_hydr, total_polar, total_cys, total_small, total_cons;
	int total_ungapped;
	int ungapped;
	int ngaps;
	int ignore;
	int totally_conserved,aa_len;
	int start,end;
	int colour;
	int boxcys;
	int motif;
	int n_helix, n_strand;
	int startrel,endrel;
        int startnonrel,endnonrel;
	int WINDOW;
	int *reliable;

	int *naa;

	char c;
	char TYPE;
	char infile[200];
	char prefix[200];
	char filename[200];
	char value[200],keyword[200];

	float CUTOFF;

	FILE *BLOC,*OUT;

	struct seqdat *bloc;
	struct stampdat *stamp;
	struct domain_loc *domain;

	ignore=0;
	colour=0;
	boxcys=0;
	motif=0;
	TYPE = 'G';
	CUTOFF = 6.0;
	WINDOW = 3;

	strcpy(&prefix[0],"alscript");
	if(argc<3) exit_error();
	aa_len=strlen(RBR_AA1);
	naa=(int*)malloc(aa_len*sizeof(int));

	 for(i=1; i<argc; ++i) {
           if(argv[i][0]!='-') exit_error();
           strcpy(keyword,&argv[i][1]);
           if((i+1)<argc) strcpy(value,argv[i+1]);
           else strcpy(value,"none");
           for(j=0; j<strlen(keyword); ++j) 
              keyword[j]=ltou(keyword[j]); /* change to upper case */
           if(strcmp(&argv[i][1],"f")==0) { 
             if((i+1)>=argc) exit_error();
             strcpy(&infile[0],argv[i+1]);
             i++;
           } else if(strcmp(&argv[i][1],"prefix")==0) {
             if((i+1)>=argc) exit_error();
             /* assignment file name */
             strcpy(&prefix[0],argv[i+1]);
             i++;
           } else if(strcmp(&argv[i][1],"ignore")==0) {
             if((i+1)>=argc) exit_error();
             /* listfile name */
             sscanf(argv[i+1],"%d",&ignore);
             i++;
	   } else if(strcmp(&argv[i][1],"colour")==0) {
             /* listfile name */
	     colour=1;
	   } else if(strcmp(&argv[i][1],"boxcys")==0) {
             /* listfile name */
             boxcys=1;
	   } else if(strcmp(&argv[i][1],"motif")==0) {
	     motif=1;
	   } else if(strcmp(&argv[i][1],"c")==0) {
             if((i+1)>=argc) exit_error();
             TYPE=argv[i+1][0];
             i++;
           } else if(strcmp(&argv[i][1],"t")==0) {
             if((i+1)>=argc) exit_error();
             sscanf(argv[i+1],"%f",&CUTOFF);
             i++;
           } else if(strcmp(&argv[i][1],"w")==0) {
             if((i+1)>=argc) exit_error();
             sscanf(argv[i+1],"%d",&WINDOW);
             i++;
           } else {
             exit_error();
           }
        }

	if((BLOC=fopen(infile,"r"))==NULL) {
	  fprintf(stderr,"error opening file %s\n",infile);
	  exit(-1);
	}
	printf("Reading block file\n");
	nbloc=0;
	while((c=getc(BLOC))!=(char)EOF) nbloc+=(c=='>');
	rewind(BLOC);
	bloc=(struct seqdat*)malloc((nbloc+1)*sizeof(struct seqdat));
	if(Agetbloc(BLOC,bloc,&nbloc)==-1) exit(-1);
	rewind(BLOC);
	bloclen=strlen(&bloc[1].seq[1]);
	printf("Searching for STAMP data\n");
	nstamp=0;
	while((c=getc(BLOC))!=(char)EOF) nstamp+=(c=='#');
	
	for(i=0; i<nbloc; ++i) {
	   for(j=0; j<strlen(bloc[i+1].id); ++j) if(bloc[i+1].id[j]=='\n') bloc[i+1].id[j]='\0';
	}
	rewind(BLOC);
	if(nstamp>0) { /* Get STAMP data if there */
	   stamp=(struct stampdat*)malloc(nstamp*sizeof(struct stampdat));
           if(getstampdat(stamp,BLOC,&nstamp,&nstampseq,&nstamppos,bloclen)==-1) exit(-1);

	   if((reliable=getstamprel(stamp,nstamp,nstamppos,TYPE,CUTOFF,WINDOW))==NULL) exit(-1); 
	}

	fclose(BLOC);
		
	printf("Block file contains %d sequences and %d stamp fields the alignment length is %d\n",
	       nbloc,nstamp,strlen(&bloc[1].seq[1]));

	sprintf(filename,"%s.als",prefix);
	if((OUT=fopen(filename,"w"))==NULL) {
	  fprintf(stderr,"error opening output file %s\n",filename);
	  exit(-1);
	}
	printf("Will now output ALSCRIPT file to %s.als\n",prefix);
	fprintf(OUT,"SILENT_MODE\n");
	fprintf(OUT,"BLOCK_FILE %s\n",infile);
	fprintf(OUT,"OUTPUT_FILE %s.ps\n",prefix);
	fprintf(OUT,"LANDSCAPE\n");
	fprintf(OUT,"POINTSIZE  8\n");
	fprintf(OUT,"IDENT_WIDTH 12\n");
	fprintf(OUT,"DEFINE_FONT 0 Helvetica      DEFAULT \n");
	fprintf(OUT,"DEFINE_FONT 1 Helvetica REL  0.75   \n");
	fprintf(OUT,"DEFINE_FONT 7 Helvetica REL 0.5\n");
	fprintf(OUT,"DEFINE_FONT 3 Helvetica-Bold DEFAULT    \n");
	fprintf(OUT,"DEFINE_FONT 4 Times-Bold     DEFAULT   \n");
	fprintf(OUT,"DEFINE_FONT 5 Helvetica-BoldOblique  DEFAULT \n");
	fprintf(OUT,"DEFINE_COLOUR 7  1 0 1\n");
	fprintf(OUT,"DEFINE_COLOUR 8  0 0 1\n");
	fprintf(OUT,"DEFINE_COLOUR 9  0 1 0\n");
	fprintf(OUT,"DEFINE_COLOUR 10 1 1 0\n");

	fprintf(OUT,"NUMBER_INT 10\n");
	fprintf(OUT,"SETUP\n");	
	fprintf(OUT,"#\n#\n");

	/* First work out how many sequences there are */
	nseq=0; nsec=0; nacc=0;
	while(nseq<nbloc && strncmp(bloc[nseq+1].id,"space",5)!=0 && strncmp(bloc[nseq+1].id,"Sec",3)!=0) nseq++;
	if(nbloc>nseq && strncmp(bloc[nseq+1].id,"space",5)==0) { space=nseq+1; }
	else { space=-1; }
	
	/* Conservation of residues */
	fprintf(OUT,"#\n#\n# Residues property conservation\n");
	fprintf(OUT,"# Totally conserved (or nearly) non-hydrophobic residues - inverse text or green\n");
	fprintf(OUT,"# Conserved hydrophobics - shaded grey or colour yellow \n");	
	fprintf(OUT,"# Conserved polar - bold\n");
	fprintf(OUT,"# Conserved small - small \n");
	if(boxcys) fprintf(OUT,"# Conserved cysteines - boxed \n");
	bloclen=strlen(&bloc[1].seq[1]);
        total_hydr=total_polar=total_cys=total_small=total_cons=0;
	for(i=0; i<bloclen; ++i) { 
	  n_hydr=0; n_polar=0; n_cys=0; n_small=0; n_interesting=0; n_neg=0; n_pos=0;
	  for(j=0; j<aa_len; ++j) naa[j]=0;
	  ungapped=1; ngaps=0;
	  for(j=0; j<nseq; ++j) {
		if(bloc[j+1].seq[i+1]==' ') { ungapped=0; ngaps++; }
		naa[(int)(bloc[j+1].seq[i+1]-'A')]++;
		if(bloc[j+1].seq[i+1]=='A' || bloc[j+1].seq[i+1]=='C' || bloc[j+1].seq[i+1]=='F' ||
		   bloc[j+1].seq[i+1]=='I' || bloc[j+1].seq[i+1]=='L' || bloc[j+1].seq[i+1]=='M' ||
		   bloc[j+1].seq[i+1]=='V' || bloc[j+1].seq[i+1]=='W' || bloc[j+1].seq[i+1]=='Y' ||
		   bloc[j+1].seq[i+1]=='G' || bloc[j+1].seq[i+1]=='H' || bloc[j+1].seq[i+1]=='S' ||
		   bloc[j+1].seq[i+1]=='T') n_hydr++;
		if(bloc[j+1].seq[i+1]=='R' || bloc[j+1].seq[i+1]=='K' || bloc[j+1].seq[i+1]=='H') n_pos++;
	        if(bloc[j+1].seq[i+1]=='E' || bloc[j+1].seq[i+1]=='D') n_neg++;
		if(bloc[j+1].seq[i+1]=='A' || bloc[j+1].seq[i+1]=='C' || bloc[j+1].seq[i+1]=='D' ||
                   bloc[j+1].seq[i+1]=='E' || bloc[j+1].seq[i+1]=='G' || bloc[j+1].seq[i+1]=='H' ||
                   bloc[j+1].seq[i+1]=='K' || bloc[j+1].seq[i+1]=='N' || bloc[j+1].seq[i+1]=='P' ||
		   bloc[j+1].seq[i+1]=='Q' || bloc[j+1].seq[i+1]=='R' || bloc[j+1].seq[i+1]=='S' ||
		   bloc[j+1].seq[i+1]=='T')  n_polar++;
	        if(bloc[j+1].seq[i+1]=='C') n_cys++;
		if(bloc[j+1].seq[i+1]=='A' || bloc[j+1].seq[i+1]=='C' || bloc[j+1].seq[i+1]=='G' ||
                    bloc[j+1].seq[i+1]=='S' || bloc[j+1].seq[i+1]=='T' || bloc[j+1].seq[i+1]=='P' ||
		    bloc[j+1].seq[i+1]=='D') n_small++;
		if(bloc[j+1].seq[i+1]=='C'  || bloc[j+1].seq[i+1]=='D' || bloc[j+1].seq[i+1]=='E' ||
                    bloc[j+1].seq[i+1]=='H' || bloc[j+1].seq[i+1]=='K' || bloc[j+1].seq[i+1]=='N' ||
                    bloc[j+1].seq[i+1]=='Q' || bloc[j+1].seq[i+1]=='R' || bloc[j+1].seq[i+1]=='S' ||
                    bloc[j+1].seq[i+1]=='T' || bloc[j+1].seq[i+1]=='W' || bloc[j+1].seq[i+1]=='Y' ||
		    bloc[j+1].seq[i+1]=='G')
					 n_interesting++;
		
	   }
	 /* Re-assess ungapped in light of 'ingnore' */
	 if((ngaps-ignore)<=0) ungapped=1;
	 if(ungapped) {
	   total_ungapped++;
	   totally_conserved=0;
	   for(j=0; j<aa_len; ++j) {
		if(naa[j]>=(nseq-ignore)) {  /* Ignore not considered here */
			totally_conserved=1;
			break;
		}
	   }
	   if(space!=-1 && motif==1) {
	     if(totally_conserved==1) {
	         fprintf(OUT,"SUB_CHARS %d %d %d %d SPACE %c\n",i+1,space,i+1,space,RBR_AA1[j]);
	     } else if(n_pos>=(nseq-ignore)) {
		 fprintf(OUT,"SUB_CHARS %d %d %d %d SPACE +\n",i+1,space,i+1,space);
	     } else if(n_neg>=(nseq-ignore)) {
		 fprintf(OUT,"SUB_CHARS %d %d %d %d SPACE -\n",i+1,space,i+1,space);
	     } else if(n_small>=(nseq-ignore)) {
		 fprintf(OUT,"SUB_CHARS %d %d %d %d SPACE s\n",i+1,space,i+1,space);
	     } else if(n_hydr>=(nseq-ignore)) {
		 fprintf(OUT,"SUB_CHARS %d %d %d %d SPACE h\n",i+1,space,i+1,space);
	     } else if(n_polar>=(nseq-ignore)) {
		 fprintf(OUT,"SUB_CHARS %d %d %d %d SPACE p\n",i+1,space,i+1,space);
	     }
	   }
		
	   
	     
	

		


	   if(!colour) {	
	     if(boxcys && n_cys>=(nseq-ignore)) { /* conserved cysteine  boxed */
	        fprintf(OUT,"BOX_REGION %d 1 %d %d \n",i+1,i+1,nseq);
		total_cys++;
  	     } else if(totally_conserved && n_interesting>ignore && n_interesting>3) { /* conserved and interesting */
	        fprintf(OUT,"INVERSE_CHARS ABCDEFGHIJKLMNOPQRSTUVWXYZ %d 1 %d %d \n",i+1,i+1,nseq);
	        fprintf(OUT,"SHADE_REGION %d 1 %d %d 0.00\n",i+1,i+1,nseq);
		total_cons++;
	     } else if(n_small>=(nseq-ignore)) { /* conserved small  */
	       fprintf(OUT,"FONT_REGION %d 1 %d %d 7 \n",i+1,i+1,nseq);
	       total_small++;
	     } else if(n_hydr>=(nseq-ignore)) { /* conserved hydrophobic shaded grey */
	        fprintf(OUT,"SHADE_CHARS ABCDEFGHIJKLMNOPQRSTUVWXYZ %d 1 %d %d 0.90\n",i+1,i+1,nseq);
		total_hydr++;
	     } else if(n_polar>=(nseq-ignore)) { /* conserved polar bold */
	        fprintf(OUT,"FONT_REGION %d 1 %d %d 3 \n",i+1,i+1,nseq);
		total_polar++;
	     }
	  } else {
	     if(boxcys && n_cys>=(nseq-ignore)) { /* conserved cysteine boxed (still) */
                fprintf(OUT,"BOX_REGION %d 1 %d %d \n",i+1,i+1,nseq);
		total_cys++;
             } else if(totally_conserved) { /* conserved residue green */
                fprintf(OUT,"COLOUR_REGION %d 1 %d %d 9 \n",i+1,i+1,nseq);
		total_cons++;
             } else if(n_small>=(nseq-ignore)) { /* conserved small (still small)  */
               fprintf(OUT,"FONT_REGION %d 1 %d %d 7 \n",i+1,i+1,nseq);
		total_small++;
             } else if(n_hydr>=(nseq-ignore) && n_polar<(nseq/2)) { /* conserved hydrophobics yellow */
                fprintf(OUT,"COLOUR_REGION %d 1 %d %d 10\n",i+1,i+1,nseq);
		total_hydr++;
             } else if(n_polar>=(nseq-ignore) && n_hydr<(nseq/2)) { /* conserved polar bold  (still) */
                fprintf(OUT,"FONT_REGION %d 1 %d %d 3 \n",i+1,i+1,nseq);
		total_polar++;
             }
	  }
	 }
	}
	fprintf(OUT,"#\n#\n# Summary of residue conservation\n");
	if(boxcys) fprintf(OUT,"# totally conserved CYS : %d\n",total_cys);
	fprintf(OUT,"# Conserved non-hydrophobic residues : %d\n",total_cons);
	fprintf(OUT,"# Conserved small residues: %d\n",total_small);
	fprintf(OUT,"# Conserved hydrophobic residues: %d\n",total_hydr);
	fprintf(OUT,"# Conserved polar residues: %d\n",total_polar);
	fprintf(OUT,"# Total number of un-gapped positions: %d \n",total_ungapped);
	fprintf(OUT,"#\n#\n");

	/* Now do the secondary structures if there are any */
	for(i=0; i<nbloc; ++i) {
	 /* check for secondary structure assignment */
	  
	 if((strstr(bloc[i+1].id,"Sec")!=NULL)     || (strstr(bloc[i+1].id,"SEC")  !=NULL) ||
	    (strstr(bloc[i+1].id,"sec")    !=NULL) || (strstr(bloc[i+1].id,"dssp") !=NULL) || 
	    (strstr(bloc[i+1].id,"DSSP")   !=NULL) || (strstr(bloc[i+1].id,"SOPMA")!=NULL) || 
	    (strstr(bloc[i+1].id,"SSPRED") !=NULL) || (strstr(bloc[i+1].id,"sopma")!=NULL) || 
	    (strstr(bloc[i+1].id,"sspred") !=NULL) || (strstr(bloc[i+1].id,"PHD")  !=NULL) || 
	    (strstr(bloc[i+1].id,"phd_sec")!=NULL) || (strstr(bloc[i+1].id,"Pred__SS")!=NULL) ||
            (strstr(bloc[i+1].id,"Clean_SS")!=NULL)|| (strstr(bloc[i+1].id,"Consensus")!=NULL)) { 
	   /* Secondary structure assignment */
	   /*  Change the secondary structure assignment to threestate */
	   nsec++;
	   threestate(&bloc[i+1].seq[1],"HG3A","EB","STIc");
	   fprintf(OUT,"#\n# Secondary structures for %s\n",bloc[i+1].id);
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d H SPACE\n",i+1,bloclen,i+1);
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d G SPACE\n",i+1,bloclen,i+1);
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d E SPACE\n",i+1,bloclen,i+1); 
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d B SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d S SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d T SPACE\n",i+1,bloclen,i+1); 
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d I SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d - SPACE\n",i+1,bloclen,i+1);
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d C SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d c SPACE\n",i+1,bloclen,i+1);


 	   n_helix=n_strand=0;
	   for(j=0; j<bloclen; ++j) {
		if(bloc[i+1].seq[j+1]=='H') {
		   n_helix++;
/*	           printf("Helix from %d - ",j+1); */
		   start=j+1;
		   while((bloc[i+1].seq[j+1]=='H' || bloc[i+1].seq[j+1]=='G') && j<(bloclen-1)) j++;
		   end=j;
/*		   printf("%d\n",j+1); */
		   fprintf(OUT,"HELIX %d %d %d\n",start,i+1,end);
		   if(colour) fprintf(OUT,"COLOUR_TEXT_REGION %d %d %d %d 8\n",start,i+1,end,i+1);
		}
	        if(bloc[i+1].seq[j+1]=='B') {
		   n_strand++;
/*		   printf("Strand from %d - ",j+1); */
		   start=j+1;
                   while((bloc[i+1].seq[j+1]=='B') && j<(bloclen-1)) j++;
		   end=j;
                   fprintf(OUT,"STRAND %d %d %d\n",start,i+1,end);
		   if(colour) fprintf(OUT,"COLOUR_TEXT_REGION %d %d %d %d 7\n",start,i+1,end,i+1);
/*		   printf("%d\n",j+1); */
                }
	   }
	   fprintf(OUT,"#  %s %d helices and %d strands\n",bloc[i+1].id,n_helix,n_strand);
	 }
         /* check for burial/accessibility data */
	 if(strstr(bloc[i+1].id,"acc")!=NULL || strstr(bloc[i+1].id,"ACC")!=NULL ||
           strstr(bloc[i+1].id,"burial")!=NULL || strstr(bloc[i+1].id,"BURIAL")!=NULL) {
	   fprintf(OUT,"#\n# Accessibility data for %s\n",bloc[i+1].id);
           threestate(&bloc[i+1].seq[1],"bB0123","uiUI456","eE789");  /* Use three-state H=buried B=half c=exposed */
	   n_buried = 0; n_half = 0; n_exposed = 0;
	   nacc++;
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d b SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d B SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d i SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d u SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d I SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d U SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d E SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d e SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 0 SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 1 SPACE\n",i+1,bloclen,i+1);
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d 2 SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 3 SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 4 SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 5 SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 6 SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 7 SPACE\n",i+1,bloclen,i+1);
           fprintf(OUT,"SUB_CHARS 1 %d %d %d 8 SPACE\n",i+1,bloclen,i+1);
	   fprintf(OUT,"SUB_CHARS 1 %d %d %d 9 SPACE\n",i+1,bloclen,i+1);
	   for(j=0; j<bloclen; ++j) {
		fprintf(OUT,"BOX_REGION %d %d %d %d\n",j+1,i+1,j+1,i+1);
                if(bloc[i+1].seq[j+1]=='H') { /* Buried, filled square */
		  fprintf(OUT,"INVERSE_CHARS ABCDEFGHIJKLMNOPQRSTUVWXYZ %d %d %d %d \n",j+1,i+1,j+1,i+1);
                  fprintf(OUT,"SHADE_REGION %d %d %d %d 0.00\n",j+1,i+1,j+1,i+1);
		  
		  n_buried++;
		} else if(bloc[i+1].seq[j+1]=='B') { /* Half-buried, shaded square */
                  fprintf(OUT,"SHADE_REGION %d %d %d %d 0.90\n",j+1,i+1,j+1,i+1);
		  n_half++;
                }  else { /* otherwise leave blank = exposed */
		  n_exposed++;
	  	}
	    }
	    fprintf(OUT,"#\n#A total of %d buried, %d half-buried and %d exposed\n#\n#\n",n_buried,n_half,n_exposed);
	  }
	}

	 while(i<=nstamppos) {
           if(reliable[i]) { /* reliable region */
              startrel=i;
              while(reliable[i] && i<=nstamppos) i++;
              endrel=i-1;
              fprintf(OUT,"BOX_REGION %d %d %d %d\n", startrel+1,1,endrel+1,nseq);
           }
	   i++;
	 }
	  

	fclose(OUT);

	printf("There are %d sequences %d sec-strucs and %d accessibilities out of %d block file entries\n",nseq,nsec,nacc,nbloc);
	exit(0);

}
void exit_error() {
	  fprintf(stderr,"format: dstamp -f <block file> -prefix <output prefix> -ignore <no. of res ignored> \n");
	  exit(-1);
}
