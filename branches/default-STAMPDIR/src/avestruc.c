#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <stamp.h>
#include <stamprel.h>

#define MAX_SEQ_LEN 10000
#define PRECISION 1000

void exit_error();

/* avestruc: reads in a stamp output file and generates an average structure file */
int main(int argc, char *argv[]) {

	int i,j,k,l,test;
	int ndomain,total,add,nbloc;
	int nstamp,nseq,nstamppos;
	int gottrans;
	int bloclen;
	int stampwindow;
	int polyA,nrel;
	int ncoords,nres;
	int all_aligned,sidechain;
	int ident,identical;
	int cons,conserved;
	int AApt;
	int npolar, nhydrophobic, naromatic, nsmall, ntiny;
        int npositive, nnegative, ncharged, naliphatic, nbranch;
	int ignore;


	int *pointer;
	int *counter;
	int *all;
	int *rel,*realrel;
	int **avecoords;
	int ***n,***c,***o,***cb;
	

	char filename[200],outfile[200];
	static char *AA3[]={
           "ALA", "ASX", "CYS", "ASP", "GLU", /* 4 */
           "PHE", "GLY", "HIS", "ILE", "UNK", /* 9 */
           "LYS", "LEU", "MET", "ASN", "UNK", /* 14 */
           "PRO", "GLN", "ARG", "SER", "THR", /* 19 */
           "UNK", "VAL", "TRP", "UNK", "TYR", /* 24 */
           "GLX", "SMA", "TIN", "POL", "HYD", /* 29  small, tiny, polar, hydrophobic */
           "POS", "NEG", "CHA", "ARO", "ALI", /* 34  positive, negative, charged, aromatic, aliphatic */
           "BRA"  			      /* 39  branched */
        };


	char stampchar;

    char *stampdir = AM_STAMPDIR;
	char *buff;
	char *aa;
	char *N,*C,*O;

	float stampthresh;
	float *average;

	FILE *TRANS,*PDB,*OUT;

	struct domain_loc *domain;
	struct seqdat *bloc;
	struct stampdat *stamp;
	struct brookn *numb;
	struct side_chain **side; 

	if(argc<3) exit_error();
	
	/* defaults */
	polyA=0; sidechain=0;
	stampchar='G';
	stampwindow=3;
	stampthresh=6.0;
	strcpy(&outfile[0],"average.pdb");
	N=(char*)malloc(5*sizeof(char));
	C=(char*)malloc(5*sizeof(char));
	O=(char*)malloc(5*sizeof(char));
	strcpy(N," N  ");
	strcpy(C," C  ");
	strcpy(O," O  ");
	
	average=(float*)malloc(3*sizeof(float));
	buff=(char*)malloc(100*sizeof(char));
	all_aligned=0;
	ident=0;
	cons=0;
	ignore=0;

	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(strcmp(&argv[i][1],"f")==0) {
	      if((i+1)>=argc) exit_error();
	      strcpy(&filename[0],argv[i+1]);
	      i++;
	   } else if(strcmp(&argv[i][1],"c")==0) {
	     if((i+1)>=argc) exit_error();
             stampchar=argv[i+1][0];
             i++;
 	   } else if(strcmp(&argv[i][1],"t")==0) {
	     if((i+1)>=argc) exit_error();
	     sscanf(argv[i+1],"%f",&stampthresh);
	     i++;
	   } else if(strcmp(&argv[i][1],"w")==0) {
	     if((i+1)>=argc) exit_error();
             sscanf(argv[i+1],"%d",&stampwindow);		
             i++;
	   } else if(strcmp(&argv[i][1],"polyA")==0) {
	     polyA=1;
	   } else if(strcmp(&argv[i][1],"ident")==0) {
             ident=1;
	   } else if(strcmp(&argv[i][1],"cons")==0) {
             cons=1;
	   } else if(strcmp(&argv[i][1],"ignore")==0) {
	     if((i+1)>=argc) exit_error();
	     sscanf(argv[i+1],"%d",&ignore);
	     i++;
	   } else if(strcmp(&argv[i][1],"o")==0) {
	     if((i+1)>=argc) exit_error();
             strcpy(&outfile[0],argv[i+1]);
             i++;
	   } else if(strcmp(&argv[i][1],"aligned")==0 || strcmp(&argv[i][1],"alignment")==0) {
	     all_aligned=1;
	   } else {
	     exit_error();
	   }
	}

    if(getenv("STAMPDIR")!=NULL) {
      /* Allow environment variable to override config setting */
      stampdir=getenv("STAMPDIR");
    }
	/* read in coordinate locations and initial transformations */
	if((TRANS = fopen(filename,"r")) == NULL) {
	   fprintf(stderr,"error: file %s does not exist\n",filename);
	   exit(-1);
	}
	if((OUT=fopen(outfile,"w"))==NULL) {
	  fprintf(stderr,"error opening file %s\n",outfile);
	  exit(-1);
	}
	printf("AVESTRUC R.B. Russell, 1995\n");
	printf(" Using STAMP file %s\n",filename);
	if(polyA) printf(" Will make an averaged poly-alanine structure\n");
	else printf(" Will make an averaged C-alpha-only structure\n");
	printf(" Averaged coordinates will be in %s\n",outfile);
	fprintf(OUT,"REMARK  Averaged coordinates \n");
	fprintf(OUT,"REMARK  Calculated using the transformations and alignments\n");
	fprintf(OUT,"REMARK   in the STAMP file %s\n",filename);
	fprintf(OUT,"REMARK   R. B. Russell, 1994\n");
	fprintf(OUT,"REMARK   \n");

	/* determine the number of domains specified */
	printf(" Reading domain descriptors...\n");
	ndomain=count_domain(TRANS);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	rewind(TRANS);
	if(getdomain(TRANS,domain,&ndomain,ndomain,&gottrans,stampdir,0,stdout)==-1) exit(-1);
	for(i=0; i<ndomain; ++i) {
	   printf("Domain %s\n",domain[i].id);
	   printdomain(stdout,domain[i],1);
	}
	pointer=(int*)malloc(ndomain*sizeof(int));
	if(polyA) {
	  n=(int***)malloc(ndomain*sizeof(int**));
	  c=(int***)malloc(ndomain*sizeof(int**));
	  o=(int***)malloc(ndomain*sizeof(int**));
	  cb=(int***)malloc(ndomain*sizeof(int**));
	  aa=(char*)malloc(MAX_SEQ_LEN*sizeof(char));
	  numb=(struct brookn*)malloc(MAX_SEQ_LEN*sizeof(struct brookn));
	}
	if(sidechain) {
	  side=(struct side_chain**)malloc(ndomain*sizeof(struct side_chain*));
	}


	for(i=0; i<ndomain; ++i) {
	  if(polyA) {
	    n[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	    o[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	    c[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	    cb[i]=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	  }
	  if(sidechain) {
	     side[i]=(struct side_chain*)malloc(MAX_SEQ_LEN*sizeof(struct side_chain));
	  }

	  for(j=0; j<domain[i].nobj; ++j) {
	     if(domain[i].start[j].cid=='_') domain[i].start[j].cid=' ';
	     if( domain[i].start[j].in=='_') domain[i].start[j].in=' ';
	     if(  domain[i].end[j].cid=='_') domain[i].end[j].cid=' '; 
	     if(   domain[i].end[j].in=='_') domain[i].end[j].in=' ';
	  }
	  fprintf(OUT,"REMARK Domain %4d Name is %s\n",i+1,domain[i].id);
	}

	/* read in the alignment */
	printf(" Reading alignment...\n");
	rewind(TRANS);
	bloc=(struct seqdat*)malloc((ndomain*2+2)*sizeof(struct seqdat));
	printf(" ");
	if(Agetbloc(TRANS,bloc,&nbloc)==-1) {
           fprintf(stderr,"error: alignment file %s appears to be incomplete\n",filename);                     
            exit(-1); 
        }
	bloclen=strlen(&bloc[1].seq[1]);
	if(nbloc!=(ndomain*2+1)) {
	   fprintf(stderr,"error: number of domains and alignment disagree\n");
	   exit(-1);
	}
	counter=(int*)malloc(bloclen*sizeof(int));
	all=(int*)malloc(bloclen*sizeof(int));

	fprintf(OUT,"REMARK   \n");
	/* read in the STAMP output and get reliable regions */
	rewind(TRANS);
	stamp=(struct stampdat*)malloc(100*sizeof(struct stampdat));
	if(getstampdat(stamp,TRANS,&nstamp,&nseq,&nstamppos,MAX_SEQ_LEN)==-1) exit(-1);
	if((rel=getstamprel(stamp,nstamp,bloclen,stampchar,stampthresh,stampwindow))==NULL) exit(-1);
	realrel=(int*)malloc(bloclen*sizeof(int));

	nrel=0;
        for(i=0; i<bloclen; ++i) {
             if(rel[i]==1) nrel++;
	     realrel[i]=rel[i];
	}

	if(all_aligned) {
	   for(i=0; i<bloclen; ++i) {
	     k=0;
	     for(j=0; j<ndomain; ++j) {
		if(bloc[j+1].seq[i+1]!=' ') k++;
	     }
	     if(k==ndomain) rel[i]=1;
	     else rel[i]=0;
	   }
	}
	fprintf(OUT,"REMARK %4d STAMP reliable or core positions found\n",nrel);
	fprintf(OUT,"REMARK  (STAMP parms '%c' %5.2f %2d) \n",stampchar,stampthresh,stampwindow);
	fprintf(OUT,"REMARK The characters to the far right of the coordinates show\n");
	fprintf(OUT,"REMARK  the residues occuring at each position within the alignment\n");
	if(all_aligned) {
	  fprintf(OUT,"REMARK  All aligned positions were used to generate the\n");
	  fprintf(OUT,"REMARK   coordinates below\n");
	}
	
	fclose(TRANS);

	printf(" Reading coordinates...\n");
	for(i=0; i<ndomain; ++i) {
	   /* output the domain descriptors again */
	   printf(" Domain %3d %s %s\n   ",i+1,domain[i].filename,domain[i].id); 
	   if((PDB=openfile(domain[i].filename,"r"))==NULL) {
	      fprintf(stderr,"error: file %s does not exist\n",domain[i].filename);
	      exit(-1);
	   }
	   domain[i].ncoords=0;
	   domain[i].coords=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	   domain[i].aa=(char*)malloc((MAX_SEQ_LEN+1)*sizeof(char)); 
	   domain[i].numb=(struct brookn*)malloc((MAX_SEQ_LEN)*sizeof(struct brookn));
	   total=0;
	   printf("    "); 
	   for(j=0; j<domain[i].nobj; ++j) {
	       if(igetca(PDB,&domain[i].coords[total],&domain[i].aa[total],&domain[i].numb[total],
		  &add,domain[i].start[j],domain[i].end[j],
		  domain[i].type[j],(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) {
		    fprintf(stderr,"Error in domain %s object %d \n",domain[i].id,j+1);
		    exit(-1);
	       }
	       if(polyA) {
		  closefile(PDB,domain[i].filename);
	          PDB=openfile(domain[i].filename,"r");
		  if(igetgen(PDB,&n[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],domain[i].type[j],N,(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
		  if(test!=add) {
			fprintf(stderr,"error: different number of N atoms in %s\n",domain[i].id);
		        exit(-1);
		  }
		  closefile(PDB,domain[i].filename);
	          PDB=openfile(domain[i].filename,"r");
		  
                  if(igetgen(PDB,&c[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],
                    domain[i].type[j],C,(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
                  if(test!=add) {
                        fprintf(stderr,"error: missing C atoms in %s\n",domain[i].id);
                        exit(-1);
                  }     
		  closefile(PDB,domain[i].filename);
	          PDB=openfile(domain[i].filename,"r");
                  if(igetgen(PDB,&o[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],
                     domain[i].type[j],O,(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
                  if(test!=add) {
                        fprintf(stderr,"error: missing O atoms in %s\n",domain[i].id);
                        exit(-1);
                  }
		  closefile(PDB,domain[i].filename);
	          PDB=openfile(domain[i].filename,"r");
                  if(igetcb(PDB,&cb[i][total],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],
                    domain[i].type[j],(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
                  if(test!=add) {
                        fprintf(stderr,"error: missing CB atoms in %s\n",domain[i].id);
                        exit(-1);
                  }
	       }
	       if(sidechain) {  
		  /* read all side chains into this sort of array of structures */
		  if(igetside(PDB,side[i],&aa[total],&numb[total],&test,domain[i].start[j],domain[i].end[j],
                    domain[i].type[j],(MAX_SEQ_LEN-total),domain[i].reverse[j],PRECISION,stdout)==-1) exit(-1);
		}


	       switch(domain[i].type[j]) {
	 	  case 1: printf(" all residues"); break; 
		  case 2: printf(" chain %c",domain[i].start[j].cid); break;
		  case 3: printf(" from %c %4d %c to %c %4d %c",
			 domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
			 domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
		} 
		if(polyA) printf("   %4d main chain/C-betas\n",add);
		else printf("   %4d C-alphas\n",add); 
	        total+=add;
		closefile(PDB,domain[i].filename);
	        PDB=openfile(domain[i].filename,"r");
	    }
	    domain[i].ncoords=total;
	    if(polyA) printf("  A total of %4d sets of main chain/C-beta atoms\n",domain[i].ncoords);
	    else printf("  A total of %4d C-alphas\n",domain[i].ncoords);
	    /* disp(domain[i],stdout); */
	    printf(" Applying transformation... \n");
/*	    printmat(domain[i].R,domain[i].V,3,stdout);
	    printf("      ...to these coordinates.\n");  */
/*	    if(polyA) {
               for(j=0; j<5; ++j) {
                  printf("Res %3d %c N : %8d %8d %8d\n",j+1,domain[i].aa[j],n[i][j][0],n[i][j][1],n[i][j][2]);
                  printf("Res %3d %c CA: %8d %8d %8d\n",j+1,domain[i].aa[j],domain[i].coords[j][0],domain[i].coords[j][1],domain[i].coords[j][2]);
                  printf("Res %3d %c C : %8d %8d %8d\n",j+1,domain[i].aa[j],c[i][j][0],c[i][j][1],c[i][j][2]);
                  printf("Res %3d %c O : %8d %8d %8d\n",j+1,domain[i].aa[j],o[i][j][0],o[i][j][1],o[i][j][2]);
                  printf("Res %3d %c CB: %8d %8d %8d\n",j+1,domain[i].aa[j],cb[i][j][0],cb[i][j][1],cb[i][j][2]);
               }
	    } */
	    matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,PRECISION);  
	    if(polyA) {
	       matmult(domain[i].R,domain[i].V,n[i],domain[i].ncoords,PRECISION);  
               matmult(domain[i].R,domain[i].V,c[i],domain[i].ncoords,PRECISION);  
               matmult(domain[i].R,domain[i].V,o[i],domain[i].ncoords,PRECISION);  
	       matmult(domain[i].R,domain[i].V,cb[i],domain[i].ncoords,PRECISION);
/*	       for(j=0; j<5; ++j) {
		  printf("Res %3d %c N : %8d %8d %8d\n",j+1,domain[i].aa[j],n[i][j][0],n[i][j][1],n[i][j][2]);
		  printf("Res %3d %c CA: %8d %8d %8d\n",j+1,domain[i].aa[j],domain[i].coords[j][0],domain[i].coords[j][1],domain[i].coords[j][2]);
 		  printf("Res %3d %c C : %8d %8d %8d\n",j+1,domain[i].aa[j],c[i][j][0],c[i][j][1],c[i][j][2]);
		  printf("Res %3d %c O : %8d %8d %8d\n",j+1,domain[i].aa[j],o[i][j][0],o[i][j][1],o[i][j][2]);
		  printf("Res %3d %c CB: %8d %8d %8d\n",j+1,domain[i].aa[j],cb[i][j][0],cb[i][j][1],cb[i][j][2]); 
	       } */
	    }
/*	    printmat(domain[i].R,domain[i].V,3,stdout); */
	    /* disp(domain[i],stdout); */
	    closefile(PDB,domain[i].filename);
	}
	
	/* Now proceed through the alignment and calculate an average coordinate using the
	 *  reliable positions */
	for(i=0; i<ndomain; ++i) pointer[i]=0;
	ncoords=0; nres=0;
	for(i=0; i<bloclen; ++i) {
	  if(rel[i]==1) {
	    AApt = 23; /* UNK by default */
	    if(cons || ident) { /* Work out if the residue is conserved/identical */
	        identical=1; 
		for(j=0; j<ndomain; ++j) {
		   if(bloc[j+1].seq[i+1]!=bloc[1].seq[i+1]) identical=0;
	        }
		if(identical) {
		     AApt = (int)(bloc[1].seq[i+1]-'A');
		     if(AApt<0 || AApt>26) { AApt = 23; } /* UNK if gone awry */
		} else {
		     AApt = 23; /* UNK if unconserved */
		}
	    }
	    if(cons && AApt == 23) { /* If not identical, look for property conservation */
		npolar=0; nhydrophobic=0; naromatic=0; nsmall=0; ntiny=0; 
		npositive=0; nnegative=0; ncharged=0; naliphatic=0; nbranch=0;
		for(j=0; j<ndomain; ++j) {
		    switch(bloc[j+1].seq[i+1]) {
			case 'R': case 'K': { npositive++; ncharged++; npolar++; } break;
			case 'H': { npositive++; ncharged++; npolar++; naromatic++; } break;
		        case 'E': case 'D': { nnegative++; ncharged++; npolar++; } break;
			case 'Q': case 'N': { npolar++; } break;
			case 'T': { nsmall++; nbranch++; npolar++; } break;
			case 'S': case 'P': case 'G': { nsmall++; ntiny++; npolar++; } break;
			case 'A': { nsmall++; ntiny++; naliphatic++; } break;
			case 'V': { nsmall++; naliphatic++; nbranch++; nhydrophobic++; } break;
			case 'I': { naliphatic++; nbranch++; nhydrophobic++; } break;
			case 'L': case 'M': { naliphatic++; nhydrophobic++; } break;
			case 'C': { nsmall++; ntiny++; nhydrophobic++; } break;
			case 'F': case 'W': { nhydrophobic++; naromatic++; } break;
			case 'Y': { nhydrophobic++; naromatic++;  npolar++; } break;
		    }
		 }
		 if(npositive>=(ndomain-ignore)) { AApt= 30; }
		 else if(nnegative>=(ndomain-ignore)) { AApt = 31; }
		 else if(ncharged>=(ndomain-ignore)) { AApt = 32; }
		 else if(npolar>=(ndomain-ignore)) { AApt = 28; }
		 else if(naromatic>=(ndomain-ignore)) { AApt = 33; }
		 else if(ntiny>=(ndomain-ignore)) { AApt = 27; }
		 else if(nbranch>=(ndomain-ignore)) { AApt = 35; }
		 else if(nsmall>=(ndomain-ignore)) { AApt = 26; }
		 else if(naliphatic>=(ndomain-ignore)) { AApt = 34; }
		 else if(nhydrophobic>=(ndomain-ignore)) { AApt = 29; }
		 else { AApt = 23; } 
	    }
	    for(j=0; j<3; ++j) average[j]=0.0;
	    if(polyA) {
		/* N atoms */
		for(j=0; j<3; ++j) average[j]=0.0;
		for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)n[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
		fprintf(OUT,"ATOM  %5d  N   %3s Z%4d    %8.3f%8.3f%8.3f  0.00 ",
                  ncoords+1,AA3[AApt],nres+1,average[0],average[1],average[2]);
		ncoords++;
		if(realrel[i]==1) fprintf(OUT," 1.00   0 \n");
		else fprintf(OUT,"99.00   0 \n");
	    }
	    /* CA atoms */	
 	    for(j=0; j<3; ++j) average[j]=0.0;
	    for(j=0; j<ndomain; ++j) {	
		for(k=0; k<3; ++k) average[k]+=((float)domain[j].coords[pointer[j]][k]/(float)ndomain)/(float)PRECISION;
	    }
	    /* Output the average coordiates */ 
	    fprintf(OUT,"ATOM  %5d  CA  %3s Z%4d    %8.3f%8.3f%8.3f  0.00 ",
		ncoords+1,AA3[AApt],nres+1,average[0],average[1],average[2]);
	    if(realrel[i]==1) fprintf(OUT," 1.00   0 ");
            else fprintf(OUT,"99.00   0 ");

	    for(j=0; j<ndomain; ++j) 
               if(j<13) fprintf(OUT,"%c",bloc[j+1].seq[i+1]); 
	    fprintf(OUT,"\n");
	    ncoords++;
	    if(polyA) {
		/* C O and CB atoms */
                for(j=0; j<3; ++j) average[j]=0.0;
                for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)c[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
                fprintf(OUT,"ATOM  %5d  C   %3s Z%4d    %8.3f%8.3f%8.3f  0.00 ",
                  ncoords+1,AA3[AApt],nres+1,average[0],average[1],average[2]);
                ncoords++;
		if(realrel[i]==1) fprintf(OUT," 1.00   0 \n");
                else fprintf(OUT,"99.00   0 \n");

                for(j=0; j<3; ++j) average[j]=0.0;
                for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)o[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
                fprintf(OUT,"ATOM  %5d  O   %3s Z%4d    %8.3f%8.3f%8.3f  0.00 ",
                  ncoords+1,AA3[AApt],nres+1,average[0],average[1],average[2]);
	        if(realrel[i]==1) fprintf(OUT," 1.00   0 \n");
                else fprintf(OUT,"99.00   0 \n");
                ncoords++;
                for(j=0; j<3; ++j) average[j]=0.0;
                for(j=0; j<ndomain; ++j) {
                   for(k=0; k<3; ++k) average[k]+=((float)cb[j][pointer[j]][k]/(float)ndomain)/(float)PRECISION;
                }
                fprintf(OUT,"ATOM  %5d  CB  %3s Z%4d    %8.3f%8.3f%8.3f  0.00 ",
                  ncoords+1,AA3[AApt],nres+1,average[0],average[1],average[2]);
	        if(realrel[i]==1) fprintf(OUT," 1.00   0 \n");
                else fprintf(OUT,"99.00   0 \n");
                ncoords++;
	    }
	    nres++; 
	  }
	  for(j=0; j<ndomain; ++j) {
	      if(bloc[j+1].seq[i+1]!=' ') pointer[j]++;
	  }
	}
	free(pointer);
	printf(" ...done.\n");
	exit(0);
}
void exit_error()
{
	 fprintf(stderr,"format:  avestruc -f <file> -c <char> -w <win> -t <thresh> -polyA -aligned\n");
	 fprintf(stderr,"          where -polyA ==> poly alanine structure (default is CAs only)\n");
 	 fprintf(stderr,"                -aligned ==> consider all aligned positions\n");
	 fprintf(stderr,"                -cons ==> label conserved AA positions in the file\n");
	 exit(-1);
}
