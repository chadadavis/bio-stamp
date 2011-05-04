#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "stamp.h"
#include "stamprel.h"

#define MAX_SEQ_LEN 2000
#define PRECISION 1000

/* gstamp: reads in a stamp output file and generates molscripts files
 *   displaying the regions of similarity */

void exit_error();

main(int argc, char *argv[]) {

	int i,j,k,l,test;
	int ndomain,total,add,nbloc;
	int nstamp,nseq,nstamppos;
	int gottrans;
	int bloclen,maxlen;
	int stampwindow;
	int nrel;
	int ncoords,nres;
	int all_aligned;
	int gen_average;
	int helix,strand,coil;
	int cons_sec;
	int colour;
	int verbose;

	int *counter;
	int *all;
	int *rel;
	int **avecoords;
	int ***n,***c,***o,***cb;
	

	char filename[200],outfile[200];
	char stampchar;
	char type;
	char *summary;

	char *env;
	char *buff;
	char *aa;
	char *sec;
	char *N,*C,*O;

	float stampthresh;
	float *average;

	FILE *TRANS,*PDB,*OUT;

	struct domain_loc *domain;
	struct seqdat *bloc;
	struct stampdat *stamp;
	struct brookn *numb;

	if(argc<3) exit_error();
	
	/* defaults */
	stampchar='G';
	stampwindow=3;
	stampthresh=6.0;
	verbose=0;
	strcpy(&outfile[0],"ave.molscript");
	
	buff=(char*)malloc(100*sizeof(char));
	all_aligned=0; gen_average=0;
	cons_sec=0; colour=0;

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
	   } else if(strcmp(&argv[i][1],"o")==0) {
	     if((i+1)>=argc) exit_error();
             strcpy(&outfile[0],argv[i+1]);
             i++;
	   } else if(strcmp(&argv[i][1],"aligned")==0) {
	     all_aligned=1;
	   } else if(strcmp(&argv[i][1],"a")==0) {
	     gen_average=1;
	   } else if(strcmp(&argv[i][1],"cons")==0) {
	     cons_sec=1;
	   } else if(strcmp(&argv[i][1],"colour")==0) {
	     colour=1;
	   } else if(strcmp(&argv[i][1],"V")==0) {
	      verbose=1;
	   } else {
	     exit_error();
	   }
	}

	printf("GSTAMP R.B. Russell, 1995\n");
        printf(" Makes input for MOLSCRIPT, to display STAMP\n");
        printf("  superimpositions in PostScript format\n");
	printf(" You must run `transform' or `avestruc' on the file:\n   %s\n",filename);
	printf("  in order to run MOLSCRIPT subsequently\n\n");

	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
      	}
	/* read in coordinate locations and initial transformations */
	if((TRANS = fopen(filename,"r")) == NULL) {
	   fprintf(stderr,"error: file %s does not exist\n",filename);
	   exit(-1);
	}

	/* determine the number of domains specified */
	printf(" Reading domain descriptors...\n");
	ndomain=count_domain(TRANS);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	rewind(TRANS);
	if(getdomain(TRANS,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	for(i=0; i<ndomain; ++i) {
	  for(j=0; j<domain[i].nobj; ++j) {
	     if(domain[i].start[j].cid=='_') domain[i].start[j].cid=' ';
	     if( domain[i].start[j].in=='_') domain[i].start[j].in=' ';
	     if(  domain[i].end[j].cid=='_') domain[i].end[j].cid=' '; 
	     if(   domain[i].end[j].in=='_') domain[i].end[j].in=' ';
	  }
	} 

	bloc=(struct seqdat*)malloc((ndomain*2+2)*sizeof(struct seqdat));
	/* read in the alignment */
	printf(" Reading alignment...\n");
	rewind(TRANS);
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
	for(i=0; i<ndomain; ++i) {
	   /* three state assignment */
	   if(threestate(&bloc[i+ndomain+2].seq[1],"HG","EB","TS-Ic")==-1) exit(-1);
	   if(smoothsec(&bloc[i+ndomain+2].seq[1],4,2)==-1) exit(-1);
	}

	/* read in the STAMP output and get reliable regions */
	rewind(TRANS);
	if(all_aligned) {
	   rel=(int*)malloc(MAX_SEQ_LEN*sizeof(int));
	   for(i=0; i<bloclen; ++i) {
	     k=0;
	     for(j=0; j<ndomain; ++j) {
		if(bloc[j+1].seq[i+1]!=' ') k++;
	     }
	     if(k==ndomain) rel[i]=1;
	     else rel[i]=0;
	    }
        } else { 
	     stamp=(struct stampdat*)malloc(100*sizeof(struct stampdat));
	     if(getstampdat(stamp,TRANS,&nstamp,&nseq,&nstamppos,MAX_SEQ_LEN)==-1) exit(-1);
	     if((rel=getstamprel(stamp,nstamp,bloclen,stampchar,stampthresh,stampwindow))==NULL) exit(-1);
	}
	nrel=0;
	   for(i=0; i<bloclen; ++i) {
             if(rel[i]==1) nrel++;
	}
	fclose(TRANS);

	/* display the secondary structure summary */
	if(verbose) {
	  j=0;
	  printf(" The alignment looks like this...\n");
	  while(j<bloclen) {
	      for(l=0; l<nbloc; ++l) if(strncmp(bloc[l+1].id,"space",5)!=0) {
		printf(" %20s ",bloc[l+1].id);
	        k=0; 
		while((k+j)<bloclen && k<60) { 
		   printf("%c",bloc[l+1].seq[j+k+1]); 
		   k++;
	        } 
	        printf("\n");
	      }
              k=0; 
	      printf("      Reliablity      ");
	      while((k+j)<bloclen && k<60) { 
		printf("%1d",rel[j+k]); 
		k++; 
	      } 
	      printf("\n\n");
	      j+=k;
	  }
	}

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

	       switch(domain[i].type[j]) {
	 	  case 1: printf(" all residues"); break; 
		  case 2: printf(" chain %c",domain[i].start[j].cid); break;
		  case 3: printf(" from %c %4d %c to %c %4d %c",
			 domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
			 domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
		} 
		printf("%4d CAs ",add); 
	        total+=add;
		closefile(PDB,domain[i].filename); PDB=openfile(domain[i].filename,"r");
	    }
	    domain[i].ncoords=total;
	    printf(" A total of %4d CAs in total\n",domain[i].ncoords);
	    /* disp(domain[i],stdout); */
	    printf(" Applying transformation... \n");
/*	    printmat(domain[i].R,domain[i].V,3,stdout);
	    printf("      ...to these coordinates.\n");  */
	    matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,PRECISION);  
	    closefile(PDB,domain[i].filename);
	}
        summary=(char*)malloc(bloclen*sizeof(char));

	printf(" Processing alignment...\n");
	for(i=0; i<bloclen; ++i) {
	   helix=strand=coil=0;
	   for(j=0; j<ndomain; ++j) {
	      /* change '_' to ' ' in brookhaven numbering */
	      if(domain[j].numb[i].cid=='_') domain[j].numb[i].cid=' ';
	      if(domain[j].numb[i].in=='_') domain[j].numb[i].in=' ';
	      if(bloc[j+ndomain+2].seq[i+1]=='H') helix++;
	      if(bloc[j+ndomain+2].seq[i+1]=='B') strand++;
	      if(bloc[j+ndomain+2].seq[i+1]=='c') coil++;
	   }
	   if(cons_sec==1) {
	     /* if we are after a consensus, then the majority secondary
	      *  structure will dominate */
	     if(helix>=(int)((float)ndomain/2)) summary[i]='H';	
             else if(strand>=(int)((float)ndomain/2)) summary[i]='B';
	     else summary[i]='c';
	   } else {
	     if(helix==ndomain) summary[i]='H';
	     else if(strand==ndomain) summary[i]='B';
	     else summary[i]='c';
	   }
	}
	printf(" Generating MOLSCRIPT file(s)...\n");
	if(gen_average) {
          if((OUT=fopen(outfile,"r"))!=NULL) {
             fprintf(stderr,"error: file %s already exists\n",outfile);
             exit(-1);
          }
	  printf(" Writing average molscript file to %s\n",outfile);
          fclose(OUT); OUT=fopen(outfile,"w");
	  fprintf(OUT,"! Molscript input file generated from STAMP output file:\n");
	  fprintf(OUT,"!  domain %s from %s\n",domain[i].id,filename);
          fprintf(OUT,"! R.B. Russell 1994\n");
          fprintf(OUT,"!\n!\n");
          fprintf(OUT,"! To use you must run AVESTRUC on the same file \n");
          fprintf(OUT,"!  (i.e. to generate the averaged PDB file)\n");
	  fprintf(OUT,"! A total of %d equivalent residues were found\n",nrel);
          fprintf(OUT,"!\n!\n");
          fprintf(OUT,"plot\n");
          fprintf(OUT,"read mol \"average.pdb\";\n");
          fprintf(OUT,"  transform atom *\n");
          fprintf(OUT,"      by centre position atom *\n");
          fprintf(OUT,"      by rotation x 0.0\n");
          fprintf(OUT,"      by rotation y 0.0\n");
	  fprintf(OUT,"      by rotation z 0.0\n");
	  fprintf(OUT,"      ;\n");
          fprintf(OUT," set shading 0.0;\n");
	  ncoords=0; nres=0;
	  j=0; nres=0;
	  while(j<bloclen) {
	     type=summary[j];
	     if(rel[j]==1) { /* reliable region */
                switch(type) {
                     case 'H':  {
				if(colour==1) fprintf(OUT,"set planecolour blue;\nset plane2colour blue;\n");
				fprintf(OUT,"helix "); break;
		     }
                     case 'B':  {
				if(colour==1) fprintf(OUT,"set planecolour red;\nset plane2colour red;\n");
			 	fprintf(OUT,"strand "); break;
		     }          
                     default:   {
			   if(colour==1) fprintf(OUT,"set planecolour white;\nset plane2colour white;\n");
			   fprintf(OUT,"coil ");
		     }
                }
	        fprintf(OUT,"from A%d to ",nres+1);
                while(rel[j]==1 && j<bloclen && summary[j]==type) {
                   j++; 
                   nres++;
                }
                fprintf(OUT,"A%d;\n",nres+1);
		j++;
             } else {
		j++;
	     }
	   }
	   fprintf(OUT,"end_plot\n");
           fclose(OUT);
	} else {
	  /* if not average, output N molscript files */
	  for(i=0; i<ndomain; ++i) {
	    sprintf(filename,"%s.molscript",domain[i].id);
	    printf(" Domain %s to file %s\n",domain[i].id,filename);
/*	    if((OUT=fopen(filename,"r"))!=NULL) {
		fprintf(stderr,"error: file %s already exists\n",filename);
		exit(-1);
	    }
	    fclose(OUT); 
*/
	    OUT=fopen(filename,"w");
	    /* Molscript preamble */
	    fprintf(OUT,"! Molscript input file generated from STAMP output file\n");
	    fprintf(OUT,"! R.B. Russell 1994\n");
	    fprintf(OUT,"!\n!\n");
	    fprintf(OUT,"! To use you must run TRANSFORM on the same file\n");
	    fprintf(OUT,"!  (i.e. to generate the PDB files)\n");
	    fprintf(OUT,"! A total of %d equivalent residues were found\n",nrel);
	    fprintf(OUT,"!\n!\n");

	    fprintf(OUT,"plot\n");
	    fprintf(OUT,"read mol \"%s.pdb\";\n",domain[i].id);
            fprintf(OUT,"  transform atom *\n");
            fprintf(OUT,"      by centre position atom *\n");
            fprintf(OUT,"      by rotation x 0.0\n");
            fprintf(OUT,"      by rotation y 0.0\n");
            fprintf(OUT,"      by rotation z 0.0\n");
            fprintf(OUT,"      ;\n");
            fprintf(OUT," set shading 0.0;\n");



	    nres=0; j=0;
	    j=0;

	    while(j<bloclen) {
	      type=bloc[i+ndomain+2].seq[j+1]; 
	      if(type!=' ') {
	        if(rel[j]==1) { /* reliable region */
		  switch(type) {
                     case 'H':  {
                                if(colour==1) fprintf(OUT,"set planecolour blue;\nset plane2colour blue;\n");
                                fprintf(OUT,"helix "); break;
                     }
                     case 'B':  {
                                if(colour==1) fprintf(OUT,"set planecolour red;\nset plane2colour red;\n");
                                fprintf(OUT,"strand "); break;
                     }
                     default:   {
                           if(colour==1) fprintf(OUT,"set planecolour white;\nset plane2colour white;\n");
                           fprintf(OUT,"coil ");
                     }
                  }
		  if(type=='c' && nres>0 && rel[j-1]==1) { /* last residue can't be trace */
	             fprintf(OUT,"from %c%d%c to ",
                        domain[i].numb[nres-1].cid,
                        domain[i].numb[nres-1].n,
                        domain[i].numb[nres-1].in);
		  } else {
		     fprintf(OUT,"from %c%d%c to ",
                        domain[i].numb[nres].cid,
                        domain[i].numb[nres].n,
                        domain[i].numb[nres].in);
                  } 
                  while(rel[j]==1 && j<bloclen && bloc[i+ndomain+2].seq[j+1]!=' ' && 
                        bloc[i+ndomain+2].seq[j+1]==type && nres<domain[i].ncoords) { 
                           j++;
                           nres++;
                  }
	 	  if(type=='c' && nres<domain[i].ncoords) {
		   fprintf(OUT,"%c%d%c;\n",
                     domain[i].numb[nres].cid,
                     domain[i].numb[nres].n,
                     domain[i].numb[nres].in);
		  } else {
		    fprintf(OUT,"%c%d%c;\n",
                     domain[i].numb[nres-1].cid,
                     domain[i].numb[nres-1].n,
                     domain[i].numb[nres-1].in);
		  }
	        } else { /* unreliable region (in trace) */
		   fprintf(OUT,"trace ");
		 if(nres>0) {
                   fprintf(OUT,"from %c%d%c to ",
                      domain[i].numb[nres-1].cid,
                      domain[i].numb[nres-1].n,
                      domain[i].numb[nres-1].in);
		  } else {
		   fprintf(OUT,"from %c%d%c to ",
                      domain[i].numb[nres].cid,
                      domain[i].numb[nres].n,
                      domain[i].numb[nres].in);
		  }
		  while(rel[j]==0 && j<bloclen) {
		    j++;
		    if(bloc[i+ndomain+2].seq[j+1]!=' ') nres++;
		  }
		  if(nres<domain[i].ncoords) {
                    fprintf(OUT,"%c%d%c;\n", 
                      domain[i].numb[nres].cid, 
                      domain[i].numb[nres].n, 
                      domain[i].numb[nres].in);
		  } else {
		     fprintf(OUT,"%c%d%c;\n",
                      domain[i].numb[nres-1].cid, 
                      domain[i].numb[nres-1].n, 
                      domain[i].numb[nres-1].in);
		  }
		}
	      } else {
	        j++;
	      } 
	      if(nres==(domain[i].ncoords-1)) break;
	    }
	    /* Molscript end */
	    fprintf(OUT,"end_plot\n");
	    fclose(OUT);
	  }
	}
	printf(" ...done.\n");
	printf(" You now must run molscript using the files created above,\n");
	printf("  to get PostScript pictures of the superimposed structures.\n");

	exit(0);
}
void exit_error()
{
	 fprintf(stderr,"format:  gstamp -f <file> -c <char> -w <win> -t <thresh> -a -aligned -cons\n");
 	 fprintf(stderr,"          -aligned ==> consider all aligned positions\n");
	 fprintf(stderr,"          -a ==> generate average structure file\n");
	 fprintf(stderr,"          -cons ==> generate consensus secondary structure summary\n");	
	 fprintf(stderr,"           (default is all or nothing, e.g. all helix all strand or coil)\n");
	 fprintf(stderr,"          -V ==> verbose output\n");
	 exit(-1);
}
