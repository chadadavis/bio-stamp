#include "alignfit.h"
#include "gjutil.h"
#include "gjnoc.h"

/* Reads an AMPS format block file containing structurally derived sequences and a
 *  file containing a description as to where the coordinates may be found */

main(int argc, char *argv[]) {

	char c;
	char ins,cid;
	char *env;
	char **ids;
	char cmnd[200];

	int number,total;
	int flag;
	int i,j,k,l,m;
	int ntofit;
	int bloclen;
	int no_gap;
	int nbloc,ndomain,ncoord;
	int maxdomain;
	int nuse,add,nclust;
	int RP,T_FLAG;

	char keyword[100];
	char value[100];
	char tmpstring[100];
	char noc_parms[200];


	int *pointers,*counter;

	float **coords1;
	float **coords2;
	float **R;
	float *V;
	float rmsd;
	double **matrix;

	FILE *IN,*OUT,*BLOC,*LIST,*PDB,*TRANS;

	struct seqdat *bloc;
	struct CA *safe;
	struct domain_loc *domain;
	struct gjnoc *gjclust;
	struct cluster *cl;
	struct parameters *parms;

	strcpy(noc_parms,"noc sim single");
/* SMJS Changed malloc to calloc to zero struct */
	parms=(struct parameters*)calloc(1,sizeof(struct parameters));
	
	if(argc<3) exit_error(); 

	/* define/get parameters etc. */
	parms[0].MAX_SEQ_LEN=1000;
	parms[0].PAIRWISE=1;
	parms[0].TREEWISE=1;
	parms[0].OLDFORMAT=0;
	RP=0;
	strcpy(&parms[0].TRANSFILE[0],"alignfit.trans");
	strcpy(&parms[0].bloc_file[0],"align_mult.out");

	/* get command line arguments */
	for(i=1; i<argc; ++i) {
           if(argv[i][0]!='-') exit_error();
           strcpy(keyword,&argv[i][1]);
           if((i+1)<argc) strcpy(value,argv[i+1]);
           else strcpy(value,"none");
           for(j=0; j<strlen(keyword); ++j) 
              keyword[j]=ltou(keyword[j]); /* change to upper case */
           T_FLAG=(value[0]=='Y' || value[0]=='y' || value[0]=='1' || 
                 value[0]=='T' || value[0]=='t' || value[0]=='o' || 
                 value[0]=='O');
           /* enables one to write '1', 'YES', 'Yes', 'yes', 
            *  'TRUE', 'True' or 'true' to 
            *  set any boolean parmsiable to one 
            * Author's note, three years later: "for some reason
	    *  I must have thought this was particularly clever".
	    *   RBR, Sept 1995 */
	   if(strcmp(&argv[i][1],"f")==0) { 
	     if((i+1)>=argc) exit_error();
             /* listfile name */
             strcpy(&parms[0].bloc_file[0],argv[i+1]);
             i++;
	   } else if(strcmp(&argv[i][1],"d")==0) {
	     if((i+1)>=argc) exit_error();
             /* listfile name */
             strcpy(&parms[0].dom_file[0],argv[i+1]);
             i++;
	   } else if(strcmp(&argv[i][1],"P")==0) {
	     if((i+1)>=argc) exit_error();
             /* listfile name */
             strcpy(&parms[0].parm_file[0],argv[i+1]);
	     RP=1;
             i++;
	   } else if(strcmp(&argv[i][1],"out")==0) {
	     if((i+1)>=argc) exit_error();
             /* listfile name */
             strcpy(&parms[0].TRANSFILE[0],argv[i+1]);
             i++;
	   } else {
	     exit_error();
	   }
	}

	if(RP) { /* read parameters from file (old method) */
	   if((IN=fopen(parms[0].parm_file,"r"))==NULL) {
	      fprintf(stderr,"error: parameter file %s does not exist\n",parms[0].parm_file);
	      exit(-1);
	   }
	   if(getpars(IN,parms)==-1) exit(-1);
	   fclose(IN);
	}

	if((env=getenv("STAMPDIR"))==NULL) {
	   fprintf(stderr,"error: you haven't set STAMPDIR to anything\n");
	   exit(-1);
	}

	/* Output file */
	if((TRANS=fopen(parms[0].TRANSFILE,"w"))==NULL) {
	   fprintf(stderr,"error opening file %s\n",parms[0].TRANSFILE);
	   exit(-1);
	} 
	fprintf(TRANS,"%% Output from the program ALIGNFIT\n");
	fprintf(TRANS,"%% ALIGNFIT was run using the files\n");
	fprintf(TRANS,"%% %s and %s\n",parms[0].bloc_file,parms[0].dom_file);
	fprintf(TRANS,"%%\n");
	fprintf(TRANS,"%% You can generate transformed coordinates by typing\n");
	fprintf(TRANS,"%% transform -f %s [-g -h ]\n",parms[0].TRANSFILE);
	fprintf(TRANS,"%% You can also get a refined StAMP alignment by\n");
	fprintf(TRANS,"%%  typing stamp -l %s [ other parameters]\n",parms[0].TRANSFILE);
	fprintf(TRANS,"%%\n");
	strcpy(&parms[0].STAMPDIR[0],env);
/*	printf("STAMPDIR is %s\n",parms[0].STAMPDIR); */


	/* read in bloc file */
	if((BLOC=fopen(parms[0].bloc_file,"r"))==NULL) {
	   fprintf(stderr,"error: block file %s does not exist\n",parms[0].bloc_file);
	   exit(-1);
	}
	printf("ALIGNFIT R.B. Russell 1995\n");
	printf(" Reading in block file...\n");
	nbloc=0;
	while((c=getc(BLOC))!=(char)EOF) nbloc+=(c=='>');
	bloc=(struct seqdat*)malloc((nbloc+1)*sizeof(struct seqdat));
	rewind(BLOC);
	printf(" ");
	if(Agetbloc(BLOC,bloc,&nbloc)==-1) exit(-1);
	/* Ok, hack to take out "P1;" if it is in the IDs */
	for(i=0; i<nbloc; ++i) {
		if(strncmp(bloc[i+1].id,"P1;",3)==0) {
		     printf(" Warning: changing ID %s in alignment file to ",bloc[i+1].id);
		     sprintf(&tmpstring[0],"%s",&bloc[i+1].id[3]);
		     strcpy(&bloc[i+1].id[0],&tmpstring[0]);
		     printf("%s\n",bloc[i+1].id);
		}
	}
	    
	counter=(int*)malloc(nbloc*sizeof(int));
	bloclen=strlen(bloc[1].seq)-1;
	fclose(BLOC);

	/* read in list of domains */
	if((LIST=fopen(parms[0].dom_file,"r"))==NULL) {
	   fprintf(stderr,"error: domain file %s does not exist\n",parms[0].dom_file);
	   exit(-1);
	}
	printf(" Reading in coordinate descriptions...\n");
	maxdomain=count_domain(LIST);
	domain=(struct domain_loc*)malloc(maxdomain*sizeof(struct domain_loc));
	rewind(LIST);
	if(getdomain(LIST,domain,&ndomain,maxdomain,&j,parms[0].STAMPDIR,0,stdout)==-1) exit(-1);
	fclose(LIST);
	matrix=GJDudarr(ndomain);
	ids=(char**)malloc(ndomain*sizeof(char*));
	for(i=0; i<ndomain; ++i) ids[i]=domain[i].id;
	
	printf(" Reading coordinates...\n");
	/* get the CA coordinates from the brookhaven files 
	 *  
	 * memory allocation */
	pointers=(int*)malloc(ndomain*sizeof(int));
	for(i=0; i<nbloc; ++i) rmsp(bloc[i+1].id);
	for(i=0; i<ndomain; ++i) {
	   pointers[i]=-1;
	   for(j=0; j<nbloc; ++j) {
	      /* finds which id in the blocfile corresponds to the id in
	       *  the domain file */
	      if(strcmp(domain[i].id,bloc[j+1].id)==0) 
		 pointers[i]=j;
	   }
	   if(pointers[i]==-1) {
	      fprintf(stderr,"error: id %s not found in block file\n",domain[i].id);
	      exit(-1);
	   }
	   fprintf(TRANS,"%%Domain %3d %s %s\n",i+1,domain[i].filename,domain[i].id);
	   if((PDB=openfile(domain[i].filename,"r"))==NULL) {
	      fprintf(stderr,"error: PDB file %s does not exist\n",domain[i].filename);
	      exit(-1);
	   }
	   domain[i].ncoords=0;
	   domain[i].coords=(float**)malloc((parms[0].MAX_SEQ_LEN+1)*sizeof(float*));
	   domain[i].aa=(char*)malloc((parms[0].MAX_SEQ_LEN+1)*sizeof(char)); 
	   domain[i].numb=(struct brookn*)malloc((parms[0].MAX_SEQ_LEN+1)*sizeof(struct brookn));
	   total=0;
	   fprintf(TRANS,"%% ");
	   for(j=0; j<domain[i].nobj; ++j) {
	      if(getca(PDB,&domain[i].coords[total],&domain[i].aa[total],
		 &domain[i].numb[total],&add,domain[i].start[j],domain[i].end[j],
		 domain[i].type[j],(parms[0].MAX_SEQ_LEN-total),0,0)==-1) exit(-1);
	      switch(domain[i].type[j]) {
		case 1: fprintf(TRANS,"all residues "); break;
		case 2: fprintf(TRANS,"chain %c ",domain[i].start[j].cid); break;
		case 3: fprintf(TRANS,"from %c %4d %c to %c %4d %c ",
		 domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
		 domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
	      }
	      fprintf(TRANS,"(%3d CAs) ",add);
	      total+=add;
	      closefile(PDB,domain[i].filename);
	      PDB=openfile(domain[i].filename,"r");
	   }
	   fprintf(TRANS,"\n");
	   domain[i].ncoords=total;
	   fprintf(TRANS,"%% %4d CAs in total\n",domain[i].ncoords);
	   closefile(PDB,domain[i].filename);
	   domain[i].use=(int*)malloc(bloclen*sizeof(int));
	   counter[i]=0;
	}


	/* some checks */
	printf(" Checking for inconsistencies...\n");
	if(ndomain!=nbloc) {
	   fprintf(stderr,"error: block file and domain location file disagree\n");
	   printf(" 	  block file has %3d, location file has %3d\n",nbloc,ndomain);
	   exit(-1);
	}
	j=0;
	for(i=0; i<ndomain; ++i) { /* now check the lengths */
	   l=0;
	   for(k=0; k<bloclen; ++k) l+=(bloc[pointers[i]+1].seq[k+1]!=' ');
	   if(domain[i].ncoords!=l) {
	      fprintf(stderr,"error: block file and domain file have a different numbers of residues\n");
	      fprintf(stderr,"        for domain %10s, block file has %4d, DOMAIN file has %4d\n",
			domain[i].id,l,domain[i].ncoords);
	      j=1;
	   }
	}
	if(j==1) exit(-1); /* checks all the lengths then quits if in error */
	/*  The arrays 'use' for each domain will be as long as the 
	 *   block file, and will consist of positions that each
	 *   sequence in the block file corresponds to in each array of
	 *   CA coordinates retrieved using getca() 
	 *  This bit will also compare the blocfile to the sequence
	 *   obtained using getca() to maintain consistency */
	nuse=0;
	for(i=0; i<bloclen; ++i) {
	  no_gap=1;
	  for(j=0; j<ndomain; ++j) 
	     no_gap*=(bloc[pointers[j]+1].seq[i+1]!=' ');
	  if(no_gap) nuse++;
          for(j=0; j<ndomain; ++j) {
	     if(bloc[pointers[j]+1].seq[i+1]!=' ') {
		domain[j].use[i]=counter[j]; 
		counter[j]++;
	     }  else {
		domain[j].use[i]=-1;
	     }
	  }
	  /* comparing the block file to the PDB sequence */
	  flag=0;
	  for(j=0; j<ndomain; ++j) {
	     if(bloc[pointers[j]+1].seq[i+1]!=' ') {
	       if(bloc[pointers[j]+1].seq[i+1]!=domain[j].aa[counter[j]-1]) flag=1;
	     }
	  }

	}

	fprintf(TRANS,"%% There are %4d non-gapped positions in block file\n",nuse);
	fprintf(TRANS,"%% These were used for fitting the coordinates in a treewise fashion\n");

	/* If necessary, calculate each pairwise RMS deviation */
	coords1=(float**)malloc(parms[0].MAX_SEQ_LEN*sizeof(float*));
	coords2=(float**)malloc(parms[0].MAX_SEQ_LEN*sizeof(float*));
	for(i=0; i<parms[0].MAX_SEQ_LEN; ++i) {
	   coords1[i]=(float*)malloc(3*sizeof(float));
	   coords2[i]=(float*)malloc(3*sizeof(float));
	}


	V=(float*)malloc(3*sizeof(float));
	R=(float**)malloc(3*sizeof(float*));
	for(i=0; i<3; ++i) 
	   R[i]=(float*)malloc(3*sizeof(float));
	
	if(parms[0].PAIRWISE) {
   	   printf(" Doing pairwise comparisons...\n");
	   l=0;
	   for(i=0; i<ndomain; ++i) {
	      for(j=i+1; j<ndomain; ++j) {
		 fprintf(TRANS,"%% Comparison: %4d, %s and %s\n",l+1,domain[i].id,domain[j].id); 
		 ntofit=0;
		 for(k=0; k<bloclen; ++k) {
		    if(bloc[pointers[i]+1].seq[k+1]!=' ' && bloc[pointers[j]+1].seq[k+1]!=' ') {
		       for(m=0; m<3; ++m) {
		          coords1[ntofit][m]=domain[i].coords[domain[i].use[k]][m];
		          coords2[ntofit][m]=domain[j].coords[domain[j].use[k]][m];
		       }
		       ntofit++;
		    } 
		 }
		 fprintf(TRANS,"%%  %4d CAs to fit, ",ntofit);
		 l++;
		 m=0;
		 rmsd=fmatfit(coords1,coords2,R,V,ntofit,m);
		 matrix[i][j-i-1]=(double)(1/rmsd);
		 fprintf(TRANS,"RMS deviation %10.4f.\n",rmsd);
	      }
	   }
	}
/*	printf(" IDs\n");
	for(i=0; i<ndomain; ++i) {
	   printf(" %s\n",ids[i]);
	}
	printf(" Upper diagonal is...\n");
	for(i=0; i<ndomain; ++i) {
	  for(j=i+1; j<ndomain; ++j) {
	     printf(" %7.5f ",matrix[i][j-i-1]);
	  }
	  printf(" \n");
	}
*/
	if(parms[0].TREEWISE) {
	   printf(" Doing treewise comparisons...\n");
	   /* Use Geoff's routine to get a tree structure */
	   cl=get_clust(matrix,ids,ndomain,noc_parms);
	   nclust=ndomain-1;
	

	   for(i=0; i<nclust; ++i) {
	      fprintf(TRANS,"%%Cluster: %4d (",i+1);
	      for(k=0; k<cl[i].a.number; ++k)
		 fprintf(TRANS,"%s ",domain[cl[i].a.member[k]].id);
	      fprintf(TRANS,") and (");
	      for(k=0; k<cl[i].b.number; ++k)
		 fprintf(TRANS,"%s ",domain[cl[i].b.member[k]].id);
	      fprintf(TRANS,")\n");
	      /* we shall generate average coordinats for each
	       *  cluster, get transformation based on these 
	       *  coordinates */
	      ntofit=0;
	      for(j=0; j<bloclen; ++j) {
		 /* determine whether to use the position or not */
		 no_gap=1;
		 for(k=0; k<cl[i].a.number; ++k) 
		    no_gap*=(bloc[pointers[cl[i].a.member[k]]+1].seq[j+1]!=' ');
		 for(k=0; k<cl[i].b.number; ++k)
		    no_gap*=(bloc[pointers[cl[i].b.member[k]]+1].seq[j+1]!=' ');
		 if(no_gap) {
		    for(l=0; l<3; ++l) coords1[ntofit][l]=coords2[ntofit][l]=0.0;
		    for(k=0; k<cl[i].a.number; ++k) {
		       m=cl[i].a.member[k];
		       for(l=0; l<3; ++l) 
			  coords1[ntofit][l]+=domain[m].coords[domain[m].use[j]][l];
		    }
		    for(k=0; k<cl[i].b.number; ++k) {
		       m=cl[i].b.member[k];
		       for(l=0; l<3; ++l)
		          coords2[ntofit][l]+=domain[m].coords[domain[m].use[j]][l];
		    }
		    for(l=0; l<3; ++l) {
		       coords1[ntofit][l]/=cl[i].a.number;
		       coords2[ntofit][l]/=cl[i].b.number;
		    }
		    ntofit++;
		  } /* end of if(no_gap... */
	      } /* end of for(j... */
	      fprintf(TRANS,"%% CAs to fit: %4d, ",ntofit);
	      rmsd=fmatfit(coords1,coords2,R,V,ntofit,1);
	      fprintf(TRANS," RMS between ave atoms: %10.5f.\n",rmsd);
	      /* now we must apply the transformation to each set of
	       *  coordinates in the 'B' cluster */
	      for(k=0; k<cl[i].b.number; ++k) {
		/* void fmatmult(float **r, float *v, float **coord, int n)  */
		 fmatmult(R,V,domain[cl[i].b.member[k]].coords,
			domain[cl[i].b.member[k]].ncoords);
		 /* void update(float **dR, float **R, float *dV, float *V) */
		 update(R,domain[cl[i].b.member[k]].R,V,domain[cl[i].b.member[k]].V);
	      } 
	   } 
	   /* outputing transformations */

	   newoutput(TRANS,domain,ndomain,1);
	} /* end of if(parms[0].TREE... */

	printf(" ALIGNFIT done.\n Look in the file %s for output and details\n",parms[0].TRANSFILE);
	/* freeing to keep purify happy */
	free(counter); free(pointers);
  	for(i=0; i<parms[0].MAX_SEQ_LEN; ++i) {
	   free(coords1[i]);
	   free(coords2[i]);
	}


	free(coords1); free(coords2);

	exit(0);
}
void exit_error()
{
	   fprintf(stderr,"format: alignfit -f <block format alignment file> -d <domain description file>\n");
	   fprintf(stderr,"               -out <output file> -P <paramter file>\n");
	   exit(-1);
}
