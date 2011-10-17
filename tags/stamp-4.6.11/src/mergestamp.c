#include <mergestamp.h>

/* MERGESTAMP
 * Given two transformation/alignment files and and optional domain ID, 
 *  this program centres all transformations/alignments such that they are
 *  expressed as a transformation of all domains onto the common ID.
 *
 * If transformations/alignments are missing from both files, then only
 *  alignments/transformations will be merged.  If these are missing from
 *  one file, then the file with them will just be echoed.
 */


int numseq(FILE *fp);

int main(int argc, char *argv[]) {
	
	int i,j,k;
	int n1, n2, n2new, which, which1, which2, nn, found, len;
	int bpt1, bpt2, pt;
	int ndomain,gotdomain;
	int ndomain2,gotdomain2;
	int got_file,got_file2;
	int got_id;
	int gottrans,gottrans2;
	int ignore;
	int nstamp1,nstamp2;
	int nstamppos1,nstamppos2;
	int nstampseq1,nstampseq2;
	int newnstamp,newnstamppos,newnstampseq;
	int gapped,allspace;
	int MAX_SEQ_LEN;
	int verbose;

	int *indx;
	int *rel1,*rel2,*newrel;
	int *neworder; /* Name of shite band */
	int *use2; /* array of booleans telling the routine which sequences within the second alignment
		    *  to use and which to ignore (ie the ones found in the first alignment) */

	float sign;
	float *negvec;
	float **invmat;
	float **R,**RI;

	char c;
	char *id;
	char *buff;
	char infile[200],infile2[200];
    char *stampdir = AM_STAMPDIR;

	FILE *DOM1,*DOM2;
	FILE *A1,*A2;

	struct domain_loc *domain,*domain2;

	struct seqdat *b1, *b2,*new;
	struct stampdat *s1,*s2,*sn;

	id=(char*)malloc(100*sizeof(char));
	buff=(char*)malloc(1000*sizeof(char));
	indx=(int*)malloc(100*sizeof(int));
/* SMJS Was sizeof(float) */
	invmat=(float**)malloc(3*sizeof(float *));
/* SMJS Was sizeof(float) */
	R=(float**)malloc(4*sizeof(float *));
/* SMJS Was sizeof(float) */
	RI=(float**)malloc(4*sizeof(float *));
	negvec=(float*)malloc(3*sizeof(float));
	for(i=0; i<4; ++i) { 
	   R[i]=(float*)malloc(4*sizeof(float));
	   RI[i]=(float*)malloc(4*sizeof(float));
	}
	for(i=0; i<3; ++i) 
	  invmat[i]=(float*)malloc(3*sizeof(float));


	if(argc<3) exit_error();

	got_file=got_file2=got_id=0;
	MAX_SEQ_LEN=100000;
	verbose=0;
	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(strcmp(&argv[i][1],"f1")==0) { 
	      if((i+1)>=argc) exit_error();
	      if((DOM1=fopen(argv[i+1],"r"))==NULL) {
		 fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
		 exit(-1);
	      }
	      got_file=1;
	      strcpy(infile,argv[i+1]);
	      i++;
	   } else if(strcmp(&argv[i][1],"f2")==0) {
              if((i+1)>=argc) exit_error();
              if((DOM2=fopen(argv[i+1],"r"))==NULL) {
                 fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
                 exit(-1);
              }
              got_file2=1;
              strcpy(infile2,argv[i+1]);
              i++;
	   } else if(argv[i][1]=='i') {
	      if((i+1)>=argc) exit_error();
	      strcpy(id,argv[i+1]);
	      got_id=1;
	      i++;
	   } else if(strcmp(&argv[i][1],"max_seq_len")==0) {
              if((i+1)>=argc) exit_error();
              sscanf(argv[i+1],"%d",&MAX_SEQ_LEN);
              i++;
	   }  else if(strcmp(&argv[i][1],"V")==0 || strcmp(&argv[i][1],"v")==0) {
		verbose=1;
	   } else exit_error();
	}

	if(!got_file || !got_file2) {
	  fprintf(stderr,"One or both files not found\n");
	  exit(-1);
	}

    if(getenv("STAMPDIR")!=NULL) {
      /* Allow environment variable to override config setting */
      stampdir=getenv("STAMPDIR");
    }

	ndomain=count_domain(DOM1);
	rewind(DOM1);
	if(ndomain>0) {
	  domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	  if(getdomain(DOM1,domain,&ndomain,ndomain,&gottrans,stampdir,0,stdout)==-1) exit(-1);
	}
	fclose(DOM1);

	ndomain2=count_domain(DOM2);
        rewind(DOM2);
	if(ndomain2>0) {
          domain2=(struct domain_loc*)malloc(ndomain2*sizeof(struct domain_loc));
          if(getdomain(DOM2,domain2,&ndomain2,ndomain2,&gottrans2,stampdir,0,stdout)==-1) exit(-1);
	}
	fclose(DOM2);

	if(ndomain>0 && ndomain2>0 && got_id) { /* User has specified ID */
  	  /* find which domain ID to use */
	  which=-1;
	  for(i=0; i<ndomain; ++i) {
	     if(strstr(domain[i].id,id)!=NULL) {
	        which=i;
	        break;
	     }
	  }
	  if(which==-1) {
	     fprintf(stderr," Error: ID %s not found in file %s\n",id,infile);
	     fprintf(stderr," The identifiers in the file are:\n"); 
	     for(i=0; i<ndomain; ++i) fprintf(stderr,"          %s\n",domain[i].id);
	     exit(-1);
	  }
	  which2=-1;
          for(i=0; i<ndomain2; ++i) {
             if(strstr(domain2[i].id,id)!=NULL) {
                which2=i;
/*              printf("which = %4d\n",which2); */
                break;
             }
          }
          if(which2==-1) {
             fprintf(stderr," Error: ID %s not found in file %s\n",id,infile2);
             fprintf(stderr," The identifiers in the file are:\n"); 
             for(i=0; i<ndomain2; ++i) fprintf(stderr,"          %s\n",domain2[i].id);
             exit(-1);
          }
	} else if(ndomain2>0) {
	  /* Just look for the first ID in common */
	  which = which2 = -1;
	  for(i=0; i<ndomain; ++i) {
		for(j=0; j<ndomain2; ++j) {
		    if(strstr(domain2[i].id,domain[j].id)!=NULL) {
			which = i; which2 = j;
			printf("%% Using domain %s to centre (%d in %s; %d in %s)\n",
			   domain[i].id,i+1,infile,j+1,infile2);
			break;
		    }
	         }
		 if(which !=-1) break;
	   }
	} else {
	   fprintf(stderr,"Warning no domains found in file %s - ignored in output\n",infile2);
	   which = 0;
	}
	if(which==-1 || which2==-1) {
	   fprintf(stderr,"Warning couldn't index domains from the two files.  Domains from %s will be ignored\n",infile2);
	   ndomain2 = 0;
	   which = 0;
	}
	if(ndomain==0) {
	    which2 =0;
	}
	if(ndomain2==0) {
	     which = 0;
	}



	/* put some comments in to tell which domain is centered */
	printf("%%\n%%\n%% MERGESTAMP R.B. Russell, 1997\n");
	printf("%%  The input files were %s and %s\n",infile,infile2);
	if(ndomain!=0 && ndomain2!=0) {
	   printf("%%  the domains in this file have been centred on domain %s\n",domain[which].id);
	} else {
	   printf("%%  the domains/transformations were not centered since they were not found or\n");
	   printf("%%  couldn't be indexed\n");
	}
		
	printf("%%  The original comments are have been removed.\n");
	printf("%%\n%%\n"); /* N.B. Comments are lost */


	if(ndomain>0) {
	  /* Do each file one at a time, then just ignore the second copy (i.e. which2) */
	  /* First domains first */
	  /* get the inverse */
	  for(i=0; i<3; ++i) {
	    for(j=0; j<3; ++j) {
	     R[i+1][j+1]=domain[which].R[i][j];
	    }
	  }
	  matinv(R,RI,&sign,indx);
	  for(i=0; i<3; ++i) { 
	   for(j=0; j<3; ++j) {
	      invmat[i][j]=RI[i+1][j+1];
	   }
	  }
	  /* store the negative of domain[which].V in negvec */
	  for(i=0; i<3; ++i) 
	    negvec[i]=-1*domain[which].V[i];

	  /* See PICKFRAME for an explanation */
	  for(i=0; i<ndomain; ++i) {
	   /* first apply the inverse to the translation */
	   for(j=0; j<3; ++j) {
	     domain[i].V[j]=domain[i].V[j]+negvec[j];
	   }
	   matvecprod(invmat,domain[i].V,domain[i].V,stdout);
	   /* now apply the inverse to the old matrix */
	   matprod(domain[i].R,invmat,domain[i].R,stdout);
	   printdomain(stdout,domain[i],1);
	  }
	} else {
	   printf("%% No domains found in file %s\n",infile);
	}

	/* Second set of domains */
	if(ndomain2>0) {
	  for(i=0; i<3; ++i) {
            for(j=0; j<3; ++j) {
               R[i+1][j+1]=domain2[which2].R[i][j];
            }
          }
          matinv(R,RI,&sign,indx);
          for(i=0; i<3; ++i) {
             for(j=0; j<3; ++j) {
                invmat[i][j]=RI[i+1][j+1];
             }
          }
          /* store the negative of domain[which].V in negvec */
          for(i=0; i<3; ++i)
            negvec[i]=-1*domain2[which2].V[i];

          /* See PICKFRAME for an explanation */
          for(i=0; i<ndomain2; ++i) {
	   ignore=0;
           for(j=0; j<ndomain; ++j) if(strcmp(domain[j].id,domain2[i].id)==0) ignore=1;
	   if(!ignore) { 
              /* first apply the inverse to the translation */
              for(j=0; j<3; ++j) {
                domain2[i].V[j]=domain2[i].V[j]+negvec[j];
              }
              matvecprod(invmat,domain2[i].V,domain2[i].V,stdout);
              /* now apply the inverse to the old matrix */
              matprod(domain2[i].R,invmat,domain2[i].R,stdout);
              printdomain(stdout,domain2[i],1);
           } else {
	      printf("%% Second copy of domain %s ignored\n",domain2[i].id);
	   }
	 }
	} else {
	  printf("%% No domains found in %s\n",infile2);
	}

	/* N.B. This was the merger of two old programs
	 *  the checks on the files below should already have been done */
	if((A1=fopen(infile,"r")) ==NULL) { exit(-1); }
	if((A2=fopen(infile2,"r"))==NULL) { exit(-1); }

	nn=0;
	i=numseq(A1);
	nn+=i;
	if(i>0) {
	  b1=(struct seqdat*)malloc((i+2)*sizeof(struct seqdat));
	  Agetbloc(A1,b1,&n1);
	  for(i=0; i<n1; ++i) {
	   rmsp(b1[i+1].id);
	   for(j=0; j<strlen(b1[i+1].id); ++j) {
		if(b1[i+1].id[j]=='\n')  b1[i+1].id[j]='\0';
	   }
	  }
	} else {
	  n1=0;
	}

	i=numseq(A2);
	if(i>0) {
	  b2=(struct seqdat*)malloc((i+2)*sizeof(struct seqdat));
	  Agetbloc(A2,b2,&n2);
	  use2=(int*)malloc(n2*sizeof(int));
	  for(i=0; i<n2; ++i) {
	    use2[i]=1;
	    rmsp(b2[i+1].id);
	    for(j=0; j<strlen(b2[i+1].id); ++j) {
		if(b2[i+1].id[j]=='\n')  b2[i+1].id[j]='\0';
	    }
	  }
	} else {
	  n2=0;
	}
	if(n1==0 && n2==0) exit(0);

	/* Get STAMP information for each file */
	rewind(A1);
	s1=(struct stampdat*)malloc(10*sizeof(struct stampdat));
        if(getstampdat(s1,A1,&nstamp1,&nstampseq1,&nstamppos1,b1[1].slen)==-1) exit(-1);
	fclose(A1);
	if(nstamp1>0) {
	  printf("%% STAMP data from first file: %4d fields\n",nstamp1);
	}

	rewind(A2);
        s2=(struct stampdat*)malloc(10*sizeof(struct stampdat));
        if(getstampdat(s2,A2,&nstamp2,&nstampseq2,&nstamppos2,b2[1].slen)==-1) exit(-1);
        fclose(A2);
	if(nstamp2>0) {
          printf("%% STAMP data from second file: %4d fields\n",nstamp2);
        }


	


        /* determine the reliable regions */
	if(nstamp1>1) {
           if((rel1=getstamprel(s1,nstamp1,nstamppos1,'G',4.0,2))==NULL) exit(-1);
	} else {
	   rel1=(int*)malloc(b1[1].slen*sizeof(int));
	   for(i=0; i<b1[1].slen; ++i) rel1[i]=1;
	}
	if(nstamp2>1) {
          if((rel2=getstamprel(s2,nstamp2,nstamppos2,'G',4.0,2))==NULL) exit(-1);
	} else {
           rel2=(int*)malloc(b2[1].slen*sizeof(int));
           for(i=0; i<b2[1].slen; ++i) rel2[i]=1;
        }

	/* now determine which sequence in thesecond bloc file are unique and
	 *  which are repeats of the first blocfile */

	/* read through each blocfile till a common member is found */
	if(got_id) { /* User has specified ID */
          /* find which domain ID to use */
          which1=-1;
	  found=0;
          for(i=1; i<=n1; ++i) {
             if(strstr(b1[i].id,id)!=NULL) {
		found=1;
		which1=i;
	     }
	  }
	  if(found==0) {
		fprintf(stderr,"error: couldn't find id %s in alignment from file %s\n",id,infile);
		exit(-1);
	  }

	  which2=-1;
          found=0;
          for(i=1; i<=n2; ++i) {
             if(strstr(b2[i].id,id)!=NULL) {
                found=1;
                which2=i;
             }
          }
	  if(found==0) {
                fprintf(stderr,"error: couldn't find id %s in alignment from file %s\n",id,infile2);
                exit(-1);
          }
	} else {
	 i=1; found=0;
	 while(i<=n1)  {
	  j=1;
	  while(j<=n2) {
/*	   printf("id1: %s, id2: %s\n",b1[i].id,b2[j].id); */
	   if(strstr(b2[j].id,b1[i].id)!=NULL) {  
	      printf("Repeated sequence: %s --- copy in file %s ignored\n",b1[i].id,argv[2]);
	      /* if we have an identity, ignore this sequence in the second blocfile */
	      use2[j-1]=0;
	      if(!found && strncmp(b1[i].id,"space",5)!=0) { 
		which1=i; 
		which2=j; 
		found=1; 
		printf(" (using this sequence to match the two alignments)\n"); 
	      }
	   }
	   j++;
	   }
	  i++;
	 }
	}
	if(!found) {
	   printf("error: no sequence matches were found between the two files\n");
	   printf("   you may have to edit them such that at least one id in each file\n");
	   printf("   is in common\n");
	   for(j=0; j<n1; ++j) { printf("%4d %s\n",j+1,b1[j+1].id); }
	   for(j=0; j<n2; ++j) { printf("%4d %s\n",j+1,b2[j+1].id); }
	   return -1;
	}
	for(i=0; i<n2; ++i) nn+=(use2[i-1]==1);

	/* Now align the two sequence creating a new bloc `new' */
	new=(struct seqdat*)malloc((n1+n2)*sizeof(struct seqdat));
	/* The maximum alignment length is simply the sum of the two initial alignment lengths */
	len=b1[1].slen+b2[1].slen;

	/* Allocate new stamp data if necessary */
	if(nstamp1>0 || nstamp2>0) {
	   newrel=(int*)malloc(MAX_SEQ_LEN*sizeof(int));
	   if(nstamp1>0) {
	       newnstamp=nstamp1; newnstamppos=nstamppos1; newnstampseq=nstampseq1;
	       sn=(struct stampdat*)malloc(10*sizeof(struct stampdat));
	       for(i=0; i<nstamp1; ++i) {
		  sn[i].what=s1[i].what;
                  sn[i].title=s1[i].title;
		  sn[i].n=(float*)malloc(MAX_SEQ_LEN*sizeof(float));
	       }
	   } else {
	      newnstamp=nstamp2; newnstamppos=nstamppos2; newnstampseq=nstampseq2;
	      sn=(struct stampdat*)malloc(10*sizeof(struct stampdat));
               for(i=0; i<nstamp2; ++i) {
                  sn[i].what=s2[i].what;
                  sn[i].title=s2[i].title;
                  sn[i].n=(float*)malloc(MAX_SEQ_LEN*sizeof(float));
               }
	   }
	} else {
	  newnstamp=0;
 	}


		

      if(n1>0 && n2>0) {
	for(i=1; i<=n1; ++i) {
	  new[i].ilen=b1[i].ilen;
	  new[i].id=b1[i].id;
	  new[i].tlen=b1[i].tlen;
	  new[i].title=b1[i].title;
	  new[i].slen=len;
	  new[i].seq=(char*)malloc(MAX_SEQ_LEN*sizeof(char));
	}
	k=0;
	for(j=1; j<=n2; ++j) if(use2[j-1]==1) {
	  i=n1+k+1;
	  new[i].ilen=b2[j].ilen;
          new[i].id=b2[j].id;
	  new[i].tlen=b2[j].tlen;
	  new[i].title=b2[j].title;
	  new[i].slen=len;
	  new[i].seq=(char*)malloc(MAX_SEQ_LEN*sizeof(char));
	  k++;
	}
	n2new=k;
	nn=n1+n2new;
	bpt1=bpt2=1;
	pt=1;
	while(bpt1<=b1[1].slen && bpt2<=b2[1].slen) {
	 /* If a gap occurs in one sequence and not the other
	  *   then we must place a gap in the other.  Otherwise
	  *   we must simply increment both pointers */
/*	 printf("%4d %4d %4d ",pt,bpt1,bpt2);
	 for(i=1; i<=n1; ++i) printf("%c",b1[i].seq[bpt1]);
	 printf(" ");
	 for(i=1; i<=n2; ++i) printf("%c",b2[i].seq[bpt2]);
	 printf("\n");
*/
	 if((nstamp1>0 || nstamp2>0)) { /* New stamp information */ 
		if((rel1[bpt1-1]==1) && (rel2[bpt2-1]==1)) {
		   newrel[pt-1]=1;
	           for(j=0; j<newnstamp; ++j) {
		      if(nstamp1>0) {
			sn[j].n[pt-1]=s1[j].n[bpt1-1];
		      } else {
			 sn[j].n[pt-1]=s2[j].n[bpt1-1];
		      }
		   }
		} else {
		  newrel[pt-1]=0;
		  for(j=0; j<newnstamp; ++j) {
		    if(sn[j].what=='T') sn[j].n[pt-1]=99;
		    else sn[j].n[pt-1]=0.0;
		  }
		}
	 }
	 if(b1[which1].seq[bpt1]==' ' && b2[which2].seq[bpt2]!=' ') { /* gap in first file */
	    for(i=1; i<=n1; ++i) 
	       new[i].seq[pt]=b1[i].seq[bpt1];

	    k=1;
	    for(i=1; i<=n2; ++i) if(use2[i-1]) {
	         new[n1+k].seq[pt]=' ';
		 k++;
	    }
	    bpt1++;
	 } else if(b1[which1].seq[bpt1]!=' ' && b2[which2].seq[bpt2]==' ') { /* gap in second file */
	    for(i=1; i<=n1; ++i) 
	       new[i].seq[pt]=' ';
	    k=1;
	    for(i=1; i<=n2; ++i) if(use2[i-1]) {
	       new[n1+k].seq[pt]=b2[i].seq[bpt2];
	       k++;
	    }
	    bpt2++;
	 } else if((b1[which1].seq[bpt1]!=' ' && b2[which2].seq[bpt2]!=' ') ||
	           (b1[which1].seq[bpt1]==' ' && b2[which2].seq[bpt2]==' ') ) { /* all gaps or no gaps */
	    for(i=1; i<=n1; ++i) {
	       new[i].seq[pt]=b1[i].seq[bpt1];
	    }
	    k=1;
	    for(i=1; i<=n2; ++i) if(use2[i-1]) {
	       new[n1+k].seq[pt]=b2[i].seq[bpt2];
	       k++;
	    }
	    bpt1++; bpt2++;
	 }
	 pt++;
	}

	/* Display the results */
	/* Get the new order -> this just puts all the secondary
	 *  structures after the sequences separted by the space */
	neworder=(int*)malloc(nn*sizeof(int));
        j=0;
	for(i=0; i<nn; ++i) { neworder[i]=-1; }
	for(i=0; i<nn; ++i) {
	   if((strstr(new[i+1].id,"_dssp")==NULL) && (strstr(new[i+1].id,"space")==NULL)) { neworder[j]=i; j++;}
	}
	for(i=0; i<nn; ++i) {
           if(strstr(new[i+1].id,"space")!=NULL) { printf("Space found at %4d\n",i); neworder[j]=i;  j++; break; } /* only count the first space */
        }
	for(i=0; i<nn; ++i) {
	   if(strstr(new[i+1].id,"_dssp")!=NULL) { neworder[j]=i; j++; }
        }
	printf("The new order is (j=%4d nn=%4d): ",j,nn);
	for(i=0; i<nn; ++i) if(neworder[i]!=-1) { printf("%4d ",neworder[i]); }
	printf("\n");

	


	   
	for(i=0; i<nn; ++i) new[i+1].slen=strlen(&new[1].seq[1])+1;

	for(i=0; i<nn; ++i) if(neworder[i]!=-1) {
          printf(">%s %s\n",new[neworder[i]+1].id,new[neworder[i]+1].title);
	}
	for(i=0; i<newnstamp; ++i) if(neworder[i]!=-1) {
	  printf("#%c %s",sn[i].what,sn[i].title);
	}
        printf("* iteration 1\n");
        for(k=1; k<new[1].slen; ++k) {
          allspace=1;
	  gapped=0;
          for(i=0; i<nn; ++i) if(neworder[i]!=-1) { 
		if(new[neworder[i]+1].seq[k]!=' ') allspace=0;
		if(new[neworder[i]+1].seq[k]==' ') gapped=1;
	  }
          if(!allspace) {
/*	      printf("%4d/%4d",k,new[1].slen); */
              for(i=0; i<nn; ++i) if(neworder[i]!=-1) {
                printf("%c",new[neworder[i]+1].seq[k]);
	      }
	      if(newnstamp>0) {
	       if(newrel[k-1]==1) { /* Output STAMP data */
		printf("  ");
		for(i=0; i<newnstamp; ++i) {
		  if(sn[i].what=='T') printf("%1.0f ",sn[i].n[k-1]);
                  else printf("%10.5f ",sn[i].n[k-1]);
		}
	       }
	      }
              printf("\n");
	  }
        }
        printf("*\n");
     }
     exit(0);
}

int numseq(FILE *fp) {
	/* note that this assumes that the vertical alignment is first */
	int n;
	char c;
	n=0;
	while((c=getc(fp))!=(char)EOF) n+=(c=='>');
	rewind(fp);
	return n;
}

void exit_error() {
 	fprintf(stderr,"format: mergestamp -f1 <file1> -f2 <file2> [options]\n");
	exit(-1);
}
