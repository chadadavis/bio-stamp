#include <stamp.h>
#include <math.h>

#define MAX_SEQ_LEN   10000
#define PRECISION  1000

#define lastmod "5 March 2003"

void exit_error();
void help_exit_error();

main(int argc, char *argv[]) {

	int i,j,k,l,test;
	int ndomain,total,add;
	int gottrans,verbose;
	int really_verbose;
        int N_within,N_within_min_dist;
	int T_FLAG;
        int already_used1[MAX_SEQ_LEN];
        int already_used2[MAX_SEQ_LEN];

	char c;
	char *env;
	char *deffile,*keyword,*value;
        char domfile[100];

        float dist,min_dist,D;

	FILE *PARMS,*TRANS,*PDB;

	struct domain_loc *domain;

	if(argc<2) exit_error();

        N_within_min_dist = 5;
        min_dist = 10.0;
        verbose=0;
        really_verbose=0;
	/* now search the command line for commands */
	keyword=(char*)malloc(1000*sizeof(char));
	value=(char*)malloc(1000*sizeof(char));
        for(i=1; i<argc; ++i) {
           if(argv[i][0]!='-') exit_error();
	   strcpy(keyword,&argv[i][1]);
	   if(i+1<argc) strcpy(value,argv[i+1]);
	   else strcpy(value,"none");
	   for(j=0; j<strlen(keyword); ++j) {
	      keyword[j]=ltou(keyword[j]); /* change to upper case */
           }
	   T_FLAG=(value[0]=='Y' || value[0]=='y' || value[0]=='1' || 
		 value[0]=='T' || value[0]=='t' || value[0]=='o' || 
		 value[0]=='O');
	   /* enables one to write '1', 'YES', 'Yes', 'yes', 'T_FLAG', 'True' or 'true' to 
	    *  set any boolean parmsiable to one */
           if((strcmp(&argv[i][1],"l")==0) || (strcmp(&argv[i][1],"f")==0) || (strcmp(&argv[i][1],"p")==0)) {
	      if(i+1>=argc) exit_error();
              /* listfile name */
	      strcpy(domfile,argv[i+1]);
	      i++;
           } else if(strcmp(&argv[i][1],"min_dist")==0) {
	      if(i+1>=argc) exit_error();
	      sscanf(argv[i+1],"%f",&min_dist);
	      i++;
           } else if(strcmp(&argv[i][1],"N")==0) {
	      if(i+1>=argc) exit_error();
	      sscanf(argv[i+1],"%d",&N_within_min_dist);
	      i++;
           } else if(strcmp(&argv[i][1],"v")==0) {
              verbose=1;
           } else if(strcmp(&argv[i][1],"rv")==0) {
              really_verbose=1; /* Reports all distances */
	   } else  {
	     exit_error();
	   }
	}
	free(keyword); 
	free(value);
	
	/* read in coordinate locations and initial transformations */
	if((TRANS = fopen(domfile,"r")) == NULL) {
	   fprintf(stderr,"error: file %s does not exist\n",domfile);
	   exit(-1);
	}
	/* determine the number of domains specified */
	ndomain=count_domain(TRANS);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	rewind(TRANS);
	if(getdomain(TRANS,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	fclose(TRANS);

	if(verbose==1) { fprintf(stdout,"Reading coordinates...\n"); }
	for(i=0; i<ndomain; ++i) {
	   if(verbose==1) { fprintf(stdout,"Domain %3d %s %s\n   ",i+1,domain[i].filename,domain[i].id); }
	   if((PDB=openfile(domain[i].filename,"r"))==NULL) {
	      fprintf(stderr,"error opening file %s\n",domain[i].filename);
	      exit(-1);
	   }
	   domain[i].ncoords=0;
	   domain[i].coords=(int**)malloc(MAX_SEQ_LEN*sizeof(int*));
	   domain[i].aa=(char*)malloc((MAX_SEQ_LEN+1)*sizeof(char)); 
	   domain[i].numb=(struct brookn*)malloc((MAX_SEQ_LEN)*sizeof(struct brookn));
	   total=0;
	   if(verbose==1) { fprintf(stdout,"    "); }
	   for(j=0; j<domain[i].nobj; ++j) {
	       if(igetca(PDB,&domain[i].coords[total],&domain[i].aa[total],&domain[i].numb[total],
		    &add,domain[i].start[j],domain[i].end[j],domain[i].type[j],(MAX_SEQ_LEN-total),
		    domain[i].reverse[j],PRECISION,stdout)==-1) {
		    fprintf(stderr,"Error in domain %s object %d \n",domain[i].id,j+1);
                    exit(-1);
	       }
               if(verbose==1) {
  	          switch(domain[i].type[j]) {
	 	     case 1: fprintf(stdout," all residues"); break;
		     case 2: fprintf(stdout," chain %c",domain[i].start[j].cid); break;
		     case 3: fprintf(stdout," from %c %4d %c to %c %4d %c",
			 domain[i].start[j].cid,domain[i].start[j].n,domain[i].start[j].in,
			 domain[i].end[j].cid,domain[i].end[j].n,domain[i].end[j].in); break;
	   	   }
		   fprintf(stdout,"%4d CAs ",add);
                }
	        total+=add;
		closefile(PDB,domain[i].filename); PDB=openfile(domain[i].filename,"r");
	    }
	    domain[i].ncoords=total;
            if(verbose==1) {
	       fprintf(stdout,"=> %4d CAs in total\n",domain[i].ncoords);
	       fprintf(stdout,"Applying the transformation... \n");
	       printmat(domain[i].R,domain[i].V,3,stdout);
	       fprintf(stdout,"      ...to these coordinates.\n");
            }
	    matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,PRECISION);
	    closefile(PDB,domain[i].filename);
	}
	if(verbose) { fprintf(stdout,"\n\n"); }
	for(i=0; i<ndomain; ++i) {
	   for(j=i+1; j<ndomain; ++j) {
               N_within=0;
               for(k=0; k<domain[i].ncoords; ++k) { already_used1[k]=0; }
               for(k=0; k<domain[i].ncoords; ++k) {
                  for(l=0; l<domain[j].ncoords; ++l) { already_used2[l]=0; }
                  for(l=0; l<domain[j].ncoords; ++l) {
                      D=idist(domain[i].coords[k],domain[j].coords[l],PRECISION);
                      if((really_verbose==1) && (D<=min_dist)) { 
                        printf("DIST_WITHIN_THRESH: %-20s %c %4d %c : %7.3f : %-20s %c %4d %c \n",
                           domain[i].id,domain[i].numb[k].cid,domain[i].numb[k].n,domain[i].numb[k].in,D,
                           domain[j].id,domain[j].numb[l].cid,domain[j].numb[l].n,domain[j].numb[l].in);
                      }
                      if((already_used1[k]==0) && (already_used2[l]==0) && (D<=min_dist)) {
                           N_within++; 
                           already_used1[k]=1;
                           already_used2[l]=1;
                      }
                  }
               }
               if(N_within>=N_within_min_dist) { 
                  if(verbose) { 
                        printf("Domains %20s %c and %20s %c ARE in contact\n",
                           domain[i].id,(char)('A'+(char)i),domain[j].id,(char)('A'+(char)j)); 
                  } else {
                        printf("%-20s %-20s %4d\n",domain[i].id,domain[j].id,N_within);
                  }
                           
               } else {
                  if(verbose) { 
                        printf("Domains %20s %c and %20s %c NOT in contact\n",
                           domain[i].id,(char)('A'+(char)i),domain[j].id,(char)('A'+(char)j)); 
                  } 
               }
           }
        }
               
         
	
	/* freeing memory to keep purify happy */
	for(i=0; i<ndomain; ++i) {
	   free(domain[i].aa);
	   free(domain[i].sec);
	   free(domain[i].v); free(domain[i].V);
	   for(j=0; j<3; ++j) {
	      free(domain[i].R[j]);
	      free(domain[i].r[j]);
	   }
	   free(domain[i].R); 
	   free(domain[i].r);
	   for(j=0; j<domain[i].ncoords; ++j) 
	      free(domain[i].coords[j]);
	   free(domain[i].coords);
	   free(domain[i].type);
	   free(domain[i].start);
	   free(domain[i].end);
	   free(domain[i].reverse);
	   free(domain[i].numb);
	}
	free(domain);

	exit(0);
}

void exit_error() {
	fprintf(stderr, 
"Usage:\n"
"  check_ints -f <domain file> [-v] [-rv] [-N <int>] [-min_dist <float>]\n"
"\n"
"  -min_dist  Residues in contact when <= min_dist (A)\n"
"  -N         Number of residues that need to be <= min_dist\n"
"  -v         verbose (show transformations applied)\n"
"  -rv        really verbose (show residue contacts)\n"
"\n"
);
	exit(-1);
}
