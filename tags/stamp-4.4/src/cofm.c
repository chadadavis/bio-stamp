#include "R.h"

/* 
 *  Reads in a domain descriptor and returns
 *   various commodities (centre of mass, radius of gyration, etc.)
 */

int exit_error() {

        fprintf(stderr,"format: cofm -f <dom file> [-v -file]\n");
        exit(-1);
} 
main(int argc, char *argv[]) {

	int i,j;
	int ndomain;
	int verbose;
	int gottrans,ct;
	int total,add;
        int add_file_name;
        int seed,n;

	float total_mass;
        float *Rg;
        float *Rm;
        float diameter;

	int **Ro;

	static char *env;
	char domfile[1000];
	char dssp_files[1000],pdb_files[1000];

	FILE *DOM,*PDB;
	struct domain_loc *domain;

	if(argc<3) exit_error();

	verbose=0;

	if((env=getenv("STAMPDIR"))==NULL) {
          fprintf(stderr,"error: environment variable STAMPDIR must be set\n");
          exit(-1);
        }
	sprintf(pdb_files,"%s/pdb.directories",env);
	sprintf(dssp_files,"%s/dssp.directories",env);

        seed = 0;
        ct = 0;
        diameter = 5.0;
        add_file_name=0;
	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(strcmp(&argv[i][1],"f")==0) {
	     if((i+1)>=argc) exit_error();
	     strcpy(domfile,argv[i+1]);
	     i++;
	   } else if(strcmp(&argv[i][1],"seed")==0) {
	     if((i+1)>=argc) exit_error();
	     sscanf(argv[i+1],"%d",&seed);
	     i++;
	   } else if(strcmp(&argv[i][1],"v")==0) {
	     verbose=1;
	   } else if(strcmp(&argv[i][1],"ct")==0) {
	     ct=1;
	   } else if(strcmp(&argv[i][1],"file")==0) {
	     add_file_name=1;
	   } else {
	     exit_error();
	   }
	}


	if((DOM=fopen(domfile,"r"))==NULL) {
	  fprintf(stderr,"error opening file %s\n",domfile);
	  exit(-1);
	}

	/* determine the number of domains specified */
	ndomain=count_domain(DOM);
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	rewind(DOM);
	if(getdomain(DOM,domain,&ndomain,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	fclose(DOM);

	if(verbose==1) { fprintf(stdout,"Reading coordinates...\n"); }
	for(i=0; i<ndomain; ++i) {
	   if(verbose==1) { fprintf(stdout,"Domain %3d %s %s\n   ",i+1,domain[i].filename,domain[i].id); }
	   if((PDB=(FILE *)openfile(domain[i].filename,"r"))==NULL) {
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
		    domain[i].reverse[j],PRECIS,stdout)==-1) {
		    fprintf(stderr,"Error in domain %s object %d \n",domain[i].id,j+1);
                    exit(-1);
	       }
           if(verbose==1) {
        	   if (domain[i].start[j].in == ' ') { domain[i].start[j].in = '_'; }
        	   if (domain[i].end[j].in == ' ')   { domain[i].end[j].in   = '_'; }
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
	    matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,PRECIS);
	    closefile(PDB,domain[i].filename);
	}
	if(verbose) { fprintf(stdout,"\n\n"); }
	

	Ro=(int**)malloc(ndomain*sizeof(int));
	Rg=(float*)malloc(ndomain*sizeof(float));
	Rm=(float*)malloc(ndomain*sizeof(float));
	for(i=0; i<ndomain; ++i) {
		/* Chain ID, or rather the chain of the first segment, when multiple */
		char chainid = domain[i].start[0].cid;
		/* Increment temperature factor, to distinguish orientations visually */
		float temp = 0;

		/* Memory leak? */
	   Ro[i]=(int*)malloc(3*sizeof(int));
	   Ro[i] = RBR_c_of_m_pep(domain[i].coords,domain[i].ncoords,domain[i].aa,1,&total_mass);
           Rg[i] = RBR_r_of_gyration_pep(Ro[i],domain[i].coords,domain[i].ncoords,domain[i].aa,1,PRECIS,2);
           Rm[i] = RBR_r_max_pep(Ro[i],domain[i].coords,domain[i].ncoords,domain[i].aa,PRECIS,2);
           printf("REMARK Domain %3d Id %10s N = %d Rg = %7.3f Rmax = %7.3f Ro   = ", 
               i+1,domain[i].id,domain[i].ncoords,Rg[i],Rm[i]);
           printf("%8.3f %8.3f %8.3f",(float)Ro[i][0]/(float)PRECIS,(float)Ro[i][1]/(float)PRECIS,(float)Ro[i][2]/(float)PRECIS);
           if(add_file_name==1) {
               printf(" %s",domain[i].filename);
           } 
           printf("\n");

                /* ATOM      2  CA  ALA A   7      25.400  -4.374  40.370  1.00 74.92           C          */
           printf("ATOM      0  CA  ALA %c   0    %8.3f%8.3f%8.3f  1.00  0.00\n",
        		  chainid,
        		((float)Ro[i][0]/(float)PRECIS),
                ((float)Ro[i][1]/(float)PRECIS),
                ((float)Ro[i][2]/(float)PRECIS));
           printf("ATOM      1  CA  ALA %c   1    %8.3f%8.3f%8.3f  1.00  4.00\n",
         		  chainid,
	        ((float)Ro[i][0]/(float)PRECIS)+diameter,
                ((float)Ro[i][1]/(float)PRECIS),
                ((float)Ro[i][2]/(float)PRECIS));
           printf("ATOM      1  CA  ALA %c   1    %8.3f%8.3f%8.3f  1.00  0.00\n",
         		  chainid,
	        ((float)Ro[i][0]/(float)PRECIS)-diameter,
                ((float)Ro[i][1]/(float)PRECIS),
                ((float)Ro[i][2]/(float)PRECIS));
           printf("ATOM      2  CA  ALA %c   2    %8.3f%8.3f%8.3f  1.00 12.00\n",
         		  chainid,
	        ((float)Ro[i][0]/(float)PRECIS),
                ((float)Ro[i][1]/(float)PRECIS)+diameter,
                ((float)Ro[i][2]/(float)PRECIS));
           printf("ATOM      2  CA  ALA %c   2    %8.3f%8.3f%8.3f  1.00  0.00\n",
         		  chainid,
	        ((float)Ro[i][0]/(float)PRECIS),
                ((float)Ro[i][1]/(float)PRECIS)-diameter,
                ((float)Ro[i][2]/(float)PRECIS));
           printf("ATOM      3  CA  ALA %c   3    %8.3f%8.3f%8.3f  1.00 20.00\n",
         		  chainid,
	        ((float)Ro[i][0]/(float)PRECIS),
                ((float)Ro[i][1]/(float)PRECIS),
                ((float)Ro[i][2]/(float)PRECIS)+diameter);
           printf("ATOM      3  CA  ALA %c   3    %8.3f%8.3f%8.3f  1.00  0.00\n",
         		  chainid,
	        ((float)Ro[i][0]/(float)PRECIS),
                ((float)Ro[i][1]/(float)PRECIS),
                ((float)Ro[i][2]/(float)PRECIS)-diameter);

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
	   for(j=0; j<domain[i].ncoords; ++j)  { free(domain[i].coords[j]); }
	   free(domain[i].coords);
	   free(domain[i].type);
	   free(domain[i].start);
	   free(domain[i].end);
	   free(domain[i].reverse);
	   free(domain[i].numb);
	}
	free(domain);
}
