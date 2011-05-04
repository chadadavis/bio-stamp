#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <stamp.h>

#define min_R_diff 1.50
#define min_V_diff 60.0
#define max_start_diff 5

/* Reads in a list of domains and their corresponding Sc and other
 *  values and sorts them as specified */

struct info_struc {
  float Sc;
  float RMS;
  float value;
  float frac, seq_id,sec_id;
  int nfit;
  int nequiv;
  int len;
  int q_len;
  int d_len;
  int n_sec;
  int ignore;
  struct brookn fitpos;
  };

void compared_sc_hsort(int n,struct domain_loc* domain);
void comparei_sc_hsort(int n,struct info_struc* info);
void compared_rms_hsort(int n,struct domain_loc* domain);
void comparei_rms_hsort(int n,struct info_struc* info);



void exit_error();

main(int argc, char *argv[]) {

	char c;
	char *buff;

	int i,j,k;
	int method,gottrans;
	int ndomain,redundant,poor,left,similar;
	int ignore_trans,finicky;
	int pt;
	int oldout;

	char *filename;

	float cutoff;

	FILE *IN;

	struct domain_loc *domain;
	struct info_struc *info;

	char *env;

	filename=(char*)malloc(1000*sizeof(char));

	if(argc<3) exit_error();
	ignore_trans=0;
	finicky=1;
	method=0; cutoff=2.0;
	for(i=1; i<argc; ++i) {
	   if(argv[i][0]!='-') exit_error();
	   if(argv[i][1]=='f') {
	      /* file name */
	      strcpy(filename,argv[i+1]);
	      if((IN=fopen(argv[i+1],"r"))==NULL) {
	        fprintf(stderr,"error: file %s does not exist\n",argv[i+1]);
	        exit(-1);
	      }
	      i++;
	   } else if(argv[i][1]=='s') { /* sort type */
	     if(strcmp(argv[i+1],"Sc")==0) method=0;
	     else if(strcmp(argv[i+1],"rms")==0) method=1;
	     else if(strcmp(argv[i+1],"nfit")==0) method=2;
	     else if(strcmp(argv[i+1],"len")==0) method=3;
	     else if(strcmp(argv[i+1],"frac")==0) method=4;
	     else if(strcmp(argv[i+1],"q_frac")==0) method=5;
	     else if(strcmp(argv[i+1],"d_frac")==0) method=6;
	     else if(strcmp(argv[i+1],"n_sec")==0) method=7;
	     else if(strcmp(argv[i+1],"n_equiv")==0) method=10;
	     else if(strcmp(argv[i+1],"seq_id")==0) method=8;
	     else if(strcmp(argv[i+1],"sec_id")==0) method=9;
	     else exit_error();
	     sscanf(argv[i+2],"%f",&cutoff);
	     i+=2;
	   } else if(argv[i][1]=='i') {
	     ignore_trans=1;
	   } else if(argv[i][1]=='n') {
	     finicky=0;
	   } else exit_error();
	}


	if((env=getenv("STAMPDIR"))==NULL) {
           fprintf(stderr,"error: you haven't set the environment parameter STAMPDIR to anything\n");
           return -1;
      	}
	/* find the number of transformations */
	ndomain=0;
	ndomain=count_domain(IN);
	rewind(IN);

	/* now allocate memory */
	info=(struct info_struc*)malloc(ndomain*sizeof(struct info_struc));
	domain=(struct domain_loc*)malloc(ndomain*sizeof(struct domain_loc));
	buff=(char*)malloc(400*sizeof(char));

	/* read in domain descriptions */
	if(getdomain(IN,domain,&i,ndomain,&gottrans,env,0,stdout)==-1) exit(-1);
	if(i!=ndomain) {
	  fprintf(stderr,"error in domain descriptors\n");
	  exit(-1);
	}

	/* read in numbers */
	rewind(IN); i=0;
	while(fgets(buff,399,IN)!=NULL) {
	  pt=0;
	  if(buff[0]=='#') {
	     pt=skiptononspace(buff,pt);  pt=skiptononspace(buff,pt);
	     sscanf(&buff[pt],"%f",&info[i].Sc);
	     pt=skiptononspace(buff,pt);  pt=skiptononspace(buff,pt);
	     sscanf(&buff[pt],"%f",&info[i].RMS);
	     pt=skiptononspace(buff,pt); pt=skiptononspace(buff,pt);
	     sscanf(&buff[pt],"%d",&info[i].len);
	     pt=skiptononspace(buff,pt); pt=skiptononspace(buff,pt);
	     sscanf(&buff[pt],"%d",&info[i].nfit);
	     pt=skiptononspace(buff,pt);
             if(strncmp(&buff[pt],"fit_pos",7)!=0) { /* check for old format style */
		 /* read new format --- additional numbers */
                 pt=skiptononspace(buff,pt);
	         sscanf(&buff[pt],"%f",&info[i].seq_id);
	         pt=skiptononspace(buff,pt); pt=skiptononspace(buff,pt);
	         sscanf(&buff[pt],"%f",&info[i].sec_id);
	         pt=skiptononspace(buff,pt); pt=skiptononspace(buff,pt);	
	         sscanf(&buff[pt],"%d",&info[i].q_len);
	         pt=skiptononspace(buff,pt); pt=skiptononspace(buff,pt);
	         sscanf(&buff[pt],"%d",&info[i].d_len);
	         pt=skiptononspace(buff,pt); pt=skiptononspace(buff,pt);
	         sscanf(&buff[pt],"%d",&info[i].n_sec);
	         pt=skiptononspace(buff,pt); pt=skiptononspace(buff,pt);
		 sscanf(&buff[pt],"%d",&info[i].nequiv);
	         pt=skiptononspace(buff,pt);
		oldout=0;
	     } else {
		oldout=1;
		if(method>4) {
		   fprintf(stderr,"error: this scan file is apparently in old STAMP format\n");
		   fprintf(stderr,"       you will have to run the scan again to use the sort\n");
		   fprintf(stderr,"       option you have selected\n");
		   exit(-1);
		}
	     }
	     pt=skiptononspace(buff,pt);
	     sscanf(&buff[pt],"%c",&info[i].fitpos.cid);
	     pt=skiptononspace(buff,pt);
	     sscanf(&buff[pt],"%d",&info[i].fitpos.n);
	     pt=skiptononspace(buff,pt);
	     sscanf(&buff[pt],"%c",&info[i].fitpos.in);
	     i++;
	     if(i>ndomain) { 
		fprintf(stderr,"error: more numbers were found than domains\n");
		exit(-1);
	     }
	  }
	}
	
	/* now sort things out accoriding to Sc or RMS as required */
	for(i=0; i<ndomain; ++i) {
	   switch(method) {
              case 0:  domain[i].value=info[i].value=info[i].Sc; break;
	      case 1:  domain[i].value=info[i].value=info[i].RMS; break;
	      case 2:  domain[i].value=info[i].value=(float)info[i].nfit; break;
	      case 3:  domain[i].value=info[i].value=(float)info[i].len; break;
	      case 4:  domain[i].value=info[i].value=(float)info[i].nfit/(float)info[i].len; break;
	      case 5:  domain[i].value=info[i].value=(float)info[i].nfit/(float)info[i].q_len; break;
	      case 6:  domain[i].value=info[i].value=(float)info[i].nfit/(float)info[i].d_len; break;
	      case 7:  domain[i].value=info[i].value=(float)info[i].n_sec; break;
	      case 8:  domain[i].value=info[i].value=info[i].seq_id; break;
	      case 9:  domain[i].value=info[i].value=info[i].sec_id; break;
	      case 10: domain[i].value=info[i].value=(float)info[i].nequiv; break;
	      default: domain[i].value=info[i].value=1.0;
	   }
	}
	if(method!=1) {
	/* from stdlib.h
	 * extern void qsort(void *, size_t, size_t, int (*)(const void *, const void *)); */

	   compared_sc_hsort(ndomain,&(domain[-1]));
	   comparei_sc_hsort(ndomain,&(info[-1]));


	} else {
	   compared_rms_hsort(ndomain,&(domain[-1]));
	   comparei_rms_hsort(ndomain,&(info[-1]));
	}
	
	/* now let us determine which transformations are worth ignoring 
	 *  because the sort has already been done, only to top scoring
	 *  of the redundant domains will be left */
	for(i=0; i<ndomain; ++i) info[i].ignore=0;
	for(i=0; i<ndomain; ++i) {
	   k=0; while(k<strlen(domain[i].id) && domain[i].id[k]!='_') k++;
	   for(j=i+1; j<ndomain; ++j) {
	      if(ignore_trans==0) similar=transcompare(domain[i].R,domain[i].V,domain[j].R,domain[j].V,3);
	      else similar=1;
/*	      printf("%s and %s => %d\n",domain[i].id,domain[j].id,similar); */

	     if(!finicky) {
	       /* if we are not finicky, all that we require is that the transformations and the file name be
		*  the same  (since similar transformations of the same file will give the same result */
	       if(strcmp(domain[i].filename,domain[j].filename)==0 &&  /* the files are the same */
		 (!gottrans || transcompare(domain[i].R,domain[i].V,domain[j].R,domain[j].V,3))) {
		    info[j].ignore=1;
	       }
	     } else {
		 /* if we are finicky, then we require an exact match */
		 if(strcmp(domain[i].filename,domain[j].filename)==0 &&  /* the files are the same */
		 strncmp(domain[i].id,domain[j].id,k)==0 &&  /* the IDs are the same */
	         domain[i].start[0].cid == domain[j].start[0].cid &&
		 abs(domain[i].start[0].n  -   domain[j].start[0].n)<max_start_diff   &&
		 domain[i].start[0].in  == domain[j].start[0].in  &&
		 domain[i].end[0].cid   == domain[j].end[0].cid &&
		 abs(domain[i].end[0].n -  domain[j].end[0].n)<max_start_diff   &&
		 domain[i].end[0].in    == domain[j].end[0].in &&  /* the start and end are the same */
		 (!gottrans || similar))  {  
		    /* the transformations are sufficiently similar */
		   info[j].ignore=1;
	         }
	     }
	   }
	}
	/* now remove things that are below (Sc) or above (RMS) the cutoff */
	poor=0; 
	for(i=0; i<ndomain; ++i) {
	   if((method==1 && info[i].value>cutoff) ||
	      (method!=1 && info[i].value<cutoff) ) {
		info[i].ignore=1;
		poor++;
	  }
	}
	   
	redundant=left=0;
	for(i=0; i<ndomain; ++i) {
	  redundant+=info[i].ignore;
	  left+=!(info[i].ignore);
	}
	/* output the results */
	printf("%% Sorted output from STAMP scan routine\n");
	printf("%% The file %s was sorted according to \n",filename);
        switch(method) {
              case 0: printf("%% Sc values\n"); break;
              case 1: printf("%% RMS deviations\n"); break;
              case 2: printf("%% number of atoms used to fit\n"); break;
	      case 3: printf("%% alignment length\n"); break;
              case 4: printf("%% fraction of nfit/len\n"); break;
              case 5: printf("%% fraction of nfit/q_len\n"); break;
              case 6: printf("%% fraction of nfit/d_len\n"); break;
              case 7: printf("%% number of equivalent sec. structures. \n"); break;
              case 8: printf("%% sequence identity\n"); break;
              case 9: printf("%% sec. structure identity\n"); break;
	      case 10: printf("%% number of STAMP equivalent residues\n"); break;
        }
	printf("%% Repeated transformations were ignored\n");
	printf("%% Out of an initial %d domains,\n",ndomain);
	printf("%% a total of %d repeated domains were found\n",redundant);
	printf("%% and %d domains were found to have \n",poor);
	switch(method) {
          case 0: printf("%%  Sc values less than %7.3f\n",cutoff); break;
	  case 1: printf("%%  RMS values greater than %7.3f\n",cutoff); break;
	  case 2: printf("%%  less than %4d atoms used in the fit\n",(int)cutoff); break;
	  case 3: printf("%%  alignment lengths less than %4d\n",(int)cutoff); break;
	  case 4: case 5: case 6: printf("%%  fraction of atoms fitted less than %7.3f\n",cutoff); break;
	  case 7: printf("%% fewer than %4d equivalent sec. structures\n",(int)cutoff); break;
	  case 8: case 9: printf("%% identities less than %5.2f %%\n",cutoff); break;
	}
	printf("%% All of these were removed\n");
	printf("%% Leaving %d domains in this file\n",left);
	
	for(i=0; i<ndomain; ++i) {
	  if(!info[i].ignore) {
	     if(oldout==0) {
	      printf("# Sc= %7.3f RMS= %7.3f len= %4d nfit= %4d seq_id= %5.2f sec_id= %5.2f q_len= %4d d_len= %4d n_sec=%4d fit_pos= %c %3d %c\n",
	       info[i].Sc,info[i].RMS,info[i].len,info[i].nfit,
	       info[i].seq_id,info[i].sec_id,info[i].q_len,info[i].d_len,info[i].n_sec,
	      info[i].fitpos.cid,info[i].fitpos.n,info[i].fitpos.in);
	     } else {
	      printf("# Sc= %7.3f RMS= %7.3f len= %4d nfit= %4d fit_pos= %c %3d %c\n",
		info[i].Sc,info[i].RMS,info[i].len,info[i].nfit,
		info[i].fitpos.cid,info[i].fitpos.n,info[i].fitpos.in);
	     }
	     printdomain(stdout,domain[i],gottrans);
	   } 
	}
	exit(0);
}

/* Compares two transformations.  Returns 1 if
 *  they are sufficiently similar to be called the
 *  same thing */

/* The parameters min_R_diff=1.50 min_V_diff=60.0 
 *  were not derived very carefully.  I simply looked at
 *  what Rdiffsq and Vdiffsq were like for transformations that I took to
 *  be the same (ie. the same domain).  These are the upper limits I saw, and 
 *  appear to work pretty well, though they will occasionally lead to two
 *  of the "same" transformations in a file.
 * These values are essentially root mean square differences for
 *  matrix and vector components respectively. */


int transcompare(float **R1, float *V1, float **R2, float *V2, int dim )
{	
	int i,j,k;
	float Rdiffsq,Rdiff;
	float Vdiffsq,Vdiff;

	Rdiffsq = Vdiffsq = 0.0;

	for(i=0; i<dim; ++i) {
	  for(j=0; j<dim; ++j) Rdiffsq+=(R1[i][j]-R2[i][j])*(R1[i][j]-R2[i][j]);
	  Vdiffsq+=(V1[i]-V2[i])*(V1[i]-V2[i]);
	}
	Rdiff = sqrt(Rdiffsq);
	Vdiff = sqrt(Vdiffsq);

	return (Rdiff<=min_R_diff && Vdiff<min_V_diff);
}
void exit_error()
{
	  fprintf(stderr,"format: sorttrans -f <filename> [-s <sort method> <value> -i -n] \n");
	  fprintf(stderr,"  sort (-s) keywords: \n");
	  fprintf(stderr,"      Sc sort by STAMP Sc value\n");
	  fprintf(stderr,"      rms   sort by RMS deviations\n");
	  fprintf(stderr,"      nfit  sort by number of atoms fitted\n");
	  fprintf(stderr,"      len   sort by alignment length\n");
	  fprintf(stderr,"      frac  sort by fraction of nfit/len\n");
	  fprintf(stderr,"      q_frac sort by fraction of nfit/(query length)\n");
	  fprintf(stderr,"      d_frac sort by fraction of nfit/(databse length)\n");
	  fprintf(stderr,"      seq_id sort by sequence identity\n");
	  fprintf(stderr,"      sec_id sort by secondary structure identity\n");
	  fprintf(stderr,"      n_sec  sort by number of equivalent secondary structures\n");
	  fprintf(stderr,"  -i only allow the best superimposition per identifier\n");
	  fprintf(stderr,"  -n ignore domain descriptors (just filename and transformation)\n");
	  exit(-1);
}
void compared_sc_hsort(int n,struct domain_loc* domain) {
	int l,j,ir,i;
	struct domain_loc rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
			rra=domain[--l];
		else {
			rra=domain[ir];
			domain[ir]=domain[1];
			if (--ir == 1) {
				domain[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && domain[j].value > domain[j+1].value) ++j;
			if (rra.value > domain[j].value){
				domain[i]=domain[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		domain[i]=rra;
	}
}

void comparei_sc_hsort(int n,struct info_struc* info) {
	int l,j,ir,i;
	struct info_struc rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
			rra=info[--l];
		else {
			rra=info[ir];
			info[ir]=info[1];
			if (--ir == 1) {
				info[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && info[j].value > info[j+1].value) ++j;
			if (rra.value > info[j].value){
				info[i]=info[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		info[i]=rra;
	}
}
void compared_rms_hsort(int n,struct domain_loc* domain) {
	int l,j,ir,i;
	struct domain_loc rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
			rra=domain[--l];
		else {
			rra=domain[ir];
			domain[ir]=domain[1];
			if (--ir == 1) {
				domain[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && domain[j].value < domain[j+1].value) ++j;
			if (rra.value < domain[j].value){
				domain[i]=domain[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		domain[i]=rra;
	}
}

void comparei_rms_hsort(int n,struct info_struc* info) {
	int l,j,ir,i;
	struct info_struc rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
			rra=info[--l];
		else {
			rra=info[ir];
			info[ir]=info[1];
			if (--ir == 1) {
				info[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && info[j].value < info[j+1].value) ++j;
			if (rra.value < info[j].value){
				info[i]=info[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		info[i]=rra;
	}
}
