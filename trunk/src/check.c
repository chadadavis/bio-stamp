/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)	            Prof. Geoffrey J. Barton (GJB)
 EMBL Heidelberg                            School of Life Sciences
 Meyerhofstrasse 1                          University of Dundee
 D-69117 Heidelberg                         Dow Street
 Germany                                    Dundee, DD1 5EH
                                          
 Tel: +49 6221 387 473                      Tel: +44 1382 345860
 FAX: +44 6221 387 517                      FAX: +44 1382 345764
 E-mail: russell@embl-heidelberg.de         E-mail geoff@compbio.dundee.ac.uk
 WWW: http://www.russell.emb-heidelberg.de  WWW: http://www.compbio.dundee.ac.uk

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "stamp.h"

#define RES "REMARK   2 RESOLUTION."
#define REF "REMARK   3"

/* This program checks pdb files for inconsitencies, and
 *  other stuff */
int main(int argc, char *argv[]) {

	int i,j,k,l,mode,yes,n_main_miss,total_main_miss;
	int ftype;
	int year;
	int found;
	char c;
	int NMR,MODEL,REFINED;
	FILE *pdb,*in;
	char keyword[7],atnum[6],atname[5],
	     resname[4],resnum[7],rest[15],
	     buff[100],oldresnum[7],code[200],
	     reftext[4000];
	char *pdbfile,*dsspfile;
	char *stampdir, *dirfile;
	int nlines;
	int reflen;
	int nchains,nres,natoms,nchainres,nchainatoms;
	char curchain;
	char *chs;
	int nspecchains;
	int start,end;
	float x,y,z,Occ,B,resolution,R_factor;

	int N,C,CA,O,ACE,FOR,HET;
	int new,n_ref_type,n_r_val;
	
	static char *ref_type[]={ 
	   "XPLOR",
	   "PROLSQ",
	   "EREF",
	   "FRODO",
	   "CORELS",
	   "NONE.",
	   "JACK",
	   "LEVITT",
	   "DIAMOND",
	   "BRUNGER",
	   "KURIYAN",
	   "KARPLUS",
	   "KONNERT",
	   "HENDRICKSON",
	   "SUSMAN"
	};

	static char *r_val[]={
	   "R VALUE IS ABOUT",
	   "R VALUE IS APPROXIMATELY",
	   "R VALUE IS",
	   "R VALUE",
           "R-VALUE IS",
	   "R-VALUES ARE",
	   "R VALUES ARE",
	   "R-VALUE",
	   "R-FACTOR",
	   "R-FACTOR IS",
	   "R-FACTOR IS APPROXIMATELY",
	   "R FACTOR",
	   "R FACTOR IS",
	   "R FACTOR IS APPROXIMATELY"
 	};

	n_ref_type=15;
	n_r_val=14;

	NMR=0; MODEL=0;  REFINED=0;

	/* PDB format in FORTRAN is
	 * FORMAT(A6,A5,2X,A4,A3,1X,A6,3X,3F8.3,2F6.2,A14)
	 */
	if((stampdir=getenv("STAMPDIR"))==NULL) {
	  fprintf(stderr,"error: environment variable STAMPDIR must be specified\n");
	  exit(-1);
	}
	if((argc!=3) || (argv[1][0] != '-')) {
	   fprintf(stderr,"format: pdbc -q/m/d/n/r <PDB code>\n");
	   fprintf(stderr,"        -q verbose mode\n");
	   fprintf(stderr,"        -m minimalist mode (just report file locations)\n");
	   fprintf(stderr,"        -d write in STAMP database format\n");
	   fprintf(stderr,"        -r write out chain, resolution, R-factor, etc.\n");
	   exit(-1);
	   }
	switch(argv[1][1]) {
	   case 'q': mode=1; break;
	   case 'd': mode=2; break;
	   case 'm': mode=4; break;
	   case 'r': mode=5; break;
	   default: {
	     fprintf(stderr,"error: mode %c not recognised\n",argv[1][1]);
	     exit(-1);
	   }
	}

	if(mode==5 && strlen(argv[2])>4) {
	  printf("Warning: chains aren't considered in this mode\n");
	}

	dirfile=(char*)malloc((strlen(stampdir)+strlen("pdb.directories")+2)*sizeof(char));
	if((in=fopen("pdb.directories","r"))!=NULL) { /* directory file found in home directory */
	   sprintf(dirfile,"pdb.directories");
	   fclose(in);
	} else { /* otherwise look in STAMP directory */
	   sprintf(dirfile,"%s/pdb.directories",stampdir);
	}
	/* try whole string first */
	pdbfile=getfile(argv[2],dirfile,strlen(argv[2]),stdout);
	ftype=1;
	if(pdbfile[0]=='\0') {
	   free(pdbfile);
	   /* now just try the first four characters */
	   strncpy(&code[0],argv[2],4); code[4]='\0';
	   for(j=0; j<strlen(code); ++j) code[j]=code[j];
	   pdbfile=getfile(code,dirfile,4,stdout);
	   ftype=0;
	}
	if(pdbfile[0]=='\0') {
	   fprintf(stderr,"error: no file found for %s\n",argv[2]);
	   exit(-1);
	}
	/* find DSSP file if necessary */
	if(mode==3 || mode==4) {
	   sprintf(dirfile,"%s/dssp.directories",stampdir);
	   dsspfile=getfile(argv[2],dirfile,strlen(argv[2]),stdout);
	   if(dsspfile[0]=='\0') {
              free(dsspfile);
              /* now just try the first four characters */
              strncpy(&code[0],argv[2],4); code[4]='\0';
              for(j=0; j<strlen(code); ++j) code[j]=code[j];
              dsspfile=getfile(code,dirfile,4,stdout);
	   }
	   if(dsspfile[0]=='\0') {
      	      strcpy(dsspfile,"Unknown");
	   }
	}
	free(dirfile);
	if(mode==4) {
	   printf("CODE %s PDB %s DSSP %s\n",argv[2],pdbfile,dsspfile);
	   exit(-1);
	  }
    nlines=0;
    curchain='-';
    nchains=0;
    if(mode==1) printf("searching file: %s ",pdbfile);
    if(ftype==0) {
       strncpy(&code[0],argv[2],4);  code[4]='\0';
       for(j=0; j<strlen(code); ++j) code[j]=code[j];
    } else {
       end=strlen(argv[2])-1; start=0;
       for(j=0; j<strlen(argv[2]); ++j) {
	  if(argv[2][j]=='/') start=j+1;
	  if(argv[2][j]=='.') end=j-1;
       }
       strncpy(&code[0],&argv[2][start],(end-start+1));
       code[end-start+1]='\0';
       for(j=0; j<strlen(code); ++j) code[j]=code[j];
    }
	

    if(ftype==0 && strlen(argv[2])>4) nspecchains=strlen(argv[2])-4;
    else nspecchains=0;
    if(nspecchains>0) {
       chs=(char*)malloc(nspecchains*sizeof(char));
       for(i=0; i<nspecchains; ++i) {
	  chs[i]=argv[2][4+i];
	  if(mode==1) printf("chain %c ",chs[i]);
       }
    }
    if(mode==1) printf("\n\n");

    new=0;
    N=C=CA=O=ACE=FOR=1;
    HET=0;
    nres=natoms=nchainres=nchainatoms=n_main_miss=total_main_miss=0;

    pdb=openfile(pdbfile,"r");

    c=' ';
    while(c!=(char)EOF) {
	i=0;
	while((c=getc(pdb))!= '\n' && c!=(char)EOF) {
	   buff[i++]=c; 
	   }
	if(c==(char)EOF) break;
	buff[i]='\0';
	if(strncmp(buff,"COMPND",6)==0 || strncmp(buff,"TITLE ",6)==0 || 
	       strncmp(buff,"HEADER",6)==0 || strncmp(buff,"SOURCE",6)==0 ||
	       strncmp(buff,"KEYWDS",6)==0 || strncmp(buff,"EXPDTA",6)==0) {
	       /* check if NMR or model */
	       if(strstr(buff,"NMR")!=NULL) NMR=1;
	       if(strstr(buff,"THEORETICAL MODEL")!=NULL) MODEL=1;
	}
	if((strncmp(buff,"HEADER",6)==0 || strncmp(buff,"COMPND",6)==0 ||
	    strncmp(buff,"AUTHOR",6)==0 || strncmp(buff,"SOURCE",6)==0) ) {
	    if(mode==1) printf("%s\n",buff);
	    if(mode==2 || mode==3) printf("%% %s\n",&buff[6]);
	    if(strncmp(buff,"HEADER",6)==0) { /* get year */
/* 012345678901234567890123456789012345678901234567890123456789 */
/* HEADER    OXYGEN STORAGE                          14-JAN-88   4MBN      4MBN   3*/
	       sscanf(&buff[57],"%d",&year);
	       /* Year 2000 compliancy (worth 300K per annum apparently).  
	        * Note that the PDB is not compliant. */
	       if(year>70) {
	         year+=1900;
	       } else {
	         year+=2000;
	       }
	    }
	}
	
	if(strncmp(buff,"ENDMDL",6)==0) {
		if(mode==1) {
		   printf("ENDMDL encountered, ignoring rest of file\n");
		} else if(mode==2 || mode==3) {
		   printf("%% ENDMDL ==> NMR structure\n");
		   NMR=1;
	 	}
		
		break;
	}
	nlines++;
	/* Copy the appropriate information to the
	 *  appropriate variables. */

	strncpy(keyword,&buff[0],6);
	keyword[6]='\0';
	if(strncmp(keyword,"ATOM  ",6)==0) {
	natoms++; nchainatoms++;
	/* check to see if the chain has changed */
	if(buff[21]!=curchain) { 
	   if(nchains>0) {
	      if(mode==1) printf("        total number of residues: %d, atoms: %d\n",nchainres,nchainatoms);
	      if(n_main_miss>0) {
		   if(mode==1) printf("        missing main chain atoms for %d residues\n",n_main_miss);
		   if(mode==2 || mode==3) {
		      yes=(nspecchains==0);
		      for(i=0; i<nspecchains; ++i) 
			 if(curchain==chs[i]) printf("%% chain %c missing main chain atoms for %d residues\n",curchain,n_main_miss+1);
	           }
	      }
	      total_main_miss+=n_main_miss;
	      nchainres=0; nchainatoms=n_main_miss=0;
	   }
	   nchains++;
	   curchain=buff[21];
	   if(mode==1) printf("  chain: `%c'\n",buff[21]);
	   if(mode==2 || mode==3) {
	      if(nspecchains>0) {
		 yes=0;
		 for(i=0; i<nspecchains; ++i) if(chs[i]==buff[21]) yes=1;
	      } else yes=1;
	      if(yes) {
	       if(mode==2) {
	        if(buff[21]==' ') printf("%s %s { ALL }\n",pdbfile,code);
	        else printf("%s %s%c { CHAIN %c }\n",pdbfile,code,(char)buff[21],(char)buff[21]);
	       } else {
                if(buff[21]==' ') printf("%s %s %s { ALL }\n",pdbfile,dsspfile,code);
                else printf("%s %s %s%c { CHAIN %c }\n",pdbfile,dsspfile,code,(char)buff[21],(char)buff[21]);
	       }
	      }
	    }
	  }
	  strncpy(atnum,&buff[6],5);
	  atnum[5]='\0';
	  strncpy(resnum,&buff[21],6);
	  resnum[6]='\0';
	  if(strcmp(resnum,oldresnum) !=0) {
	     nres++;  nchainres++;
	     if(!(N*CA*C*O)) {
		if(mode==1 && n_main_miss<5) printf("    ****missing main chain atoms in residue: %s\n",oldresnum);
		n_main_miss++;
		if(n_main_miss==5 && mode==1) 
		   printf("    ****missing many main chain atoms, not reporting further\n");
	     }
	     N=CA=C=O=ACE=FOR=0;
	     if(strncmp(&buff[17],"ACE",3)==0) { if(mode==1) printf("    ****acetylation \n"); N=CA=C=O=ACE=1;}
	     if(strncmp(&buff[17],"FOR",3)==0) { if(mode==1) printf("    ****formylation \n"); N=CA=C=O=FOR=1;; }
	     } /* End of if(strncmp(res... */
	  strncpy(atname,&buff[13],4);
	  atname[4]='\0';
	  if(!N) N=(!strcmp(atname,"N   "));
	  if(!CA) CA=(!strcmp(atname,"CA  "));
	  if(!C) C=(!strcmp(atname,"C   "));
	  if(!O) O=(!strcmp(atname,"O   "));
	  strncpy(resname,&buff[17],3);
	  resname[3]='\0';

	  sscanf(&buff[27],"%f%f%f%f%f",&x,&y,&z,&Occ,&B);
	  strcpy(rest,&buff[66]);

	  strcpy(oldresnum,resnum);
	} else if(!HET && strcmp(keyword,"HETATM")==0) {
	   if(mode==1) printf("    ****contains heteroatoms (`HETATM')\n");
	   HET=1;
	} /* end of if(strncmp... */

    } /* End of while(c... */
    /* check the last residue */
    if(!(N*CA*C*O)) {
       if(mode==1 && n_main_miss<5) printf("    ****missing main chain atoms in residue: %s\n",oldresnum);
           n_main_miss++;
           if(n_main_miss==5 && mode==1)
              printf("    ****missing many main chain atoms, not reporting further\n");
    }

    if(nchains>0) {
        if(mode==1) printf("        total number of residues: %d, atoms: %d\n",nchainres,nchainatoms);
	nchainres=0; nchainatoms=0;
    } 
    if(nchains==1) total_main_miss=n_main_miss; 
    else total_main_miss+=n_main_miss; 

    if(mode==1) {
      printf("\n  summary: total number of chains: %d\n      number of atoms: %d\n      number of residues: %d\n",nchains,natoms,nres);
    }
    if(total_main_miss>0) {
      if(mode==1) printf("      total of %d residues are missing main chain atoms\n",total_main_miss);
      if((mode==2 || mode==3) && total_main_miss>10) printf("%% warning: main chain atoms are missing for at least 10 residues\n");
    }
     
    /* determine resolution  */
    closefile(pdb,pdbfile); 
    pdb=openfile(pdbfile,"r");
    i=0; found=0;
    while((c=getc(pdb))!=(char)EOF) {
       if(c==RES[i]) i++;
       else i=0;
       if(i==22) {
	  fscanf(pdb,"%f",&resolution); 
	  if(resolution>0 && resolution <1000) {
	    if(mode==1) printf("%s; resolution: %6.3f angstroms\n",pdbfile,resolution); 
	    if(mode==2 || mode==3) printf("%% resolution: %6.3f angstroms\n",resolution);
	    found=1;
	  }
       }
     }
     if(found==0) {
	if(mode==2 || mode==3) printf("%% No resolution found!\n");
	if(NMR) resolution=-1;
	else resolution=-2;	
     }

    /* davis: trying to flush before whatever is causing the next crash */   
    fflush(stdout);
 
    /* get the refinement details */
    closefile(pdb,pdbfile); 
    pdb=openfile(pdbfile,"r");

    reftext[0]='\0';
    found=0;

    while(fgets(buff,99,pdb)!=NULL) {
       /* copy all REMARK  3 lines into a string */
       if(strncmp(REF,buff,10)==0) {
	  buff[72]='\0';
	  sprintf(&reftext[strlen(reftext)],"%s",&buff[10]);
	  reflen+=strlen(&buff[10]);
	  found=1;
       }
    }
    if(found==1) {
       /* remove the double spaces */
       found=0;
       for(i=0; i<strlen(reftext)-1; ++i) {
          if(reftext[i]==' ' && reftext[i+1]==' ') {
   	     sprintf(&reftext[i+1],"%s",&reftext[i+2]);
	     i--;
          }
       }
/*     printf("%s\n",reftext);  */
       /* search through the text for refinement details */
       if(mode==1) printf("Refinement keywords: \n");
       if(mode==2 || mode==3) printf("%% refinement: ");
       for(i=0; i<strlen(reftext); ++i) {
         for(j=0; j<n_ref_type; ++j) {
   	    if(strncmp(&reftext[i],ref_type[j],strlen(ref_type[j]))==0) {
	      if(mode==1 && mode==2 && mode==3) printf("%s ",ref_type[j]);
	      if(mode==1) printf("\n");
	    }
         }
       }
       /* search through the text for the R factor */
       if(mode==1) printf("R value =  ");
       if(mode==2 || mode==3) printf(" R = ");
       for(i=0; i<strlen(reftext); ++i) {
         for(j=0; j<n_r_val; ++j) {
   	   if(strncmp(&reftext[i],r_val[j],strlen(r_val[j]))==0) {
/* SMJS Added if */
	     if (sscanf(&reftext[i+strlen(r_val[j])],"%f",&R_factor))
             {
	       if(mode>=1 && mode<=3) printf("%8.5f",R_factor);
	       found=1;
	       REFINED=1;
	       break;
             }
	   }
         }
       }
    } 
    if(found==0) {
	if((mode==2 || mode==3)) printf("%% R factor not found");
	if(NMR) R_factor=-1;
	else R_factor=-2;
    }
    if(NMR) R_factor=resolution=-1;

    if(mode==2 || mode==3) printf("\n%%\n%%\n");
    if(mode==1) printf("\n");
    if(mode==5) {
	for(i=0; i<4; ++i) printf("%c",code[i]);
	printf(" %8.5f %8.5f %4d %1d %1d %1d %5d %5d %5d\n",resolution,R_factor,year,REFINED,NMR,MODEL,nchains,nres,total_main_miss);
    }
    closefile(pdb,pdbfile);
    free(pdbfile);
    return 0;
}
