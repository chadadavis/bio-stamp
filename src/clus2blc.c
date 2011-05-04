/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 Appendix A of the file LICENSE.  Conditions and restrictions for use of
 this package are also in this file.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Biomolecular Modelling Laboratory   Laboratory of Molecular Biophysics
 Imperial Cancer Research Fund       The Rex Richards Building
 Lincoln's Inn Fields, P.O. Box 123  South Parks Road
 London, WC2A 3PX, U.K.              Oxford, OX1 3PG, U.K.
 Tel: +44 171 269 3583               Tel: +44 865 275368
 FAX: +44 171 269 3417               FAX: 44 865 510454
 e-mail: russell@icrf.icnet.uk       e-mail gjb@bioch.ox.ac.uk
 WWW: http://bonsai.lif.icnet.uk/    WWW: http://geoff.biop.ox.ac.uk/

 The WORK is Copyright (1992,1993,1995,1996) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <gjutil.h>
#include <array.h>
#include <defaults.h>

#define TOKENS " \t\n"


main(int argc,char *argv[])
{
	struct seqdat *seqs;
	FILE *fp,*fout;
	int nseq;
	int found;
	int i,j;
        char *token,*sbit;
        char *line;
        extern FILE *std_err,*std_in,*std_out;
        char *msffile;
        char *blocfile;
        int quiet;
        char c;
	int allen;  /* total alignment length */
        
        std_err = stderr;
        std_in = stdin;
        std_out = stdout;
        
        line = GJstrcreate(MAX_INLEN," ");
        msffile = GJstrcreate(MAX_INLEN,NULL);
        blocfile = GJstrcreate(MAX_INLEN,NULL);

        nseq = 0;
        found = 0;
        quiet = 0;
	allen = 0;

        if(argc > 1){
	  if(strcmp(argv[1],"-q")==0){
            /* Quiet mode - read .MSF file from stdin and output block file to stdout */
            quiet = 1;
            fp = std_in;
            fout = std_out;
	  }
        }else{
          /* Verbose mode - prompt for all filenames */
          fprintf(std_out,"\n\n");
          fprintf(std_out,"CLUSTAL NBRF-PIR format to AMPS Blockfile conversion\n");
          fprintf(std_out,"Copyright: University of Oxford (1992)\n");
          fprintf(std_out,"Author: G. J. Barton (1992)\n\n");
          fprintf(std_out,"Max number/length of alignment - Defined by System\n");
          fprintf(std_out,"If you get a malloc error message - see manual\n\n");
          fprintf(std_out,"Enter CLUSTAL NBRF-PIR alignment filename: ");
          
          fscanf(std_in,"%s",msffile);
          fprintf(std_out,"Opening: %s\n",msffile);
          fp = GJfopen(msffile,"r",1);
          
          fprintf(std_out,"Enter Block filename: ");
          fscanf(std_in,"%s",blocfile);
          fprintf(std_out,"Opening: %s\n",blocfile);
          fout = GJfopen(blocfile,"w",1);
        }
	
	fprintf(fout,"\n");
	fprintf(fout,"Conversion of CLUSTAL NBRF-PIR file to AMPS BLOCKFILE format\n");
	fprintf(fout,"clus2blc:  Geoffrey J. Barton (1992)\n\n");

        seqs = (struct seqdat *) GJmalloc(sizeof(struct seqdat));

       	if(!quiet)fprintf(std_out,"Reading .pir file\n");
       	nseq = 0;
        while(fgets(line,MAX_INLEN,fp) != NULL){
	  if(line[0] == '>'){
	    /* found an identifier */
	    token = strtok(&line[1]," \n");
            if(token != NULL){
              seqs = (struct seqdat *) GJrealloc(seqs,sizeof(struct seqdat) * (nseq + 1));
              seqs[nseq].id = GJstrdup(token);
              if(fgets(line,MAX_INLEN,fp) != NULL){
                /* read the title line */
                seqs[nseq].title = GJstrdup(line);
                seqs[nseq].seq = GJstrcreate(MAX_SEQ_LEN,NULL);
                seqs[nseq].slen = 0;
                seqs[nseq].seq = (char *) GJmalloc(sizeof(char));
                i=0;
                while((c = fgetc(fp)) != '*'){
                    /* read characters until * */
                    if(isalpha(c) || c == '-' || c == '.'){
                        seqs[nseq].seq = (char *) GJrealloc(seqs[nseq].seq,sizeof(char) * (i+1));
                        seqs[nseq].seq[i] = c;
                        ++i;
                    }else if(c == EOF){
                        break;
                    }
		}
	      }
              seqs[nseq].slen = i;
	      if(i > allen) allen = i;
              ++nseq;
	    }
	  }
	}

        if(!quiet)fprintf(std_out,"All %d sequences read in\n",nseq);
        if(!quiet)fprintf(std_out,"Writing .blc file\n");
        
        for(i=0;i<nseq;++i){
            fprintf(fout,">%s %s",seqs[i].id,seqs[i].title);
        }
        fprintf(fout,"* iteration 1\n");
        for(i=0;i<allen;++i){
            for(j=0;j<nseq;++j){
	        if(seqs[j].slen <= i){
		  fprintf(fout,"%c",' ');
		}else{
		  fprintf(fout,"%c",seqs[j].seq[i]);
		}
            }
            fprintf(fout,"\n");
        }
        fprintf(fout,"*\n");
        if(!quiet)fprintf(std_out,"All done\n");
        
        for(i=0;i<nseq;++i){
	  GJfree(seqs[i].seq);
	  GJfree(seqs[i].id);
  	  GJfree(seqs[i].title);
	}
	GJfree(seqs);
	GJfree(line);
	GJfree(blocfile);
	GJfree(msffile);

}	
