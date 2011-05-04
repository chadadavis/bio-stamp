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

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
/***************************************************************************

msf2blc:  A program to convert a GCG .MSF file into an AMPS blockfile.


   Copyright:  University of Oxford (1992)

TERMS OF USE:

The computer software and associated documentation called ASSP hereinafter
referred to as the WORK is more particularly identified and described in 
Appendix A.

The WORK was written and developed by: Geoffrey J. Barton

Laboratory of Molecular Biophysics
University of Oxford
Rex Richards Building
South Parks Road
Oxford OX1 3QU
U.K.

Tel:  (+44) 865-275368
Fax:  (+44) 865-510454

Internet: gjb@bioch.ox.ac.uk
Janet:    gjb@uk.ac.ox.bioch

The WORK is Copyright (1992) University of Oxford

Administrative Offices
Wellington Square
Oxford OX1 2JD
U.K.

CONDITIONS:

The WORK is made available for educational and non-commercial research 
purposes.

For commercial use, a commercial licence is required - contact the author
at the above address for details.

The WORK may be copied and redistributed, provided that no charge is
made for the WORK, that the WORK is not incorporated into a program or
package for sale and that this text, all copyright notices and the
authors' address is left unchanged on every copy.

The WORK may be modified, however this text, all copyright notices and
the authors' address must be left unchanged on every copy, and all
changes must be documented after this notice.  A copy of the
modified WORK must be supplied to the author.

All use of the WORK must cite:  
R.B. Russell and G.J. Barton, "The limits of secondary structure prediction
accuracy from multiple sequence alignment", J. Mol. Biol., 234, 951-957, 1993.

See APPENDIX A in the file LICENCE

****************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "gjutil.h"
#include "array.h"
#include "defaults.h"

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
        
        std_err = stderr;
        std_in = stdin;
        std_out = stdout;
        
        line = GJstrcreate(MAX_INLEN," ");
        msffile = GJstrcreate(MAX_INLEN,NULL);
        blocfile = GJstrcreate(MAX_INLEN,NULL);

        nseq = 0;
        found = 0;
        quiet = 0;

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
          fprintf(std_out,"GCG .MSF to AMPS Blockfile conversion\n");
          fprintf(std_out,"Copyright: University of Oxford (1992)\n");
          fprintf(std_out,"Author: G. J. Barton (1992)\n\n");
          fprintf(std_out,"Max number/length of alignment - Defined by System\n");
          fprintf(std_out,"If you get a malloc error message - see manual\n\n");
          fprintf(std_out,"Enter MSF filename: ");
          
          fscanf(std_in,"%s",msffile);
          fprintf(std_out,"Opening: %s\n",msffile);
          fp = GJfopen(msffile,"r",1);
          
          fprintf(std_out,"Enter Block filename: ");
          fscanf(std_in,"%s",blocfile);
          fprintf(std_out,"Opening: %s\n",blocfile);
          fout = GJfopen(blocfile,"w",1);
        }
	
	fprintf(fout,"\n");
	fprintf(fout,"Conversion of GCG .MSF file to AMPS BLOCKFILE format\n");
	fprintf(fout,"msf2blc:  Geoffrey J. Barton (1992)\n\n");

        seqs = (struct seqdat *) GJmalloc(sizeof(struct seqdat));

       	if(!quiet)fprintf(std_out,"Reading .blc file\n");
        while(fgets(line,MAX_INLEN,fp) != NULL){
	  if(line[0] != '\n'){
             token = strtok(line,TOKENS);
             if(token != NULL){
               if(strcmp(token,"Name:") == 0){
                 /* This is a seq id name */
                  token = strtok(NULL,TOKENS);
                  seqs = (struct seqdat *) GJrealloc(seqs,sizeof(struct seqdat) * (nseq +1));
                  seqs[nseq].id = GJstrdup(token);
                  seqs[nseq].title = GJstrdup(line);
                  seqs[nseq].slen = 0;
                  seqs[nseq].seq = (char *) GJmalloc(sizeof(char));
                  ++nseq;
                  if(!quiet)fprintf(std_out,"%s\n",seqs[nseq-1].id);
	       }else if((strcmp(token,"//") == 0) || found){
                  /* this signals the end of identifiers so process sequences*/
                  found = 1;
                  if(token != NULL){
                    /* find out which seq this is */
                    i=0;
		    for(i=0;i<nseq;++i){
	               if(strcmp(token,seqs[i].id) == 0){
		         break;
		       }
		     }
                     /* read in the sequence */
                     if(i < nseq){
                       token = strtok(NULL,"\n");
                       if(token == NULL){
                         GJerror("Cannot find sequence in line");
                         fprintf(std_err,"%s",line);
                         exit(1);
		       }
                       j=0;
                       while(token[j] != '\0'){
                         if(isalpha(token[j]) || token[j] == '.'){
                           seqs[i].seq = (char *) GJrealloc(seqs[i].seq,sizeof(char) * (seqs[i].slen +1));
                           seqs[i].seq[seqs[i].slen] = token[j];
                           ++seqs[i].slen;
			 }
                         ++j;
		       }
		     }
		  }
		}else{
                  /* this is a comment line - just echo */
                  fprintf(fout,"%s\n",line);
		}
	     }
	   }
	}
        if(!quiet)fprintf(std_out,"All %d sequences read in\n",nseq);
        if(!quiet)fprintf(std_out,"Writing .blc file\n");
        
        for(i=0;i<nseq;++i){
            fprintf(fout,">%s %s\n",seqs[i].id,seqs[i].title);
        }
        fprintf(fout,"* iteration 1\n");
        for(i=0;i<seqs[0].slen;++i){
            for(j=0;j<nseq;++j){
                fprintf(fout,"%c",seqs[j].seq[i]);
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

