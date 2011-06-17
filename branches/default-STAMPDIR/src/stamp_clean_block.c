#include <stdio.h>

#include <stamp.h>
#include <stamprel.h>

/* as for clean_block, only this cleans up stamp information as well */

int stamp_clean_block(struct seqdat *bloc, int nbloc, int window, float min_frac,
	struct stampdat *stamp, int nstamp) {

	int i,j,k,l;
	int N,C;
	int neighbors;
	int firstgap;
	int halfgap;
	int bloclen;
	int nstart,nend;
	int count;
	int ngap;
	int nungapped;

        int min_n;
	int nskip;
	int *keep,*start,*end;
	int *max,*min;
	int *skip;
	char *temp;

	bloclen=strlen(&bloc[1].seq[1]);
	keep=(int*)malloc(bloclen*sizeof(int));
	start=(int*)malloc(bloclen*sizeof(int));
	end=(int*)malloc(bloclen*sizeof(int));
	min=(int*)malloc(bloclen*sizeof(int));
	max=(int*)malloc(bloclen*sizeof(int));
	temp=(char*)malloc(bloclen*sizeof(char));
	skip=(int*)malloc(nbloc*sizeof(int));

	/* if spaces are encountered, we must ignore them */
	nskip=0;
	for(i=0; i<nbloc; ++i) {
	   rmsp(bloc[i+1].id);
	   skip[i]=(strncmp(bloc[i+1].id,"space",5)==0);
	   nskip+=(skip[i]);
	}
        min_n = (int)(min_frac * (nbloc-nskip)+0.5);
        printf("Min N of sequence to be ungapped is %d\n",min_n);

	/* first thing we need to to is to find the bits to keep intact 
	 *  first we simply assign keep to '1' if the position is
	 *  devoid of gaps */
	for(i=0; i<bloclen; ++i) {
	  keep[i]=1;
          nungapped=0;
	  for(j=0; j<nbloc; ++j)  {
	     if(!skip[j]) {
                if(bloc[j+1].seq[i+1]!=' ') { nungapped++; }
             }
          }
          if(nungapped<min_n) {
		keep[i]=0;
          }
	}

	/* now we must smooth the array keep out */
        printf("before smoothing\n");
        for(i=0; i<bloclen; ++i) printf("%1d",keep[i]); printf("\n");   

        /* now smooth the array out according to the window */
        for(i=0; i<bloclen; ++i) {
         if(keep[i]==1) {
           neighbors=0;
           for(j=1; j<=window; ++j) {
              if(keep[i+j]==0 || i+j>(bloclen-1)) break;
              else neighbors++;
           }
           for(j=1; j<=window; ++j) {
              if(keep[i-j]==0 || (i-j)<0) break;
              else neighbors++;
           }
           keep[i]=(neighbors>=(window-1)); 
         }
        }
        printf("after smoothing\n");  
     	for(i=0; i<bloclen; ++i) printf("%1d",keep[i]); printf("\n");   

	/* ignore all STAMP values outside the smoothed region */
	for(i=0; i<bloclen; ++i) if(!keep[i]) stamp[0].n[i]=-1.0;

	/* now lets find the maximum and minimum length of gaps */

	/* now we must clean up the gaps.  This will be done by defining
	 *  several start and end points in the alignment and working out
	 *  from these */
	nstart=1; nend=1;
	start[0]=end[0]=-1;
	if(keep[0]) start[nstart++]=0;
	for(i=0; i<bloclen; ++i) {
	   if(i>0 && keep[i] && !keep[i-1]) start[nstart++]=i;
	   if(i<(bloclen-1) && keep[i] && !keep[i+1]) end[nend++]=i;
	}
	start[nstart++]=bloclen;
	start[nstart+1]=bloclen;
	end[nend++]=bloclen;
/*	printf("nstart: %d, nend: %d\n",nstart,nend); 
	for(i=0; i<nstart; ++i) 
           printf("i=%d, start[i]=%d, end[i]=%d\n",i,start[i],end[i]);  */
	/* now find the maximum and minimum length of each gap */
	for(i=1; i<nstart; ++i) {
	   max[i]=0; min[i]=(start[i]-end[i-1]-1);
	   for(j=0; j<nbloc; ++j) if(!skip[j]) {
	      l=0;
	      for(k=end[i-1]+1; k<start[i]; ++k) 
		 l+=(bloc[j+1].seq[k+1]!=' ');
	      if(l<min[i]) min[i]=l;
	      if(l>max[i]) max[i]=l;
	   }
/*	   printf("i=%d, end[i-1]=%d, start[i]=%d, min[i]=%d, max[i]=%d\n",i,end[i-1],start[i],min[i],max[i]); */
	}

	/* now go through each gap for each sequence  and compress it to fit into the maximum in
	 *  a sensible way */
	for(i=0; i<nbloc; ++i) if(!skip[i]) {
	   for(j=1; j<nstart; ++j) {
	      l=0;
	      for(k=end[j-1]+1; k<start[j]; ++k) {
		 if(bloc[i+1].seq[k+1]!=' ') temp[l++]=bloc[i+1].seq[k+1];
		 bloc[i+1].seq[k+1]=' ';
	      }
	      temp[l]='\0';
/*	      printf("sequence %d, between %d and %d, sequence is %s; ",i,end[j-1]+1,start[j],temp); */
	      /* now put spaces in the middle of the sequence according to the maximum length
	       *  of gap allowed */
	      count=end[j-1]+2;
/*	      printf("spreading to : "); */
	      if(strlen(temp)==0) {
		 N=C=0;
	      } else if(strlen(temp)%2==0) {
		 N=C=(strlen(temp)/2);
	      } else {
		 N=(strlen(temp)/2);
		 C=((strlen(temp)/2)+1);
    	      }
	      if(j==1) {
		 N=0; 
		 C=strlen(temp);
	      }
	      if(j==nstart-1) {
		 N=strlen(temp);
		 C=0;
	      }
/*	      printf("temp is %d characters long, N=%d,C=%d\n",strlen(temp),N,C); */
	      for(k=0; k<N; ++k) {
		 bloc[i+1].seq[count]=temp[k];
/*		 printf("%c",bloc[i+1].seq[count]); */
		 count++;
	      }
	      for(k=0; k<max[j]-N-C; ++k) {
		 bloc[i+1].seq[count++]=' ';
/*		 printf(" "); */
	      }
	      for(k=0; k<C; ++k) {
		 bloc[i+1].seq[count]=temp[strlen(temp)-C+k];
/*		 printf("%c",bloc[i+1].seq[count]); */
		 count++;
	      }
/*	      printf("\n"); */
	      /* that is it */
	   }
	}
	/* now we must remove blank spaces */
	for(i=0; i<bloclen; ++i) {
	  ngap=0;
	  for(j=0; j<nbloc; ++j) if(!skip[j]) ngap+=(bloc[j+1].seq[i+1]==' ');
	  if(ngap==(nbloc-nskip)) {
	     bloclen--;
	     for(k=i; k<bloclen; ++k) {
	         for(j=0; j<nbloc; ++j) if(!skip[j]) 
		   bloc[j+1].seq[k+1]=bloc[j+1].seq[k+2]; 
	         for(j=0; j<nstamp; ++j) 
		     stamp[j].n[k]=stamp[j].n[k+1];
	     }
	     /* move next position into the blank */
	     i--; /* repeat this position */
	  }
	} 
	for(j=0; j<nbloc; ++j) bloc[j+1].seq[bloclen+1]='\0';
	/* output the results */
/*	printf("after cleaning \n"); 
	for(i=0; i<bloclen; ++i) {
	  printf("%3d: ",i);
	  for(j=0; j<nbloc; ++j) printf("%c",bloc[j+1].seq[i+1]);
	  printf("  %1d\n",keep[i]);
	}
*/
	free(keep); free(start); free(end); free(min); free(max); free(temp); free(skip);
	return 0;
}

