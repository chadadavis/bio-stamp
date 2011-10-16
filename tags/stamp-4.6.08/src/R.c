#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rbr_aadat.h>
#include <rbr_vector.h>

/* returns the power distance between a pair of integer coordinates */
float RBR_ipowdist(int *x1, int *x2, int PRECISION, int power) {
	int i,j;

	float A,B,D;

	D=0.0;
	for(i=0; i<3; ++i) {
	   A=fabs((float)(x1[i]-x2[i]))/(float)PRECISION; 
	   B=1.0;
	   for(j=0; j<power; ++j) B*=A;
	   D+=B;
	}

	return D;
}
/* calculate center of mass for a set of integer protein CA coordinates */
int *RBR_c_of_m_pep(int **coords, int ncoords, char *seq, int w, float *total_mass) {
	int i,j;
	int *R;

	R=(int*)malloc(3*sizeof(int));
	for(i=0; i<3; ++i) R[i]=0;
	
	if(w) {  /* weighting for molar weight has been requested */
	  (*total_mass)=0.0;
	  for(i=0; i<ncoords; ++i) (*total_mass)+=(float)RBR_AA_MASS[(int)(seq[i]-'A')];
	} else (*total_mass)=(float)ncoords;

	for(i=0; i<ncoords; ++i) {
	  for(j=0; j<3; ++j) {
	    if(w) R[j]+=(int)( (RBR_AA_MASS[seq[i]-'A']/(*total_mass))*(float)coords[i][j]);
	    else  R[j]+=coords[i][j]/(*total_mass);
	  }
	}

	return R;
}

/* calculate the radius of gyration from a set of integer CA coordinates and the Ro */
float RBR_r_of_gyration_pep(int *Ro, int **coords, int ncoords, char *seq, int w, int PRECISION, int DIST_POWER) {
	int i,j;

	float R,total_mass,distsq;
	 
	R=0.0;

	if(w) {  /* weighting for molar weight has been requested */
	  total_mass=0.0; 
	  for(i=0; i<ncoords; ++i) total_mass+=RBR_AA_MASS[seq[i]-'A'];
	} else total_mass=(float)ncoords;

/*	printf("total mass: %f\n",total_mass);
	printf("centre of mass: %8.5f %8.5f %8.5f\n", 
	(float)Ro[0]/(float)PRECISION,(float)Ro[1]/(float)PRECISION,(float)Ro[2]/(float)PRECISION);
	printf("distances\n");
*/	
	for(i=0; i<ncoords; ++i) {
	   distsq=RBR_ipowdist(coords[i],Ro,PRECISION,DIST_POWER);
/*	   printf("%4d: %f\n",i+1,sqrt(distsq)); */
	   if(w) R+=(RBR_AA_MASS[seq[i]-'A']/total_mass)*(float)distsq;
	   else  R+=distsq/total_mass;
	}

	R=(float)sqrt((double)R);
	return R;
}
float RBR_r_max_pep(int *Ro, int **coords, int ncoords, char *seq, int PRECISION, int DIST_POWER) {
	int i,j;
	float R,distsq,dist;

/*	printf("centre of mass: %8.5f %8.5f %8.5f\n", 
	(float)Ro[0]/(float)PRECISION,(float)Ro[1]/(float)PRECISION,(float)Ro[2]/(float)PRECISION);
	printf("distances\n");
 */

	for(i=0; i<ncoords; ++i) {
	   distsq=RBR_ipowdist(coords[i],Ro,PRECISION,DIST_POWER);
           dist = sqrt(distsq);
/*	   printf("d for : %8.5f %8.5f %8.5f is %8.5f\n",(float)coords[i][0]/(float)PRECISION,(float)coords[i][1]/(float)PRECISION,(float)coords[i][2]/(float)PRECISION,sqrt(distsq)); */
           if(dist>R) { 
              R = dist;
           }
	}
       /* printf("Max is %8.5f\n",R); */
	return R;
}

int RBR_dist_stats(int **coords, int ncoords, float *mean, float *sd, float *longest, float *shortest, 
    int PRECISION, int DIST_POWER) {
	int i,j,k;
	int total;
	float sum,sumpow;
	float powdist,dist;

	sum=0.0; sumpow=0.0; (*longest)=0.0; (*shortest)=100.0; total=0.0;

	for(j=0; j<ncoords; ++j) {
	     for(k=j+1; k<ncoords; ++k) {
		total++;
		powdist=RBR_ipowdist(coords[j],coords[k],PRECISION,DIST_POWER);
		if(DIST_POWER==2) dist=sqrt((double)powdist);
		else dist=pow((double)powdist,(1/(float)DIST_POWER));
		sum+=dist;
		sumpow+=powdist;
		if(dist>(*longest)) (*longest)=dist;
		if(dist<(*shortest)) (*shortest)=dist;
	     }
	}
	(*mean)=(sum/(float)total);
	(*sd)=sqrt((double)((sumpow-(sum*sum)/(float)total)/(float)total));

	return 0;
}


/* calculate center of mass for a set of integer coordinates
 *  with a given array of integer point masses */

int *RBR_c_of_m(int **coords, int ncoords, int *masses, int w) {
	int i,j;
	int *R;

	float total_mass;

	R=(int*)malloc(3*sizeof(int));
	for(i=0; i<3; ++i) R[i]=0;
	
	if(w) {  /* weighting for molar weight has been requested */
	  total_mass=0.0;
	  for(i=0; i<ncoords; ++i) total_mass+=(float)masses[i];
	} else total_mass=(float)ncoords;

	for(i=0; i<ncoords; ++i) {
	  for(j=0; j<3; ++j) {
	    if(w) R[j]+=(int)(((float)masses[i]/total_mass)*(float)coords[i][j]);
	    else  R[j]+=coords[i][j]/total_mass;
	  }
	}

	return R;
}

/* calculate the radius of gyration from a set of integer coordinates,
 * the Ro, and an array of point masses */
float RBR_r_of_gyration(int *Ro, int **coords, int ncoords, int *masses, int w,   int PRECISION, int DIST_POWER) {
	int i,j;

	float R,total_mass,distsq;
	 
	R=0.0;

	if(w) {  /* weighting for molar weight has been requested */
	  total_mass=0.0; 
	  for(i=0; i<ncoords; ++i) total_mass+=(float)masses[i];
	} else total_mass=(float)ncoords;

/*	printf("total mass: %f\n",total_mass);
	printf("centre of mass: %8.5f %8.5f %8.5f\n", 
	(float)Ro[0]/(float)PRECISION,(float)Ro[1]/(float)PRECISION,(float)Ro[2]/(float)PRECISION);
	printf("distances\n");
*/	
	for(i=0; i<ncoords; ++i) {
	   distsq=RBR_ipowdist(coords[i],Ro,PRECISION,DIST_POWER);
/*	   printf("%4d: %f\n",i+1,sqrt(distsq)); */
	   if(w) R+=((float)masses[i]/total_mass)*(float)distsq;
	   else  R+=distsq/total_mass;
	}

	R=(float)sqrt((double)R);
	return R;
}
