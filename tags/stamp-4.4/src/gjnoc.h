/*
Copyright (1997,1998,1999,2010) Robert B Russell & Geoffrey J Barton

This file is part of STAMP

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE  See the
GNU General Public License for more details A copy of the license
can be found in the LICENSE file in the STAMP installation directory

STAMP was developed by Robert B Russell and Geoffrey J Barton of
current addresses:

 Prof Robert B Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
 Cell Networks, University of Heidelberg            College of Life Sciences
 Room 564, Bioquant                                 University of Dundee
 Im Neuenheimer Feld 267                            Dow Street
 69120 Heidelberg                                   Dundee DD1 5EH
 Germany                                            UK
                                                
 Tel: +49 6221 54 513 62                            Tel: +44 1382 385860
 Fax: +49 6221 54 514 86                            FAX: +44 1382 385764
 Email: robert.russell@bioquant.uni-heidelberg.de   E-mail g.j.barton@dundee.ac.uk
 WWW: http://www.russell.embl-heidelberg.de         WWW: http://www.compbio.dundee.ac.uk

 All use of STAMP must cite: 

 RB Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992)
*/
#ifndef GJ_NOC_H
#define GJ_NOC_H
/* Include file of definitions for program oc (oc) */

#define val_A_B(arr,n,i,j) arr[i][j-i-1]

/* structure to hold cluster details
Note:  	Normally, parentA and parentB point to pre-existing
	clusters stored above them in the array of clusters
	However, 
       	If the cluster is a single entity, then parentA and parentB
       	are set to NULL
       	If one or more of the parents are single entites, then
       	the storage for the entites is malloc'd independently
       	of the master cluster array
*/

struct sc{
	double score;	/* score for this cluster */
	int n;		/* number of entitites in this cluster */
	int *members;	/* entities - this is redundant but speeds code */
	struct sc *parentA; /* one parent of this cluster */
	struct sc *parentB; /* the other parent */
	char lab;           /* label set this to 1 normally  0 if the cluster has been processed */
};

double **read_up_diag(FILE *infile, int n);
char   **read_idents(FILE *infile, int n);
void write_up_diag(FILE *infile, double **ret_val,int n);
void show_entity(struct sc *entity,FILE *out);
void show_inx_entity(struct sc *entity,int *inx,FILE *out);
struct sc *make_entity(int i,double base_val);
/* double val_A_B(double **arr,int n,int i,int j);*/
double val_Grp_B(
	double **arr,	/* upper diagonal array */
	int n,		/* side of array */
	int *grpA,	/* array containing entities of group */
	int ng,		/* number of members in grpA */
	int j,		/* entity to compare to grpA */
	int choice	/* see above */
);
double val_Grp_Grp(
	double **arr,	/* upper diagonal array */
	int n,		/* side of array */
	int *grpA,	/* array containing entities of group */
	int ngA,	/* number of members in grpA */
	int *grpB,	/* array containing entities of group B*/
	int ngB,	/* number of members in grpB */
	int choice	/* see above */
);
void remove_unclust(int *unclust,int *nunclust,int val);
void iprintarr(int *arr,int n,FILE *out);
int no_share(
	    int *grpA,
	    int ngA,
	    int *grpB,
	    int ngb
);
void add_notparent(int *np,int *n,int A);
void remove_notparent(int *np,int *n,int A);
void sub_notparent(int *np,int *n,int A,int B);    	
double **GJDudarr(int n);
void draw_dendrogram(FILE *fout,
                     struct sc *clust,
		     char **idents,
                     int *inx,
                     int n,
		     int sim);

double get_mean(int *arr,int *inx,int n);
void draw_arch(FILE *fout,
               float x1,
               float x2,
               float x3,
               float y1,
               float y2);

double trans_y(double X,
	       double Xmin,
	       double Xmax,
	       double Xrange,
	       float Twidth,
	       float Xoffset);

double trans_x(double x,
	       int sim,
	       double Xmin,
	       double Xmax,
	       double Xrange,
	       float Twidth,
	       float Xoffset);
void PSPutText(float x,float y,char *text,FILE *outf);
void show_id_entity(struct sc *entity,char **id,FILE *out);

void mark_parents(struct sc *entity);
void write_unclustered(struct sc *entity,char **idents,FILE *fout);
void up_diag_log(double **arr,int n);
void print_amps_cluster(struct sc *entity,int *inx,FILE *fp);

/* Rob's cluster structures for returning values to STAMP */
#if !defined(CLUST_STRUCT)
struct indclust {
        int number;
        int *member;
        };
struct cluster {
        struct indclust a;
        struct indclust b;
        };
#define CLUST_STRUCT
#endif

struct gjnoc {
        struct cluster *clusters;
        int *order;
};

struct gjnoc *GJnoc(
		    double **arr,           /* upper diagonal array */
		    char **idents,          /* identifiers for entities */
		    int n,                  /* number of entities in array */
		    char *parms             /* string of parameters separated by spaces */
		    );


#endif /* GJ_NOC_H */






