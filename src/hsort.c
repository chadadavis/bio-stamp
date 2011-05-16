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
/* modified heap sort from nr routines */

void hsort(n,sortarr)
int n;
/*float ra[];*/
struct path *sortarr;
{
	int l,j,ir,i;
/*	float rra;*/
	struct path rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
/*			rra=ra[--l];*/
			rra=sortarr[--l];
		else {
/*			rra=ra[ir];*/
			rra=sortarr[ir];
/*			ra[ir]=ra[1];*/
			sortarr[ir]=sortarr[1];
			if (--ir == 1) {
/*				ra[1]=rra;*/
				sortarr[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
/*			if (j < ir && ra[j] < ra[j+1]) ++j;*/
			if (j < ir && sortarr[j]score > sortarr[j+1]score) ++j;
/*			if (rra < ra[j]) {*/
			if (rrascore > sortarr[j]score){
/*				ra[i]=ra[j];*/
				sortarr[i]=sortarr[j];
				j += (i=j);
			}
			else j=ir+1;
		}
/*		ra[i]=rra;*/
		sortarr[i]=rra;
	}
}
