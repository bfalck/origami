//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "params.h"

#define BIGFLT 1e30 /* Biggest possible floating-point number */

typedef struct Particle {
  float vol;
  int nadj;
  char tag;
  int *adj;
} PARTICLE;


int main(int argc,char **argv) {

  FILE *adj, *vol, *tag, *hn;
  PARTICLE *p;
  char *paramfile;
  char adjfile[80], volfile[80], tagfile[80], halolistfile[80];
  int i, j,k,l, h, h2,hl,n,np, np2, nhaloes, nhl, nhlcount, nhl2;
  int *halonum, *halolist, *halolist2, *lonecore;
  int nchanges,adjj,nin;
  int change,minhalonum,numcut,ncp;
  int **hm,*numinh,*haloid,nh,*nph,haloidtoadd;

  if (argc < 2) {
    printf("Please specify a parameter file\n");
    exit(0);
  } else paramfile = argv[1];
  read_parameters(paramfile);

  sprintf(adjfile,"%s%sadj.dat",outdir,taglabel);
  sprintf(volfile,"%s%svol.dat",outdir,taglabel);
  sprintf(tagfile,"%s%stag.dat",outdir,taglabel);
  sprintf(halolistfile,"%s%shns.dat",outdir,halolabel);

  adj = fopen(adjfile, "r");
  if (adj == NULL) {
    printf("Unable to open %s\n",adjfile);
    exit(0);
  }
  fread(&np,1, sizeof(int),adj);
  
  p = (PARTICLE *)malloc(np*sizeof(PARTICLE));
  /* Adjacencies*/
  for (i=0;i<np;i++) {
    fread(&p[i].nadj,1,sizeof(int),adj); 
    /* The number of adjacencies per particle */
    if (p[i].nadj > 0)
      p[i].adj = (int *)malloc(p[i].nadj*sizeof(int));
    p[i].nadj = 0; /* Temporarily, it's an adj counter */
  }
  for (i=0;i<np;i++) {
    fread(&nin,1,sizeof(int),adj);
    if (nin > 0)
      for (k=0;k<nin;k++) {
	fread(&j,1,sizeof(int),adj);
	
	/* Set both halves of the pair */
	p[i].adj[p[i].nadj] = j;
	p[j].adj[p[j].nadj] = i;
	p[i].nadj++; p[j].nadj++;
      }
  }
  fclose(adj);

  /* Check that we got all the pairs */
  adj = fopen(adjfile, "r");
  fread(&np,1, sizeof(int),adj);
  for (i=0;i<np;i++) {
    fread(&nin,1,sizeof(int),adj); /* actually nadj */
    if (nin != p[i].nadj) {
      printf("We didn't get all of %d's adj's; %d != %d.\n",i,nin,p[i].nadj);
      /*exit(0);*/
    }
  }
  fclose(adj);

  /* Volumes */
  vol = fopen(volfile, "r");
  if (vol == NULL) {
    printf("Unable to open volume file %s.\n\n",volfile);
    exit(0);
  }
  fread(&np2,1, sizeof(int),adj);
  if (np != np2) {
    printf("Number of particles doesn't match! %d != %d\n",np,np2);
    exit(0);
  }
  for (i=0;i<np;i++) {
    fread(&p[i].vol,1,sizeof(float),vol);
    if ((p[i].vol < 1.1e-30) || (p[i].vol > 9.9e29)) {
      printf("Whacked-out volume found, of particle %d: %f\n",i,p[i].vol);
      p[i].vol = 1.;
    }
  }
  fclose(vol);

  tag = fopen(tagfile,"r");
  if (tag == NULL) {
    printf("Unable to open tag file %s.\n\n",tagfile);
    exit(0);
  }
  fread(&np,1, sizeof(int),tag);
  for (i=0;i<np;i++) {
    fread(&p[i].tag,1,sizeof(char),tag);
  }
  fclose(tag);

  halonum = (int *)malloc(np*sizeof(int));
  numinh = (int *)malloc(np*sizeof(int));
  lonecore = (int *)malloc(np*sizeof(int));
  numcut = 0;
  for (i=0;i<np;i++){
    numinh[i] = 0;
    lonecore[i] = 0;
    if (p[i].tag == 3) {
      halonum[i] = i;
      numcut ++;
    } else{
      halonum[i] = -1;
    }
  }

  printf("\n%d particles make the cut.\n",numcut);
  printf("\nFinding groups ...\n");
  /* First find halo cores */
  nchanges=1;
  while (nchanges > 0) {
    nchanges = 0;
    for (i = 0; i < np; i++) {
      if ((p[i].tag == 3)&&(p[i].vol<volcut)) {
	change=0;
	/* find minimum halonum among particle & adjacencies */
	minhalonum = halonum[i];
	for (j=0; j< p[i].nadj; j++) {
	  adjj = p[i].adj[j];
	  if ((p[adjj].tag == 3)&&(p[adjj].vol<volcut)) {
	    /*if (((p[adjj].tag == 3)&&(p[adjj].vol<volcut)) || (p[adjj].vol < p[i].vol)) {*/
	    if (halonum[adjj] < minhalonum) {
	      change=1;
	      minhalonum = halonum[adjj];
	    }
	  }
	}
	if (change == 1){
	  /* assign that minhalonum to particle & adjacencies */
	  halonum[i] = minhalonum;
	  for (j=0; j< p[i].nadj; j++) {
	    adjj = p[i].adj[j];
	    if ((p[adjj].tag == 3)&&(p[adjj].vol<volcut)) {
	      halonum[adjj] = minhalonum;
	    }
	  }
	}
	nchanges += change;
      }
    }
    printf("nchanges=%d\n",nchanges);
  
  }
  printf ("Halo cores found. Removing spurious cores.\n\n");

  /* Group halos */
  for (i=0;i<np;i++){
    if (halonum[i] >= 0) { // true for ungrouped halo particles + cores
      numinh[halonum[i]] ++; // will be > 1 only for cores
    }
  }

  /* remove cores with < NCP particles */
  nh = 0;
  ncp = 3;
  for (i=0;i<np;i++){
    if (numinh[halonum[i]] >= ncp) {
      if (halonum[i] == i) {
	nh ++;
      }
    } else if (p[i].tag == 3) {
      numinh[halonum[i]] = 0;
      halonum[i] = i;
      if (p[i].vol < volcut) lonecore[i] = 1;
    }
  }
  printf("%d halo cores\n",nh);
  for (i=0;i<np;i++) {
    numinh[i] = 0;
  }


  printf ("Associating lower-density outskirts.\n");

  nchanges = 0;
  for (i = 0; i < np; i++) {
    if ((p[i].tag == 3)&&(p[i].vol < volcut)&&(lonecore[i]==0)) { /* part of a halo core */
      
      /* associate adjacencies with that */
      for (j=0; j< p[i].nadj; j++) {
	adjj = p[i].adj[j];
	if (p[adjj].tag == 3) {
	  halonum[adjj] = halonum[i];
	  nchanges ++;
	}
      }
    }
  }
  printf("nchanges=%d\n",nchanges);
  nchanges=1;

  /* Now do friends of halo core particles */
  while (nchanges > 0) {
    nchanges = 0;
    for (i = 0; i < np; i++) {
      if ((p[i].tag == 3) && ((p[i].vol >= volcut)||(lonecore[i]==1))) {
	if ((p[halonum[i]].vol < volcut)&&!(halonum[i]==i)) { /* Attached to a halo core */
	  for (j=0; j< p[i].nadj; j++) {
	    adjj = p[i].adj[j];
	    if (((p[adjj].vol >= volcut)||(lonecore[adjj]==1)) && (halonum[adjj] == adjj) && (p[adjj].tag == 3)) { /* if p[i].tag == 3 and not previously assigned */
	      halonum[adjj] = halonum[i];
	      nchanges ++;
	    }
	  }
	}
      }
    } 
    printf("nchanges=%d\n",nchanges);
  }

  /* Now group uncored haloes */
  printf ("'Uncored' haloes now being grouped.\n");
  nchanges=1;
  while (nchanges > 0) {
    nchanges = 0;
    for (i = 0; i < np; i++) {
      if ((p[i].tag == 3)&&((p[i].vol >= volcut)||(lonecore[i]==1))) {
	if (p[halonum[i]].vol >= volcut) {
	  change=0;
	  minhalonum = halonum[i];
	  for (j=0; j< p[i].nadj; j++) {
	    adjj = p[i].adj[j];
	    if ((p[adjj].tag == 3)) {
	      if (halonum[adjj] < minhalonum) {
		if ((p[halonum[adjj]].vol >= volcut)||(lonecore[halonum[adjj]]==1)) {
		  change=1;
		  minhalonum = halonum[adjj];
		}
	      }
	    }
	  }
	  if (change == 1){
	    halonum[i] = minhalonum;
	    for (j=0; j< p[i].nadj; j++) {
	      adjj = p[i].adj[j];
	      if ((p[adjj].tag == 3)) {
		if ((p[halonum[adjj]].vol >= volcut)||(lonecore[halonum[adjj]]==1)) {
		  halonum[adjj] = minhalonum;
		}
	      }
	    }
	  }
	}
	nchanges += change;
      }
    }
    printf("nchanges=%d\n",nchanges);
  }

  /* Group haloes */
  nh=0;
  for (i=0;i<np;i++){
    if (halonum[i] >= 0) {
      numinh[halonum[i]] ++;
    }
    if (halonum[i] == i) {
      nh ++;
    }
  }
  printf("nh = %d\n",nh);fflush(stdout);

  /* remove halos with < npmin particles */
  nh = 0;
  for (i=0;i<np;i++){
    if (numinh[halonum[i]] >= npmin) {
      if (halonum[i] == i) {
	nh ++;
      }
    } else {
      numinh[halonum[i]] = 0;
      halonum[i] = -1;
    }
  }
  printf("nh >= %d = %d\n",npmin,nh);fflush(stdout);

  haloid = (int *)malloc(np*sizeof(int *));
  /* halo membership array of lists */
  hm = (int **)malloc(nh*sizeof(int **));
  nph = (int *)malloc(nh*sizeof(int *));
  h = 0;
  for (i=0;i<np;i++){
    if (numinh[i] > 0) {
      haloid[i] = h;
      hm[h] = (int *)malloc(numinh[i]*sizeof(int));
      h++;
    }
  }
  for (h=0;h<nh;h++){
    nph[h] = 0;
  }
  /* assign hm[h] arrays */
  for (i=0; i<np; i++){
    if (halonum[i] >= 0) {
      haloidtoadd = haloid[halonum[i]]; 
      hm[haloidtoadd][nph[haloidtoadd]] = i;
      nph[haloidtoadd]++;
    }
  }      
  
  /* Output */
  
  hn = fopen(halolistfile,"w");
  if (hn == NULL) {
    printf("Unable to open %s\n",halolistfile);
    exit(0);
  }
  fwrite(&nh,1, sizeof(int),hn);
  for (h=0;h<nh;h++){
    fwrite(&(nph[h]),1, sizeof(int),hn);
    fwrite(hm[h],nph[h], sizeof(int),hn);
  }
  fclose(hn);
  
} // main
