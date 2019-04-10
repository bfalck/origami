//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "voz.h"
#include "params.h"

int main(int argc, char *argv[]) {

  FILE *part, *adj, *vol;
  char *paramfile, partfile[80], adjfile[80], volfile[80];
  float *vols, volstemp;
  
  PARTADJ *adjs;

  int np,np2,na;

  int i,j,k,p,nout;
  int nvp,npnotdone,nvpmax, nvpsum, *orig;
  double avgnadj, avgvol;
  
  if (argc < 2) {
    printf("Please specify a parameter file\n");
    exit(0);
  } else paramfile = argv[1];
  read_parameters(paramfile);

  if (numdiv < 2) {
    printf("Cannot have a number of divisions less than 2.  Resetting to 2:\n");
    numdiv = 2;
  }
  
  np = -1; nvpmax = -1; nvpsum = 0;

  for (i = 0; i < numdiv; i++) {
   for (j = 0; j < numdiv; j++) {
    for (k = 0; k < numdiv; k++) {
      sprintf(partfile,"part.%s.%02d.%02d.%02d",taglabel,i,j,k);
      part = fopen(partfile,"r");
      if (part == NULL) {
	printf("Unable to open file %s.\n\n",partfile);
	exit(0);
      }
      fread(&np2,1,sizeof(int),part);
      fread(&nvp,1,sizeof(int),part);
      if (np == -1)
	np = np2;
      else 
	if (np2 != np) {
	  printf("Incompatible total particle numbers: %d,%d\n\n",np,np2);
	  exit(0);
	}
      if (nvp > nvpmax) nvpmax = nvp;
      fclose(part);
    }
   }
  }
  printf("We have %d particles to tie together.\n",np); fflush(stdout);
  printf("The maximum number of particles in a file is %d.\n",nvpmax);

  adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
  if (adjs == NULL) printf("Couldn't allocate adjs.\n");
  vols = (float *)malloc(np*sizeof(float));
  if (vols == NULL) printf("Couldn't allocate vols.\n");
  orig = (int *)malloc(nvpmax*sizeof(int));
  if (orig == NULL) printf("Couldn't allocate orig.\n");
  if ((vols == NULL) || (orig == NULL) || (adjs == NULL)) {
    printf("Not enough memory to allocate. Exiting.\n");
    exit(0);
  }
  for (p=0;p<np;p++)
    vols[p] = -1.;

  for (i = 0; i < numdiv; i++) {
   for (j = 0; j < numdiv; j++) {
    for (k = 0; k < numdiv; k++) {
      sprintf(partfile,"part.%s.%02d.%02d.%02d",taglabel,i,j,k);
      part = fopen(partfile,"r");
      if (part == NULL) {
	printf("Unable to open file %s.\n\n",partfile);
	exit(0);
      }
      fread(&np2,1,sizeof(int),part);
      fread(&nvp,1,sizeof(int),part);
      /*printf("nvp = %d\n",nvp);fflush(stdout);*/

      nvpsum += nvp;

      fread(orig,nvp,sizeof(int),part);
      for (p=0;p<nvp;p++) {
	fread(&volstemp,1,sizeof(float),part);
	if (vols[orig[p]] > -1.)
	  if (vols[orig[p]] != volstemp) {
	    printf("Inconsistent volumes for p. %d: (%g,%g)!\n",
		   orig[p],vols[orig[p]],volstemp);
	    /*exit(0);*/
	  }
	vols[orig[p]] = volstemp;
      }
      
      for (p=0;p<nvp;p++) {
	fread(&na,1,sizeof(int),part);
	if (na > 0) {
	  adjs[orig[p]].nadj = na;
	  adjs[orig[p]].adj = (int *)malloc(na*sizeof(int));
	  if (adjs[orig[p]].adj == NULL) {
	    printf("Couldn't allocate adjs[orig[%d]].adj.\n",p);
	    exit(0);
	  }
	  fread(adjs[orig[p]].adj,na,sizeof(int),part);
	} else {
	  printf("0"); fflush(stdout);
	}
      }
      fclose(part);
      printf("%d ",k);
    }
   }
  }
  printf("\n");
  npnotdone = 0; avgnadj = 0.; avgvol = 0.;
  for (p=0;p<np;p++) {
    if (vols[p] == -1.) npnotdone++;
    avgnadj += (double)(adjs[p].nadj);
    avgvol += (double)(vols[p]);
  }
  if (npnotdone > 0)
    printf("%d particles not done!\n");
  printf("%d particles done more than once.\n",nvpsum-np);
  avgnadj /= (double)np;
  avgvol /= (double)np;
  printf("Average # adjacencies = %lf (%f for Poisson)\n",avgnadj,
	 48.*3.141593*3.141593/35.+2.);
  printf("Average volume = %lf\n",avgvol);
    
  /* Now the output! */

  sprintf(adjfile,"%s%sadj.dat",outdir,taglabel);
  sprintf(volfile,"%s%svol.dat",outdir,taglabel);

  printf("Outputting to %s, %s\n\n",adjfile,volfile);

  adj = fopen(adjfile,"w");
  if (adj == NULL) {
    printf("Unable to open %s\n",adjfile);
    exit(0);
  }
  fwrite(&np,1, sizeof(int),adj);
  /* Adjacencies: first the numbers of adjacencies, 
     and the number we're actually going to write per particle */
  for (i=0;i<np;i++)
    fwrite(&adjs[i].nadj,1,sizeof(int),adj);
    
  /* Now the lists of adjacencies (without double counting) */
  for (i=0;i<np;i++)
    if (adjs[i].nadj > 0) {
      nout = 0;
      for (j=0;j<adjs[i].nadj; j++) if (adjs[i].adj[j] > i) nout++;
      fwrite(&nout,1,sizeof(int),adj);      
      for (j=0;j<adjs[i].nadj; j++) 
	if (adjs[i].adj[j] > i) 
	  fwrite(&(adjs[i].adj[j]),1,sizeof(int),adj);
    }

  fclose(adj);
  
  /* Volumes */
  vol = fopen(volfile,"w");
  if (vol == NULL) {
    printf("Unable to open %s\n",volfile);
    exit(0);
  }
  fwrite(&np,1, sizeof(int),vol);
  fwrite(vols,sizeof(float),np,vol);

  fclose(vol);

  return(0);
}
