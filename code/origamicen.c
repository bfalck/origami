//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "voz.h"
#include "params.h"

#define BF 1e30

typedef struct Halo {
  int nmem;
  float hcent[3];
  int *mem;
} HALO;

int posread(char *posfile, float ***p, float fact);


int main(int argc, char *argv[]) {

  FILE *pos,*vol,*hid,*cat,*nhalo;
  HALO *h;
  char volfile[80],hidfile[80],catfile[80],*paramfile;

  int np,nv,*haloid,nh,thisn;
  float **p, *v, wtot;

  float b2,negb2,temp[3],thisp;
  int i,j,k,l,l0;

  if (argc < 2) {
    printf("Please specify a parameter file\n");
    exit(0);
  } else paramfile = argv[1];
  read_parameters(paramfile);

  sprintf(volfile,"%s%svol.dat",outdir,taglabel);
  sprintf(hidfile,"%s%shns.dat",outdir,halolabel);
  sprintf(catfile,"%s%scen.dat",outdir,halolabel);
  b2 = boxsize/2.;
  negb2 = -boxsize/2.;

  /* Read files: */
  if (numfiles > 0) {
    np = readgadget(posfile,&p);
  } else np = posread(posfile,&p,1.);
  printf("%d particles\n",np);fflush(stdout);


  vol = fopen(volfile,"r");
  if (vol == NULL) {
    printf("Unable to open volume file %s.\n\n",volfile);
    exit(0);
  }
  fread(&nv,1, sizeof(int),vol);
  if (np != nv) {
    printf("Number of particles doesn't match! %d != %d\n",np,nv);
    exit(0);
  }
  v = (float *)malloc(np*sizeof(float));
  fread(v,np,sizeof(float),vol);
  fclose(vol);

  hid = fopen(hidfile,"r");
  if (hid == NULL) {
    printf("Unable to open hn file %s.\n\n",hidfile);
    exit(0);
  }
  fread(&nh,1,sizeof(int),hid);
  h = (HALO *)malloc(nh*sizeof(HALO));
  for (i=0; i<nh; i++) {
    fread(&h[i].nmem,1,sizeof(int),hid);
    h[i].mem = (int *)malloc(h[i].nmem*sizeof(int));
    k = h[i].nmem;
    for (j=0;j<k;j++) {
      fread(&l,1,sizeof(int),hid);
      h[i].mem[j] = l;
    }
  }
  fclose(hid);

  printf("%d halos\n",nh);

  printf("Calculating halo centers etc.\n");
  for (i=0; i<nh; i++) {
    temp[0] = 0.; temp[1] = 0.; temp[2] = 0.;
    l0 = h[i].mem[0];
    thisn = h[i].nmem;
    wtot = 0.;
    for (j=0; j<thisn; j++) {
      l = h[i].mem[j];
      for (k=0; k<3; k++) {
	thisp = p[l][k] - p[l0][k];
	thisp += (thisp <= negb2)*boxsize;
	thisp -= (thisp > b2)*boxsize;
	//	temp[k] += thisp;
	temp[k] += thisp/v[l];
      }
      wtot += 1./v[l];
    }
    thisp = (float )h[i].nmem;
    for (k=0;k<3;k++) {
      //      h[i].hcent[k] = temp[k]/thisp + p[l0][k];
      h[i].hcent[k] = temp[k]/wtot + p[l0][k];
      h[i].hcent[k] += (h[i].hcent[k] <= 0)*boxsize;
      h[i].hcent[k] -= (h[i].hcent[k] > boxsize)*boxsize;
    }
  }

//  printf("i = %d\n",i);

  /* Output */
  printf("Writing centers.\n");
  cat = fopen(catfile,"w");
  fwrite(&nh,1,sizeof(int),cat);
  //  for (i=0; i<nh; i++) fwrite(&h[i].nmem,1,sizeof(int),cat);
  for (i=0; i<nh; i++) fwrite(&h[i].hcent[0],1,sizeof(float),cat);
  for (i=0; i<nh; i++) fwrite(&h[i].hcent[1],1,sizeof(float),cat);
  for (i=0; i<nh; i++) fwrite(&h[i].hcent[2],1,sizeof(float),cat);
  fclose(cat);




}
