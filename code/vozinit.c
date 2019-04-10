//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

#include "qhull_a.h"
#include "voz.h"
#include "params.h"

#define DL for (d=0;d<3;d++)
#define BF 1e30

void read_parameters(char *filename);
int posread(char *posfile, float ***p, float fact);
int readgadget(char *filename, float ***p);

int main(int argc, char *argv[]) {
  int i, np;
  float **rfloat, rtemp[3];
  FILE *pos, *scr;
  char scrfile[80], systemstr[90],*paramfile;
  float xmin,xmax,ymin,ymax,zmin,zmax;

  int isitinbuf;
  char isitinmain, d;
  int nvp, nvpall, nvpbuf, nvpmin, nvpmax, nvpbufmin, nvpbufmax; /* yes, the insurance */
  float width, width2, totwidth, totwidth2, bf, s, g;
  float c[3];
  int b[3];

  if (argc < 2) {
    printf("Please specify a parameter file\n");
    exit(0);
  } else paramfile = argv[1];
  read_parameters(paramfile);

  if (numdiv < 2) {
    printf("Cannot have a number of divisions less than 2. Resetting to 2:\n");
    numdiv = 2;
  }

  /* Read the position file */
  if (numfiles > 0) {
    np = readgadget(posfile,&rfloat);
    for (i=0; i<np; i++) DL rfloat[i][d] /= boxsize;
  } else np = posread(posfile,&rfloat,1./boxsize);
  /* Boxsize should be the range in r, yielding a range 0-1 */

  width = 1./(float)numdiv;
  width2 = 0.5*width;
  if (buffer > 0.) bf = buffer;
  else bf = 0.1;

  /* In units of 0-1, the thickness of each subregion's buffer*/
  totwidth = width+2.*bf;
  totwidth2 = width2 + bf;

  s = width/(float)NGUARD;
  if ((bf*bf - 2.*s*s) < 0.) {
    printf("Not enough guard points for given border.\nIncrease guards to >= %f\n.",
	   totwidth/sqrt(0.5*bf*bf));
    printf("bf = %f\n",bf);
    exit(0);
  }
  g = (bf / 2.)*(1. + sqrt(1 - 2.*s*s/(bf*bf)));
  printf("s = %f, bf = %f, g = %f.\n",s,bf,g);
  
  nvpmax = 0; nvpbufmax = 0; nvpmin = np; nvpbufmin = np;
  
  for (b[0] = 0; b[0] < numdiv; b[0]++) {
   c[0] = ((float)b[0]+0.5)*width;
   for (b[1] = 0; b[1] < numdiv; b[1]++) {
    c[1] = ((float)b[1]+0.5)*width;
    for (b[2] = 0; b[2] < numdiv; b[2]++) {
      c[2] = ((float)b[2]+0.5)*width;

      nvp = 0; /* Number of particles excluding buffer */
      nvpbuf = 0; /* Number of particles to tesselate, including
		     buffer */
      xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
      for (i=0; i<np; i++) {
	isitinbuf = 1; isitinmain = 1;
	for (d=0; d<3; d++) {
	  rtemp[d] = rfloat[i][d] - c[d];
	  if (rtemp[d] > 0.5) rtemp[d] --;
	  if (rtemp[d] < -0.5) rtemp[d] ++;
	  isitinbuf = isitinbuf && (fabs(rtemp[d]) < totwidth2);
	  isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
	}
	if (isitinbuf) {
	  nvpbuf++;
	}
	if (isitinmain) {
	  nvp++;
	  if (rtemp[0] < xmin) xmin = rtemp[0];
	  if (rtemp[0] > xmax) xmax = rtemp[0];
	  if (rtemp[1] < ymin) ymin = rtemp[1];
	  if (rtemp[1] > ymax) ymax = rtemp[1];
	  if (rtemp[2] < zmin) zmin = rtemp[2];
	  if (rtemp[2] > zmax) zmax = rtemp[2];
	}
      }
      if (nvp > nvpmax) nvpmax = nvp;
      if (nvpbuf > nvpbufmax) nvpbufmax = nvpbuf;
      if (nvp < nvpmin) nvpmin = nvp;
      if (nvpbuf < nvpbufmin) nvpbufmin = nvpbuf;

      printf("b=(%d,%d,%d), c=(%f,%f,%f), nvp=%d, nvpbuf=%d\n",
	     b[0],b[1],b[2],c[0],c[1],c[2],nvp,nvpbuf);
    }
   }
  }
  printf("Nvp range: %d,%d\n",nvpmin,nvpmax);
  printf("Nvpbuf range: %d,%d\n",nvpbufmin,nvpbufmax);

  /* Output script file */
  sprintf(scrfile,"%sscr",taglabel);
  printf("Writing script file to %s.\n",scrfile);fflush(stdout);
  scr = fopen(scrfile,"w");
  if (scr == NULL) {
    printf("Problem opening script file.\n");
    fflush(stdout);
    exit(0);
  }
  fprintf(scr,"#!/bin/bash -f\n");
  for (b[0]=0;b[0]<numdiv; b[0]++) {
   for (b[1] = 0; b[1] < numdiv; b[1]++) {
    for (b[2] = 0; b[2] < numdiv; b[2]++) {
      fprintf(scr,"./voz1b1 %s %d %d %d\n",
	     paramfile,b[0],b[1],b[2]);
    }
   }
  }
  fprintf(scr,"./voztie %s\n",paramfile);

  /* Delete temporary 'part' files: */
  for (b[0]=0;b[0]<numdiv; b[0]++) {
   for (b[1] = 0; b[1] < numdiv; b[1]++) {
    for (b[2] = 0; b[2] < numdiv; b[2]++) {
      fprintf(scr,"rm -f part.%s.%02d.%02d.%02d\n",taglabel,b[0],b[1],b[2]);
    }
   }
  }


  fclose(scr);

  sprintf(systemstr,"chmod u+x %s",scrfile);
  system(systemstr);

  return(0);
}
