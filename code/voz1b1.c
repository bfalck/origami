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

int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs);
int vorvol (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, float *vol);
int posread(char *posfile, float ***p, float fact);
int readgadget(char *filename, float ***p);
void read_parameters(char *filename);


int main(int argc, char *argv[]) {
  int exitcode;
  int i, j, np;
  float **r;
  coordT rtemp[3], *parts;
  coordT deladjs[3*MAXVERVER], points[3*MAXVERVER];
  pointT intpoints[3*MAXVERVER];
  FILE *pos, *out;
  char *paramfile, outfile[80];
  PARTADJ *adjs;
  float *vols;
  float predict, xmin,xmax,ymin,ymax,zmin,zmax;
  int *orig;
  
  int isitinbuf;
  char isitinmain, d;
  int nvp, nvpall, nvpbuf;
  float width, width2, totwidth, totwidth2, bf, s, g;
  float c[3];
  int b[3];
  double totalvol;

  if (argc != 5) {
    printf("Wrong number of arguments.\n");
    printf("arg1: parameter file\n");
    printf("arg2-4: b[0-2]\n\n");
    exit(0);
  }
  paramfile = argv[1];
  read_parameters(paramfile);

  if (numdiv == 1) {
    printf("Only using one division; should only use for an isolated segment.\n");
  }
  if (numdiv < 1) {
    printf("Cannot have a number of divisions less than 1. Resetting to 1.\n");
    numdiv = 1;
  }
  if (sscanf(argv[2],"%d",&b[0]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  if (sscanf(argv[3],"%d",&b[1]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  if (sscanf(argv[4],"%d",&b[2]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  
  /* Read the position file */
  if (numfiles > 0) {
    np = readgadget(posfile,&r);
    for (i=0; i<np; i++) DL r[i][d] /= boxsize;
  } else np = posread(posfile,&r,1./boxsize);
  /* Boxsize should be the range in r, yielding a range 0-1 */
  printf("%d particles\n",np);fflush(stdout);
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
    if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
    if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  width = 1./(float)numdiv;
  width2 = 0.5*width;
  if (buffer > 0.) bf = buffer;
  else bf = 0.1;
      /* In units of 0-1, the thickness of each subregion's buffer*/
  totwidth = width+2.*bf;
  totwidth2 = width2 + bf;
  
  s = width/(float)NGUARD;
  if ((bf*bf - 2.*s*s) < 0.) {
    printf("bf = %f, s = %f.\n",bf,s);
    printf("Not enough guard points for given border.\nIncrease guards to >= %f\n.",
	   sqrt(2.)*width/bf);
    exit(0);
  }
  g = (bf / 2.)*(1. + sqrt(1 - 2.*s*s/(bf*bf)));
  printf("s = %f, bf = %f, g = %f.\n",s,bf,g);
  
  fflush(stdout);

  adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
  if (adjs == NULL) {
    printf("Unable to allocate adjs\n");
    exit(0);
  }
  
  DL c[d] = ((float)b[d]+0.5)*width;
  printf("c: %f,%f,%f\n",c[0],c[1],c[2]);
  /* Assign temporary array*/
  nvpbuf = 0; /* Number of particles to tesselate, including
		 buffer */
  nvp = 0; /* Without the buffer */
  for (i=0; i<np; i++) {
    isitinbuf = 1;
    isitinmain = 1;
    DL {
      rtemp[d] = (double)r[i][d] - (double)c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinbuf = isitinbuf && (fabs(rtemp[d]) < totwidth2);
      isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
    }
  
    if (isitinbuf) nvpbuf++;
    if (isitinmain) nvp++;
  }
  
  nvpbuf += 6*(NGUARD+1)*(NGUARD+1); /* number of guard
					points */

  parts = (coordT *)malloc(3*nvpbuf*sizeof(coordT));
  orig = (int *)malloc(nvpbuf*sizeof(int));

  if (parts == NULL) {
    printf("Unable to allocate parts\n");
    fflush(stdout);
  }
  if (orig == NULL) {
    printf("Unable to allocate orig\n");
    fflush(stdout);
  }

  nvp = 0; nvpall = 0; /* nvp = number of particles without buffer */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np; i++) {
    isitinmain = 1;
    DL {
      rtemp[d] = r[i][d] - c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
    }
    if (isitinmain) {
      parts[3*nvp] = rtemp[0];
      parts[3*nvp+1] = rtemp[1];
      parts[3*nvp+2] = rtemp[2];
      orig[nvp] = i;
      nvp++;
      if (rtemp[0] < xmin) xmin = rtemp[0];
      if (rtemp[0] > xmax) xmax = rtemp[0];
      if (rtemp[1] < ymin) ymin = rtemp[1];
      if (rtemp[1] > ymax) ymax = rtemp[1];
      if (rtemp[2] < zmin) zmin = rtemp[2];
      if (rtemp[2] > zmax) zmax = rtemp[2];
    }
  }
  printf("nvp = %d\n",nvp);
  printf("x: %f,%f; y: %f,%f; z:%f,%f\n",xmin,xmax,ymin,ymax,zmin,zmax);
  nvpbuf = nvp;
  for (i=0; i<np; i++) {
    isitinbuf = 1;
    DL {
      rtemp[d] = r[i][d] - c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinbuf = isitinbuf && (fabs(rtemp[d])<totwidth2);
    }
    if ((isitinbuf > 0) &&
	((fabs(rtemp[0])>width2)||(fabs(rtemp[1])>width2)||(fabs(rtemp[2])>width2))) {
      
      /*printf("%3.3f ",sqrt(rtemp[0]*rtemp[0] + rtemp[1]*rtemp[1] +
	rtemp[2]*rtemp[2]));
	printf("|%2.2f,%2.2f,%2.2f,%f,%f",r[i][0],r[i][1],r[i][2],width2,totwidth2);*/
      parts[3*nvpbuf] = rtemp[0];
      parts[3*nvpbuf+1] = rtemp[1];
      parts[3*nvpbuf+2] = rtemp[2];
      orig[nvpbuf] = i;

      nvpbuf++;
      if (rtemp[0] < xmin) xmin = rtemp[0];
      if (rtemp[0] > xmax) xmax = rtemp[0];
      if (rtemp[1] < ymin) ymin = rtemp[1];
      if (rtemp[1] > ymax) ymax = rtemp[1];
      if (rtemp[2] < zmin) zmin = rtemp[2];
      if (rtemp[2] > zmax) zmax = rtemp[2];
    }
  }
  printf("nvpbuf = %d\n",nvpbuf);
  printf("x: %f,%f; y: %f,%f; z:%f,%f\n",xmin,xmax,ymin,ymax,zmin,zmax);
  nvpall = nvpbuf;
  predict = pow(width+2.*bf,3)*(float)np;
  printf("There should be ~ %f points; there are %d\n",predict,nvpbuf);

  for (i=0;i<np;i++) free(r[i]);
  free(r);
  
  /* Add guard points */
  for (i=0; i<NGUARD+1; i++) {
    for (j=0; j<NGUARD+1; j++) {
      /* Bottom */
      parts[3*nvpall]   = -width2 + (float)i * s;
      parts[3*nvpall+1] = -width2 + (float)j * s;
      parts[3*nvpall+2] = -width2 - g;
      nvpall++;
      /* Top */
      parts[3*nvpall]   = -width2 + (float)i * s;
      parts[3*nvpall+1] = -width2 + (float)j * s;
      parts[3*nvpall+2] = width2 + g;
      nvpall++;
    }
  }
  for (i=0; i<NGUARD+1; i++) { /* Don't want to overdo the corners*/
    for (j=0; j<NGUARD+1; j++) {
      parts[3*nvpall]   = -width2 + (float)i * s;
      parts[3*nvpall+1] = -width2 - g;
      parts[3*nvpall+2] = -width2 + (float)j * s;
      nvpall++;
      
      parts[3*nvpall]   = -width2 + (float)i * s;
      parts[3*nvpall+1] = width2 + g;
      parts[3*nvpall+2] = -width2 + (float)j * s;
      nvpall++;
    }
  }
  for (i=0; i<NGUARD+1; i++) {
    for (j=0; j<NGUARD+1; j++) {
      parts[3*nvpall]   = -width2 - g;
      parts[3*nvpall+1] = -width2 + (float)i * s;
      parts[3*nvpall+2] = -width2 + (float)j * s;
      nvpall++;
      
      parts[3*nvpall]   = width2 + g;
      parts[3*nvpall+1] = -width2 + (float)i * s;
      parts[3*nvpall+2] = -width2 + (float)j * s;
      nvpall++;
    }
  }
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=nvpbuf;i<nvpall;i++) {
    if (parts[3*i] < xmin) xmin = parts[3*i];
    if (parts[3*i] > xmax) xmax = parts[3*i];
    if (parts[3*i+1] < ymin) ymin = parts[3*i+1];
    if (parts[3*i+1] > ymax) ymax = parts[3*i+1];
    if (parts[3*i+2] < zmin) zmin = parts[3*i+2];
    if (parts[3*i+2] > zmax) zmax = parts[3*i+2];
  }
  
  printf("Added guard points to total %d points (should be %d)\n",nvpall,
	 nvpbuf + 6*(NGUARD+1)*(NGUARD+1));
  printf("x: %f,%f; y: %f,%f; z:%f,%f\n",xmin,xmax,ymin,ymax,zmin,zmax);
  
  /* Do tesselation*/
  printf("File read.  Tessellating ...\n"); fflush(stdout);
  exitcode = delaunadj(parts, nvp, nvpbuf, nvpall, &adjs);
  
  /* Now calculate volumes*/
  printf("Now finding volumes ...\n"); fflush(stdout);
  vols = (float *)malloc(nvp*sizeof(float));
  
  for (i=0; i<nvp; i++) { /* Just the original particles
			     Assign adjacency coordinate array*/
    /* Volumes */
    for (j = 0; j < adjs[i].nadj; j++)
      DL {
	deladjs[3*j + d] = parts[3*adjs[i].adj[j]+d] - parts[3*i+d];
	if (deladjs[3*j+d] < -0.5) deladjs[3*j+d]++;
	if (deladjs[3*j+d] > 0.5) deladjs[3*j+d]--;
      }
    
    exitcode = vorvol(deladjs, points, intpoints, adjs[i].nadj, &(vols[i]));
    vols[i] *= (float)np;
   /* if (i % 1000 == 0)
      printf("%d: %d, %f\n",i,adjs[i].nadj,vols[i]); */
  }

  /* Get the adjacencies back to their original values */

  for (i=0; i<nvp; i++)
    for (j = 0; j < adjs[i].nadj; j++)
      adjs[i].adj[j] = orig[adjs[i].adj[j]];
  
  totalvol = 0.;
  for (i=0;i<nvp; i++) {
    totalvol += (double)vols[i];
  }
  printf("Average volume = %g\n",totalvol/(float)nvp);
  
  /* Now the output!
     First number of particles */
  sprintf(outfile,"part.%s.%02d.%02d.%02d",taglabel,b[0],b[1],b[2]);

  printf("Output to %s\n\n",outfile);
  out = fopen(outfile,"w");
  if (out == NULL) {
    printf("Unable to open %s\n",outfile);
    exit(0);
  }
  fwrite(&np,1, sizeof(int),out);
  fwrite(&nvp,1, sizeof(int),out);
  printf("nvp = %d\n",nvp);

  /* Tell us where the original particles were */
  fwrite(orig,sizeof(int),nvp,out);
  /* Volumes*/
  fwrite(vols,sizeof(float),nvp,out);
  /* Adjacencies */
  for (i=0;i<nvp;i++) {
    fwrite(&(adjs[i].nadj),1,sizeof(int),out);
    if (adjs[i].nadj > 0)
      fwrite(adjs[i].adj,adjs[i].nadj,sizeof(int),out);
    else printf("0");
  }
  fclose(out);
  
  return(0);
}
