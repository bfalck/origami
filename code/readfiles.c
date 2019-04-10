//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "params.h"

#define DL for (d=0;d<3;d++) /* Dimension loop */
#define BF 1e30

/* Gadget variables: */
struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

int     NumPart, Ngas;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P;

long long *Id;

double  Time, Redshift;

/* Global parameters */
int maxline = 80;
char *posfile;
char *outdir;
char *taglabel;
float boxsize;
int np1d;
int nsplit;
int numfiles;
float buffer;
int numdiv;
float volcut;
int npmin;
char *halolabel;

void read_parameters(char *filename) {

  FILE *fd;
  char inbuf[maxline];

  posfile = (char *)malloc(maxline*sizeof(char));
  outdir = (char *)malloc(maxline*sizeof(char));
  taglabel = (char *)malloc(maxline*sizeof(char));
  halolabel = (char *)malloc(maxline*sizeof(char));

  fd = fopen(filename,"r");
  if (fd == NULL) {
    printf("Unable to open parameter file %s\n\n",filename);
    exit(0);
  }
  printf("Reading parameter file %s\n",filename);

  fgets(inbuf,maxline,fd);

  fscanf(fd,"%s %s",inbuf,posfile);
  fscanf(fd,"%s %s",inbuf,outdir);
  fscanf(fd,"%s %s",inbuf,taglabel);
  fscanf(fd,"%s %f",inbuf,&boxsize);
  fscanf(fd,"%s %d",inbuf,&np1d);
  fscanf(fd,"%s %d",inbuf,&nsplit);
  fscanf(fd,"%s %d",inbuf,&numfiles);

  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);

  fscanf(fd,"%s %f",inbuf,&buffer);
  fscanf(fd,"%s %d",inbuf,&numdiv);

  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);

  fscanf(fd,"%s %f",inbuf,&volcut);
  fscanf(fd,"%s %d",inbuf,&npmin);
  fscanf(fd,"%s %s",inbuf,halolabel);

  fclose(fd);


}


/* Positions */
/* Returns number of particles read */
int posread(char *posfile, float ***p, float fact) {

  FILE *pos;
  int np,dum,d,i;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *ptemp;

  //printf("fact=%f\n",fact);
  fflush(stdout);

  pos = fopen(posfile, "r");
  if (pos == NULL) {
    printf("Unable to open position file %s\n\n",posfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  /*fread(&dum,1,4,pos); */
  fread(&np,1, sizeof(int),pos); 
  /*fread(&dum,1,4,pos);*/

  /* Allocate the arrays */
  (*p) = (float **)malloc(np*sizeof(float *));
  ptemp = (float *)malloc(np*sizeof(float));

  //printf("np = %d\n",np);

  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) {
    (*p)[i] = (float *)malloc(3*sizeof(float));
    if ((*p)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
    (*p)[i][0] = ptemp[i];
  }
  /*fread(&dum,1,4,pos); 
    fread(&dum,1,4,pos); */
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) (*p)[i][1] = ptemp[i];
  /*fread(&dum,1,4,pos);
    fread(&dum,1,4,pos); */
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) (*p)[i][2] = ptemp[i];
  /*fread(&dum,1,4,pos); */

  fclose(pos);
  free(ptemp);

  /* Get into physical units */

  for (i=0; i<np; i++) DL (*p)[i][d] *= fact;


  /* Test range -- can comment out */
  /*
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
    if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
    if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
  */
  return(np);
}

/* Velocities */
/* Returns number of particles read */
int velread(char *velfile, float ***v, float fact) {

  FILE *vel;
  int np,dum,d,i;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *vtemp;

  vel = fopen(velfile, "r");
  if (vel == NULL) {
    printf("Unable to open velocity file %s\n\n",velfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  /*fread(&dum,1,4,vel);*/ fread(&np,1, sizeof(int),vel); /*fread(&dum,1,4,vel);*/
  printf("np=%d\n",np); fflush(stdout);

  /* Allocate the arrays */
  (*v) = (float **)malloc(np*sizeof(float *));
  vtemp = (float *)malloc(np*sizeof(float));

  printf("about to read velocity array\n"); fflush(stdout);

  /* Fill the arrays */
  /*fread(&dum,1,4,vel);*/
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) {
    (*v)[i] = (float *)malloc(3*sizeof(float));
    if ((*v)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
    (*v)[i][0] = vtemp[i];
  }
  /*fread(&dum,1,4,vel); 
    fread(&dum,1,4,vel); */
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) (*v)[i][1] = vtemp[i];
  /*fread(&dum,1,4,vel);
    fread(&dum,1,4,vel); */
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) (*v)[i][2] = vtemp[i];
  /*fread(&dum,1,4,vel); */
  fclose(vel);

  /* Convert from code units into physical units (km/sec) */
  
  //printf("np=%d\n",np); fflush(stdout);
  for (i=0; i<np; i++) DL (*v)[i][d] *= fact;

  /* Test range -- can comment out */
  /*
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*v)[i][0] < xmin) xmin = (*v)[i][0]; if ((*v)[i][0] > xmax) xmax = (*v)[i][0];
    if ((*v)[i][1] < ymin) ymin = (*v)[i][1]; if ((*v)[i][1] > ymax) ymax = (*v)[i][1];
    if ((*v)[i][2] < zmin) zmin = (*v)[i][2]; if ((*v)[i][2] > zmax) zmax = (*v)[i][2];
  }
  printf("vx: %f,%f; vy: %f,%f; vz: %f,%f\n",xmin,xmax, ymin,ymax, zmin,zmax);fflush(stdout);
  */
  return(np);
}


int readgadget(char *filename, float ***pos) {

  int i,d;

  load_snapshot(filename,numfiles);
  printf("NumPart = %d\n",NumPart);

  reordering();  /* call this routine only if your IDs are set properly */

  // Define p from *P:

  /* Allocate the arrays */
  (*pos) = (float **)malloc(NumPart*sizeof(float *));
//  printf("Allocated memory for pos in readgadget\n");

  for (i=0; i<NumPart; i++) {
    (*pos)[i] = (float *)malloc(3*sizeof(float));
    if ((*pos)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
  }
  printf("Allocated memory for pos in readgadget\n");

  P++; /* start with offset 0 */
  for (i=0; i<NumPart; i++) {
    (*pos)[i] = P[i].Pos;
  }


  return(NumPart);

}



/* (From Volker Springel's Gadget snapshot reader)
 * this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.)
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
    

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&Id[pc_new], sizeof(long long), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses>0)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;

	      if(header1.mass[k]==0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
      

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }

  Time= header1.time;
  Redshift= header1.time;
}

/* (From Volker Springel's Gadget snapshot reader)
 * this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(long long))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
}

/* (From Volker Springel's Gadget snapshot reader)
 * This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;   
  free(Id);

  printf("space for particle ID freed\n");
}

