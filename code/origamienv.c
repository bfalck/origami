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
#include "params.h"

#define BIGFLT 1e30 /* Biggest possible floating-point number */

typedef struct Particle {
  float vol;
  int nadj;
  char tag;
  int *adj;
} PARTICLE;

typedef struct Halo {
  int nmem;
  float hcent[3];
  int *mem;
} HALO;

int main(int argc,char **argv) {

  FILE *adj, *tag, *hid, *env, *nbr;
  PARTICLE *p;
  HALO *h;
  char adjfile[80],tagfile[80],hidfile[80],envfile[80],nbrfile[80];
  char *paramfile;
  int i,j,k,l,ll, np,nin, adjj, nh;
  int *haloid, *hnbrnum, hnbr[8000], hnum;
  float *hfrac, *ffrac, *wfrac, *vfrac, thisn;
  int inArray;
  unsigned char *htag;
  float havg, favg, wavg, vavg;

  if (argc < 2) {
    printf("Please specify a parameter file\n");
    exit(0);
  } else paramfile = argv[1];
  read_parameters(paramfile);

  sprintf(adjfile,"%s%sadj.dat",outdir,taglabel);
  sprintf(tagfile,"%s%stag.dat",outdir,taglabel);
  sprintf(hidfile,"%s%shns.dat",outdir,halolabel);
  sprintf(envfile,"%s%senv.dat",outdir,halolabel);
  sprintf(nbrfile,"%s%shnbr.dat",outdir,halolabel);

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

  hid = fopen(hidfile,"r");
  if (hid == NULL) {
    printf("Unable to open hid file %s.\n\n",hidfile);
    exit(0);
  }
  fread(&nh,1,sizeof(int),hid);
  h = (HALO *)malloc(nh*sizeof(HALO));
  for (i=0; i<nh; i++) {
    fread(&h[i].nmem,1,sizeof(int),hid);
    h[i].mem = (int *)malloc(h[i].nmem*sizeof(int));
    k = h[i].nmem;
    for (j=0; j<k; j++) {
      fread(&l,1,sizeof(int),hid);
      h[i].mem[j] = l;
    }
  }
  fclose(hid);


  /* Calculate fraction of Delaunay neighbors with each M for every 
     halo particle, and exclude neighbors that are part of same halo */

  hfrac = (float *)malloc(np*sizeof(float));
  ffrac = (float *)malloc(np*sizeof(float));
  wfrac = (float *)malloc(np*sizeof(float));
  vfrac = (float *)malloc(np*sizeof(float));

  haloid = (int *)malloc(np*sizeof(int));
  for (i=0; i<np; i++) {
    haloid[i] = 0;
    vfrac[i] = 0;
    wfrac[i] = 0;
    ffrac[i] = 0;
    hfrac[i] = 0;
  }
  /* Get halo id of each particle (0 if not in halo) */
  for (i=0; i<nh; i++) {
    k = h[i].nmem;
    for (j=0; j<k; j++) haloid[h[i].mem[j]] = i+1;
  }

  hnbrnum = (int *)malloc(nh*sizeof(int));
  for (i=0; i<nh; i++) {
    /* re-set ids of neighboring halos to 0 */
    for (ll=0; ll<80; ll++) hnbr[ll] = -1;
    hnbrnum[i] = 0;
    hnum = 0; // running index of hnbr array
    for (j=0; j<h[i].nmem; j++) {
      l = h[i].mem[j];
      thisn = 0;
      for (k=0; k<p[l].nadj; k++) {
	adjj = p[l].adj[k];
	if (haloid[adjj] != (i+1)) {
	  thisn++;
	  if (p[adjj].tag == 3) hfrac[l]++;
	  if (p[adjj].tag == 2) ffrac[l]++;
	  if (p[adjj].tag == 1) wfrac[l]++;
	  if (p[adjj].tag == 0) vfrac[l]++;
	}
	// If halo neighbor is new, update count 
	if ((haloid[adjj] != (i+1))&&(p[adjj].tag == 3)&&(haloid[adjj]!=0)) {
	  if (hnum == 0) {
	    hnbr[0] = haloid[adjj];
	    //	    printf("first hnbr = %d",hnbr[hnbrnum[i]]);
	    //	    hnbrnum[i]++;
	    hnum++;
	  } else {
	    inArray = 0;
	    for (ll=0; ll<hnum; ll++) {
	      if (hnbr[ll] == haloid[adjj]) {
		inArray = 1;
		break;
	      }
	    }
	    if (inArray == 0) {
	      hnbr[hnum] = haloid[adjj];
	      hnum++;
	    }
	  } 
	}
      } // kth neighbor of particle l
      if (thisn != 0) {
	hfrac[l] /= thisn;
	ffrac[l] /= thisn;
	wfrac[l] /= thisn;
	vfrac[l] /= thisn;
      }
    } // particle l in halo i+1
    hnbrnum[i] = hnum;
    //    printf("hnbrnum = %d",hnbrnum[i]);
  } // halo i+1


  htag = (unsigned char *)malloc(nh*sizeof(unsigned char));
  for (i=0; i<nh; i++) {
    // If adjacent to 3 or more halos, halo is in a cluster:
    if (hnbrnum[i] > 2) {
      htag[i] = 3;
    } else {
      // calculate average of v, w, f, and h fracs for all particles in halo:
      thisn = 0;
      havg = 0;
      favg = 0;
      wavg = 0;
      vavg = 0;
      for (j=0; j<h[i].nmem; j++) {
	l = h[i].mem[j];
	// exclude center, where all fracs are 0:
	if (!((hfrac[l] == 0) && (ffrac[l] == 0) && (wfrac[l] == 0) && (vfrac[l] == 0))) {
	  thisn++;
	  havg += hfrac[l];
	  favg += ffrac[l];
	  wavg += wfrac[l];
	  vavg += vfrac[l];
	  //  vavg += (1. - hfrac[l] - ffrac[l] - wfrac[l]);
	}
      } // particle l in halo i+1
      // Halo env't corresponds to max of avg. fracs:
      if (thisn != 0) {
	havg /= thisn;
	favg /= thisn;
	wavg /= thisn;
	vavg /= thisn;
      }
      if ((havg > favg) && (havg > wavg) && (havg > vavg)) htag[i] = 3;
      if ((favg > havg) && (favg > wavg) && (favg > vavg)) htag[i] = 2;
      if ((wavg > havg) && (wavg > favg) && (wavg > vavg)) htag[i] = 1;
      if ((vavg > havg) && (vavg > favg) && (vavg > wavg)) htag[i] = 0;
    }
  } // halo i+1



  /* Write output: */
  printf("Writing environment file %s\n",envfile);

  env = fopen(envfile,"w");
  fwrite(&nh,1,sizeof(int),env);
  fwrite(htag,nh,sizeof(unsigned char),env);
  //  fwrite(&np,1,sizeof(int),env);
  //  fwrite(hfrac,np,sizeof(float),env);
  //  fwrite(ffrac,np,sizeof(float),env);
  //  fwrite(wfrac,np,sizeof(float),env);
  fclose(env);

  //  nbr = fopen(nbrfile,"w");
  //  fwrite(&nh,1,sizeof(int),nbr);
  //  fwrite(hnbrnum,nh,sizeof(int),nbr);
  //  fclose(nbr);

} // main
