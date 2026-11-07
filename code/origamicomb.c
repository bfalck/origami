#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

  FILE **sub, *tag;
  char *label, *outdir, tagsubfile[80], tagoutfile[80];

  int num_part, div_factor;
  int i,j,k,ijk,l,o,n, np3,df2,df3, n1d1,n1d3,n2d2,n3d3,n3d4,n3d5, irem,idiv;
  int *indices;
  unsigned char *m1, *m2;

  label = argv[1];
  outdir = argv[2];
  sscanf(argv[3],"%d",&num_part);
  sscanf(argv[4],"%d",&div_factor);

  np3 = num_part*num_part*num_part;
  df2 = div_factor*div_factor;
  df3 = df2*div_factor;

  if (num_part % df3 != 0) {
    printf("Particle grid size must be divisible by number of subdivisions.\n");
    exit(0);
  }

  n1d1 = num_part/div_factor;
  n1d3 = num_part/df3;
  n2d2 = n1d1*n1d1;
  n3d3 = n2d2*n1d1;
  n3d4 = n3d3/div_factor;
  n3d5 = n3d4/div_factor;

  indices = (int *)malloc(n3d3*sizeof(int));
  i = 0;
  o = 0;
  while (i<n3d5) {
    j = 0;
    while (j<n3d3) {
      k = 0;
      while (k<n2d2) {
        l = 0;
	    ijk = i+j+k;
	    while (l<n3d4) {
	      for (n=0; n<n1d1; n++) {
            indices[o] = ijk+l+n;
	        ++o;
	      }
	      l += n3d5;
	    }
	    k += n1d1;
      }
      j += n3d4;
    }
    i += n2d2;
  }

  sub = (FILE **)malloc(df3*sizeof(FILE *));
  for (i=0; i<div_factor; i++) {
    for (j=0; j<div_factor; j++) {
      for (k=0; k<div_factor; k++) {
        sprintf(tagsubfile,"%s%s_tag_%d.dat",outdir,label,i+(div_factor*j)+(df2*k)); //fixed filename mismatch with oriagmidiv
	    ijk = (df2*i)+(div_factor*j)+k;
	    sub[ijk] = fopen(tagsubfile,"r");
		//Safety check that gives clearer error message than 'Segmentation fault'
		if (sub[ijk] == NULL) {			
			printf("Unable to open sub-cube file %s\n",tagsubfile);
			exit(0);
		}
      }
    }
  }

  m1 = (unsigned char *)malloc(n3d3*sizeof(unsigned char));
  m2 = (unsigned char *)malloc(n3d5*sizeof(unsigned char));
  sprintf(tagoutfile,"%s%stag.dat",outdir,label);
  tag = fopen(tagoutfile,"w"); // program now overwrites toytag.dat (used to append to previous toytag.dat)
  fwrite(&np3,1,sizeof(int),tag);  //previously was 'fwrite(&np3,1,sizeof(int));
  for (i=0; i<df3; i++) {
    irem = i%df2;
    idiv = i/df2;
    o = 0;
    for (j=0; j<div_factor; j++) {
      for (k=0; k<div_factor; k++) {
	    ijk = (df2*k)+(div_factor*j)+idiv;
	    fread(m2,sizeof(unsigned char),n3d5,sub[ijk]);
	    for (l=0; l<n3d5; l++) {
	      m1[o] = m2[l];
	      ++o;
 	    }
      }
    }
    for (o=0; o<n3d3; o++) {
      fwrite(&m1[indices[o]],1,sizeof(unsigned char),tag);
    }
  }

  return(1);
}
