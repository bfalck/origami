//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

/* ===========================================================================
 * origamidiv.c    
 * 
 * WHAT IT DOES
 *   A version of origamitag.c that computes ORIGAMI morphology (crossing
 *   number: 0 void, 1 wall, 2 filament, 3 halo) for ONE sub-cube of the
 *   simulation instead of the whole box.
 *
 *   The box's initial particle grid (np1d per side) is cut into
 *   ndiv x ndiv x ndiv sub-cubes, each npd = np1d/ndiv particles per side.
 *   This program loads one sub-cube PLUS a "buffer" layer that is nbuf
 *   thick on every side (so np2 = npd + 2*nbuf per side), runs the usual
 *   crossing test, then throws the buffer away and saves only the
 *   sub-cube's own particles. The buffer is there so particles near the
 *   sub-cube's edge can still be compared with neighbours that belong to
 *   the next sub-cube over.
 *
 *   Run it once per sub-cube (divid = 0 .. ndiv^3 - 1); the runs are
 *   independent and can go in parallel. origamicomb.c then stitches the
 *   pieces back into one full-box tag file.
 *
 * TO RUN CALL: 
 *   ./origamidiv parameterfile ndiv divid
 * IN TERMINAL. Where:
 *     ndiv   number of sub-cubes per side (must divide np1d evenly)
 *     divid  which sub-cube to process. 
 *   Parameter file: same as upstream plus the new "nbuf" line (see
 *   readfiles.c). Note ndiv here comes from the command line, NOT from the
 *   "ndiv" line in the parameter file (that one is for the VOBOZ tessellation).
 *
 * INPUT   posfile (numfiles must be 0; see KNOWN ISSUES). The file holding particle 
 *         positions post dark matter simulation -- will be ref. in the parameter file.
 *         Format read by divread() in readfiles.c: 8-byte particle count, then all x,
 *         then all y, then all z as 4-byte floats, sorted by particle ID.
 * OUTPUT  outdir/<taglabel>_tag_<divid>.dat : npd^3 bytes, one morphology
 *         value (0-3) per particle, no header. Order: x fastest, then y,
 *         then z, within the sub-cube. 
 *
 * KNOWN ISSUES (found while annotating; marked [ISSUE] below)
 *   - Output name "<label>_tag_<divid>.dat" does not match what
 *     origamicomb.c reads ("<label>tag_<x><y><z>.dat"); files must be renamed.
 *   - The Gadget branch (numfiles > 0) loads the whole box with the
 *     non-subdivided reader, which does not fit the sub-cube layout used here.
 *   - With nsplit > 1, np2 must be divisible by nsplit or some rows are skipped.
 *
 * ======================================================================== */

#include <stdio.h>     /* file reading/writing and printf */
#include <stdlib.h>    /* malloc (memory allocation), exit, abs */
#include "params.h"    /* names of the global parameters filled by read_parameters() */

/* Shorthand macros (text substitutions done before compiling):
 *   DL          loop d over the 3 dimensions (0=x, 1=y, 2=z)
 *   BF          a huge number, used as a starting min/max
 *   max(A,B)    the larger of A and B
 *   goodmod(A,B) wrap index A back into the range 0..B-1 (periodic wrap),
 *               assuming A is at most one period outside that range */

#define DL for (d=0;d<3;d++)
#define BF 1e30
#define max(A,B) (((A)>(B)) ? (A):(B))
#define goodmod(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))

/* Returns 1 if h is negative, else 0 */
int isneg(int h) {
    return (int)(h < 0);
}

/* Converts 3-D grid coordinates (i,j,k) on a grid with ng points
 * per side into a single array position. x changes fastest, then y, then z. */
int par(int i, int j, int k, int ng) {
    return i + (j + k*ng)*ng;
}

/* Declarations of functions defined in readfiles.c, so the compiler
 * knows their argument types. divread() is Veyrat's sub-cube reader. */
void read_parameters(char *filename);

/* Reads the particle final pos file, asks for ***p -- where final positions are held -- 
 * total cube size, number of particles, number of particles per sub-cube, number of divisions,
 * sub-cube we're running on (divid), and buffer.
 * This allows code to determine initial positions which determine if particle is in subcube, 
 * then final positions of particles in the desired subcube are returned by divread in initial pos order */ 

long divread(char *filename, float ***p, float fact, float boxsize, int np1d, int npd, int ndiv, int divid, int nbuf);  


int readgadget(char *filename, float ***p);


/* divid is the subcube id this run is for
* if ndiv = 4 there are four subcubes per dimension so divid=0 is subcube (0,0,0), divid=1 is (1,0,0)
* divid = 3 is (3, 0, 0), divid = 4 is (0, 1, 0), divid = 5 is (1, 1, 0), etc. */

int main(int argc, char *argv[]) {
    /* variable declarations */ 
    
    long np;            /* particle count reported by the reader (total in the file) */ 
    float **r;          /* particle positions: r[i][0..2] = x,y,z of particle i */

    FILE *tag;          /* output file */ 
    char *paramfile,tagoutfile[80];    /* param file argv[1] points to and outputfile placeholder */ 
    float xmin,xmax,ymin,ymax,zmin,zmax; /* unused? */ 

    int s;              /* block # in OpenMP loop */ 
    float negb2,b2;     /* -boxsize/2 and +boxsize/2, for periodic distance wrapping */
    int ng2,ng4, np2,np3, h, i,i2, j, x,y,z,nhalo,nhalo0,nhalo1,nhalo2,nhaloany, ndiv,divid,npd;  
    
    unsigned char *m,*m0,*m1,*m2,*M, mn,m0n,m1n,m2n; /*Morphology tag */ 
    /*   m, m0, m1, m2: holds CN for each particle in cube (including buffer), one
     *   array per set of orthogonal axes (see the crossing section).
     * M: final CN array for the sub-cube's own particles only (buffer removed). */
   
    int zstart, ystart, xstart, zend, yend, xend;
    
    float dx,d1,d2; /* Difference in pos. for cartesian test and two diag. tests */ 

    /* END variable declarations */ 

    /*
     * argv[] holds the arguments typed after the program name:
     * argv[1] = parameter file, argv[2] = ndiv, argv[3] = divid. 
     * Quinn: Changed from argc < 2 to argc < 4 b/c 3 args are actually needed (param, ndiv, divid)
     * also swapped exit(0) for exit(1)
    */
    
    if (argc < 4) {
      printf("Not enough arguments. Need parameter file, ndiv, and divid. \n");
      exit(1);
    } else paramfile = argv[1];
    read_parameters(paramfile);   /* Fills posfile, outdir, boxsize, np1d, nbuf, ... */
    
    /* sscanf converts text to a number; it returns 1 on success. 
     * If the argument provided for ndiv or divid is not a number it returns 0 & triggers error msg
    */ 
    if (sscanf(argv[2],"%d",&ndiv) != 1) { 
        printf("That's no ndiv; try again.\n");  
        exit(1);  
    }  
    if (sscanf(argv[3],"%d",&divid) != 1) {  
        printf("That's no divid; try again.\n"); 
        exit(1);  
    }  
    /* 
    Grid of particles must split evenly among the subcubes aka np1d/ndiv has 0 remainder.
    */
    if (np1d%ndiv != 0) {  
        printf("ndiv must be factor of np1d.\n");  
	exit(1);  
    }  

            
    npd = np1d/ndiv;   /* Particles per side of one sub-cube (no buffer) */

    b2 = boxsize/2.; /* half sizes of entire box, not subcube */ 
    negb2 = -boxsize/2.;

    /* Builds the output filename as outdir/taglabel_tag_divid.dat
     * Requires outdir in params.txt to be trailed by a "/" 
     * EXAMPLE: If outdir = '/data/' and  taglabel = 'n256', output file is: "/data/n256_tag_5.dat".
    */
    sprintf(tagoutfile,"%s%s_tag_%d.dat",outdir,taglabel,divid); 

    np2 = npd+nbuf+nbuf;   /* Particles per side INCLUDING the buffer on both sides */  
    np3 = np2*np2*np2;     /* Total particles loaded (sub-cube + buffer) */ 

    /* ng4 is the maximum lattice distance searched for a crossing.
     * due to sub-cube logic we max search at the distance of the buffer
     * (a particle at edge of non-buffer cube can't cross with a particle past the buffer) 
    */
    ng4=nbuf; 

    /* ---- Load positions ----
     * divread() fills pointer r with the np3 particles of this sub-cube + buffer,
     * laid out on the np2^3 grid via par(x,y,z,np2).
     * divread() returns total particle count which is defined as np here 
     * [ISSUE] The readgadget() branch loads the full box in the full-box
     * layout, which does not match the np2 grid used below; use numfiles = 0. */
    if (numfiles > 0) {
      np = readgadget(posfile,&r);
    } else { 
        np = divread(posfile,&r,1.,boxsize,np1d,npd,ndiv,divid,nbuf); 
    } 
    
    /* Prints total particles in the file */
    printf("%ld particles\n",np);fflush(stdout); 

    /* ---- Allocate and initialise the tag arrays ----
     * malloc(n) reserves n bytes of memory and returns its address.
     * n is np3 (total particles (sub-cube + buff) * 1 byte (sizeof(unsigned char))
     * Every tag starts at 1; detecting a crossing multiplies it by a prime
     * (2, 3 or 5), one prime per axis, so later we can tell which axes
     * crossed by checking divisibility. */
    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
    m = (unsigned char *)malloc(np3*sizeof(unsigned char));  /* crossings in cartesian */ 
    m0 = (unsigned char *)malloc(np3*sizeof(unsigned char)); /* crossing w diag. in y-z plane */
    m1 = (unsigned char *)malloc(np3*sizeof(unsigned char));  /* crossing w diags. in x-z */
    m2 = (unsigned char *)malloc(np3*sizeof(unsigned char));  /* crossings w diags. in x-y */
    M = (unsigned char *)malloc(npd*npd*npd*sizeof(unsigned char));  

    /* == NULL means malloc failed (out of memory) 
    this check was previously mistakenly after the for loop, 
    also switched exit(0) to exit(1) */ 
    
    if (m == NULL || m0 == NULL || m1 == NULL || m2 == NULL || M == NULL) {
    printf("Morphology arrays cannot be allocated.\n");
    exit(1);
    }
    
    for (i=0; i<np3;i++) {  
        /* Track the position range (only used by a commented-out printf). */
        if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
        if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
        if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];

        m[i] = 1;
        m0[i] = 1;
        m1[i] = 1;
        m2[i] = 1;

	if (i<(npd*npd*npd)) { 
	    M[i] = 1;  
	} 
    }

    //  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
    printf("Calculating ORIGAMI morphology.\n");

    /* ==== Crossing detection ====
     * "#pragma omp parallel for" tells the compiler (with -fopenmp) to run the
     * iterations of the next loop on several CPU cores at once. The loaded
     * np2^3 grid is split into nsplit^3 blocks, one loop iteration per block.
     * shared(...) = variables all threads see; private(...) = each thread
     * gets its own copy.
     * [ISSUE] np2/nsplit is whole-number division, so if np2 is not a
     * multiple of nsplit the last few rows are never processed. If those
     * rows are inside the sub-cube (not buffer), their tags stay 0 (void). */
#pragma omp parallel for default(none) shared(np1d, np2, ng4, r, boxsize, negb2, b2, m, m0, m1, m2, nsplit) private (x, y, z, h, dx, i, i2, d1, d2, xstart, xend, ystart, yend, zstart, zend)  
    for (s = 0; s < nsplit*nsplit*nsplit; s++) {
        /* Turn block number s into block coordinates, then into
         * the grid range [start, end) that this block covers.
         * In C, "/" between whole numbers drops the remainder. */
        zstart = s / (nsplit*nsplit);
        zend = zstart + 1;
        ystart = (s - zstart * nsplit * nsplit)/nsplit;
        yend = ystart + 1;
        xstart = (s - zstart * nsplit * nsplit - ystart * nsplit);
        xend = xstart + 1;
        printf("%d\n",xend);      /* leftover debugging output */  
        zstart *= (np2/nsplit);  
        ystart *= (np2/nsplit);  
        xstart *= (np2/nsplit);  
        zend *= (np2/nsplit);  
        yend *= (np2/nsplit);  
        xend *= (np2/nsplit);  
	printf("%d\n",xend);          /* leftover debugging output */ 

        /* Visit every particle (x,y,z = its position on the
         * INITIAL grid) in this block. Buffer particles are processed too,
         * but their results are discarded at the end. */
        for (x=xstart; x<xend; x++){
            //    printf("%d\n",x);fflush(stdout);
            for (y=ystart; y<yend; y++) {
                for (z=zstart; z<zend; z++) {
                    i = par(x,y,z,np2);   /* array index of this particle */  

                    /* ---- How every test below works ----
                     * Walk along one line of the initial grid, checking the
                     * neighbour h steps away, with h = 1, -1, 2, -2, ...
                     * (h = -h + isneg(h) produces that sequence) until
                     * |h| reaches ng4. A crossing means the neighbour is now
                     * on the WRONG side: it started at +h (ahead) but its
                     * current position is behind (or vice versa). That is
                     * the test "dx*h < 0" (dx and h have opposite signs).
                     * The "while" lines fold dx into -boxsize/2..+boxsize/2
                     * so particles on opposite sides of the periodic box
                     * are measured the short way round.
                     * "break" stops the search at the first crossing found. */

                    /* First just along the Cartesian axes */
                    /* x-direction */
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),y,z,np2);  
                        dx = r[i2][0]-r[i][0];
                        while (dx < negb2) dx += boxsize;  
                        while (dx > b2) dx -= boxsize;  
                        if (dx*h < 0.) {
                            /* x crossing: multiply by 2 (only once) */
                            if (m[i] % 2 > 0) {
                                m[i] *= 2;
                            }
                            break;
                        }
                    }
                    /* y-direction: same test on the y coordinate; prime 3 */
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h, np2),z,np2);  
                        dx = r[i2][1]-r[i][1];
                        while (dx < negb2) dx += boxsize;  
                        while (dx > b2) dx -= boxsize;  
                        if (dx*h < 0.) {
                            /*printf("y:%d %d %d %d %f\n",x,y,z,h,dx);*/
                            if (m[i] % 3 > 0) {
                                m[i] *= 3;
                            }
                            break;
                        }
                    }
                    /* z-direction: prime 5 */
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(x,y,goodmod(z+h, np2),np2);  
                        dx = r[i2][2]-r[i][2];
                        while (dx < negb2) dx += boxsize;  
                        while (dx > b2) dx -= boxsize; 
                        if (dx*h < 0.) {
                            /*printf("z:%d %d %d %d %f\n",x,y,z,h,dx);*/
                            if (m[i] % 5 > 0) {
                                m[i] *= 5;
                            }
                            break;
                        }
                    }

                    /* ---- 45-degree diagonal axes ----
                     * Distance along a
                     * diagonal is the sum (d1+d2) or difference (d1-d2) of
                     * the two coordinate differences.
                     * m0 = diagonals in the y-z plane (pairs with x),
                     * m1 = x-z plane (pairs with y), m2 = x-y plane (pairs with z).
                     * Prime 2 = the "+,+" diagonal, prime 3 = the "+,-" diagonal. */

                    
                    /* x triplet: y-z plane, direction (0,+1,+1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h,np2),goodmod(z+h,np2),np2); 
                        d1 = r[i2][1]-r[i][1];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize;  
                        while (d1 > b2) d1 -= boxsize;  
                        while (d2 < negb2) d2 += boxsize; 
                        while (d2 > b2) d2 -= boxsize;  
                        if ((d1 + d2)*h < 0.) {
                            m0[i] *= 2;
                            break;
                        }
                    }
                    /* x triplet: y-z plane, direction (0,+1,-1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h,np2),goodmod(z-h,np2),np2);  
                        d1 = r[i2][1]-r[i][1];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize;  
                        while (d1 > b2) d1 -= boxsize; 
                        while (d2 < negb2) d2 += boxsize; 
                        while (d2 > b2) d2 -= boxsize; 
                        if ((d1 - d2)*h < 0.) {
                            m0[i] *= 3;
                            break;
                        }
                    }
                 
                    /* y triplet: x-z plane, direction (+1,0,+1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),y,goodmod(z+h,np2),np2);  
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize; 
                        while (d1 > b2) d1 -= boxsize; 
                        while (d2 < negb2) d2 += boxsize;
                        while (d2 > b2) d2 -= boxsize;  
                        if ((d1 + d2)*h < 0.) {
                            m1[i] *= 2;
                            break;
                        }
                    }
                    /* y triplet: x-z plane, direction (+1,0,-1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),y,goodmod(z-h,np2),np2); 
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize;  
                        while (d1 > b2) d1 -= boxsize;  
                        while (d2 < negb2) d2 += boxsize;  
                        while (d2 > b2) d2 -= boxsize;  
                        if ((d1 - d2)*h < 0.) {
                            m1[i] *= 3;
                            break;
                        }
                    }
                    
                    /* z triplet: x-y plane, direction (+1,+1,0) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),goodmod(y+h,np2),z,np2);   
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][1]-r[i][1];
                        while (d1 < negb2) d1 += boxsize; 
                        while (d1 > b2) d1 -= boxsize; 
                        while (d2 < negb2) d2 += boxsize; 
                        while (d2 > b2) d2 -= boxsize;  
                        if ((d1 + d2)*h < 0.) {
                            m2[i] *=2;
                            break;
                        }
                    }
                    /* z triplet: x-y plane, direction (+1,-1,0) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),goodmod(y-h,np2),z,np2);  
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][1]-r[i][1];
                        while (d1 < negb2) d1 += boxsize;  
                        while (d1 > b2) d1 -= boxsize; 
                        while (d2 < negb2) d2 += boxsize;  
                        while (d2 > b2) d2 -= boxsize;  
                        if ((d1 - d2)*h < 0.) {
                            m2[i] *= 3;
                            break;
                        }
                    }
                }
            }
        }
    }

    /* ==== Count crossed axes and keep only the sub-cube ====
     * For each triplet, count how many of its 3 axes crossed:
     *   (a % p == 0) is 1 if a is divisible by prime p (that axis crossed), else 0.
     *   mn  = Cartesian x,y,z           (m: primes 2,3,5)
     *   m0n = x + the two y-z diagonals (m prime 2, m0 primes 2,3)
     *   m1n = y + the two x-z diagonals (m prime 3, m1 primes 2,3)
     *   m2n = z + the two x-y diagonals (m prime 5, m2 primes 2,3)
     * The morphology is the largest of the four counts (0 void .. 3 halo).
     * The nhalo* counters count halo particles per triplet (just printed).
     * Note: they include buffer particles. */
    nhalo = 0;
    nhalo0 = 0;
    nhalo1 = 0;
    nhalo2 = 0;
    nhaloany = 0;
    j = 0;   /* next free slot in the output array M */  
    for (i=0;i<np3;i++){
        mn = (m[i]%2 == 0) + (m[i]%3 == 0) + (m[i]%5 == 0);
        m0n = (unsigned char)(m[i]%2 == 0) + (unsigned char)(m0[i]%2 == 0) + (unsigned char)(m0[i]%3 == 0);
        m1n = (unsigned char)(m[i]%3 == 0) + (unsigned char)(m1[i]%2 == 0) + (unsigned char)(m1[i]%3 == 0);
        m2n = (unsigned char)(m[i]%5 == 0) + (unsigned char)(m2[i]%2 == 0) + (unsigned char)(m2[i]%3 == 0);
        m[i] = max(mn,max(m0n,max(m1n,m2n)));
        if (mn == 3) nhalo++;
        if (m0n == 3) nhalo0++;
        if (m1n == 3) nhalo1++;
        if (m2n == 3) nhalo2++;
        if (m[i] == 3) {
            nhaloany++;
        }
	/* Recover this particle's grid coordinates (inverse of par())
	 * and shift by nbuf so the sub-cube itself runs 0..npd-1. If any
	 * coordinate is outside that range the particle is in the buffer:
	 * "continue" skips to the next particle without saving it. */
	x = (i%np2)-nbuf; 
	if (x < 0 || x >= npd) {  
	    continue;  
	}  
	y = ((i/np2)%np2)-nbuf;  
	if (y < 0 || y >= npd) {  
	    continue;  
	}  
	z = (i/(np2*np2))-nbuf;  
	if (z < 0 || z >= npd) {  
	    continue; 
	}  
	/* Sub-cube particle: save it. Because i runs x-fastest,
	 * M ends up in x-fastest, then y, then z order within the sub-cube. */
	M[j] = m[i]; 
	++j;  
    }
    printf("nhalo=%d,%d,%d,%d,%d\n",nhalo,nhalo0,nhalo1,nhalo2,nhaloany);

    /* Output m */
    /* ==== Write the sub-cube's tags ====
     * fopen(..., "w") creates/overwrites the file; fwrite copies npd^3
     * bytes from M into it. Unlike upstream origamitag, no particle-count
     * header is written (that line is commented out); origamicomb.c adds
     * one to the combined file. */
    printf("Writing morphology file to %s.\n",tagoutfile);fflush(stdout);
    tag = fopen(tagoutfile,"w");
    if (tag == NULL) {
        printf("Unable to open %s\n",tagoutfile);
        exit(0);
    }
    //fwrite(&np,1, sizeof(int),tag);
    fwrite(M,npd*npd*npd,sizeof(unsigned char),tag); 

    return(1);
}
