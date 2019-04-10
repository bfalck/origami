# ORIGAMI

---

ORIGAMI is a structure-finding routine that identifies halo, filament, 
wall, and void particles in an *N*-body simulation by counting the 
number of orthogonal axes along which a particle has undergone 
shell-crossing. Halo particles are then grouped by connecting 
neighbors on a Delaunay tessellation, and the web environment of 
a halo can be found by checking the morphology of its Delaunay 
neighbors.

This document provides instructions for the installation and 
use of ORIGAMI. For details of the algorithm, 
please see the [code paper](https://ui.adsabs.harvard.edu/abs/2012ApJ...754..126F/abstract) 
[(arXiv:1201.2353)](https://arxiv.org/abs/1201.2353). The ORIGAMI code consists of a few 
simple C routines and makes use of the Delaunay tessellation 
calculation routines from the [VOBOZ package](http://skysrv.pha.jhu.edu/~neyrinck/zobov), 
which itself relies on the [Qhull package](http://www.qhull.org). 


### Authors
ORIGAMI was developed by Bridget Falck & Mark Neyrinck and parallelized by Nuala McCullagh.

### Acknowledgments
ORIGAMI is released under the Modified BSD license; for details see LICENSE. 
If you use ORIGAMI, please cite [Falck, Neyrinck, & Szalay 2012, ApJ, 754, 126](https://ui.adsabs.harvard.edu/abs/2012ApJ...754..126F/abstract).


---
## Installing and Compiling
---

The ORIGAMI tarball contains the necessary routines from 
VOBOZ to run the tessellation, as well as version 2012.1 of 
the Qhull package required by VOBOZ. The main routine, 
`origamitag.c`, that determines the cosmic web morphology of dark 
matter particles, does not require VOBOZ or Qhull; these are 
required by `origamihalo.c` to group the halo particles into halos.

1. Unpack the ORIGAMI tarball.
2. Unpack and install Qhull into its own directory. See its 
README for details on installation. 
3. `cd` to `code/` and edit the QLIBS and QINC lines of the 
ORIGAMI `Makefile` to point to your Qhull installation. 
Remember  to update your `LD_LIBRARY_PATH` for Qhull to work 
properly.
4. Type `make` to compile everything, `make origami` for just the ORIGAMI code, 
or `make voboz` for just the VOBOZ tessellation code.


---
## Input File Formats
---

### Parameter File

ORIGAMI consists of several stand-alone C routines, explained 
below in more detail: VOBOZ, `origamitag`, `origamihalo`, `origamienv`, 
and `origamicen`. All input variables needed for these are read 
from a parameter file: an example, `parameters_example.txt`, is 
included. These parameters are listed below, including which of 
the above five routines use them. The parameter file  format is 
not flexible; all parameters need to be set to something, 
regardless of whether the current routine you are running needs 
them.

|Parameter | Description | Used by |
| :--- | :--- | :--- |
|`posfile` | Full location of particle position file | voboz, tag, cen
|`outdir` | Directory where output will be written | voboz, tag, halo, env, cen
|`taglabel` | Label NAME describing this run | voboz, tag, halo, env, cen
|`boxsize` | Length of simulation cube in same units as positions | voboz, tag, cen
|`np1d` | 1-D number of particles, i.e. cube root of total | tag
|`nsplit` | Split volume into `nsplit`^3 sub-volumes | tag
|`numfiles` | For Gadget input, # of particles per snapshot; 0 otherwise | voboz, tag, cen
|`buffer` | Buffer size around each sub-cube (defaults to 0.1) | voboz
|`ndiv` | Number of divisions (defaults to 2, giving 8 sub-cubes) | voboz
|`volcut` | Volume (i.e. 1/density) cut which defines a halo "core" | halo
|`npmin` | Minimum number of particles per halo | halo
|`halolabel` | Label HNAME describing these halos | halo, env, cen


### Particle Positions

ORIGAMI currently only works with simulations that have grid-like 
initial conditions (see the paper for more details). ORIGAMI 
uses the particle ID to determine the particle's initial 
position (i.e., ID = *x* + *yN* + *zN*^2, where *x*, *y*, and *z* go from 0 
to *N*-1 and *N* is the one-dimensional number of particles), so the code requires 
particles to be sorted by their particle ID. 

As of ORIGAMI version 2, the positions can be read directly from 
Gadget snapshot files, which are then sorted according to their IDs.
This option is chosen by setting `numfiles` to be greater than 0, i.e. 
if Gadget was run on 64 cores this number should be 64. In this 
case, the posfile excludes the final '.i', where *i* goes from 0 to 
`numfiles - 1`. Each array is sorted by the particle IDs, which are 
read as 8 byte long integers. The extra `readfiles_int.c` can be 
substituted for `readfiles.c` if your Gadget IDs are the standard 
4 byte integers.


If `numfiles` is set to 0, the positions must be binary files 
containing the number of particles as a 4 byte integer, followed by 
*x*, *y*, and *z* as single precision float arrays, which are read 
separately. 

The file `readfiles.c` can be modified to read positions (and 
velocities) in a different format, for both the tessellation and 
morphology identification, but particles must be in the correct 
order for ORIGAMI to work correctly. The particle velocities are 
only required for the post-processing step but must also be in this 
format and sorted by the above particle ID.

---
## Calculating the Tessellation
---

The VOBOZ documentation contains many details of how to set 
the parameters and what they mean, so here I will only 
describe in brief how to calculate the tessellation with the 
small changes that make the VOBOZ code fit with the rest of 
the ORIGAMI routines. The tessellation is not required for 
the morphology identification step (see next section) so can 
be run separately, but both are required to group the 
particles into halos.

The first step is to run `vozinit`, which creates a script 
that will run the actual tessellation on sub-cubes (or you 
can run each sub-cube in parallel) with `voz1b1`, and then 
combine the results with `voztie`. `vozinit` can be run via:
`> ./vozinit parameterfile`.

Running `vozinit` will create a script called `NAMEscr` in the 
`code/` directory; running `NAMEscr` will temporarily create 
files starting with 'part' in the current directory then 
delete them when finished, and will output the 
adjacency and volume files to `outdir/NAMEadj.dat` and 
`outdir/NAMEvol.dat`. The adjacency file gives you the 
Delaunay neighbors of each particle, and the volume file 
gives you the volume of each particle normalized to unity.


---
## Running ORIGAMI
---

### Morphology Identification

This step is done by origamitag and does not require the 
tessellation:
`> ./origamitag parameterfile`. To run this in parallel, set `nsplit`
to be greater than 1, and it will split the calculation into `nsplit`^3 
sub-cubes, using shared-memory OpenMP parallelization.

`origamitag` will calculate the ORIGAMI morphology of every 
particle and save this as a byte array in `outdir/NAMEtag.dat`. 
This morphology index is 0 for voids, 1 for walls, 2 for 
filaments, and 3 for halo particles, corresponding to the 
maximum number of orthogonal axes along which shell-crossing 
has been detected. This should run much faster than the 
tessellation calculation.

### Grouping Halos

Once you have the tessellation and morphology files, 
`origamihalo` groups the identified halo particles into 
distinct halos:
`> ./origamihalo parameterfile`. 
Note that the `outdir` parameter is both where output will be 
written and where previous `adj`, `vol`, and `tag` files must be 
contained.

Running `origamihalo` produces the file `outdir/HNAMEhns.dat` 
containing the particle IDs in each halo. This binary file 
first contains the (4 byte integer) number of halos, then the 
number of particles in the first halo, then the particle IDs 
in the first halo, etc. Refer to the IDL routine 
`origamicat.pro` for an example of how to read this file. See 
the paper for details on the effect of the volume cut 
parameter. This step should run the fastest.

### Halo Web Environment

The particle morphologies and tessellation files can be used to 
determine the morphological environment of individual halos. This 
routine, `origamienv`, was developed for [Falck, Koyama, Zhao, & Li 
2014, JCAP, 7, 58](https://ui.adsabs.harvard.edu/#abs/2014JCAP...07..058F/abstract) 
to determine the effect of the Vainshtein 
screening mechanism on halos in different web environments. This 
step requires the tessellation, morphology, and halo files 
`outdir/NAMEadj.dat`, `outdir/NAMEtag.dat`, and `outdir/HNAMEhns.dat`. 
We first calculate the number of halo neighbors of each halo; if 
a halo is connected to three or more different halos on the 
tessellation, then it is in a cluster. If not, we calculate the 
fraction of halo, filament, wall, and void particles that the 
halo is connected to and assign the halo to the morphology with 
the highest fraction of connected particles. 

Again this is run via:
`> ./origamienv parameterfile`.
This produces the binary file `outdir/HNAMEenv.dat` containing 
first the (long int) number of halos, then a byte array of the 
halo morphology indices; as with the particle morphologies, this 
is 0 for a void halo, 1 for wall, 2 for filament, and 3 for a 
halo in a cluster.

### Halo Catalog

If you are used to using a different halo-finder, it might 
be a good idea to run your favorite post-processing code 
on the ORIGAMI halo file, since different methods of 
calculating halo properties can lead to quite different 
catalogs (see [Knebe et al. 2011](https://ui.adsabs.harvard.edu/abs/2011MNRAS.415.2293K/abstract)). 
In any case, for completeness (and convenience) I have also 
included the code I used to calculate halo properties in 
the paper. First, the centers of the halos are calculated 
with `origamicen`, which averages the particle positions 
weighted by each particle's VTFE density. 

The centers are calculated via:
`> ./origamicen parameterfile`.
Note that the outdir parameter must contain the output of previous
steps. This produces the file `outdir/HNAMEcen.dat` containing the *x*, 
*y*, and *z* centers of each halo. The IDL routine `origamicat.pro` 
calculates the rest of the halo properties and writes them to 
the file `outdir/HNAMEcat.dat` (see `origamicat.pro` for more details 
on how to run that routine). Note: if you've changed `readfiles.c` 
to read the position and velocity files in a different format, you 
must update `origamicat.pro` as well.


---
## Version History
---
- 2.0: This version, 2017
    -  the main routine, `origamitag.c`, is now OpenMP parallel,
	set by the parameter `nsplit`
    -  all routines now read parameters from a file instead
	of the command line
    -  the routine `readfiles.c` has been updated so that particle
	positions can be read directly from Gadget snapshot files,
	if the parameter `numfiles` is larger than 0.
- 1.2: April 2014
    -  bugs fixed in `origamitag.c` and `origamicat.pro`
    -  addition of `origamienv.c` to determine halos' web 
	environment
    -  added new input to `origamihalo.c` and `origamicen.c` to 
	specify filename of halo catalogs, so that they can have a 
	different label than `tag`, `adj`, and `vol` files
    -  require that each halo "core" contains at least 3 
	particles in `origamihalo.c`, to reduce the number of spurious 
	halos caused by density spikes
- 1.1: Public release, October 2012
- 1.0: Initial development, 2010 - 2012
---
