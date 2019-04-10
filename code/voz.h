//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

#define MAXVERVER 2000
#define NGUARD 42 /*Actually, the number of SPACES between guard points
		    in each dim */

//#define BF REALmax /* needed by jozov for Qhull ?*/

typedef struct Partadj {
  int nadj;
  int *adj;
} PARTADJ;
