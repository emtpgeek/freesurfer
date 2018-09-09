/**
 * @file  topology_test.cpp
 * @brief mrisurf_topology tests
 *
 */
/*
 * Original Author: B. R. Brett
 *
 * Copyright © 2018 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <iostream>
#include <math.h>

extern "C"
{
#include "mrisurf.h"
#include "utils.h"
const char *Progname = "topology_test";
}

using namespace std;

int main(int argc, char *argv[])
{

  int const max_vertices = 4, max_faces = 4, nvertices = 4, nfaces = 4;

  MRIS mrisObject;
  MRIS* mris = &mrisObject;
  
  MRISctr(mris, max_vertices, max_faces, nvertices, nfaces);

  mrisAddEdge(mris, 0, 1);
  mrisAddEdge(mris, 0, 2);
  mrisAddEdge(mris, 0, 3);
  mrisAddEdge(mris, 1, 2);
  mrisAddEdge(mris, 2, 3);
  mrisAddEdge(mris, 3, 1);

  size_t     dst_nvertices;
  int const* dst_mapToDstVno;
  size_t     dst_pnfaces;
  int const* dst_mapToDstFno;
  
  MRIScreateSimilarTopologyMapsForNonripped(
    mris,
    &dst_nvertices,
    &dst_mapToDstVno,
    &dst_nfaces,
    &dst_mapToDstFno);
  
  MRIS* dst = 
    MRIScreateWithSimilarTopologyAsSubset(
      mris,
      dst_nvertices,
      dst_mapToDstVno,
      dst_nfaces,
      dst_mapToDstFno);
    
  MRISFree(dst);
  
  MRISdtr(mris);
  
  if (fails) return 77;
  return 0;
}
