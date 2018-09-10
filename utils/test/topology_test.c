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

#include "mrisurf.h"
#include "utils.h"
#include <stdio.h>

#define CHECK(COND) if (!(COND)) { printf("%s failed at line %d\n", #COND, __LINE__); fails++; }

const char *Progname = "topology_test";

int main(int argc, char *argv[])
{
  int fails = 0;

  // Make the src
  //
  int const max_vertices = 4, max_faces = 2, 
               nvertices = 4,    nfaces = 2;

  MRIS mrisObject;
  MRIS* src = &mrisObject;
  
  MRISctr(src, max_vertices, max_faces, nvertices, nfaces);
  CHECK(src->nvertices == 4);
  CHECK(src->nfaces == 2);
   
  mrisAddEdge(src, 0, 1);                   // will get removed
  mrisAddEdge(src, 0, 2);
  mrisAddEdge(src, 0, 3);
  mrisAddEdge(src, 1, 2);
  mrisAddEdge(src, 2, 3);
  mrisAddEdge(src, 3, 1);

  mrisAttachFaceToVertices(src, 0, 0,1,2);  // will get removed
  mrisAttachFaceToVertices(src, 1, 1,2,3);  // will survive

  CHECK(src->nvertices == 4);
  CHECK(src->nfaces == 2);
   
  // Copy the src
  //
  size_t     dst_nvertices;
  int const* dst_mapToDstVno;
  size_t     dst_nfaces;
  int const* dst_mapToDstFno;

  MRIS* dst;

  MRIS_HASH srcHash, dstHash;
  
  MRIScreateSimilarTopologyMapsForNonripped(
    src,
    &dst_nvertices,
    &dst_mapToDstVno,
    &dst_nfaces,
    &dst_mapToDstFno);
  CHECK(dst_nvertices == 4);
  CHECK(dst_nfaces == 2);
  
  dst =
    MRIScreateWithSimilarTopologyAsSubset(
      src,
      dst_nvertices,
      dst_mapToDstVno,
      dst_nfaces,
      dst_mapToDstFno);
  freeAndNULL(dst_mapToDstVno);
  freeAndNULL(dst_mapToDstFno);
  
  CHECK(dst->nvertices == 4);
  CHECK(dst->nfaces == 2);
  
  mris_print_diff(stdout, src, dst);
  mris_hash_init(&srcHash, src);
  mris_hash_init(&dstHash, dst);
  CHECK(srcHash.hash == dstHash.hash);

  MRISfree(&dst);

  // Subset the src
  //
  src->vertices[0].ripflag = 1;
  src->faces   [0].ripflag = 1;
  
  MRIScreateSimilarTopologyMapsForNonripped(
    src,
    &dst_nvertices,
    &dst_mapToDstVno,
    &dst_nfaces,
    &dst_mapToDstFno);
  CHECK(dst_nvertices == 3);
  CHECK(dst_nfaces == 1);
      
  dst =
    MRIScreateWithSimilarTopologyAsSubset(
      src,
      dst_nvertices,
      dst_mapToDstVno,
      dst_nfaces,
      dst_mapToDstFno);
  freeAndNULL(dst_mapToDstVno);
  freeAndNULL(dst_mapToDstFno);

  CHECK(dst->nvertices == 3);
  CHECK(dst->nfaces == 1);

  // TODO verify the shape
  //

  // TODO neighbourhood check
  //

  // Done
  //  
  MRISfree(&dst);
  MRISdtr(src);
  
  return fails;
}
