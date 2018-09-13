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

int test1()
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

  // Done
  //  
  MRISfree(&dst);
  MRISdtr(src);

  return fails;
}

int bitCount(int x) {
  int count = 0;
  int i = 1;
  while (x) {
    if (x&i) { x ^= i; count++; }
    i *= 2;
  }
  return count;
}

int test2()
{
  if (bitCount(5) != 2) {
    printf("bitCount(5) != 2\n");
    return 1;
  }
  
  int fails = 0;

  // Make the src
  //
  int const verticesLog2 = 8;
  int const max_vertices = 1<<verticesLog2, max_faces = 0, 
               nvertices = max_vertices,    nfaces    = 0;

  MRIS mrisObject;
  MRIS* src = &mrisObject;
  
  MRISctr(src, max_vertices, max_faces, nvertices, nfaces);

  int i,j;
  for (i = 0; i < nvertices; i++) 
  for (j = 0; j < i; j++)
    if (bitCount(i^j) == 1) 
      mrisAddEdge(src, i, j);
  
  MRISfindNeighborsAtVertex(src, vno1);
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno1];
  
  int nlinks;
  for (nlinks = 1; nlinks < 4; nlinks++) {
    size_t b,e;
    
    MRISgetNeighborsBeginEnd(
        mris, 
        vno1, 
        nlinks, 
        nlinks, 
        &b,
        &e);

    size_t i;
    for (i = b; i < e; i++) {
      int vno2  = vt->v[i];
      mris->vertices[vno2]->marked = links;
    }
  }

  for (i = 0; i < nvertices; i++) {
    int expected = (bitCount(i^j);
    if (expected > 3) expected = 0;
    if (mris->vertices[vno2]->marked != expected) {
      fails++;
      if (fails < 10) printf("mris->vertices[vno2:%d]->marked:%d != expected:%d\n", vno2, mris->vertices[vno2]->marked, expected);
    }
  }
  
  // Done
  //  
  MRISdtr(src);

  return fails;
}


int main(int argc, char *argv[])
{
    return 
      //test1() ? 1 : 
        test2() ? 1 :
        0;
}
