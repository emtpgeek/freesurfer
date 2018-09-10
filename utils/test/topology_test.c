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
  
  int* vlist = (int*)malloc(nvertices*sizeof(int));
  
  int nlinks;
  for (nlinks = 0; nlinks < 4; nlinks++) {
    int acquiredMarked = MRIS_acquireTemp(src, MRIS_TempAssigned_Vertex_marked);                           
    MRISclearMarks(src);
    
    int const vno1 = 0x0;
    int const vlistSize = MRISfindNeighborsAtVertex(src, acquiredMarked, vno1, nlinks, vlist);
    
    printf("Returned vlistSize:%d at nlinks %d\n", vlistSize, nlinks);

    // Check the vno1 marked is correct
    //
    if (src->vertices[vno1].marked != -1) {
        printf("src->vertices[vno1:%p].marked:%d should be -1\n", 
            (void*)(long)vno1, src->vertices[vno1].marked);
            fails++;
    }
    src->vertices[vno1].marked = 0;

    // Check the other marked are correct
    //
    int expected_vlistSize = 0;
    int vno2;    
    for (vno2 = 0; vno2 < nvertices; vno2++) {
      int marked = bitCount(vno1^vno2);
      if (marked > nlinks) marked = 0;
      if (marked != src->vertices[vno2].marked) {
        printf("Marked %p with the wrong distance %d\n", (void*)(long)vno2, src->vertices[vno2].marked);
        fails++;
      }
      if (marked > 0) expected_vlistSize++;
    }

    // Check the count is right
    //
    if (expected_vlistSize != vlistSize) {
      printf("Returned the wrong number of elements %d, instead of %d\n", 
        expected_vlistSize, vlistSize);
      fails++;
    }
    
    // Check that the list is in the correct order
    //
    int furtherestSoFar = 1;    // should not find a 0
    int distancesFound  = 0;
    for (i = 0; i < vlistSize; i++) {
      int vno2 = vlist[i];
      int distance = src->vertices[vno2].marked;
      if (distance < furtherestSoFar) {
        printf("Distances in wrong order %d before %d\n", furtherestSoFar, distance);
        fails++;
      }
      int distancesFoundMask = 1 << distance;
      if (!(distancesFound & distancesFoundMask)) {
        printf("Listed some at distance %d\n", distance);
        distancesFound |= distancesFoundMask;
      } 
    }
    
    MRIS_releaseTemp(src, MRIS_TempAssigned_Vertex_marked, acquiredMarked);
  }
  
  freeAndNULL(vlist);
  
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
