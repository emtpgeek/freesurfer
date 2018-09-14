#include "stdio.h"

#include "mrisurf.h"
#include "../mrisurf_topology.h"

int main() {

  int fails = 0;
  
  int const max_vertices = 4, max_faces = 2, 
               nvertices = 4,    nfaces = 0;

  MRIS* src = MRISoverAlloc(max_vertices, max_faces, nvertices, nfaces);
  
  mrisAddEdge(src, 0, 1);
  mrisAddEdge(src, 0, 2);
  mrisAddEdge(src, 0, 3);
  mrisAddEdge(src, 1, 2);
  mrisAddEdge(src, 2, 3);
  mrisAddEdge(src, 3, 1);

  mrisAddFace(src, 0,1,2);
  mrisAddFace(src, 1,2,3);

  const char* fnm = "./test_TriangleFile_readWrite.tmp";
  printf("MRISwrite %s returned %d\n",
    fnm, MRISwrite(src, fnm));
    
  MRIS* dst = MRISread(fnm) ;
  printf("MRISread %s returned mris with dims nvertices:%d nfaces:%d\n",
    fnm, dst->nvertices, dst->nfaces);

  if (src->nvertices != dst->nvertices) {
    fails++;
    printf("FAIL src->nvertices:%d != dst->nvertices:%d\n",
      dst->nvertices, dst->nfaces);
  }
  
  if (src->nfaces != dst->nfaces) {
    fails++;
    printf("FAIL src->nfaces:%d != dst->nfaces:%d\n",
      dst->nfaces, dst->nfaces);
  }
  
  return fails ? 1 : 0;
}
