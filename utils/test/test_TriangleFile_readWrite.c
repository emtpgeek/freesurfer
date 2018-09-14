#include "stdio.h"

#include "mrisurf.h"

int main() {

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

  const char* fnm = "./test_TriangleFile_readWrite.tmp";
  printf("MRISwrite %s returned %d\n",
    MRISwrite(src, fnm));
    
  MRIS* dst = MRISread(fnm) ;
  printf("MRISread %s returned mris with dims nvertices:%d nfaces:%d\n",
    dst->nvertices, dst->nfaces);

  return 0;
}
