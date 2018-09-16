#include "stdio.h"

#include "mrisurf.h"
#include "../mrisurf_topology.h"

static const char* fnm_base = "./test_TriangleFile_readWrite.tmp";

static const char* extensions[] = { 
    "ANNOT",
    "Any_other_means_MRIS_BINARY_QUADRANGLE_FILE",
    "GEO",
    "ICO",      // "TRI" is equivalent
    "VTK",
    "STL",
    "GII",
    "MGH",
    "ASC",
    NULL};

static char fnm[1024];

int main() {
  
  int fails = 0;

  int const max_vertices = 4, max_faces = 2, 
               nvertices = 4,    nfaces = 0;

  MRIS* src = MRISoverAlloc(max_vertices, max_faces, nvertices, nfaces);

  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    // The STL format uses location to distinquish vertices
    VERTEX* v = &src->vertices[vno];
    v->x = vno &  1;
    v->y = vno &  2;
    v->z = vno & ~3;
  }

  mrisAddEdge(src, 0, 1);
  mrisAddEdge(src, 0, 2);
  mrisAddEdge(src, 0, 3);
  mrisAddEdge(src, 1, 2);
  mrisAddEdge(src, 2, 3);
  mrisAddEdge(src, 3, 1);

  mrisAddFace(src, 0,1,2);
  mrisAddFace(src, 1,2,3);

  const char* const * pext;
  for (pext = &extensions[0]; *pext; pext++) {
    const char* const ext = *pext;

    sprintf(fnm, "%s.%s", fnm_base, ext);

    printf("Trying %s\n", fnm);
      
    printf("MRISwrite %s returned %d\n",
      fnm, MRISwrite(src, fnm));

    if (ext[0] == 'M') {
      printf("MRISread does not support %s files\n",
        ext);
      continue;
    }
    
    MRIS* dst = MRISread(fnm) ;
    if (!dst) {
      fails++;
      printf("FAIL could not read %s\n", fnm);
      continue;
    }
    
    printf("MRISread %s returned mris with dims nvertices:%d nfaces:%d\n",
      fnm, dst->nvertices, dst->nfaces);

    if (src->nvertices != dst->nvertices) {
      fails++;
      printf("FAIL src->nvertices:%d != dst->nvertices:%d\n",
        dst->nvertices, dst->nvertices);
    }

    if (src->nfaces != dst->nfaces) {
      fails++;
      printf("FAIL src->nfaces:%d != dst->nfaces:%d\n",
        src->nfaces, dst->nfaces);
    }
    
    MRISfree(&dst);
  }
  
  MRISfree(&src);
  
  return fails ? 1 : 0;
}
