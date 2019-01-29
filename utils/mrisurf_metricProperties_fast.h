// Not a header file - is included into mrisurf_metricProperties.c
//
// By having mrisurf_metricProperties_{fast,slow}.h it is easier to keep the two versions in sync


// WEIRD DISCOVERY - v_area is written twice!
//      once  by MRISMP_computeNormals 
//      later by MRISMP_computeTriangleProperties
// and I don't think it is used in between
// so the first calculation is a waste of time, except that it also writes origarea which the later one seems not to...
//
// There is a lot of overlap between these two functions that needs to be rationalized


// These should become invalid when XYZ are changed - NYI
//
static void MRISMP_computeSurfaceDimensions(MRIS_MP* mris)
{
  float xlo, ylo, zlo, xhi, yhi, zhi;

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    float const x = mris->v_x[vno];
    float const y = mris->v_y[vno];
    float const z = mris->v_z[vno];
    if (x > xhi) xhi = x;
    if (x < xlo) xlo = x;
    if (y > yhi) yhi = y;
    if (y < ylo) ylo = y;
    if (z > zhi) zhi = z;
    if (z < zlo) zlo = z;
  }
  
  mris->xlo = xlo;
  mris->xhi = xhi;
  mris->ylo = ylo;
  mris->yhi = yhi;
  mris->zlo = zlo;
  mris->zhi = zhi;
  
  mris->xctr = 0.5f * (float)((double)xlo + (double)xhi);
  mris->yctr = 0.5f * (float)((double)ylo + (double)yhi);
  mris->zctr = 0.5f * (float)((double)zlo + (double)zhi);
}


static float mrismp_TriangleArea(MRIS_MP *mris, int fno, int n)
{
  FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];

  int const n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  int const n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;

  float v0[3], v1[3];

  int const vno0 = ft->v[n0];
  int const vno1 = ft->v[n1];
  int const vno2 = ft->v[n ];
  
  v0[0] = mris->v_x[vno2] - mris->v_x[vno0];
  v0[1] = mris->v_y[vno2] - mris->v_y[vno0];
  v0[2] = mris->v_z[vno2] - mris->v_z[vno0];
  v1[0] = mris->v_x[vno1] - mris->v_x[vno2];
  v1[1] = mris->v_y[vno1] - mris->v_y[vno2];
  v1[2] = mris->v_z[vno1] - mris->v_z[vno2];
  
  float d1 = -v1[1] * v0[2] + v0[1] * v1[2];
  float d2 =  v1[0] * v0[2] - v0[0] * v1[2];
  float d3 = -v1[0] * v0[1] + v0[0] * v1[1];
  return sqrt(d1 * d1 + d2 * d2 + d3 * d3) / 2;
}

static void mrismp_setFaceNorm(MRIS_MP* mris, int fno, float nx, float ny, float nz) {
  mris->f_normSet[fno] = true;
  mris->f_norm   [fno].x = nx;
  mris->f_norm   [fno].y = ny;
  mris->f_norm   [fno].z = nz;
}

static void mrismp_setFaceAreaNormal(MRIS_MP* mris, int fno)
{
  // The code seems to not care about the length of the face normal
  // so use the cross product of any two edges
  //
  FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];
  int  const * const pv = ft->v;
  
  int  const vno0 = pv[0];
  int  const vno1 = pv[1];
  int  const vno2 = pv[2];
  
  float const x0 = mris->v_x[vno0];
  float const y0 = mris->v_y[vno0];
  float const z0 = mris->v_z[vno0];
  float const x1 = mris->v_x[vno1];
  float const y1 = mris->v_y[vno1];
  float const z1 = mris->v_z[vno1];
  float const x2 = mris->v_x[vno2];
  float const y2 = mris->v_y[vno2];
  float const z2 = mris->v_z[vno2];
  
  float v20[3], v12[3];

  // Any two sides of the triangle will give the same answer
  v20[0] = x2 - x0;
  v20[1] = y2 - y0;
  v20[2] = z2 - z0;
  v12[0] = x1 - x2;
  v12[1] = y1 - y2;
  v12[2] = z1 - z2;

  // compute cross product  
  float nx = -v12[1]*v20[2] + v20[1]*v12[2];
  float ny =  v12[0]*v20[2] - v20[0]*v12[2];
  float nz = -v12[0]*v20[1] + v20[0]*v12[1];

  mrismp_setFaceNorm(mris, fno, nx, ny, nz);
}


/*!
  \fn int mrismp_NormalFace(MRIS *mris, int fac, int n,float norm[])
  \brief Computes the normal to a triangle face. The normal will not
  have a unit length unless global variable UnitizeNormalFace=1. fac
  is the face index in mris->faces. n is the nth (0,1,2) vertex 
  
  The definition here results in the length being the sin of the angle at the vertex
  and is used to bias the sum of the face normal vectors used to compute a vertex normal vector
  
 */ 
static void mrismp_NormalFace(MRIS_MP* mris, int fno, int n, float norm[])
{
  FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];
  int  const * const pv = ft->v;
  
  int const n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  int const n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
  int const n2 = n;
  
  int  const vno0 = pv[n0];
  int  const vno1 = pv[n1];
  int  const vno2 = pv[n2];
  
  float const x0 = mris->v_x[vno0];
  float const y0 = mris->v_y[vno0];
  float const z0 = mris->v_z[vno0];
  float const x1 = mris->v_x[vno1];
  float const y1 = mris->v_y[vno1];
  float const z1 = mris->v_z[vno1];
  float const x2 = mris->v_x[vno2];
  float const y2 = mris->v_y[vno2];
  float const z2 = mris->v_z[vno2];

  float v0[3], v1[3];
  v0[0] = x2 - x0;
  v0[1] = y2 - y0;
  v0[2] = z2 - z0;
  v1[0] = x1 - x2;
  v1[1] = y1 - y2;
  v1[2] = z1 - z2;

  mrisNormalize(v0);
  mrisNormalize(v1);

  // compute cross product
  norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
  norm[1] =  v1[0]*v0[2] - v0[0]*v1[2];
  norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];

  // Note: cross product is not a unit vector even if inputs
  // are. Inputs do not need to be unit.  Until Oct 2017, this
  // function always returned a non-unitized vector. UnitizeNormalFace is a
  // global variable defined in mrisurf.h that can turn unitization on
  // and off to test its effect. Note: skull stripping uses this
  // function.
  if (UnitizeNormalFace) mrisNormalize(norm);
}


static void MRISMP_computeNormals(MRIS_MP* mris, bool check)
{
  static const double RAN = 0.001; /* one thousandth of a millimeter */

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
    // mris_print_hash(stdout, mris, "mris ", "\n");
  }

  int k;

  // For every face, 
  // if   it is ripped then mark its vertices as .border
  // else compute its norm so we don't have to compute it later
  //
  ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(shown_reproducible)
#endif
  for (k = 0; k < mris->nfaces; k++) {
    ROMP_PFLB_begin

    FACE_TOPOLOGY const * const ft = &mris->faces_topology[k];
    
    if (mris->f_ripflag[k]) {
      int n;
      for (n = 0; n < VERTICES_PER_FACE; n++) {
#ifdef HAVE_OPENMP
        #pragma omp critical
#endif
        {
          mris->v_border[ft->v[n]] = TRUE;
        }
      }
    } else {
      
      // The old code only adjusts the face norm
      // if the face is adjacent to a non-ripped vertex.
      //
      // Mimic that behavior here, for no good reason other than compatibility
      //
      int vi;
      for (vi = 0; vi < VERTICES_PER_FACE; vi++) {
          if (!mris->v_ripflag[ft->v[vi]]) break;
      }
      
      if (vi < VERTICES_PER_FACE) {
        mrismp_setFaceAreaNormal(mris, k);
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  // Build the initial pending list
  //
  int  pendingCapacity = mris->nvertices;
  int* pending         = (int*)malloc(pendingCapacity * sizeof(int));
  int* nextPending     = (int*)malloc(pendingCapacity * sizeof(int));
  int  pendingSize     = 0;

  for (k = 0; k < mris->nvertices; k++) {
    if (mris->v_ripflag[k]) continue;
    pending[pendingSize++] = k;
  }
     
  // Process the pending vertices, keeping those that need more work on the nextPending list
  // Try only a few times, because the problem might be insoluable
  //
  int trial = 0;
  for (trial = 0; (trial < 5) && (pendingSize > 0); trial++) {

    int nextPendingSize = 0;
  
    if (debugNonDeterminism) {
      fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
      // mris_print_hash(stdout, mris, "mris ", "\n");
    }

    int p;
    ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(shown_reproducible)
#endif
    for (p = 0; p < pendingSize; p++) {
      ROMP_PFLB_begin

      int const     k = pending[p];
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      //RTEX                * const v  = &mris->vertices         [k];

      // calculate the vertex area (sum of the face areas)
      // and       the average of the face normals
      //
      float snorm[3];
      snorm[0] = snorm[1] = snorm[2] = 0;

      float area = 0;

      int count = 0;

      int n;
      for (n = 0; n < vt->num; n++) {
        int const fno = vt->f[n];
        //FACE* face = &mris->faces[fno];
        if (mris->f_ripflag[fno]) continue;
        
        count++;
        
        float norm[3];
        mrismp_NormalFace(mris, vt->f[n], (int)vt->n[n], norm);
            // The normal is NOT unit length OR area length
            // Instead it's length is the sin of the angle of the vertex
            // The vertex normal is biased towards being perpendicular to 90degree contributors...

        snorm[0] += norm[0];
        snorm[1] += norm[1];
        snorm[2] += norm[2];

        area += mrismp_TriangleArea(mris, vt->f[n], (int)vt->n[n]);
      }
      
      if (!count || mrisNormalize(snorm) > 0.0) {       // Success?

        if (fix_vertex_area)
          area = area / 3.0;                            // Since each face is added to three vertices...
        else
          area = area / 2.0;

        mris->v_area[k] = area;
        
        if (mris->v_origarea[k] < 0)                   // has never been set
            mris->v_origarea[k] = area;

        mris->v_nx[k] = snorm[0];
        mris->v_ny[k] = snorm[1];
        mris->v_nz[k] = snorm[2];
        
        ROMP_PFLB_continue;
      }
      
        
#ifdef HAVE_OPENMP
      #pragma omp critical                              // Retry after various adjustments below
#endif
      {
        nextPending[nextPendingSize++] = k;
      }
      
      ROMP_PFLB_end
    }
    ROMP_PF_end

    // The test I was using ALWAYS took this path!
    //
    if (nextPendingSize == 0) {
        pendingSize = 0;
        break;
    }
    
    // THIS HAS BEEN OBSERVED
    
    // Sort the nextPending list because the above appends are not in order
    //
    qsort(nextPending, nextPendingSize, sizeof(int), int_compare);
    
    // Randomly move nextPending vertices and their neighbors
    //
    // This can not be done easily in parallel because they share faces
    // If this is a performance problem, can be fixed, but I doubt it will be
    //
    for (p = 0; p < nextPendingSize; p++) {
      int     const k = nextPending[p];
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      //RTEX                * const v  = &mris->vertices         [k];

      // Warn
      //      
      if ((mris->status == MRIS_SPHERICAL_PATCH      || 
           mris->status == MRIS_PARAMETERIZED_SPHERE ||
           mris->status == MRIS_SPHERE) &&
          DIAG_VERBOSE_ON) {
        fprintf(stderr, "vertex %d: degenerate normal\n", k);
      }
    
      // randomly move x and y
      // When written, assumed being done in parallel, hence fnv_hash giving reproducible movement
      //
      int random_counter = 0;

      mris->v_x[k] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
      mris->v_y[k] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
      
      if (mris->status == MRIS_PLANE || mris->status == MRIS_CUT) {
        // nothing
      } else {
        mris->v_z[k] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
        
        // I don't know why this was not done for the MRIS_PLANE and _CUT
        //
        int n;
        for (n = 0; n < vt->vnum; n++) { /* if (!mris->faces[v->f[n]].ripflag) */
          int const vnon = vt->v[n];
          //RTEX* vn = &mris->vertices[vnon];
          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
            printf("   k=%5d %d nbr = %5d / %d\n", k, n, vt->v[n], vt->vnum);
          mris->v_x[vnon] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
          mris->v_y[vnon] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
          mris->v_z[vnon] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
        }
        
        // Recompute the face norms for the affected faces
        // Note: this will recompute some, but I suspect the number is too small to be relevant
        //
        for (n = 0; n < vt->num; n++) {
          int const fno = vt->f[n];
          if (mris->f_ripflag[fno]) continue;
          mrismp_setFaceAreaNormal(mris, fno);
        }
      }
    }
    
    // Try the moved ones again - although not their moved neighbors, weirdly enough
    //
    int* tempPending = pending; pending = nextPending; nextPending = tempPending;
    pendingSize = nextPendingSize;

  } // trials

#if 0  
  if (pendingSize > 0) {
    fprintf(stderr, "%s:%d MRISMP_computeNormals could not do all vertices after %d attempts, %d remain\n",
      __FILE__, __LINE__, trial, pendingSize);
  }
#endif

  freeAndNULL(nextPending);
  freeAndNULL(pending);
}


#define COMPILING_MRIS_MP


#define FUNCTION_NAME MRISMP_computeVertexDistancesWkr
#define INPUT_X v_x
#define INPUT_Y v_y
#define INPUT_Z v_z
#define OUTPUT_DIST dist
#define OUTPUT_MAKER MRISMP_makeDist
#include "mrisComputeVertexDistancesWkr_extracted.h"

static void MRISMP_computeVertexDistances(MRIS_MP *mris) {
  MRISMP_computeVertexDistancesWkr(mris, mris->nsize, false);
  mris->dist_nsize = mris->nsize;
}

#define FUNCTION_NAME MRISMP_computeTriangleProperties
#include "MRIScomputeTriangleProperties_extracted.h"


/*-------------------------------------------------------------
  MRIScomputeAvgInterVertexDist() - computes the average and stddev of
  the distance between neighboring vertices. If StdDev is NULL,
  it is ignored. Requires that mrisComputeVertexDistances()
  have been run in order to compute vertex->dist[n].
  -------------------------------------------------------------*/
static void MRISMP_computeAvgInterVertexDist(MRIS_MP *mris, double *StdDev)
{
  bool const showHashs = false || debugNonDeterminism;

  if (showHashs) {
    fprintf(stdout, "%s:%d MRISMP_computeAvgInterVertexDist starting ",__FILE__,__LINE__);
    //mris_print_hash(stdout, mris, "mris ", "\n");
  }
  
  double Sum = 0, Sum2 = 0;

  double N = 0.0;

  #define ROMP_VARIABLE       vno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  Sum
  #define ROMP_SUMREDUCTION1  Sum2
  #define ROMP_SUMREDUCTION2  N
    
  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define Sum  ROMP_PARTIALSUM(0)
    #define Sum2 ROMP_PARTIALSUM(1)
    #define N    ROMP_PARTIALSUM(2)

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];

    if (mris->v_ripflag[vno]) {
      continue;
    }
    int const vnum = vt->vnum;
    int m;
    for (m = 0; m < vnum; m++) {
      int const vno2 = vt->v[m];
      
      if (mris->v_ripflag[vno2]) {
        continue;
      }
      float* dist = mris->v_dist[vno];
      double d = dist[m];
      Sum  += d;
      Sum2 += (d * d);
      N    += 1;
    }

    #undef Sum
    #undef Sum2
    #undef N

  #include "romp_for_end.h"
  
  // NOTE - This is a poor algorithm for computing the std dev because of how the floating point errors accumulate
  // but it seems to work for us because the double has enough accuracy to sum the few hundred thousand small but not
  // too small floats that we have
  double Avg = Sum / N;
  if (StdDev != NULL) {
    *StdDev = sqrt(N * (Sum2 / N - Avg * Avg) / (N - 1));
  }

  // printf("\n\nN = %ld, Sum = %g, Sum2 = %g, Avg=%g, Std = %g\n\n",
  // N,Sum,Sum2,Avg,*StdDev);

  mris->avg_vertex_dist = Avg;
}


static void mrismp_OrientEllipsoid(MRIS_MP *mris)
{
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f;

  double total_area = 0.0, neg_area = 0.0, neg_orig_area = 0.0;

  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces

  #define ROMP_SUMREDUCTION0  total_area
  #define ROMP_SUMREDUCTION1  neg_area
  #define ROMP_SUMREDUCTION2  neg_orig_area

  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible

#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define total_area      ROMP_PARTIALSUM(0)
    #define neg_area        ROMP_PARTIALSUM(1)
    #define neg_orig_area   ROMP_PARTIALSUM(2)
    
    FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];
    
    if (mris->f_ripflag[fno]) {
      ROMP_PFLB_continue;
    }

    FloatXYZ const * const fNorm = &mris->f_norm[fno];
    float const nx = fNorm->x;
    float const ny = fNorm->y;
    float const nz = fNorm->z;

    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
    */
    int const vno0 = ft->v[0];
    int const vno1 = ft->v[1];
    int const vno2 = ft->v[2];
    
    float x0 = mris->v_x[vno0];
    float y0 = mris->v_y[vno0];
    float z0 = mris->v_z[vno0];
    float x1 = mris->v_x[vno1];
    float y1 = mris->v_y[vno1];
    float z1 = mris->v_z[vno1];
    float x2 = mris->v_x[vno2];
    float y2 = mris->v_y[vno2];
    float z2 = mris->v_z[vno2];

    float   const xc = (x0 + x1 + x2) /* / 3 */;   // These divides by three are a waste of time
    float   const yc = (y0 + y1 + y2) /* / 3 */;   // since we only use the magnitude of the dot product
    float   const zc = (z0 + z1 + z2) /* / 3 */;

    float   const dot = xc * nx + yc * ny + zc * nz;
    
    if (dot < 0.0f) /* not in same direction, area < 0 and reverse n */
    {
      mris->f_area[fno]   *= -1.0f;

      mris->f_norm[fno].x *= -1.0f;
      mris->f_norm[fno].y *= -1.0f;
      mris->f_norm[fno].z *= -1.0f;

      angles_per_triangle_t* angle = &mris->f_angle[fno];
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        (*angle)[ano]      *= -1.0f;
      }
    }
    
    if (mris->f_area[fno] >= 0.0f) {
      total_area += mris->f_area[fno];
    } else {
      neg_area      += -mris->f_area[fno];
      neg_orig_area +=  mris->f_norm_orig_area[fno];
    }

    #undef total_area
    #undef neg_area
    #undef neg_orig_area

  #include "romp_for_end.h"

  mris->total_area    = total_area;
  mris->neg_area      = neg_area;
  mris->neg_orig_area = neg_orig_area;
}



void MRISMP_updateEllipsoidSurface(MRIS_MP* mris)
{
  if (mris->status != MRIS_UNORIENTED_SPHERE) {
    mrismp_OrientEllipsoid(mris); /* orient the normals and angles */
  }
}



static void mrismp_OrientPlane(MRIS_MP* mris)
{
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f;

  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {

    if (mris->f_ripflag[fno]) continue;

    // give the area an orientation: if the unit normal is pointing
    // downwards in the plane then the area should be negative.
    //
    FloatXYZ* const fNorm = &mris->f_norm[fno];
    float const nx = fNorm->x;
    float const ny = fNorm->y;
    float const nz = fNorm->z;
    
    if (nz < 0.0f) {
      /* not in same direction, area < 0 and reverse n */
      mris->f_area[fno]   *= -1.0f;
      
      fNorm->x = -nx;
      fNorm->y = -ny;
      fNorm->z = -nz;

      angles_per_triangle_t* angle = &mris->f_angle[fno];
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        (*angle)[ano]      *= -1.0f;
      }
    }
    
    if (mris->f_area[fno] >= 0.0f) {
      mris->total_area += mris->f_area[fno];
    } else {
      mris->neg_area      += -mris->f_area[fno];
      mris->neg_orig_area +=  mris->f_norm_orig_area[fno];
    }
    
  }

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    
    if (mris->v_ripflag[vno]) continue;

    if (mris->v_nz[vno] < 0) {
      mris->v_nz [vno] *= -1.0f;
      mris->v_neg[vno] = 1;
    } else {
      mris->v_neg[vno] = 0;
    }
    float v_area = 0.0f;
    int fn;
    for (fn = 0; fn < vt->num; fn++) {
      int fno = vt->f[fn];
      v_area += mris->f_area[fno];
    }
    if (fix_vertex_area) {
      v_area /= 3.0;
    } else {
      v_area /= 2.0;
    }
    mris->v_area[vno] = v_area;
  }
}


int mrismp_OrientSurface(MRIS_MP *mris)
{
  switch (mris->status) {
    case MRIS_RIGID_BODY:
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
    case MRIS_ELLIPSOID:
    case MRIS_SPHERICAL_PATCH:
      MRISMP_updateEllipsoidSurface(mris);
      break;
    case MRIS_PLANE:
      mrismp_OrientPlane(mris);
      break;
    default:
      /*    MRISupdateSurface(mris) ;*/
      break;
  }
  return (NO_ERROR);
}


void MRISMP_computeMetricProperties(MRIS_MP* mris) {
  // fprintf(stdout,"%s:%d %s\n",__FILE__,__LINE__,__MYFUNCTION__);

  MRISMP_computeNormals(mris, false);           // changes XYZ
  MRISMP_computeSurfaceDimensions(mris);
  MRISMP_computeVertexDistances(mris);
  MRISMP_computeTriangleProperties(mris);       // compute face areas and normals
  
  mris->avg_vertex_area = mris->total_area / mris->nvertices;

  MRISMP_computeAvgInterVertexDist(mris, &mris->std_vertex_dist);

  mrismp_OrientSurface(mris);
  // See also MRISrescaleMetricProperties()
  
  if (mris->status == MRIS_PARAMETERIZED_SPHERE || mris->status == MRIS_RIGID_BODY || mris->status == MRIS_SPHERE) {
    mris->total_area = M_PI * mris->radius * mris->radius * 4.0;
  }
}


float mrismp_SampleMinimizationEnergy(
    MRIS_MP *mris, int const vno, INTEGRATION_PARMS *parms, float cx, float cy, float cz)
{
  project_point_onto_sphere(cx, cy, cz, mris->radius, &cx, &cy, &cz);
  float const xw = mris->v_whitex[vno], yw = mris->v_whitey[vno], zw = mris->v_whitez[vno];

  float xp, yp, zp;
  MRISMP_sampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_VERTICES, &xp, &yp, &zp);

  float dx = xp - xw;
  float dy = yp - yw;
  float dz = zp - zw;
  float thick_sq = dx * dx + dy * dy + dz * dz;

  return (thick_sq);
}


int MRISMP_sampleFaceCoordsCanonical(
    MHT *mht, MRIS_MP *mris, float x, float y, float z, int which, float *px, float *py, float *pz)
{
  double norm = sqrt(x * x + y * y + z * z);
  if (!FEQUAL(norm, mris->radius))  // project point onto sphere
  {
    DiagBreak();
    project_point_onto_sphere(x, y, z, mris->radius, &x, &y, &z);
  }

  double lambda[3], fdist;
  int   fno;
  MHTfindClosestFaceGeneric2(mht, MRISBaseConstCtr(mris,NULL), x, y, z, 8, 8, 1, &fno, &fdist);
  if (fno < 0) {
    DiagBreak();
    MHTfindClosestFaceGeneric2(mht, MRISBaseConstCtr(mris,NULL), x, y, z, 1000, -1, -1, &fno, &fdist);
    lambda[0] = lambda[1] = lambda[2] = 1.0 / 3.0;
  }
  else {
    face_barycentric_coords2(MRISBaseConstCtr(mris,NULL), fno, CANONICAL_VERTICES, x, y, z, &lambda[0], &lambda[1], &lambda[2]);
  }

  FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];

  *px = *py = *pz = 0;
  int n;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    float xv, yv, zv;
    MRISBase_getWhichXYZ(MRISBaseConstCtr(mris,NULL), ft->v[n], which, &xv, &yv, &zv);
    *px += lambda[n] * xv;
    *py += lambda[n] * yv;
    *pz += lambda[n] * zv;
  }

  return NO_ERROR;
}
#undef COMPILING_MRIS_MP
