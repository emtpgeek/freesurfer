#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
#define COMPILING_MRISURF_METRIC_PROPERTIES_FRIEND
/*
 * @file utilities operating on Original
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2018 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_project.h"

/* project onto the sphere of radius DEFAULT_RADIUS */
void mrisSphericalProjectXYZ(float xs, float ys, float zs, float *xd, float *yd, float *zd)
{
  double dist, lambda;

  dist = sqrt(SQR(xs) + SQR(ys) + SQR(zs));
  lambda = DEFAULT_RADIUS / dist;

  /* making sure things are stable : double projection */
  *xd = xs * lambda;
  *yd = ys * lambda;
  *zd = zs * lambda;

  xs = *xd;
  ys = *yd;
  zs = *zd;
  dist = sqrt(SQR(xs) + SQR(ys) + SQR(zs));
  lambda = DEFAULT_RADIUS / dist;

  *xd = xs * lambda;
  *yd = ys * lambda;
  *zd = zs * lambda;
}


void mrisSphericalProjection(MRIS *mris)
{
  int n;
  VERTEX *v;
  //    fprintf(stderr,"spherical projection\n");
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }

    /*
      if(n == 88 )
      fprintf(stderr,"bf sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 89 )
      fprintf(stderr,"bf sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 209 )
      fprintf(stderr,"bf sp: nvertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
    */

    mrisSphericalProjectXYZ(v->x, v->y, v->z, &v->x, &v->y, &v->z);
    v->cx = v->x;
    v->cy = v->y;
    v->cz = v->z;

    /*
      if(n == 88 )
      fprintf(stderr,"af sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 89 )
      fprintf(stderr,"af sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 209 )
      fprintf(stderr,"af sp: nvertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);

      mrisSphericalProjectXYZ(v->x,v->y,v->z,&v->x,&v->y,&v->z);

      if(n == 88 )
      fprintf(stderr,"af 2 sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 89 )
      fprintf(stderr,"af 2 sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 209 )
      fprintf(stderr,"af 2 sp: nvertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
    */
  }
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Perform a projection onto an sphere moving each
  point on the cortical surface to the closest spherical
  coordinate.
  ------------------------------------------------------*/
  
#define FUNCTION_NAME MRISprojectOntoSphereWkr
#include "mrisurf_project_projectOntoSphereWkr.h"

#define COMPILING_MRIS_MP
#define FUNCTION_NAME MRISMP_projectOntoSphereWkr
#include "mrisurf_project_projectOntoSphereWkr.h"
#undef  COMPILING_MRIS_MP


MRIS* MRISprojectOntoSphere(MRIS* mris_src, MRIS* mris_dst, double r)
{
  if (!mris_dst) {
    mris_dst = MRISclone(mris_src);
    mris_dst->status = mris_src->status;    // added this, think it right - Bevin
  }

  if ((mris_dst->status != MRIS_SPHERE) && (mris_dst->status != MRIS_PARAMETERIZED_SPHERE)) {
    MRIScenter(mris_dst, mris_dst);
  }

  MRISfreeDistsButNotOrig(mris_dst);

  MRISprojectOntoSphereWkr(mris_dst, r);
  
  return (mris_dst);
}


void MRISMP_projectOntoSphere(MRIS_MP* mris, double r) 
{
  if ((mris->status != MRIS_SPHERE) && (mris->status != MRIS_PARAMETERIZED_SPHERE)) {
    cheapAssert(false);
  }

  MRISMP_projectOntoSphereWkr(mris, r);
}


void mrisAssignFaces(MRIS *mris, MHT *mht, int which_vertices)
{
  int vno;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin 
    
    VERTEX *v = &mris->vertices[vno];
    if (v->ripflag) continue;

    if (vno == Gdiag_no) DiagBreak();

    project_point_onto_sphere(v->x, v->y, v->z, mris->radius, &v->x, &v->y, &v->z);

    int fno;
    FACE *face;
    double fdist;

    MHTfindClosestFaceGeneric(mht, mris, v->x, v->y, v->z, 8, 8, 1, &face, &fno, &fdist);
    if (fno < 0) MHTfindClosestFaceGeneric(mht, mris, v->x, v->y, v->z, 1000, -1, -1, &face, &fno, &fdist);

    v->fno = fno;
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
}
