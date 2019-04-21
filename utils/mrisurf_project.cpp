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
#include "mrisurf_mp.h"
#include "mrisurf_MRISBase.h"

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


void mrisAssignFaces(MRISBase mrisBase, MHT *mht, int which_vertices)
{
  MRISBaseConst mris = MRISBaseToConst(mrisBase);
  
  int const nvertices = MRISBase_nvertices(mris);
  int vno;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < nvertices; vno++) {
    ROMP_PFLB_begin 
    
    if (MRISBase_v_ripflag(mris, vno)) continue;
    
    if (vno == Gdiag_no) DiagBreak();

    float x,y,z;
    MRISBase_v_xyz(mris,vno, &x,&y,&z);
    project_point_onto_sphere(x, y, z, MRISBase_radius(mris), &x, &y, &z);
    MRISBase_set_v_xyz(mrisBase,vno, x,y,z);

    int fno;
    double fdist;
    MHTfindClosestFaceGeneric2(mht, mris, x, y, z, 8, 8, 1, &fno, &fdist);
    if (fno < 0) MHTfindClosestFaceGeneric2(mht, mris, x, y, z, 1000, -1, -1, &fno, &fdist);

    MRISBase_set_v_assigned_fno(mrisBase,vno, fno);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
}
