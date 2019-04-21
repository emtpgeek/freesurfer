#pragma once
/*
 * @file  dense representation of data needed to compute the SSE
 *
 */
/*
 * surfaces Author: Bevin Brett
 *
 * $ © copyright-2018,2019 The General Hospital Corporation (Boston, MA) "MGH"
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
// The following format provides a abstract representation of just the data needed to compute the SSE
// so that code that must apply to both MRIS_MP and to MRIS does not need to be duplicated
//
#include "mrisurf_mp.h"

static void MRISBase_getWhichXYZ(MRISBaseConst mrisBase, int vno, int which, float* x, float* y, float* z);

int face_barycentric_coords2(MRISBaseConst mrisBase, int fno, int which_vertices,
                            double cx, double cy, double cz, double *pl1, double *pl2, double *pl3) ;

static double MRISBase_radius           (MRISBaseConst mrisBase);
static int    MRISBase_nvertices        (MRISBaseConst mrisBase);
static int    MRISBase_nfaces           (MRISBaseConst mrisBase);

static bool MRISBase_v_ripflag          (MRISBaseConst mrisBase, int vno);
static void MRISBase_v_xyz              (MRISBaseConst mrisBase, int vno, float* x, float* y, float* z);
static void MRISBase_set_v_xyz          (MRISBase      mrisBase, int vno, float  x, float  y, float  z);
static void MRISBase_set_v_assigned_fno (MRISBase      mrisBase, int vno, int to);


// implementations
//
static void MRISBase_getWhichXYZ(MRISBaseConst mrisBase, int vno, int which, float* x, float* y, float* z) {
    cheapAssert(!mrisBase.mris_mp);
    MRISvertexCoord2XYZ_float(&mrisBase.mris->vertices[vno], which, x, y, z);
}

static double MRISBase_radius   (MRISBaseConst mrisBase) {
  return mrisBase.mris_mp ? mrisBase.mris_mp->radius : mrisBase.mris->radius;
}

static int  MRISBase_nvertices  (MRISBaseConst mrisBase) {
  return mrisBase.mris_mp ? mrisBase.mris_mp->nvertices : mrisBase.mris->nvertices;
}

static int  MRISBase_nfaces     (MRISBaseConst mrisBase) {
  return mrisBase.mris_mp ? mrisBase.mris_mp->nfaces : mrisBase.mris->nfaces;
}

static bool MRISBase_v_ripflag         (MRISBaseConst mrisBase, int vno) {
  return mrisBase.mris_mp ? mrisBase.mris_mp->v_ripflag[vno] : mrisBase.mris->vertices[vno].ripflag;
}

static void MRISBase_v_xyz             (MRISBaseConst mrisBase, int vno, float* x, float* y, float* z) {
  if (mrisBase.mris_mp) {
    *x = mrisBase.mris_mp->v_x[vno];
    *y = mrisBase.mris_mp->v_y[vno];
    *z = mrisBase.mris_mp->v_z[vno];
  } else {
    VERTEX const * const v = &mrisBase.mris->vertices[vno];
    *x = v->x; *y = v->y; *z = v->z;
  }
}

static void MRISBase_set_v_xyz         (MRISBase      mrisBase, int vno, float  x, float  y, float  z) {
  if (mrisBase.mris_mp) {
    mrisBase.mris_mp->v_x[vno] = x;
    mrisBase.mris_mp->v_y[vno] = y;
    mrisBase.mris_mp->v_z[vno] = z;
  } else {
    MRISsetXYZ(mrisBase.mris, vno, x,y,z);
  }
}

static void MRISBase_set_v_assigned_fno(MRISBase      mrisBase, int vno, int to) {
  if (mrisBase.mris_mp) mrisBase.mris_mp->v_assigned_fno[vno] = to; else mrisBase.mris->vertices[vno].assigned_fno = to;
}
