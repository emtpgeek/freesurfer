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
#include "mrisurf.h"

static void MRISBase_getWhichXYZ(MRISBaseConst mrisBase, int vno, int which, float* x, float* y, float* z);

int face_barycentric_coords2(MRISBaseConst mrisBase, int fno, int which_vertices,
                            double cx, double cy, double cz, double *pl1, double *pl2, double *pl3) ;



static void MRISBase_getWhichXYZ(MRISBaseConst mrisBase, int vno, int which, float* x, float* y, float* z) {
    cheapAssert(!mrisBase.mris_mp);
    MRISvertexCoord2XYZ_float(&mrisBase.mris->vertices[vno], which, x, y, z);
}


static int face_barycentric_coords(MRIS const *mris, int p1, int p2,
                            double p3, double p4, double p5, double *p6, double *p7, double *p8) 
{
  return face_barycentric_coords2(MRISBaseConstCtr(NULL,mris),p1,p2,p3,p4,p5,p6,p7,p8);
}


