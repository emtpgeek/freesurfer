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

