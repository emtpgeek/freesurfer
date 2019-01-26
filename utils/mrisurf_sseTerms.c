/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2019 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_sseTerms.h"
#include "mrisurf_project.h"


// The SSE terms are either computed by iterating over the vertices or the faces
// Ideally there would be one pass over each, to 
//      a) minimize the scanning logic
//      b) minimize the cache traffic
// Furthermore it would be more efficient to do these using the dense MRIS_MP format rather than the MRIS format
// We are working towards this ideal...
//
// For now, only the terms used in recon-all bert have been moved here
// The rest will follow



// Utilities
//



//====================================================================================
// VERTEX SSE TERMS
//
static double mrisComputeCorrelationErrorTerm(MRIS *mris, int vno, INTEGRATION_PARMS *parms, int use_stds, bool vertexTrace) {
  VERTEX const * const v  = &mris->vertices[vno];
  
  float const x = v->x;
  float const y = v->y;
  float const z = v->z;

#if 0
  double const src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0, vertexTrace) ;                                      // YUCK
#else
  double const src = v->curv;                                                                                             // YUCK - this is written by mrisComputeThicknessMinimizationEnergyTerm
#endif

  double const target = MRISPfunctionValTraceable(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no, vertexTrace);    // YUCK

#define DEFAULT_STD 4.0f
#define DISABLE_STDS 0
#if DISABLE_STDS
  double const std = 1.0f;
#else
  double std = MRISPfunctionValTraceable(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no + 1, vertexTrace);         // YUCK
  std = sqrt(std);
  if (FZERO(std)) {
    std = DEFAULT_STD /*FSMALL*/;
  }
  if (!use_stds) {
    std = 1.0f;
  }
#endif

  double const delta = (src - target) / std;
  if (!isfinite(target) || !isfinite(delta)) {
    DiagBreak();
  }

  if (parms->geometry_error) {
    parms->geometry_error[vno] = (delta * delta);                                                                         // YUCK
  }

  if (parms->abs_norm) {
    return fabs(delta);
  } else {
    return delta * delta;
  }
}


static double mrisComputeRepulsiveRatioEnergyTerm(MRIS *mris, int vno) {
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];
  
  float const x = v->x;
  float const y = v->y;
  float const z = v->z;
  float const cx = v->cx;
  float const cy = v->cy;
  float const cz = v->cz;
  
  double v_sse = 0.0;

  int n;
  for (n = 0; n < vt->vnum; n++) {
    VERTEX const * const vn = &mris->vertices[vt->v[n]];

    if (vn->ripflag) continue;

    double dx = x - vn->x;
    double dy = y - vn->y;
    double dz = z - vn->z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    
    double cdx = cx - vn->cx;
    double cdy = cy - vn->cy;
    double cdz = cz - vn->cz;
    double canon_dist = sqrt(cdx * cdx + cdy * cdy + cdz * cdz) + REPULSE_E;
    
    dist /= canon_dist;
    dist += REPULSE_E;

#if 0
    v_sse += REPULSE_K / (dist*dist*dist*dist) ;
#else
    v_sse += REPULSE_K / (dist * dist);
#endif
  }

  return v_sse;
}


static double mrisComputeSpringEnergyTerm(MRIS *mris, int vno, double area_scale) {
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];
  
  double v_sse = 0.0;
  int n;
  for (n = 0; n < vt->vnum; n++) {
    v_sse += (v->dist[n] * v->dist[n]);
  }

  return area_scale * v_sse;
}


static double mrisComputeThicknessMinimizationEnergyTerm(MRIS *mris, int vno, INTEGRATION_PARMS *parms) {
  VERTEX /* const */ * const v = &mris->vertices[vno];                                      // YUCK

  float thick_sq = mrisSampleMinimizationEnergy(mris, v, parms, v->x, v->y, v->z);
  v->curv = sqrt(thick_sq);                                                                 // YUCK

  return thick_sq;
}


static double mrisComputeThicknessNormalEnergyTerm(MRIS *mris, int vno, INTEGRATION_PARMS *parms) {
  VERTEX const * const v = &mris->vertices[vno];
  return mrisSampleNormalEnergy(mris, v, parms, v->x, v->y, v->z);
}


static double mrisComputeThicknessSpringEnergyTerm(MRIS *mris, int vno, INTEGRATION_PARMS *parms) {
  VERTEX const * const v = &mris->vertices[vno];
  return mrisSampleSpringEnergy(mris, vno, v->x, v->y, v->z, parms);
}


static double mrisComputeThicknessParallelEnergyTerm(MRIS *mris, int vno, INTEGRATION_PARMS *parms) {
  VERTEX * const v = &mris->vertices[vno];
  return mrisSampleParallelEnergy(mris, vno, parms, v->x, v->y, v->z);
}


static double mrisComputeThicknessSmoothnessEnergyTerm(MRIS *mris, int vno, INTEGRATION_PARMS *parms) {
  VERTEX * const v = &mris->vertices[vno];
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];

  float xp,yp,zp;
  MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);

  double d0 = SQR(xp - v->whitex) + SQR(yp - v->whitey) + SQR(zp - v->whitez);
  double v_sse = 0.0;
  
  int n;
  for (n = 0; n < vt->vnum; n++) {
    VERTEX const * const vn = &mris->vertices[vt->v[n]];
    if (vn->ripflag) continue;

    float xp, yp, zp;
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, vn->x, vn->y, vn->z, PIAL_VERTICES, &xp, &yp, &zp);

    double dx = xp - vn->whitex;
    double dy = yp - vn->whitey;
    double dz = zp - vn->whitez;
    double dn = (dx * dx + dy * dy + dz * dz);
    v_sse += (dn - d0) * (dn - d0);
  }

  return v_sse;
}


static double mrisComputeNonlinearAreaSSETerm(MRIS *mris, int fno, double area_scale) {
  FACE *face = &mris->faces[fno];
  
  double ratio;
#define SCALE_NONLINEAR_AREA 0
#if SCALE_NONLINEAR_AREA
  if (!FZERO(face->orig_area)) {
    ratio = area_scale * face->area / face->orig_area;
  } else {
    ratio = 0.0f;
  }
#else
  ratio = area_scale * face->area;
#endif

  if (ratio > MAX_NEG_RATIO) {
    ratio = MAX_NEG_RATIO;
  } else if (ratio < -MAX_NEG_RATIO) {
    ratio = -MAX_NEG_RATIO;
  }
  
#if 0
  double error = (1.0 / NEG_AREA_K) * log(1.0+exp(-NEG_AREA_K*ratio)) ;
#else
  double error = (log(1.0 + exp(NEG_AREA_K * ratio)) / NEG_AREA_K) - ratio;
#endif

#undef SCALE_NONLINEAR_AREA

  return error;
}


// WALK ALL THE VERTICES
//
double mrisComputeCorrelationError(MRIS *mris, INTEGRATION_PARMS *parms, int use_stds, bool trace)
{
  float const l_corr = parms->l_corr + parms->l_pcorr; /* only one will be nonzero */

  if (FZERO(l_corr)) return (0.0);

  double sse = 0.0;
  
  #define ROMP_VARIABLE       vno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse  ROMP_PARTIALSUM(0)
    
    bool const vertexTrace = trace && (vno == 0);
    
    VERTEX *v = &mris->vertices[vno];

    if (v->ripflag) {
      ROMP_PF_continue;
    }

    sse += mrisComputeCorrelationErrorTerm(mris, vno, parms, use_stds, vertexTrace);

    #undef sse
  #include "romp_for_end.h"

  return (sse);
}


double mrisComputeRepulsiveRatioEnergy(MRIS *mris, double l_repulse)
{
  int vno;
  double sse_repulse;

  if (FZERO(l_repulse))
    return (0.0);

  for (sse_repulse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v  = &mris->vertices[vno];

    if (v->ripflag) 
      continue;

    sse_repulse += mrisComputeRepulsiveRatioEnergyTerm(mris, vno);
  }
  
  return (l_repulse * sse_repulse);
}


double mrisComputeSpringEnergy(MRIS *mris)
{
  double area_scale;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  double sse_spring = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    if (v->ripflag) continue;

    sse_spring += mrisComputeSpringEnergyTerm(mris, vno, area_scale);
  }
  
  return sse_spring;
}


double mrisComputeThicknessMinimizationEnergy(MRIS *mris, double l_thick_min, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_thick_min)) {
    return (0.0);
  }

  double sse_tmin = 0.0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX* v = &mris->vertices[vno];
    if (v->ripflag) continue;
    
    sse_tmin += mrisComputeThicknessMinimizationEnergyTerm(mris, vno, parms);
  }

  sse_tmin /= 2;

  return (sse_tmin);
}


double mrisComputeThicknessNormalEnergy(MRIS *mris, double l_thick_normal, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_thick_normal)) return (0.0);

  double sse_tnormal = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    if (v->ripflag) continue;

    sse_tnormal += mrisComputeThicknessNormalEnergyTerm(mris, l_thick_normal, parms);
  }

  sse_tnormal /= 2;

  return (sse_tnormal);
}


double mrisComputeThicknessSpringEnergy(MRIS *mris, double l_thick_spring, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_thick_spring)) {
    return (0.0);
  }

  double sse_spring = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) continue;

    sse_spring += mrisComputeThicknessSpringEnergyTerm(mris,vno,parms);
  }

  sse_spring /= 2;
  return (sse_spring);
}


double mrisComputeThicknessParallelEnergy(MRIS *mris, double l_thick_parallel, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_thick_parallel)) {
    return (0.0);
  }

  mrisAssignFaces(mris, (MHT *)(parms->mht), CANONICAL_VERTICES);  // don't look it up every time

  double sse_tparallel = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) continue;

    sse_tparallel += mrisComputeThicknessParallelEnergyTerm(mris, vno, parms);
  }
  
  sse_tparallel /= 2;
  return (sse_tparallel);
}


double mrisComputeThicknessSmoothnessEnergy(MRIS *mris, double l_tsmooth, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_tsmooth)) {
    return (0.0);
  }

  double sse_tsmooth = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;

    sse_tsmooth += mrisComputeThicknessSmoothnessEnergyTerm(mris,vno,parms);
  }
  
  return (sse_tsmooth);
}


//====================================================================================
// FACE SSE TERMS
//
double mrisComputeNonlinearAreaSSE(MRIS *mris)
{
  double area_scale;
#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  double sse = 0;
  
  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    #define sse  ROMP_PARTIALSUM(0)

    FACE *face = &mris->faces[fno];
    if (face->ripflag) ROMP_PF_continue;

    sse += mrisComputeNonlinearAreaSSETerm(mris, fno, area_scale);
    
    #undef sse
  #include "romp_for_end.h"
  
  return (sse);
}
