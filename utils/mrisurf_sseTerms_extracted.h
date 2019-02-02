#ifdef COMPILING_MRIS_MP
    #define MRIS_INFO MRIS_MP
    #define GET_MRISBase                                \
        MRISBase mrisBase = MRISBaseCtr(mris,NULL);
    #define FUNCTION_NAME(MRIS_NAME,MRISMP_NAME) MRISMP_NAME
    #define GET_V
    #define GET_V_RIPFLAG                               \
        bool const ripflag = mris->v_ripflag[vno];
    #define GET_V_XYZ                                   \
        float const x = mris->v_x[vno];                 \
        float const y = mris->v_y[vno];                 \
        float const z = mris->v_z[vno];                 \
        // end of macro
    #define GET_V2_XYZ_Ripflag2(VNO2)                   \
        int const vno2 = VNO2;                          \
        bool  const ripflag2 = mris->v_ripflag[vno2];   \
        float const x2 = mris->v_x[vno2];               \
        float const y2 = mris->v_y[vno2];               \
        float const z2 = mris->v_z[vno2];               \
        // end of macro
    #define GET_CXYZ                                    \
        float const cx = mris->v_cx[vno];               \
        float const cy = mris->v_cy[vno];               \
        float const cz = mris->v_cz[vno];               \
        // end of macro
    #define GET_CXYZ2                                   \
        float const cx2 = mris->v_cx[vno2];             \
        float const cy2 = mris->v_cy[vno2];             \
        float const cz2 = mris->v_cz[vno2];             \
        // end of macro
    #define GET_V_DIST                                  \
        float const * v_dist = mris->v_dist[vno];       \
        // end of macro
    #define GET_CURV                                    \
        float const curv = mris->v_curv[vno];
    #define GET_PCURV                                   \
        float* const pcurv = &mris->v_curv[vno];
    #define GET_WHITEXYZ                                \
        float const whitex = mris->v_whitex[vno];       \
        float const whitey = mris->v_whitey[vno];       \
        float const whitez = mris->v_whitez[vno];       \
        // end of macro
    #define GET_WHITEXYZ2                               \
        float const whitex2 = mris->v_whitex[vno2];     \
        float const whitey2 = mris->v_whitey[vno2];     \
        float const whitez2 = mris->v_whitez[vno2];     \
        // end of macro
    #define GET_F                                       \
        // end of macro
    #define GET_F_RIPFLAG                               \
        bool const ripflag = mris->f_ripflag[fno];      \
        // end of macro
    #define GET_F_AREA                                  \
        float const area = mris->f_area[fno];           \
        // end of macro
    #define GET_F_AREA_ANGLE_ORIGAREA_ORIGANGLE         \
        float const area      = mris->f_area     [fno]; \
        float const orig_area = mris->f_norm_orig_area[fno]; \
        float const* const angle      = mris->f_angle     [fno];  \
        float const* const orig_angle = mris->f_orig_angle[fno];  \
        // end of macro
#else
    #define MRIS_INFO MRIS
    #define GET_MRISBase                                \
        MRISBase mrisBase = MRISBaseCtr(NULL,mris);
    #define FUNCTION_NAME(MRIS_NAME,MRISMP_NAME) MRIS_NAME
    #define GET_V                                       \
        VERTEX const * const v  = &mris->vertices[vno];
    #define GET_V_RIPFLAG                               \
        GET_V                                           \
        bool const ripflag = v->ripflag;
    #define GET_V_XYZ                                   \
        GET_V                                           \
        float const x = v->x;                           \
        float const y = v->y;                           \
        float const z = v->z;                           \
        // end of macro
    #define GET_V2_XYZ_Ripflag2(VNO2)                   \
        int const vno2 = VNO2;                          \
        VERTEX const * const v2 = &mris->vertices[vno2];\
        bool  const ripflag2 = v2->ripflag;             \
        float const x2 = v2->x;                         \
        float const y2 = v2->y;                         \
        float const z2 = v2->z;                         \
        // end of macro
    #define GET_CXYZ                                    \
        float const cx = v->cx;                         \
        float const cy = v->cy;                         \
        float const cz = v->cz;                         \
        // end of macro
    #define GET_CXYZ2                                   \
        float const cx2 = v2->cx;                       \
        float const cy2 = v2->cy;                       \
        float const cz2 = v2->cz;                       \
        // end of macro
    #define GET_V_DIST                                  \
        GET_V                                           \
        float const * v_dist = v->dist;                 \
        // end of macro
    #define GET_CURV                                    \
        float const curv = v->curv;
    #define GET_PCURV                                   \
        float* const pcurv = &mris->vertices[vno].curv;
    #define GET_WHITEXYZ                                \
        float const whitex = v->whitex;                 \
        float const whitey = v->whitey;                 \
        float const whitez = v->whitez;                 \
        // end of macro
    #define GET_WHITEXYZ2                               \
        float const whitex2 = v2->whitex;               \
        float const whitey2 = v2->whitey;               \
        float const whitez2 = v2->whitez;               \
        // end of macro
    #define GET_F_RIPFLAG                               \
        bool const ripflag = mris->faces[fno].ripflag;  \
        // end of macro
    #define GET_F_AREA                                  \
        FACE const * const face = &mris->faces[fno];    \
        float const area      = face->area;             \
        // end of macro
    #define GET_F_AREA_ANGLE_ORIGAREA_ORIGANGLE         \
        FACE const * const face = &mris->faces[fno];    \
        FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno); \
        float const area      = face->area;             \
        float const orig_area = fNorm->orig_area;       \
        float const* const angle      = face->angle;        \
        float const* const orig_angle = face->orig_angle;   \
        // end of macro
        
#endif


void FUNCTION_NAME(mrisComputeFaceRelevantAngleAndArea,mrismp_ComputeFaceRelevantAngleAndArea) (
  MRIS_INFO *mris, INTEGRATION_PARMS *parms, double* p_relevant_angle, double* p_computed_neg_area, double* p_computed_area)
{
  double const area_scale =
#if METRIC_SCALE
    (mris->patch || mris->noscale) ? 1.0 : mris->orig_area / mris->total_area;
#else
    1.0;
#endif

  double relevant_angle = 0.0, computed_neg_area = 0.0, computed_area = 0.0;

  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces
    
  #define ROMP_SUMREDUCTION0  relevant_angle
  #define ROMP_SUMREDUCTION1  computed_neg_area
  #define ROMP_SUMREDUCTION2  computed_area
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define relevant_angle    ROMP_PARTIALSUM(0)
    #define computed_neg_area ROMP_PARTIALSUM(1)
    #define computed_area     ROMP_PARTIALSUM(2)

      GET_F_RIPFLAG
      if (ripflag) ROMP_PF_continue;
      
      GET_F_AREA_ANGLE_ORIGAREA_ORIGANGLE
      
      {
        double const delta = (double)(area_scale * area - orig_area);
#if ONLY_NEG_AREA_TERM
        if (area < 0.0f) computed_neg_area += delta * delta;
#endif
        computed_area += delta * delta;
      }
      
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        double delta = deltaAngle(angle[ano], orig_angle[ano]);
#if ONLY_NEG_AREA_TERM
        if (angle[ano] >= 0.0f) delta = 0.0f;

#endif
        relevant_angle += delta * delta;
      }
      
      if (!isfinite(computed_area) || !isfinite(relevant_angle)) {
        ErrorExit(ERROR_BADPARM, "sse not finite at face %d!\n", fno);
      }

    #undef relevant_angle
    #undef computed_neg_area
    #undef computed_area

  #include "romp_for_end.h"

   
  *p_relevant_angle    = relevant_angle;
  *p_computed_neg_area = computed_neg_area;
  *p_computed_area     = computed_area;
}


//====================================================================================
// VERTEX SSE TERMS
//
static double FUNCTION_NAME(mrisComputeCorrelationErrorTerm,mrismp_ComputeCorrelationErrorTerm) (
  MRIS_INFO *mris, int vno, INTEGRATION_PARMS *parms, int use_stds, bool vertexTrace) 
{
  GET_V_XYZ
  GET_CURV
  
#if 0
  double const src = MRISPfunctionVal(parms->mrisp, mris->radius, x, y, z, 0, vertexTrace) ;
#else
  double const src = curv;
#endif

  double const target = MRISPfunctionValTraceable(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no, vertexTrace);

#define DEFAULT_STD 4.0f
#define DISABLE_STDS 0
#if DISABLE_STDS
  double const std = 1.0f;
#else
  double std = MRISPfunctionValTraceable(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no + 1, vertexTrace);
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


static double FUNCTION_NAME(mrisComputeRepulsiveRatioEnergyTerm, mrismp_ComputeRepulsiveRatioEnergyTerm) (MRIS_INFO *mris, int vno) 
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  GET_V_XYZ
  GET_CXYZ
  
  double v_sse = 0.0;

  int n;
  for (n = 0; n < vt->vnum; n++) {
    GET_V2_XYZ_Ripflag2(vt->v[n])
    GET_CXYZ2
    
    if (ripflag2) continue;

    double dx = x - x2;
    double dy = y - y2;
    double dz = z - z2;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    
    double cdx = cx - cx2;
    double cdy = cy - cy2;
    double cdz = cz - cz2;
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


static double FUNCTION_NAME(mrisComputeSpringEnergyTerm, mrismp_ComputeSpringEnergyTerm) (MRIS_INFO *mris, int vno, double area_scale) 
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  GET_V_DIST
  
  double v_sse = 0.0;
  int n;
  for (n = 0; n < vt->vnum; n++) {
    v_sse += SQR(v_dist[n]);
  }

  return area_scale * v_sse;
}


static double FUNCTION_NAME(mrisComputeThicknessMinimizationEnergyTerm,mrismp_ComputeThicknessMinimizationEnergyTerm) (
  MRIS_INFO *mris, int vno, INTEGRATION_PARMS *parms) 
{
  GET_V_XYZ
  GET_PCURV                                                                                // YUCK

  float thick_sq = FUNCTION_NAME(mrisSampleMinimizationEnergy,mrismp_SampleMinimizationEnergy)(mris, vno, parms, x, y, z);
  *pcurv = sqrt(thick_sq);                                                                 // YUCK

  return thick_sq;
}


static double FUNCTION_NAME(mrisComputeThicknessNormalEnergyTerm,mrismp_ComputeThicknessNormalEnergyTerm) (
  MRIS_INFO *mris, int vno, INTEGRATION_PARMS *parms) 
{
  GET_V_XYZ
  return FUNCTION_NAME(mrisSampleNormalEnergy,mrismp_SampleNormalEnergy)(mris, vno, parms, x, y, z);
}


static double FUNCTION_NAME(mrisComputeThicknessSpringEnergyTerm, mrismp_ComputeThicknessSpringEnergyTerm) (
  MRIS_INFO *mris, int vno, INTEGRATION_PARMS *parms)
{
  GET_V_XYZ
  return FUNCTION_NAME(mrisSampleSpringEnergy,mrismp_SampleSpringEnergy)(mris, vno, x, y, z, parms);
}


static double FUNCTION_NAME(mrisComputeThicknessParallelEnergyTerm, mrismp_ComputeThicknessParallelEnergyTerm) (
  MRIS_INFO *mris, int vno, INTEGRATION_PARMS *parms)
{
  GET_V_XYZ
  return FUNCTION_NAME(mrisSampleParallelEnergy,mrismp_SampleParallelEnergy)(mris, vno, parms, x, y, z);
}


static double FUNCTION_NAME(mrisComputeThicknessSmoothnessEnergyTerm,mrismp_ComputeThicknessSmoothnessEnergyTerm) (
  MRIS_INFO *mris, int vno, INTEGRATION_PARMS *parms) 
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  GET_V_XYZ
  GET_WHITEXYZ

  float xp,yp,zp;
  FUNCTION_NAME(MRISsampleFaceCoordsCanonical,MRISMP_sampleFaceCoordsCanonical)((MHT *)(parms->mht), mris, x, y, z, PIAL_VERTICES, &xp, &yp, &zp);

  double d0 = SQR(xp - whitex) + SQR(yp - whitey) + SQR(zp - whitez);
  double v_sse = 0.0;
  
  int n;
  for (n = 0; n < vt->vnum; n++) {
    GET_V2_XYZ_Ripflag2(vt->v[n]);
    if (ripflag2) continue;
    
    float xp, yp, zp;
    FUNCTION_NAME(MRISsampleFaceCoordsCanonical,MRISMP_sampleFaceCoordsCanonical)((MHT *)(parms->mht), mris, x2, y2, z2, PIAL_VERTICES, &xp, &yp, &zp);

    GET_WHITEXYZ2

    double dx = xp - whitex2;
    double dy = yp - whitey2;
    double dz = zp - whitez2;
    double dn = (dx * dx + dy * dy + dz * dz);
    v_sse += (dn - d0) * (dn - d0);
  }

  return v_sse;
}


static double FUNCTION_NAME(mrisComputeNonlinearAreaSSETerm,mrismp_ComputeNonlinearAreaSSETerm) (
  MRIS_INFO *mris, int fno, double area_scale) 
{
  GET_F_AREA
  
  double ratio;
#define SCALE_NONLINEAR_AREA 0
#if SCALE_NONLINEAR_AREA
  if (!FZERO(face->orig_area)) {
    ratio = area_scale * area / face->orig_area;
  } else {
    ratio = 0.0f;
  }
#else
  ratio = area_scale * area;
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
double FUNCTION_NAME(mrisComputeCorrelationError, mrismp_ComputeCorrelationError) (
  MRIS_INFO *mris, INTEGRATION_PARMS *parms, int use_stds, bool trace)
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
    
    GET_V_RIPFLAG

    if (ripflag) {
      ROMP_PF_continue;
    }

    sse += FUNCTION_NAME(mrisComputeCorrelationErrorTerm,mrismp_ComputeCorrelationErrorTerm)(mris, vno, parms, use_stds, vertexTrace);

    #undef sse
  #include "romp_for_end.h"

  return (sse);
}


double FUNCTION_NAME(mrisComputeRepulsiveRatioEnergy,mrismp_ComputeRepulsiveRatioEnergy) (
  MRIS_INFO *mris, double l_repulse)
{
  int vno;
  double sse_repulse;

  if (FZERO(l_repulse))
    return (0.0);

  for (sse_repulse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    GET_V_RIPFLAG

    if (ripflag) 
      continue;

    sse_repulse += FUNCTION_NAME(mrisComputeRepulsiveRatioEnergyTerm,mrismp_ComputeRepulsiveRatioEnergyTerm)(mris, vno);
  }
  
  return (l_repulse * sse_repulse);
}


double FUNCTION_NAME(mrisComputeSpringEnergy,mrismp_ComputeSpringEnergy) (
  MRIS_INFO *mris)
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
    GET_V_RIPFLAG
    
    if (ripflag) continue;

    sse_spring += FUNCTION_NAME(mrisComputeSpringEnergyTerm,mrismp_ComputeSpringEnergyTerm)(mris, vno, area_scale);
  }
  
  return sse_spring;
}


double FUNCTION_NAME(mrisComputeThicknessMinimizationEnergy,mrismp_ComputeThicknessMinimizationEnergy) (
  MRIS_INFO *mris, double l_thick_min, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_thick_min)) {
    return (0.0);
  }

  double sse_tmin = 0.0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    GET_V_RIPFLAG
    
    if (ripflag) continue;
    
    sse_tmin += FUNCTION_NAME(mrisComputeThicknessMinimizationEnergyTerm,mrismp_ComputeThicknessMinimizationEnergyTerm)(mris, vno, parms);
  }

  sse_tmin /= 2;

  return (sse_tmin);
}


double FUNCTION_NAME(mrisComputeThicknessNormalEnergy,mrismp_ComputeThicknessNormalEnergy) (
  MRIS_INFO *mris, double l_thick_normal, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_thick_normal)) return (0.0);

  double sse_tnormal = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    GET_V_RIPFLAG
    
    if (ripflag) continue;

    sse_tnormal += FUNCTION_NAME(mrisComputeThicknessNormalEnergyTerm,mrismp_ComputeThicknessNormalEnergyTerm)(mris, l_thick_normal, parms);
  }

  sse_tnormal /= 2;

  return (sse_tnormal);
}


double FUNCTION_NAME(mrisComputeThicknessSpringEnergy,mrismp_ComputeThicknessSpringEnergy) (
  MRIS_INFO *mris, double l_thick_spring, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_thick_spring)) {
    return (0.0);
  }

  double sse_spring = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    GET_V_RIPFLAG
    
    if (ripflag) continue;

    sse_spring += FUNCTION_NAME(mrisComputeThicknessSpringEnergyTerm,mrismp_ComputeThicknessSpringEnergyTerm)(mris,vno,parms);
  }

  sse_spring /= 2;
  return (sse_spring);
}


double FUNCTION_NAME(mrisComputeThicknessParallelEnergy,mrismp_ComputeThicknessParallelEnergy) (
  MRIS_INFO *mris, double l_thick_parallel, INTEGRATION_PARMS *parms)
{
  GET_MRISBase
  
  if (FZERO(l_thick_parallel)) {
    return (0.0);
  }

  mrisAssignFaces(mrisBase, (MHT *)(parms->mht), CANONICAL_VERTICES);  // don't look it up every time

  double sse_tparallel = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    GET_V_RIPFLAG
    
    if (ripflag) continue;

    sse_tparallel += FUNCTION_NAME(mrisComputeThicknessParallelEnergyTerm,mrismp_ComputeThicknessParallelEnergyTerm)(mris, vno, parms);
  }
  
  sse_tparallel /= 2;
  return (sse_tparallel);
}


double FUNCTION_NAME(mrisComputeThicknessSmoothnessEnergy,mrismp_ComputeThicknessSmoothnessEnergy) (
  MRIS_INFO *mris, double l_tsmooth, INTEGRATION_PARMS *parms)
{
  if (FZERO(l_tsmooth)) {
    return (0.0);
  }

  double sse_tsmooth = 0.0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    GET_V_RIPFLAG
    
    if (ripflag) continue;

    sse_tsmooth += FUNCTION_NAME(mrisComputeThicknessSmoothnessEnergyTerm,mrismp_ComputeThicknessSmoothnessEnergyTerm)(mris,vno,parms);
  }
  
  return (sse_tsmooth);
}


//====================================================================================
// FACE SSE TERMS
//
double FUNCTION_NAME(mrisComputeNonlinearAreaSSE,mrismp_ComputeNonlinearAreaSSE) (
  MRIS_INFO *mris)
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

    GET_F_RIPFLAG
    if (ripflag) ROMP_PF_continue;

    sse += FUNCTION_NAME(mrisComputeNonlinearAreaSSETerm,mrismp_ComputeNonlinearAreaSSETerm)(mris, fno, area_scale);
    
    #undef sse
  #include "romp_for_end.h"
  
  return (sse);
}

#undef GET_F_AREA_ANGLE_ORIGAREA_ORIGANGLE
#undef GET_F_AREA
#undef GET_F_RIPFLAG
#undef GET_WHITEXYZ2
#undef GET_WHITEXYZ
#undef GET_PCURV
#undef GET_CURV
#undef GET_V_DIST
#undef GET_CXYZ2
#undef GET_CXYZ
#undef GET_V2_XYZ_Ripflag2
#undef GET_V_XYZ
#undef GET_V_RIPFLAG
#undef GET_V
#undef MRIS_INFO
#undef GET_MRISBase
#undef FUNCTION_NAME
