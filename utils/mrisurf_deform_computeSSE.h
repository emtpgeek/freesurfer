// Do not include, intended only to be included into mrisurf_deform.c

// These are in the order the original code computed them, so that side effects are not reordered
// In older code the ashburner_triangle is computed but not used , here it is not computed at all
//
#ifndef MRIS_COMPUTESSE_CUDA
    #define COMPUTE_DISTANCE_ERROR mrisComputeDistanceError(mris, parms)
#else
    #define COMPUTE_DISTANCE_ERROR mrisComputeDistanceErrorCUDA(mris, mrisc, parms)
#endif

// ELT  are only doable using MRIS
// ELTM are also doable using MRIS_MP
//
#ifndef COMPILING_MRIS_MP
#define ELTM(NAME,MULTIPLIER,COND,EXPR,EXPRM) ELT(NAME, MULTIPLIER, COND, EXPR)
#endif

#define SSE_TERMS \
      ELT(sse_area                  , parms->l_parea,                            true,    computed_area                                                                   ) SEP \
      ELT(sse_neg_area              , parms->l_area,                             true,    computed_neg_area                                                               ) SEP \
      ELT(sse_repulse               , 1.0,                     (parms->l_repulse > 0),    mrisComputeRepulsiveEnergy(mris, parms->l_repulse, mht_v_current, mht_f_current)) SEP \
      ELT(sse_repulsive_ratio       , 1.0,                                       true,    mrisComputeRepulsiveRatioEnergy(mris, parms->l_repulse_ratio)                   ) SEP \
      ELT(sse_tsmooth               , 1.0,                                       true,    mrisComputeThicknessSmoothnessEnergy(mris, parms->l_tsmooth, parms)             ) SEP \
      ELT(sse_thick_min             , parms->l_thick_min,                        true,    mrisComputeThicknessMinimizationEnergy(mris, parms->l_thick_min, parms)         ) SEP \
      ELT(sse_ashburner_triangle    , parms->l_ashburner_triangle,               false,   mrisComputeAshburnerTriangleEnergy(mris, parms->l_ashburner_triangle, parms)    ) SEP \
      ELT(sse_thick_parallel        , parms->l_thick_parallel,                   true,    mrisComputeThicknessParallelEnergy(mris, parms->l_thick_parallel, parms)        ) SEP \
      ELT(sse_thick_normal          , parms->l_thick_normal,                     true,    mrisComputeThicknessNormalEnergy(mris, parms->l_thick_normal, parms)            ) SEP \
      ELT(sse_thick_spring          , parms->l_thick_spring,                     true,    mrisComputeThicknessSpringEnergy(mris, parms->l_thick_spring, parms)            ) SEP \
      ELT(sse_nl_area               , parms->l_nlarea,        !FZERO(parms->l_nlarea),    mrisComputeNonlinearAreaSSE(mris)                                               ) SEP \
      ELT(sse_nl_dist               , parms->l_nldist,        !DZERO(parms->l_nldist),    mrisComputeNonlinearDistanceSSE(mris)                                           ) SEP \
      ELT(sse_dist                  , parms->l_dist,          !DZERO(parms->l_dist),      COMPUTE_DISTANCE_ERROR                                                          ) SEP \
      ELT(sse_spring                , parms->l_spring,        !DZERO(parms->l_spring),    mrisComputeSpringEnergy(mris)                                                   ) SEP \
      ELT(sse_lap                   , parms->l_lap,           !DZERO(parms->l_lap),       mrisComputeLaplacianEnergy(mris)                                                ) SEP \
      ELT(sse_tspring               , parms->l_tspring,       !DZERO(parms->l_tspring),   mrisComputeTangentialSpringEnergy(mris)                                         ) SEP \
      ELT(sse_nlspring              , parms->l_nlspring,      !DZERO(parms->l_nlspring),  mrisComputeNonlinearSpringEnergy(mris, parms)                                   ) SEP \
      ELT(sse_curv                  , l_curv_scaled,          !DZERO(parms->l_curv),      mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv)                           ) SEP \
      ELT(sse_corr                  , l_corr,                 !DZERO(l_corr),             mrisComputeCorrelationError(mris, parms, 1, false)                                     ) SEP \
      ELT(sse_val                   , parms->l_intensity,     !DZERO(parms->l_intensity), mrisComputeIntensityError(mris, parms)                                          ) SEP \
      ELT(sse_loc                   , parms->l_location,      !DZERO(parms->l_location),  mrisComputeTargetLocationError(mris, parms)                                     ) SEP \
      ELT(sse_dura                  , parms->l_dura,          !DZERO(parms->l_dura),      mrisComputeDuraError(mris, parms)                                               ) SEP \
      ELT(sse_histo                 , parms->l_histo,         !DZERO(parms->l_histo),     mrisComputeHistoNegativeLikelihood(mris, parms)                                 ) SEP \
      ELT(sse_map                   , parms->l_map,           !DZERO(parms->l_map),       mrisComputeNegativeLogPosterior(mris, parms, NULL)                              ) SEP \
      ELT(sse_map2d                 , parms->l_map2d,         !DZERO(parms->l_map2d),     mrisComputeNegativeLogPosterior2D(mris, parms, NULL)                            ) SEP \
      ELT(sse_grad                  , parms->l_grad,          !DZERO(parms->l_grad),      mrisComputeIntensityGradientError(mris, parms)                                  ) SEP \
      ELT(sse_sphere                , parms->l_sphere,        !DZERO(parms->l_sphere),    mrisComputeSphereError(mris, parms->l_sphere, parms->a)                         ) SEP \
      ELT(sse_shrinkwrap            , parms->l_shrinkwrap,    !DZERO(parms->l_shrinkwrap),mrisComputeShrinkwrapError(mris, parms->mri_brain, parms->l_shrinkwrap)         ) SEP \
      ELT(sse_expandwrap            , parms->l_expandwrap,    !DZERO(parms->l_expandwrap),mrisComputeExpandwrapError(mris, parms->mri_brain, parms->l_expandwrap, parms->target_radius)) SEP \
      ELT(sse_vectorCorrelationError, 1.0,                    use_multiframes,            mrisComputeVectorCorrelationError(mris, parms, 1)                               )     \
      // end of list

#if defined(MRIS_COMPUTESSE_CUDA)
static double MRIScomputeSSE_CUDA(MRIS *mris, MRI_CUDA_SURFACE *mrisc, INTEGRATION_PARMS *parms)

#elif defined(COMPILING_MRIS_MP)
bool MRISMP_computeSSE_canDo(INTEGRATION_PARMS *parms)
{
  bool   const use_multiframes  = !!(parms->flags & IP_USE_MULTIFRAMES);
  double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);

  bool result = true;
#define SEP
#define ELT(NAME, MULTIPLIER, COND, EXPR) \
  if (COND) { static bool reported = false; \
    if (!reported) { reported = true; fprintf(stdout, "%s:%d can't do %s %s\n", __FILE__,__LINE__,#NAME,#EXPR); } \
    result = false; \
  }
  SSE_TERMS
#undef ELT
#undef SEP
  return result;
}


double MRISMP_computeSSE(MRIS_MP* mris, INTEGRATION_PARMS *parms)

#else
double MRIScomputeSSE(MRIS* mris, INTEGRATION_PARMS *parms)
#endif

#if defined(COMPILING_MRIS_MP)
{
  cheapAssert(false);
  return 0.0;
}
#else
{
  bool   const use_multiframes  = !!(parms->flags & IP_USE_MULTIFRAMES);
  double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);
  double const l_curv_scaled    = (double)parms->l_curv * CURV_SCALE;
  double const area_scale =
#if METRIC_SCALE
    (mris->patch || mris->noscale) ? 1.0 : mris->orig_area / mris->total_area;
#else
    1.0;
#endif

  double relevant_angle = 0, computed_neg_area = 0, computed_area = 0;

  if (!FZERO(parms->l_angle) || !FZERO(parms->l_area) || (!FZERO(parms->l_parea))) {

    relevant_angle = 0; computed_neg_area = 0; computed_area = 0;

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

      FACE const * const face = &mris->faces[fno];
      if (face->ripflag) ROMP_PF_continue;
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);

      {
        double const delta = (double)(area_scale * face->area - fNorm->orig_area);
#if ONLY_NEG_AREA_TERM
        if (face->area < 0.0f) computed_neg_area += delta * delta;
#endif
        computed_area += delta * delta;
      }
      
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        double delta = deltaAngle(face->angle[ano], face->orig_angle[ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[ano] >= 0.0f) delta = 0.0f;

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

  }

  MHT* mht_v_current = NULL;
  MHT* mht_f_current = NULL;
  if (!FZERO(parms->l_repulse)) {
#ifdef COMPILING_MRIS_MP
    cheapAssert(false);
#else
    double vmean, vsigma;
    vmean = MRIScomputeTotalVertexSpacingStats(mris, &vsigma, NULL, NULL, NULL, NULL);
    mht_v_current = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, vmean);
    mht_f_current = MHTcreateFaceTable_Resolution  (mris, CURRENT_VERTICES, vmean);
#endif
  }


#define SEP
#define ELT(NAME, MULTIPLIER, COND, EXPR) double const NAME = (COND) ? (EXPR) : 0.0;
    SSE_TERMS
#undef ELT
#undef SEP

  if (parms->l_thick_spring > 0 || parms->l_thick_min > 0 || parms->l_thick_parallel > 0 /* && DIAG_VERBOSE_ON*/)
    printf("min=%2.3f, parallel=%2.4f, normal=%2.4f, spring=%2.4f, ashburner=%2.3f, tsmooth=%2.3f\n",
           sse_thick_min            / (float)mris->nvertices,
           sse_thick_parallel       / (float)mris->nvertices,
           sse_thick_normal         / (float)mris->nvertices,
           sse_thick_spring         / (float)mris->nvertices,
           sse_ashburner_triangle   / (float)mris->nvertices,
           sse_tsmooth              / (float)mris->nvertices);
           
  double sse_init = 0;

  if (gMRISexternalSSE) {
    sse_init = (*gMRISexternalSSE)(mris, parms);
  }
  
  double sse = sse_init +
#define SEP +
#define ELT(NAME, MULTIPLIER, COND, EXPR) (MULTIPLIER) * (NAME)
    SSE_TERMS ;
#undef ELT
#undef SEP

  if (false) {
    #define SEP
    double sum = 0;
    #define ELT(NAME, MULTIPLIER, COND, EXPR) fprintf(stdout, "new %s : %f \n", #NAME, (MULTIPLIER) * (NAME));  sum += (MULTIPLIER) * (NAME);
    ELT(sse_init, 1, true, sse_init)
    SSE_TERMS
    fprintf(stdout, "new sum = %f \n", sum);
    #undef ELT
    #undef SEP
  }
  

  if (mht_v_current) MHTfree(&mht_v_current);
  if (mht_f_current) MHTfree(&mht_f_current);

  if (!devFinite(sse)) {
    DiagBreak();
  }

#undef COMPUTE_DISTANCE_ERROR

  return sse;
}
#endif

