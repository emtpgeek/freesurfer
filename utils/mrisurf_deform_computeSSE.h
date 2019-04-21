// Do not include, intended only to be included into mrisurf_deform.c

// These are in the order the original code computed them, so that side effects are not reordered
// In older code the ashburner_triangle is computed but not used , here it is not computed at all
//
#if defined(COMPILING_MRIS_MP)
    #define COMPUTE_DISTANCE_ERROR mrismp_ComputeDistanceError(mris, parms)
#else
    #define COMPUTE_DISTANCE_ERROR mrisComputeDistanceError(mris, parms)
#endif

// ELT  are only doable using MRIS
// ELTM are also doable using MRIS_MP
//
// Note: mrisComputeCorrelationError reads v->curv, which is written by mrisComputeThicknessMinimizationEnergy!
//      but the definition of curv here seems totally unrelated to that used elsewhere in the code
//      so I suspect it is being used as a convenient temporary... 
//
#define SSE_TERMS \
      ELTM(sse_area                 , parms->l_parea,                            true,    computed_area,computed_area                                                     ) SEP \
      ELTM(sse_neg_area             , parms->l_area,                             true,    computed_neg_area,computed_neg_area                                             ) SEP \
      ELT(sse_repulse               , 1.0,                     (parms->l_repulse > 0),    mrisComputeRepulsiveEnergy(mris, parms->l_repulse, mht_v_current, mht_f_current)) SEP \
      \
      ELTM(sse_repulsive_ratio      , 1.0,                                       true,    mrisComputeRepulsiveRatioEnergy   (mris, parms->l_repulse_ratio),                     \
                                                                                          mrismp_ComputeRepulsiveRatioEnergy(mris, parms->l_repulse_ratio)                ) SEP \
      ELTM(sse_tsmooth              , 1.0,                                       true,    mrisComputeThicknessSmoothnessEnergy(mris, parms->l_tsmooth, parms),                  \
                                                                                          mrismp_ComputeThicknessSmoothnessEnergy(mris, parms->l_tsmooth, parms)          ) SEP \
      ELTM(sse_thick_min            , parms->l_thick_min,                        true,    mrisComputeThicknessMinimizationEnergy(mris, parms->l_thick_min, parms),              \
                                                                                          mrismp_ComputeThicknessMinimizationEnergy(mris, parms->l_thick_min, parms)      ) SEP \
      \
      ELT(sse_ashburner_triangle    , parms->l_ashburner_triangle,               false,   mrisComputeAshburnerTriangleEnergy(mris, parms->l_ashburner_triangle, parms)    ) SEP \
      \
      ELTM(sse_thick_parallel       , parms->l_thick_parallel,                   true,    mrisComputeThicknessParallelEnergy(mris, parms->l_thick_parallel, parms),             \
                                                                                          mrismp_ComputeThicknessParallelEnergy(mris, parms->l_thick_parallel, parms)     ) SEP \
      ELTM(sse_thick_normal         , parms->l_thick_normal,                     true,    mrisComputeThicknessNormalEnergy(mris, parms->l_thick_normal, parms),                 \
                                                                                          mrismp_ComputeThicknessNormalEnergy(mris, parms->l_thick_normal, parms)         ) SEP \
      ELTM(sse_thick_spring         , parms->l_thick_spring,                     true,    mrisComputeThicknessSpringEnergy(mris, parms->l_thick_spring, parms),                 \
                                                                                          mrismp_ComputeThicknessSpringEnergy(mris, parms->l_thick_spring, parms)         ) SEP \
      \
      ELTM(sse_nl_area              , parms->l_nlarea,        !FZERO(parms->l_nlarea),    mrisComputeNonlinearAreaSSE(mris),                                                    \
                                                                                          mrismp_ComputeNonlinearAreaSSE(mris)                                            ) SEP \
      \
      ELT(sse_nl_dist               , parms->l_nldist,        !DZERO(parms->l_nldist),    mrisComputeNonlinearDistanceSSE(mris)                                           ) SEP \
      ELTM(sse_dist                 , parms->l_dist,          !DZERO(parms->l_dist),      COMPUTE_DISTANCE_ERROR,COMPUTE_DISTANCE_ERROR                                   ) SEP \
      \
      ELTM(sse_spring               , parms->l_spring,        !DZERO(parms->l_spring),    mrisComputeSpringEnergy(mris),                                                        \
                                                                                          mrismp_ComputeSpringEnergy(mris))                                                 SEP \
      \
      ELT(sse_lap                   , parms->l_lap,           !DZERO(parms->l_lap),       mrisComputeLaplacianEnergy(mris)                                                ) SEP \
      ELT(sse_tspring               , parms->l_tspring,       !DZERO(parms->l_tspring),   mrisComputeTangentialSpringEnergy(mris)                                         ) SEP \
      ELT(sse_nlspring              , parms->l_nlspring,      !DZERO(parms->l_nlspring),  mrisComputeNonlinearSpringEnergy(mris, parms)                                   ) SEP \
      ELT(sse_curv                  , l_curv_scaled,          !DZERO(parms->l_curv),      mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv)                           ) SEP \
      \
      ELTM(sse_corr                 , l_corr,                 !DZERO(l_corr),             mrisComputeCorrelationError(mris, parms, 1, false),                                   \
                                                                                          mrismp_ComputeCorrelationError(mris, parms, 1, false))                            SEP \
      \
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

#if defined(COMPILING_MRIS_MP)
bool MRISMP_computeSSE_canDo(INTEGRATION_PARMS *parms)
{
  bool   const use_multiframes  = !!(parms->flags & IP_USE_MULTIFRAMES);
  // double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);

  bool result = true;
#define SEP
#define ELTM(NAME,MULTIPLIER,COND,EXPR,EXPRM)
#define ELT(NAME, MULTIPLIER, COND, EXPR) \
  if (COND) { static bool reported = false; \
    if (!reported) { reported = true; fprintf(stdout, "%s:%d can't do %s %s\n", __FILE__,__LINE__,#NAME,#EXPR); } \
    result = false; \
  }
  SSE_TERMS
  ELT(sse_init,1.0,gMRISexternalSSE,)
#undef ELT
#undef ELTM
#undef SEP
  return result;
}
#endif

#if defined(COMPILING_MRIS_MP)
double MRISMP_computeSSE(MRIS_MP* mris, INTEGRATION_PARMS *parms)
#else
double MRIScomputeSSE(MRIS* mris, INTEGRATION_PARMS *parms)
#endif
{
  static const bool debug = false;
  
  bool   const use_multiframes  = !!(parms->flags & IP_USE_MULTIFRAMES);
  double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);
#ifndef COMPILING_MRIS_MP
  // not used so causes a warning message
  double const l_curv_scaled    = (double)parms->l_curv * CURV_SCALE;
#endif
  
  double relevant_angle = 0, computed_neg_area = 0, computed_area = 0;

  if (!FZERO(parms->l_angle) || !FZERO(parms->l_area) || (!FZERO(parms->l_parea))) {
#if defined(COMPILING_MRIS_MP)
    mrismp_ComputeFaceRelevantAngleAndArea(mris, parms, &relevant_angle, &computed_neg_area, &computed_area);
#else
    mrisComputeFaceRelevantAngleAndArea(mris, parms, &relevant_angle, &computed_neg_area, &computed_area);
#endif
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

#ifdef COMPILING_MRIS_MP
#define ELT(NAME, MULTIPLIER, COND, EXPR)     cheapAssert(!(COND))
#define ELTM(NAME,MULTIPLIER,COND,EXPR,EXPRM) double const NAME = (COND) ? (EXPRM) : 0.0;
#else
#define ELT(NAME, MULTIPLIER, COND, EXPR)     double const NAME = (COND) ? (EXPR) : 0.0;
#define ELTM(NAME,MULTIPLIER,COND,EXPR,EXPRM) ELT(NAME, MULTIPLIER, COND, EXPR)
#endif

    SSE_TERMS

#undef ELT
#undef ELTM
#undef SEP

  if (parms->l_thick_spring > 0 || parms->l_thick_min > 0 || parms->l_thick_parallel > 0 /* && DIAG_VERBOSE_ON*/)
    printf("min=%2.3f, parallel=%2.4f, normal=%2.4f, spring=%2.4f, ashburner=%2.3f, tsmooth=%2.3f\n",
           sse_thick_min            / (float)mris->nvertices,
           sse_thick_parallel       / (float)mris->nvertices,
           sse_thick_normal         / (float)mris->nvertices,
           sse_thick_spring         / (float)mris->nvertices,
#ifdef COMPILING_MRIS_MP
           -666.666,
#else
           sse_ashburner_triangle   / (float)mris->nvertices,
#endif
           sse_tsmooth              / (float)mris->nvertices);
           
  double sse_init = 0;

  if (gMRISexternalSSE) {
#ifdef COMPILING_MRIS_MP
    cheapAssert(false);
#else
    sse_init = (*gMRISexternalSSE)(mris, parms);
#endif
  }
  
  double sse = sse_init +
#define SEP
#define ELTM(NAME,MULTIPLIER,COND,EXPR,EXPRM) (MULTIPLIER) * (NAME) +

#ifdef COMPILING_MRIS_MP
#define ELT(NAME, MULTIPLIER,COND,EXPR)
#else
#define ELT(NAME, MULTIPLIER,COND,EXPR)       ELTM(NAME, MULTIPLIER,COND,NotUsed,NotUsed)
#endif

    SSE_TERMS 0.0 ;

#undef ELT
#undef ELTM
#undef SEP

  static int logSSECount, logSSE;
  if (!logSSECount) { logSSE = !!getenv("FREESURFER_logSSE"); }
  logSSECount++;
  
  if (debug || logSSE) {
    double sum = 0;

#define SEP
#define ELTM(NAME, MULTIPLIER, COND, EXPR, EXPRM) fprintf(stdout, "new %s : %f \n", #NAME, (MULTIPLIER) * (NAME));  sum += (MULTIPLIER) * (NAME);

#ifdef COMPILING_MRIS_MP
#define ELT(NAME, MULTIPLIER, COND, EXPR)
#else
#define ELT(NAME, MULTIPLIER, COND, EXPR) ELTM(NAME, MULTIPLIER, NotUsed, NotUsed, NotUsed)
#endif
    ELT(sse_init, 1, true, sse_init)
    SSE_TERMS
    fprintf(stdout, "new sum = %f \n", sum);

#undef ELT
#undef ELTM
#undef SEP
  }
  
  // This code matches code Bevin added to the previous good code to compare old and new runs
  //
  if (false || logSSE) {
    fprintf(stdout, "logSSE:%d \n", logSSECount);
    
    if (parms->l_dist) {
      bool dist_avail  = 
#ifdef COMPILING_MRIS_MP
        !!mris->v_dist[0];
#else
        !!(mris->dist_alloced_flags & 1);
#endif
      #define ELT(X) fprintf(stdout, " %s:%f\n", #X, (float)(X));
      ELT(dist_avail)
      if (dist_avail) {
        VERTEX_TOPOLOGY const * const vt     = &mris->vertices_topology[0];
#ifndef COMPILING_MRIS_MP
        VERTEX          const * const v      = &mris->vertices         [0];
        float           const * const v_dist      = v->dist;
        float           const * const v_dist_orig = v->dist_orig;
#else
        float           const * const v_dist      = mris->v_dist       [0];
        float           const * const v_dist_orig = mris->v_dist_orig  [0];
#endif
        int n;
        for (n = 0; n < vt->vtotal; n++) {
          float const dist_n      = !v_dist      ? 0.0 : v_dist     [n];
          float const dist_orig_n = !v_dist_orig ? 0.0 : v_dist_orig[n];
          ELT(dist_n);
          ELT(dist_orig_n);
        }
      }
      ELT(mris->patch)
      ELT(mris->status)
      ELT(mris->orig_area)
      ELT(mris->total_area)
      ELT(mris->neg_area)
#undef ELT
    }

#define SEP
#define ELTM(NAME, MULTIPLIER, COND, EXPR, EXPRM) \
    { double term = (MULTIPLIER) * (NAME); \
      if (term != 0.0) { fprintf(stdout, "new %s : %f \n", #NAME, term);  } \
    }

#ifdef COMPILING_MRIS_MP
#define ELT(NAME, MULTIPLIER, COND, EXPR)
#else
#define ELT(NAME, MULTIPLIER, COND, EXPR) ELTM(NAME, MULTIPLIER, NotUsed, NotUsed, NotUsed)
#endif
    ELT(sse_init, 1, true, sse_init)
    SSE_TERMS
    fprintf(stdout, "new sum = %f \n", sse);

#undef ELT
#undef ELTM
#undef SEP
  }
  //
  // end of Bevin added

  if (mht_v_current) MHTfree(&mht_v_current);
  if (mht_f_current) MHTfree(&mht_f_current);

  if (!devFinite(sse)) {
    DiagBreak();
  }

#undef COMPUTE_DISTANCE_ERROR

  return sse;
}

#undef COMPUTE_DISTANCE_ERROR
#undef SSE_TERMS
