// included into mrisurf_project.c as a template
  
static void FUNCTION_NAME(
#ifdef COMPILING_MRIS_MP
    MRIS_MP* mris, 
#else
    MRIS*    mris, 
#endif
    double r) {

  if (FZERO(r)) {
    r = DEFAULT_RADIUS;
  }

  mris->radius = r;

  bool const totalDistNeeded = (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON;
  double total_dist = 0.0;
  
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {

#ifndef COMPILING_MRIS_MP
    VERTEX *v = &mris->vertices[vno];
    if (v->ripflag) continue;
    double x = (double)v->x;
    double y = (double)v->y;
    double z = (double)v->z;
#else
    if (mris->v_ripflag[vno]) continue;
    double x = (double)mris->v_x[vno];
    double y = (double)mris->v_y[vno];
    double z = (double)mris->v_z[vno];
#endif

    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;
    double dist = sqrt(x2 + y2 + z2);
    
    double d;
    if (FZERO(dist)) {
      d = 0;
    }
    else {
      d = 1 - r / dist;
    }
    
    double dx = d*x;
    double dy = d*y;
    double dz = d*z;
    
    x = x - dx;
    y = y - dy;
    z = z - dz;

    if (!isfinite(x) || !isfinite(y) || !isfinite(z)) {
      DiagBreak();
    }

#ifndef COMPILING_MRIS_MP
    v->x = x;
    v->y = y;
    v->z = z;
#else
    mris->v_x[vno] = x;
    mris->v_y[vno] = y;
    mris->v_z[vno] = z;
#endif

    if (!totalDistNeeded) continue;

    dist = sqrt(dx*dx + dy*dy + dz*dz);
    total_dist += dist;
  }

  if (totalDistNeeded) {
    fprintf(stdout, "sphere_project: total dist = %f\n", total_dist);
  }
  
#ifndef COMPILING_MRIS_MP
  MRISupdateEllipsoidSurface(mris);
#else
  MRISMP_updateEllipsoidSurface(mris);
#endif
  
  mris->status = (mris->status == MRIS_PARAMETERIZED_SPHERE) ? MRIS_PARAMETERIZED_SPHERE : MRIS_SPHERE;
}


#undef FUNCTION_NAME
