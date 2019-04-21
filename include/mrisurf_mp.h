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
// The following format provides a denser representation of just the data needed to compute the SSE
//
#include "mrisurf_aaa.h"

struct MRIS_MP {

  // Avoid the 'uninitialized const member' errors in g++ > 4
  MRIS_MP() :
    nvertices(0),
    nfaces(0),
    avg_nbrs(0),
    nsize(0),
    patch(0),
    noscale(0),
    orig_area(0),
    vertices_topology(nullptr),
    faces_topology(nullptr)
    {};

#define SEP
#define ELTX(C,T,N) ELT(C,T,N)

    MRIS* underlyingMRIS;   // allows access in a few rare cases where there is a real benefit 

    MRIS_MP* in_src;        // since the in are not written, they can be shared by copies
    int      in_ref_count;  // check the src doesn't go away
    
  // MRIS
  //
  // In
  //
#define ELT(C,T,N) T C N;

  #define MRIS_MP__LIST_MRIS_IN \
    ELT(const,  int,    nvertices   ) SEP \
    ELT(const,  int,    nfaces      ) SEP \
    ELT(const,  float,  avg_nbrs    ) SEP \
    ELT(const,  int,    nsize       ) SEP \
    ELT(const,  int,    patch       ) SEP \
    ELT(const,  int,    noscale     ) SEP \
    ELT(const,  float,  orig_area   ) SEP \
    ELT(const,  VERTEX_TOPOLOGY const *, vertices_topology) \
    ELTX(const, FACE_TOPOLOGY   const *, faces_topology)

    MRIS_MP__LIST_MRIS_IN
      
  // In out
  #define MRIS_MP__LIST_MRIS_IN_OUT \
    ELT(,       MRIS_Status, status ) SEP \
    ELT(,       double, radius      ) SEP \
    ELT(,       int,    dist_nsize  ) \

    MRIS_MP__LIST_MRIS_IN_OUT
  
  // Out
  #define MRIS_MP__LIST_MRIS_OUT \
    ELT(,       float,   xlo        ) SEP   \
    ELT(,       float,   xhi        ) SEP   \
    ELT(,       float,   ylo        ) SEP   \
    ELT(,       float,   yhi        ) SEP   \
    ELT(,       float,   zlo        ) SEP   \
    ELT(,       float,   zhi        ) SEP   \
    ELT(,       float,   xctr       ) SEP   \
    ELT(,       float,   yctr       ) SEP   \
    ELT(,       float,   zctr       ) SEP   \
    ELT(,       float,   total_area ) SEP   \
    ELT(,       double,  avg_vertex_area ) SEP   \
    ELTX(,      double,  avg_vertex_dist ) SEP   \
    ELT(,       double,  std_vertex_dist ) SEP   \
    ELT(,       float,   neg_orig_area   ) SEP   \
    ELT(,       float,   neg_area        ) 

    MRIS_MP__LIST_MRIS_OUT

#undef ELT

  // Vertices
  //
#define ELT(C,T,N) T C * v_##N;

  // In
  #define MRIS_MP__LIST_V_IN                \
    ELT(const,  char,   ripflag     ) SEP   \
    ELTX(const, int,    VSize       ) SEP   \
    /* following needed for SSE */          \
    ELT(const,  float,  whitex      ) SEP   \
    ELT(const,  float,  whitey      ) SEP   \
    ELT(const,  float,  whitez      ) SEP   \
    ELT(const,  float,  pialx       ) SEP   \
    ELT(const,  float,  pialy       ) SEP   \
    ELT(const,  float,  pialz       ) SEP   \
    ELT(const,  float,  wnx         ) SEP   \
    ELT(const,  float,  wny         ) SEP   \
    ELT(const,  float,  wnz         ) SEP   \
    ELTX(const,  ptr_to_const_float, dist_orig)      /* note: these keep pointing to the original ones in the MRIS - change if code wants to change these values */
    
    MRIS_MP__LIST_V_IN

  // In out
  #define MRIS_MP__LIST_V_IN_OUT_XYZ        \
    ELT(,       float,  x           ) SEP   \
    ELT(,       float,  y           ) SEP   \
    ELT(,       float,  z           )

  #define MRIS_MP__LIST_V_IN_OUT_NOXYZ      \
    ELTX(,      int,    dist_capacity) SEP  \
    ELT(,       char,   border      ) SEP   \
    ELT(,       float,  cx          ) SEP   \
    ELT(,       float,  cy          ) SEP   \
    ELT(,       float,  cz          ) SEP   \
    ELT(,       float,  curv        ) SEP   \
    ELT(,       float,  origarea    ) SEP   \
    ELT(,       int,    assigned_fno) 

  #define MRIS_MP__LIST_V_IN_OUT            \
    MRIS_MP__LIST_V_IN_OUT_XYZ SEP MRIS_MP__LIST_V_IN_OUT_NOXYZ
  
    MRIS_MP__LIST_V_IN_OUT
    
  // Out
  #define MRIS_MP__LIST_V_OUT               \
    ELT(,       float,  area        ) SEP   \
    ELT(,       float,  nx          ) SEP   \
    ELT(,       float,  ny          ) SEP   \
    ELT(,       float,  nz          ) SEP   \
    ELTX(,      char,   neg         ) SEP   \
    ELTX(,      float*, dist        )

    MRIS_MP__LIST_V_OUT

    float**     v_dist_buffer;

#undef ELT

  // Faces
  //
#define ELT(C,T,N) C T * f_##N;
  #define MRIS_MP__LIST_F_IN                                    \
    ELT(const,  char,                   ripflag         ) SEP   \
    ELTX(const, float,                  norm_orig_area  ) SEP   \
    ELTX(const, angles_per_triangle_t,  orig_angle      )
    
    MRIS_MP__LIST_F_IN

  #define MRIS_MP__LIST_F_OUT                                   \
    ELT(,       float,                  area            ) SEP   \
    ELTX(,      char,                   normSet         ) SEP   \
    ELTX(,      FloatXYZ,               norm            ) SEP   \
    ELTX(,      angles_per_triangle_t,  angle           )

    MRIS_MP__LIST_F_OUT
    
#undef ELT

#undef ELTX
#undef SEP

};

void MRISMP_ctr(MRIS_MP* mp);
void MRISMP_dtr(MRIS_MP* mp);
void MRISMP_copy(MRIS_MP* dst, MRIS_MP* src, 
    bool only_inputs,           // if true, copy the in and in_out fields only
    bool ignore_xyz);           // if true, don't copy the v_x[*] v_y[*] v_z[*]     NYI
    
void MRISMP_load(MRIS_MP* mp, MRIS* mris,
  bool loadOutputs,
  float * dx_or_NULL, float * dy_or_NULL, float * dz_or_NULL);      // loaded if not NULL, the dx,dy,dz for ripped set to zero


