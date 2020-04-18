/*
 * Original Author: Rudolph Pienaar / Christian Haselgrove
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef __SURFACE_H__
#define __SURFACE_H__

#include "mri.h"
#include "mrisurf.h"

void  mark_geodesic(MRIS *surf, int vno_i, int vno_f, int mark);

#endif
