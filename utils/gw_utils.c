/**
 * @file  gw_utils.c
 * @brief miscellaneous utility functions contributed by Graham Wideman
 *
 */
/*
 * Original Author: Graham Wideman
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.5 $
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

#ifndef Darwin
#include <values.h>  // MAXSHORT
#endif
#ifndef MAXSHORT
#define MAXSHORT 32767
#endif
#include <time.h>

#include "error.h"
#include "gw_utils.h"
#include "mrishash_internals.h"

/*-----------------------------------------------
  GWU_make_surface_from_lists
  Adapted from ic642_make_surface()

  Assumes:
  -- Zero-based vertex and face numbering

  -------------------------------------------------*/
MRI_SURFACE *GWU_make_surface_from_lists(GWUTILS_VERTEX *vertices, int vertexcount, GWUTILS_FACE *faces, int facecount)
{
  MRI_SURFACE *mris = MRISalloc(vertexcount, facecount);

  //-----------------------------------------
  // Read vertex data into mris
  //-----------------------------------------
  int vno;
  for (vno = 0; vno < vertexcount; vno++) {
    VERTEX* const v = &mris->vertices[vno];

    v->x = vertices[vno].x;
    v->y = vertices[vno].y;
    v->z = vertices[vno].z;
  }

  //-----------------------------------------
  // Build the faces
  //-----------------------------------------
  int fno;
  for (fno = 0; fno < facecount; fno++) {
    mrisAttachFaceToVertices(mris, fno, 
        faces[fno].vno[0],
        faces[fno].vno[1],
        faces[fno].vno[2]);
  }

  //---------------------------------------------
  // Final housekeeping
  //---------------------------------------------
  MRIScomputeMetricProperties(mris);

  mris->type = MRIS_ICO_SURFACE;
  MRISsetNeighborhoodSizeAndDist(mris, 1);
  return (mris);
}

/*-------------------------------------------------------
  Translate MHT data to an MRI volume in various different ways so
  that it can be visualized.
  Does not attempt to be a proper spatial transform to MRI space,
  (because what to do when mht->vres is large?)
  just a way to read out the MHT data in 3D
  ---------------------------------------------------------*/
MRI *MRIFromMHTandMRIS(MHT *mht, MRIS *mris, MFMM_Option_t mfmm_option)
{
  MRI *amri;
  int mhtvx, mhtvy, mhtvz;  // MHT "voxels"
  int mrix, mriy, mriz;     // MRI voxels
  int binnum;
  // int fno_usage; //
  int outval;

#define HALFMHTFOV 200
#define HALFMRIFOV 128

  // fno_usage = mht->fno_usage;

  amri = MRIalloc(256, 256, 256, MRI_SHORT);
  for (mriz = 0; mriz < 255; mriz++) {
    for (mriy = 0; mriy < 255; mriy++) {
      for (mrix = 0; mrix < 255; mrix++) {
        MRISvox(amri, mrix, mriy, mriz) = 100;
      }
    }
  }
  //  goto done;

  for (mriz = 0; mriz < 255; mriz++) {
    mhtvz = HALFMHTFOV - HALFMRIFOV + mriz;
    for (mriy = 0; mriy < 255; mriy++) {
      mhtvy = HALFMHTFOV - HALFMRIFOV + mriy;
      for (mrix = 0; mrix < 255; mrix++) {
        mhtvx = HALFMHTFOV - HALFMRIFOV + mrix;

        MHBT* bucket = MHTacqBucketAtVoxIx(mht, mhtvx, mhtvy, mhtvz);

        outval = 0;

        if (bucket) {
          if (MFMM_None == mfmm_option) {
            outval = 1;
            goto outval_done;
          }
          for (binnum = 0; binnum < bucket->nused; binnum++) {
            MRIS_HASH_BIN* bin = &(bucket->bins[binnum]);

            switch (mfmm_option) {
              case MFMM_None:
                break;
              case MFMM_Num:
                outval = bin->fno;
                goto outval_done;
              case MFMM_NumDiv16:
                outval = bin->fno >> 4;
                goto outval_done;
              case MFMM_Count:
                outval++;
                break;
            }
          }
      outval_done:
          MHTrelBucket(&bucket);
        }
        if (outval > MAXSHORT) outval = MAXSHORT;

        // MRI?vox is a type-specific macro!
        MRISvox(amri, mrix, mriy, mriz) = outval;

      }  // for mrix
    }    // for mriy
  }      // for mriz
         // done:
  return amri;
}

// Some simple log functions for test programs.

static char local_Progname[500] = "uninitialized";
static char local_Progversion[100] = "uninitialized";
static char local_Logfilepath[1000] = "uninitialized";

//----------------------------------------
int gw_log_init(char *AProgname, char *AProgversion, char *ALogfilepath, int newfile)
{  // 0 for OK
  //----------------------------------------
  FILE *afile;
  int rslt = 0;
  strcpy(local_Progname, AProgname);
  strcpy(local_Progversion, AProgversion);
  strcpy(local_Logfilepath, ALogfilepath);
  if (newfile) {
    afile = fopen(local_Logfilepath, "w");
  }
  else {
    afile = fopen(local_Logfilepath, "a");
  }

  if (afile) {
    fclose(afile);
  }
  else {
    rslt = 1;
  }
  return rslt;
}

//------------------------------
void gw_log_message(char *msg)
{
  //------------------------------
  FILE *afile;
  afile = fopen(local_Logfilepath, "a");
  fprintf(afile, "%s\n", msg);
  fclose(afile);
}

//------------------------------
static void nowstr(char *buf)
{
  //------------------------------
  time_t tim;
  struct tm *tmr;
  // int rslt;

  time(&tim);
  tmr = localtime(&tim);
  // rslt =
  strftime(buf, 100, "%Y-%m-%d %H:%M:%S", tmr);
}

//------------------------------
void gw_log_timestamp(char *label)
{
  //------------------------------
  char datestr[100];
  char msg[200];

  nowstr(datestr);
  sprintf(msg, "---[%s]--- %s version %s at %s", label, local_Progname, local_Progversion, datestr);
  gw_log_message(msg);
}

//------------------------------
void gw_log_begin(void)
{
  //------------------------------
  gw_log_timestamp("Begin");
}

//------------------------------
void gw_log_end(void)
{
  //------------------------------
  gw_log_timestamp("End");
}
