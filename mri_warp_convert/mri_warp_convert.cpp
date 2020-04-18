/**
 * @brief A program to convert non-linear deformation field file formats
 *
 */

/*
 * Original Author: Oliver Hinds
 *
 * Copyright © 2016 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <string>
#include <iostream>
#include <fstream>

#include "error.h"
#include "gcamorph.h"
#include "macros.h"
#include "mri.h"
#include "mri_circulars.h"
#include "version.h"

using namespace std;

namespace filetypes {
enum FileType { UNKNOWN, M3Z, FSL, ITK, VOX, RAS };
}

struct Parameters
{
  string in_warp;
  string out_warp;
  string in_src_geom;
  filetypes::FileType in_type;
  filetypes::FileType out_type;
  bool downsample;
};

static struct Parameters P =
{ "", "", "", filetypes::UNKNOWN, filetypes::UNKNOWN, false};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);

static char vcid[] =
    "$Id: mri_warp_convert.cpp,v 1.1 2016/06/16 19:57:06 ohinds Exp $";
const char *Progname = NULL;

GCAM* readM3Z(const string& warp_file)
// Read an m3z file. Just calls down to GCAMread
{
  GCAM* gcam = GCAMread(warp_file.c_str());
  if (gcam == NULL)
  {
    cerr << "ERROR readM3Z: cannot read " << warp_file << endl;
    exit(1);
  }

  return gcam;
}

GCAM* readFSL(const string& warp_file)
// Read in an FSL warp. This is the code that used to reside in
// mri_warp_convert.c.
{
  MRI* mri = MRIread(warp_file.c_str()) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read warp volume %s\n",
              Progname, warp_file.c_str()) ;

  MATRIX* m = MRIgetVoxelToRasXform(mri) ;

  // NOTE: this assumes a standard siemens image orientation in which
  // case a neurological orientation means that the first frame is
  // flipped

  if ( MatrixDeterminant(m) > 0 )
  {
    fprintf(stdout, "non-negative Jacobian determinant -- converting to radiological ordering\n");
  }
  {
    // 2012/feb/08: tested with anisotropic voxel sizes

    MRI *mri2 = NULL ;
    int c=0,r=0,s=0;
    float v;

    mri2 = MRIcopy(mri,NULL);
    for(c=0; c < mri->width; c++)
    {
      for(r=0; r < mri->height; r++)
      {
        for(s=0; s < mri->depth; s++)
        {
          // only flip first frame (by negating relative shifts)
          v = MRIgetVoxVal(mri, c,r,s,0) / mri->xsize;
          if ( MatrixDeterminant(m) > 0 )
            MRIsetVoxVal(    mri2,c,r,s,0,-v);
          else
            MRIsetVoxVal(    mri2,c,r,s,0, v);

          v = MRIgetVoxVal(mri, c,r,s,1) / mri->ysize;
          MRIsetVoxVal(    mri2,c,r,s,1, v);

          v = MRIgetVoxVal(mri, c,r,s,2) / mri->zsize;
          MRIsetVoxVal(    mri2,c,r,s,2, v);

        }
      }
    }
    MRIfree(&mri);
    mri = mri2;

  }
  MatrixFree(&m) ;


  // this does all the work! (gcamorph.c)
  GCA_MORPH* gcam = GCAMalloc(mri->width, mri->height, mri->depth) ;
  GCAMinitVolGeom(gcam, mri, mri) ;

  // not sure if removing singularities is ever a bad thing
#if 1
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri) ;
#else
  GCAMreadWarpFromMRI(gcam, mri) ;
#endif

  return gcam;
}

// Read a warp file containing displacements in RAS or LPS space.
GCAM* read_world(const string& warp_file, const string& src_geom,
    bool is_lps=false)
{
  MRI* in = MRIread( warp_file.c_str() );
  if (in == NULL) {
    cerr << "ERROR: couldn't read input warp from " << warp_file << endl;
    return NULL;
  }
  MRI* src = MRIread( src_geom.c_str() );
  if (src == NULL) {
	cerr << "ERROR: couldn't read source geometry from " << src_geom << endl;
    return NULL;
  }

  GCA_MORPH* out = GCAMalloc(in->width, in->height, in->depth) ;
  GCAMinitVolGeom(out, src, in) ;
  out->type = GCAM_VOX;

  MATRIX* dst_vox2mm = MRIgetVoxelToRasXform(in);
  MATRIX* src_vox2mm = MRIgetVoxelToRasXform(src);
  if (is_lps) {
      MATRIX *ras2lps = MatrixIdentity(4, NULL);
      ras2lps->rptr[1][1] = -1;
      ras2lps->rptr[2][2] = -1;
      dst_vox2mm = MatrixMultiplyD(ras2lps, dst_vox2mm, dst_vox2mm);
      src_vox2mm = MatrixMultiplyD(ras2lps, src_vox2mm, src_vox2mm);
      MatrixFree(&ras2lps);
  }
  MATRIX* src_mm2vox = MatrixInverse(src_vox2mm, NULL);

  VECTOR* src_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR* src_mm = VectorAlloc(4, MATRIX_REAL);
  VECTOR* dst_mm = VectorAlloc(4, MATRIX_REAL);
  VECTOR* wrp_mm = VectorAlloc(4, MATRIX_REAL);

  VECTOR_ELT(wrp_mm, 4) = 0;
  VECTOR_ELT(dst_vox, 4) = 1;
  for(int s=0; s < in->depth; s++) {
    for(int c=0; c < in->width; c++) {
      for(int r=0; r < in->height; r++) {
        GCA_MORPH_NODE* node = &out->nodes[c][r][s];
        node->origx = c;
        node->origy = r;
        node->origz = s;
        node->xn = c;
        node->yn = r;
        node->zn = s;

        VECTOR3_LOAD(dst_vox, c, r, s);
        dst_mm = MatrixMultiplyD(dst_vox2mm, dst_vox, dst_mm);

        VECTOR3_LOAD(wrp_mm, MRIgetVoxVal(in, c, r, s, 0),
            MRIgetVoxVal(in, c, r, s, 1),
            MRIgetVoxVal(in, c, r, s, 2));
        src_mm = VectorAdd(dst_mm, wrp_mm, src_mm);
        src_vox = MatrixMultiplyD(src_mm2vox, src_mm, src_vox);

        node->x = VECTOR_ELT(src_vox, 1);
        node->y = VECTOR_ELT(src_vox, 2);
        node->z = VECTOR_ELT(src_vox, 3);
      }
    }
  }

  MRIfree(&in);
  MRIfree(&src);
  VectorFree(&src_vox);
  VectorFree(&src_mm);
  VectorFree(&wrp_mm);
  VectorFree(&dst_mm);
  VectorFree(&dst_vox);
  MatrixFree(&dst_vox2mm);
  MatrixFree(&src_vox2mm);
  MatrixFree(&src_mm2vox);
  return out;
}

// Read a warp file as displacements in source-voxel space.
GCAM* read_voxel(const string& warp_file, const string& src_geom)
{
  MRI* in = MRIread( warp_file.c_str() );
  if (in == NULL) {
    cerr << "ERROR: couldn't read input warp from " << warp_file << endl;
    return NULL;
  }
  MRI* src = MRIread( src_geom.c_str() );
  if (src == NULL) {
	cerr << "ERROR: couldn't read source geometry from " << src_geom << endl;
    return NULL;
  }

  GCA_MORPH* out = GCAMalloc(in->width, in->height, in->depth) ;
  GCAMinitVolGeom(out, src, in) ;
  out->type = GCAM_VOX;

  MATRIX* dst_vox2ras = MRIgetVoxelToRasXform(in);
  MATRIX* src_vox2ras = MRIgetVoxelToRasXform(src);
  MATRIX* src_ras2vox = MatrixInverse(src_vox2ras, NULL);
  MATRIX* dst2src_vox = MatrixMultiplyD(src_ras2vox, dst_vox2ras, NULL);
  MATRIX* src_vox = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(dst_vox, 4) = 1;

  for(int s=0; s < in->depth; s++) {
    for(int c=0; c < in->width; c++) {
      for(int r=0; r < in->height; r++) {
        GCA_MORPH_NODE* node = &out->nodes[c][r][s];
        node->origx = c;
        node->origy = r;
        node->origz = s;
        node->xn = c;
        node->yn = r;
        node->zn = s;

        VECTOR3_LOAD(dst_vox, c, r, s);
        src_vox = MatrixMultiplyD(dst2src_vox, dst_vox, src_vox);
        node->x = MRIgetVoxVal(in, c, r, s, 0) + VECTOR_ELT(src_vox, 1);
        node->y = MRIgetVoxVal(in, c, r, s, 1) + VECTOR_ELT(src_vox, 2);
        node->z = MRIgetVoxVal(in, c, r, s, 2) + VECTOR_ELT(src_vox, 3);
      }
    }
  }

  MRIfree(&in);
  MRIfree(&src);
  VectorFree(&src_vox);
  VectorFree(&dst_vox);
  MatrixFree(&dst_vox2ras);
  MatrixFree(&src_vox2ras);
  MatrixFree(&src_ras2vox);
  MatrixFree(&dst2src_vox);
  return out;
}

void writeM3Z(const string& fname, GCAM *gcam, bool downsample=false)
// Write an m3z file. Just calls down to GCAMwrite
{
  GCA_MORPH* out = downsample ? GCAMdownsample2(gcam) : gcam;
  GCAMwrite(out, fname.c_str());
  if (downsample) {
      GCAMfree(&out);
  }
}

void writeFSL(const string& fname, const GCAM *gcam)
// Write an FSL warp file.
// NOT IMPLEMENTED
{
  cerr << "ERROR writeFSL is not implemented, sorry!" << endl;
  exit(1);
}

// Write warp as displacements in RAS or LPS space.
void write_world(const string& fname, GCAM* gcam, bool is_lps=false)
{
  MATRIX* dst_vox2mm = VGgetVoxelToRasXform(&gcam->atlas, NULL, 0);
  MATRIX* src_vox2mm = VGgetVoxelToRasXform(&gcam->image, NULL, 0);

  MRI* out = MRIallocSequence(gcam->atlas.width, gcam->atlas.height,
      gcam->atlas.depth, MRI_FLOAT, 3);
  MRIsetResolution(out, gcam->atlas.xsize, gcam->atlas.ysize,
      gcam->atlas.zsize);
  MRIsetVox2RASFromMatrix(out, dst_vox2mm);
  MRIcopyVolGeomToMRI(out, &gcam->atlas);

  if (is_lps) {
      MATRIX *ras2lps = MatrixIdentity(4, NULL);
      ras2lps->rptr[1][1] = -1;
      ras2lps->rptr[2][2] = -1;
      dst_vox2mm = MatrixMultiplyD(ras2lps, dst_vox2mm, dst_vox2mm);
      src_vox2mm = MatrixMultiplyD(ras2lps, src_vox2mm, src_vox2mm);
      MatrixFree(&ras2lps);
  }

  int x, y, z;
  float xw, yw, zw;
  MATRIX* src_vox = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(src_vox, 4) = 1;
  VECTOR_ELT(dst_vox, 4) = 1;
  MATRIX* src_mm = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_mm = VectorAlloc(4, MATRIX_REAL);
  const bool is_same_size = out->width==gcam->width && out->height==gcam->height
      && out->depth==gcam->depth;
  for (x = 0; x < out->width; x++) {
    for (y = 0; y < out->height; y++) {
      for (z = 0; z < out->depth; z++) {
        if (is_same_size) {
            GCA_MORPH_NODE* node = &gcam->nodes[x][y][z];
            xw = node->x;
            yw = node->y;
            zw = node->z;
        }
        else {
            GCAMsampleMorph(gcam, x, y, z, &xw, &yw, &zw);
        }
        VECTOR3_LOAD(dst_vox, x, y, z);
        MatrixMultiplyD(dst_vox2mm, dst_vox, dst_mm);
        VECTOR3_LOAD(src_vox, xw, yw, zw);
        MatrixMultiplyD(src_vox2mm, src_vox, src_mm);

        MRIsetVoxVal(out, x, y, z, 0, VECTOR_ELT(src_mm,1)-VECTOR_ELT(dst_mm,1));
        MRIsetVoxVal(out, x, y, z, 1, VECTOR_ELT(src_mm,2)-VECTOR_ELT(dst_mm,2));
        MRIsetVoxVal(out, x, y, z, 2, VECTOR_ELT(src_mm,3)-VECTOR_ELT(dst_mm,3));
      }
    }
  }

  if (MRIwriteType(out, fname.c_str(), ITK_MORPH) != 0) {
    cerr << "Error writing warp to " << fname << endl;
  }
  MRIfree(&out);
  MatrixFree(&dst_vox2mm);
  MatrixFree(&src_vox2mm);
  MatrixFree(&src_vox);
  MatrixFree(&dst_vox);
  MatrixFree(&src_mm);
  MatrixFree(&dst_mm);
}

// Write a warp file as displacements in source-voxel space.
void write_voxel(const string& fname, GCAM* gcam)
{
  MATRIX* dst_vox2ras = VGgetVoxelToRasXform(&gcam->atlas, NULL, 0);
  MATRIX* src_vox2ras = VGgetVoxelToRasXform(&gcam->image, NULL, 0);
  MATRIX* src_ras2vox = MatrixInverse(src_vox2ras, NULL);
  MATRIX* dst2src_vox = MatrixMultiplyD(src_ras2vox, dst_vox2ras, NULL);

  MRI* out = MRIallocSequence(gcam->atlas.width, gcam->atlas.height,
      gcam->atlas.depth, MRI_FLOAT, 3);
  MRIsetResolution(out, gcam->atlas.xsize, gcam->atlas.ysize,
      gcam->atlas.zsize);
  MRIsetVox2RASFromMatrix(out, dst_vox2ras);
  MRIcopyVolGeomToMRI(out, &gcam->atlas);

  MATRIX* src_vox = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(dst_vox, 4) = 1;
  const bool is_same_size = out->width==gcam->width && out->height==gcam->height
      && out->depth==gcam->depth;
  float x, y, z;
  for (int c = 0; c < out->width; c++) {
    for (int r = 0; r < out->height; r++) {
      for (int s = 0; s < out->depth; s++) {
        if (is_same_size) {
			GCA_MORPH_NODE* node = &gcam->nodes[c][r][s];
			x = node->x;
			y = node->y;
			z = node->z;
		}
		else {
		  GCAMsampleMorph(gcam, c, r, s, &x, &y, &z);
		}
        VECTOR3_LOAD(dst_vox, c, r, s);
        MatrixMultiplyD(dst2src_vox, dst_vox, src_vox);
        MRIsetVoxVal(out, c, r, s, 0, x - VECTOR_ELT(src_vox, 1));
        MRIsetVoxVal(out, c, r, s, 1, y - VECTOR_ELT(src_vox, 2));
        MRIsetVoxVal(out, c, r, s, 2, z - VECTOR_ELT(src_vox, 3));
      }
    }
  }
  if (MRIwrite(out, fname.c_str()) != 0)
  {
    cerr << "Error writing VOX warp to " << fname << endl;
  }
  MRIfree(&out);
  MatrixFree(&dst_vox2ras);
  MatrixFree(&src_vox2ras);
  MatrixFree(&src_ras2vox);
  MatrixFree(&dst2src_vox);
}

int main(int argc, char *argv[])
{
  cout << vcid << endl << endl;

  // Default initialization
  int nargs = handleVersionOption(argc, argv, "mri_warp_convert");
  if (nargs && argc - nargs == 1)
  {
    exit(0);
  }
  argc -= nargs;
  Progname = argv[0];
  argc--;
  argv++;
  ErrorInit(NULL, NULL, NULL);

  // Parse command line
  if (!parseCommandLine(argc, argv, P))
  {
    //printUsage();
    exit(1);
  }

  GCA_MORPH* gcam = NULL;
  bool is_lps = false;
  switch (P.in_type) {
    case filetypes::M3Z:
      gcam = readM3Z(P.in_warp.c_str());
      break;
    case filetypes::FSL:
      gcam = readFSL(P.in_warp.c_str());
      break;
    case filetypes::ITK:
      is_lps = true;
      gcam = read_world(P.in_warp.c_str(), P.in_src_geom, is_lps);
      break;
    case filetypes::VOX:
      gcam = read_voxel(P.in_warp.c_str(), P.in_src_geom);
      break;
    case filetypes::RAS:
      is_lps = false;
      gcam = read_world(P.in_warp.c_str(), P.in_src_geom, is_lps);
      break;
    default:
      ErrorExit(ERROR_BADFILE, "%s: Unknown input type for %s",
                Progname, P.in_warp.c_str());
  }

  if (!gcam)
  {
    ErrorExit(ERROR_BADFILE, "%s: can't read input file %s",
              Progname, P.in_warp.c_str());
  }

  switch (P.out_type) {
    case filetypes::M3Z:
      writeM3Z(P.out_warp.c_str(), gcam, P.downsample);
      break;
    case filetypes::FSL:
      writeFSL(P.out_warp.c_str(), gcam);
      break;
    case filetypes::ITK:
      is_lps = true;
      write_world(P.out_warp.c_str(), gcam, is_lps);
      break;
    case filetypes::VOX:
      write_voxel(P.out_warp.c_str(), gcam);
      break;
    case filetypes::RAS:
      is_lps = false;
      write_world(P.out_warp.c_str(), gcam, is_lps);
      break;
    default:
      ErrorExit(ERROR_BADFILE, "%s: Unknown output type for %s",
                Progname, P.out_warp.c_str());
  }

  GCAMfree(&gcam);
  printf("%s successful.\n", Progname);
  return (0);
}

#include "mri_warp_convert.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(mri_warp_convert_help_xml, mri_warp_convert_help_xml_len);
}

/*!
 \fn int parseNextCommand(int argc, char **argv)
 \brief Parses the command-line for next command
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       number of used arguments for this command
 */
static int parseNextCommand(int argc, char *argv[], Parameters & P)
{
  bool have_input = false;
  bool have_output = false;

  int nargs = 0;
  char *option;

  option = argv[0] + 1;                     // remove '-'
  if (option[0] == '-')
  {
    option = option + 1;  // remove second '-'
  }
  StrUpper(option);

  if (!strcmp(option, "INM3Z") )
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::M3Z;
    nargs = 1;
    cout << "--inm3z: " << P.in_warp << " input M3Z warp." << endl;
  }
  else if (!strcmp(option, "INFSL"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::FSL;
    nargs = 1;
    cout << "--infsl: " << P.in_warp << " input FSL warp." << endl;
  }
  else if (!strcmp(option, "INITK") || !strcmp(option, "INLPS"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::ITK;
    nargs = 1;
    cout << "--inlps: " << P.in_warp << " input LPS warp." << endl;
  }
  else if (!strcmp(option, "INVOX"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::VOX;
    nargs = 1;
    cout << "--invox: " << P.in_warp << " input VOX warp." << endl;
  }
  else if (!strcmp(option, "INRAS"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::RAS;
    nargs = 1;
    cout << "--inras: " << P.in_warp << " input RAS warp." << endl;
  }
  else if (!strcmp(option, "OUTM3Z") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::M3Z;
    nargs = 1;
    cout << "--outm3z: " << P.out_warp << " output M3Z." << endl;
  }
  else if (!strcmp(option, "OUTFSL") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::FSL;
    nargs = 1;
    cout << "--outfsl: " << P.out_warp << " output FSL warp." << endl;
  }
  else if (!strcmp(option, "OUTITK") || !strcmp(option, "OUTLPS"))
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::ITK;
    nargs = 1;
    cout << "--outlps: " << P.out_warp << " output LPS warp." << endl;
  }
  else if (!strcmp(option, "OUTVOX") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::VOX;
    nargs = 1;
    cout << "--outvox: " << P.out_warp << " output VOX warp." << endl;
  }
  else if (!strcmp(option, "OUTRAS") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::RAS;
    nargs = 1;
    cout << "--outras: " << P.out_warp << " output RAS warp." << endl;
  }
  else if (!strcmp(option, "INSRCGEOM") || !strcmp(option, "G"))
  {
    P.in_src_geom = string(argv[1]);
    nargs = 1;
    cout << "--insrcgeom: " << P.in_src_geom
         << " source image (geometry)." << endl;
  }
  else if (!strcmp(option, "DOWNSAMPLE") || !strcmp(option, "D"))
  {
    if (!P.downsample)
        cout << "--downsample: save M3Z at half resolution." << endl;
    P.downsample = true;
    nargs = 0;
  }
  else if (!strcmp(option, "HELP") )
  {
    printUsage();
    exit(1);
  }
  else
  {
    cerr << endl << endl << "ERROR: Option: " << argv[0]
         << " unknown (see --help) !! " << endl << endl;
    exit(1);
  }

  fflush(stdout);

  return (nargs);
}

/*!
 \fn int parseCommandLine(int argc, char **argv)
 \brief Parses the command-line
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       if all necessary parameters were set
 */
static bool parseCommandLine(int argc, char *argv[], Parameters & P)
{
  int nargs;
  int inputargs = argc;
  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv, P);
    argc -= nargs;
    argv += nargs;
  }

  if (inputargs == 0)
  {
    printUsage();
    exit(1);
  }

  bool need_geom;
  switch (P.in_type) {
    case filetypes::ITK:
    case filetypes::RAS:
    case filetypes::VOX:
        need_geom = true;
        break;
    default:
        need_geom = false;
  }
  if (P.in_src_geom.empty() && need_geom) {
    cerr << endl << endl << "ERROR: specified input warp requires --insrcgeom"
         << endl << endl;
    return false;
  }

  if (P.out_type != filetypes::M3Z && P.downsample) {
    cerr << endl << endl;
    cerr << "ERROR: --downsample flag only valid for output type M3Z"
         << endl << endl;
    return false;
  }

  return true;
}
