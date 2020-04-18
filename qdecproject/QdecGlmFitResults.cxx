/**
 * @brief Contains the results of a GLM fit run.
 *
 * The bulk of the result data is stored in files on disk, and this object
 * contains paths to that data, as well as some local results.
 */
/*
 * Original Author: Nick Schmansky
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

#include "QdecGlmFitResults.h"


// Constructors/Destructors
//

QdecGlmFitResults::QdecGlmFitResults 
( QdecGlmDesign* iGlmDesign,
  vector< string > iContrastSigFiles,  /* /<contrast>/sig.mgh */
  string iConcatContrastSigFile,       /* contrast.sig.mgh */
  string ifnResidualErrorStdDevFile,   /* rstd.mgh */
  string ifnRegressionCoefficientsFile,/* beta.mgh */
  string ifnFsgdFile                   /* y.fsgd */ )
{
  assert( iGlmDesign );
  assert( iContrastSigFiles.size() );

  this->mGlmDesign = iGlmDesign;
  this->mfnContrastSigFiles = iContrastSigFiles;
  this->mfnConcatContrastSigFile = iConcatContrastSigFile;
  this->mfnResidualErrorStdDevFile = ifnResidualErrorStdDevFile;
  this->mfnRegressionCoefficientsFile = ifnRegressionCoefficientsFile;
  this->mfnFsgdFile = ifnFsgdFile;
}

QdecGlmFitResults::~QdecGlmFitResults ( )
{ }

//
// Methods
//


/**
 * Returns the design object used as input to the GLM fitter
 * @return QdecGlmDesign*
 */
QdecGlmDesign*  QdecGlmFitResults::GetGlmDesign ( )
{
  return this->mGlmDesign;
}


/**
 * Returns the names given to the contrast results produced by glmfit.
 * Example of one of the possible names: "Avg-thickness-Age-Cor"
 * @return vector< string >
 */
vector< string > QdecGlmFitResults::GetContrastNames ( )
{
  return this->mGlmDesign->GetContrastNames();
}


/**
 * Returns the human-readable questions associated with each contrast.
 * Example of one question:
 * "Does the correlation between thickness and Age differ from zero?".
 * @return vector< string >
 */
vector< string > QdecGlmFitResults::GetContrastQuestions ( )
{
  return this->mGlmDesign->GetContrastQuestions();
}


/**
 * Returns pathname to the concatenated contrast significance file, 
 * ie sig.mgh for all contrasts.
 * @return string
 */
string QdecGlmFitResults::GetConcatContrastSigFile ( )
{
  return this->mfnConcatContrastSigFile;
}

/**
 * Returns pathnames to the contrast significance file, ie sig.mgh for that
 * contrast.
 * @return vector< string >
 */
vector< string > QdecGlmFitResults::GetContrastSigFiles ( )
{
  return this->mfnContrastSigFiles;
}


/**
 * Returns pathnames to the contrast gamma file, ie gamma.mgh for
 * that contrast.
 * @return vector< string >
 */
vector< string > QdecGlmFitResults::GetContrastGammaFiles ( )
{
  vector<string> tmp;
  return tmp; //TODO
}


/**
 * Returns pathnames to the contrast F-test file, ie F.mgh for that contrast.
 * @return vector< string >
 */
vector< string > QdecGlmFitResults::GetContrast_F_Files ( )
{
  vector<string> tmp;
  return tmp; //TODO
}


/**
 * Returns pathname to the beta.mgh file.
 * @return string
 */
string QdecGlmFitResults::GetRegressionCoefficientsFile ( )
{
  return this->mfnRegressionCoefficientsFile;
}


/**
 * Returns pathname to eres.mgh
 * @return string
 */
string QdecGlmFitResults::GetResidualErrorFile ( )
{
  return this->mfnResidualErrorFile;
}


/**
 * Returns pathname to rstd.mgh
 * @return string
 */
string QdecGlmFitResults::GetResidualErrorStdDevFile ( )
{
  return this->mfnResidualErrorStdDevFile;
}


/**
 * Returns pathname to y.fsgd
 * @return string
 */
string QdecGlmFitResults::GetFsgdFile ( )
{
  return this->mfnFsgdFile;
}

