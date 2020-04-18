/**
 * @brief Interactor for volume cropping in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
 *
 */

#ifndef Interactor2DVolumeCrop_h
#define Interactor2DVolumeCrop_h

#include "Interactor2D.h"

class Region2D;

class Interactor2DVolumeCrop : public Interactor2D
{
  Q_OBJECT
public:
  Interactor2DVolumeCrop(QObject* parent);
  virtual ~Interactor2DVolumeCrop();

  // return true if to have parent Interactor2D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view );

protected:
  void UpdateCursor( QEvent* event, QWidget* wnd );

  bool        m_bSelected;
};

#endif


