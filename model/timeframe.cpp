/***************************************************************************
                          timeframe.cpp  -  description
                             -------------------
    begin                : Tue May 9 2000
    copyright            : (C) 2000 by Michael Peeters
    email                : Michael.Peeters@vub.ac.be
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "timeframe.h"

namespace MODEL{

  /** Define a global timeframe */
  TimeFrame	Universal;

  TimeFrame::TimeFrame(time resolution, time t_init) 
	:  t(t_init),dt(resolution) 
  {
  }

  void	
  TimeFrame::reset(void) 
  {
	t=0.;
  }

  const time&	
  TimeFrame::step() 
  {
	time now(t);
	execute_list(); // Calls cycle for each cylcer
	t=now+dt; // Unecessary precaution. Otherwise we might end up too far.
	return t;
  }

}

/*********************************************************************
$Id: timeframe.cpp,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
