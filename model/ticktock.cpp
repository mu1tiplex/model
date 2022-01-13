/***************************************************************************
                          ticktock.cpp  -  description
                             -------------------
    begin                : Tue Aug 12 2000
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

#include "ticktock.h"

namespace MODEL 
{
  TickTock::TickTock(TimeFrame& T,time resolution) :
	TimeFrame(resolution),Cycler(T),resmatch(0) 
  {
	// This recalibrates the resolution to a whole fraction of the timeframe
	time big_dt=T.get_dt();
	resmatch=counter(big_dt/dt);
	dt=(big_dt/number(resmatch));
  }

  TimeFrame& 
  TickTock::get_timeframe(void)
  {
	return *dynamic_cast<TimeFrame*>(get_boss());
  }

  void 
  TickTock::execute(void)
  {
	// Ooooooh you dirty ... !
	// But of course it is the ticktock which is responsible,
	// as it can only be added to a TimeFrame. The TimeFrame,
	// on the other hand, does not know what kind of objects
	// is it calling
	t = dynamic_cast<TimeFrame*>(get_boss())->get_time();

	for(counter repeat=get_resmatch();repeat>0;--repeat)
	  {
		tick();
		// Update time - don't worry:local time only
		// t+=dt;
		/** \bug We are programming by coincidence here */
	  }
  }
  
}

/*********************************************************************
$Id: ticktock.cpp,v 1.2 2001-07-26 12:02:24 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
