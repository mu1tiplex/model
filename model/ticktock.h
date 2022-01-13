/***************************************************************************
                          ticktock.h  -  description
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

#ifndef TICKTOCK_H
#define TICKTOCK_H 

#include "numerictypes.h"
#include "numerictraits.h"
#include "vectorfunction.h"
#include "timeframe.h"
#include "cycler.h"
#include "parameter.h"
#include "cowner.h"
#include <string>
#include <stdexcept>

namespace MODEL {
  
  /** Base mix-in class for ticktocks: anything that should be called when the time
	  changes by one tick: odesystem, probes, etc...
	  \todo Re-evaluate the timing model: there are problems when
	  TickTock get nested and multiple ticktocks are attached to the
	  same system. Each one is then evaluated sequentially for the
	  number of steps needed inside the master timeframe!
 */
  class TickTock : public TimeFrame, public Cycler {
  public:
	TickTock(TimeFrame& T,time resolution);
	
	/** called to make things work, for each one. Points to the step()
		function for the timeframe if not overloaded.
		\todo We have to make sure that each class which is derived
		from this one actually calls step() at some time, otherwise
		you cannot attach other things to probes and stuff (as if you
		would want to.)
	*/

	virtual void tick(){step();}
	
	/** Return the number of times the resolution fits inside the
		masters TimeFrame's */
	const counter& get_resmatch(void){return resmatch;}

	TimeFrame& get_timeframe(void);

  private:
	/** Called by Cowner(In this case, TimeFrame) to make things work */
	void execute(void);
	counter resmatch;
  };
  
} // end namespace
#endif

/*********************************************************************
$Id: ticktock.h,v 1.3 2002-07-30 19:57:19 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2.2.1  2002/07/22 21:24:03  mpeeters
Defined a default implementation for tick(), so the TickTock class can be used as a way to define a timeframe hierarchy.

Revision 1.2  2001/07/26 12:02:24  mpeeters
Removed nasty bug where time would run backwards in certain cases (see todo list).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
