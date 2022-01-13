/***************************************************************************
                          timeframe.h  -  description
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

#ifndef TIMEFRAME_H
#define TIMEFRAME_H

#include "numerictypes.h"
#include "cowner.h"
#include "utility.h"

namespace MODEL {

  /**This class handle timeframes and provides a general way
	 to keep track of the time dimension in problems (which
	 is quite often privileged and separate from the others)

	 \bug Sometimes the system jumps back, writing out the wrong
	 timestamp. Probably this is related to multiple cyclers doing the
	 same thing. 
	 */
  class TickTock;

  class TimeFrame : public COwner
  {
	friend class TickTock;
  public:
	/** resolutions */
	static const	time 	default_dt=1.E-3; 
	
	TimeFrame(time resolution=default_dt, time t_init=0.);
	
	void	reset(void);
	
	time	get_time(void) {return t;}
	void	set_time(time newt) {t=newt;}
		
	time	get_dt(void) {return dt;}
	void	set_dt(time newdt) {dt=newdt;}
		
	/** Increase time with a full step.
		Take one of the cyclers and run
		it thru its paces, calling cycle every time. Then
		reset the time, and do the same for the next one, until all
		are happy.
	 */
	const time&	step();

	/** Increase time with a full step */
	const time&	operator++() {return step();}

	/** to make so we can use it transparently
		time is a real is a double */
	operator time() {return t;}

	/** Remove one of the cyclers */
	void stop(Cycler& t){remove(t);}

	/** TimeFrame should never be copied */
	NO_COPY(TimeFrame);

  protected:
	time	t;
	time	dt;

  };

  /** The global timeframe */
  extern TimeFrame	Universal;	

} // MODEL
#endif

/*********************************************************************
$Id: timeframe.h,v 1.2 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
