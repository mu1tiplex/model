/***************************************************************************
                          modulator.h  -  description
                             -------------------
    begin                : Mon Aug 11 2000
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

#ifndef MODULATOR_H
#define MODULATOR_H 

#include "numerictypes.h"
#include "numerictraits.h"
#include "vectorfunction.h"
#include "timeframe.h"
#include "parameter.h"
#include "cycler.h"
#include "ticktock.h"
#include <string>

namespace MODEL {

  class TimeFrame;
  
  /** Base class for modulators.
	  These are objects which update a certain Parameter as a function of
	  time. You can add these to ODESystem, if you want them to be
	  executed for each integrator step, or to another TimeFrame if
	  the granularity isn't that important 
	  @param T the TimeFrame to be added to.
	  @param v the VectorFunction that will supply the Parameter to be
	  modulated.
	  @param param A string describing the Parameter.
	  \todo add a check to see if the Parameter exists.
  */
  class Modulator : public TickTock {
  public:
	template<integer dims, typename nelem, class NT  >
	Modulator(TimeFrame&
			  T,VectorFunction<dims,nelem,NT> & v, string param)
	  : TickTock(T,T.get_dt()) // a modulator has the same resolution !
	{
	  Parameter p=v.get_parameter(param);
	  _p=&p;
	}

	/** called by master to make thing work, for each modulator */
	virtual void tick() {(*_p)=modulate();}

	/** to be written for each modulator */
	virtual number modulate(void)=0;
	
  private:
	number* _p;

  };

  //--------------------------------------------------------------------------------

  /** A linear modulator: ramps up and then down */
  class LinearMod : public Modulator{
  public:
	template<integer dims, typename nelem, class NT  >
	LinearMod( TimeFrame& T,
			   VectorFunction<dims,nelem,NT>& v, 
			   string param,
			   number start,number end, 
			   time from, time until ) : 
	  Modulator(T,v,param), 
	  st(start), 
	  en(end), d(en-st), f (from), u(until/2.),du(u-f) {}

	virtual number modulate();
	  
  private:
	number st;
	number en;
	number d;
	time f;
	time u;
	time du;
  };

  //--------------------------------------------------------------------------------
  /** A block modulator: generates a block wave with period P,
	  starting in the low state*/
  class BlockMod : public Modulator{
  public:
	template<integer dims, typename nelem, class NT  >
	BlockMod( TimeFrame& T,
			   VectorFunction<dims,nelem,NT>& v, 
			   string param,
			   number low,number high, 
			   time period) : 
	  Modulator(T,v,param), 
	  l(low), h(high), p(period){}

	virtual number modulate(void);
	  
  private:
	number l;
	number h;
	time p;
  };

  //-------------------------------------------------------------------------------
  /** A step at a certain time */

  class StepMod : public Modulator {
  public:
	template<integer dims, typename nelem, class NT  >
	StepMod( TimeFrame& T,
			 VectorFunction<dims,nelem,NT>& v, 
			 string param,
			 number start,number end, 
			 time at) 
	  : Modulator(T,v,param), st(start), en(end), a (at) {}

	virtual number modulate();

  private:
	number st;
	number en;
	time a;
  };

  // ------------------------------------------------------------

  /** This class implements a fading triangular modulator. It stays
	  for a number of periods at one amplitude, and then reduces it by
	  2*deltapar (one for the top, one for the bottom boundary ),
	  until we are at 0 */
  class TriangleFadeMod : public Modulator{
  public:
	template<integer dims, typename nelem, class NT  >
	TriangleFadeMod( TimeFrame& T,
			   VectorFunction<dims,nelem,NT>& v, 
			   string param,
			   number start,number end, 
			   time period, number deltapar, integer stay) : 
	  Modulator(T,v,param), 
	  st(start), 
	  en(end), d(en-st), p(period), dp (deltapar), s(stay),cs(1) 
	{ last_t=get_time(); in_t=0; }

	virtual number modulate();
	  
  private:
	number st;
	number en;
	number d;
	time p ;
	number dp ;
	integer s;

	integer cs;
	time last_t;
	time in_t;
	
  };


} // end namespace
#endif

/*********************************************************************
$Id: modulator.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
