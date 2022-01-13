/***************************************************************************
			modulator.cpp
			-----------
                             
    begin                : Fri Sep 15 12:22:12 CEST 2000
    author               : (C) 2000 by Michael Peeters
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

#include "modulator.h"

namespace MODEL 
{
  number 
  LinearMod::modulate()
  {
	time t=get_time();
	  
	if (t<f) return st;
	if ((t>f)&(t<u)) {
	  time pos=(t-f)/du; // where in the period from 0 to 1 am I ?
	  return st+d*pos;
	}
	if (t>u) {
	  time pos=(t-f-u)/du;
	  return en-d*pos;
	}
	
	return en;
  }

  number 
  TriangleFadeMod::modulate()
  {
	// We have faded out
	if (st>en) return 0.;
	
	time t=get_time();
	time dt=t-last_t;
	last_t=t;
	
	in_t += dt;
	
	// Not an easy task: first check if we crossed a timeperiod boundary
	// We assume this never happens twice in one timestep
	if ( in_t>= (s*p) )
	  {
		in_t-=s*p;
		st +=dp;
		en -=dp;
		d=en-st;
		if (st>en) return 0.;
	  }
	
	time du=p/2.;
	time in_p = in_t-integer(in_t/p)*p;

	if (in_p <= du) {
	  time pos=in_p/du; // where in the period from 0 to 1 am I ?
	  return st+d*pos;
	}
	if (in_p > du) {
	  time pos=in_p/du - 1.;
	  return en-d*pos;
	}
	
	return 0.;
  }

  number 
  BlockMod::modulate(void)
  {
	time t=get_time();

	time v=integer(t/p);
	time now=(t-p*v)/p;
	  
	if (now<0.5) return l; 
	else return h;
  }

  number StepMod::modulate()
  {
	time t=get_time();
	if (t<a) return st; else return en;
  }

} // end namespace

/*********************************************************************
$Id: modulator.cpp,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
