/***************************************************************************
                          switchprobe.h  -  description
                             -------------------
    begin                : Thu 22 Feb 2001
    copyright            : (C) 2001 by Michael Peeters
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

#ifndef SWPROBE_H
#define SWPROBE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "odesystem.h"
#include "numerictraits.h"

namespace MODEL 
{
  
  /** A class of probes, very fast, which just look if a variable is
	  greater than a given threshold. If so, they return true and keep
	  the time at which is happenenen */
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class SwitchProbe : public TickTock
  {
  public:	
	typedef ODESystem<dims,nelem,NT> system;
	SwitchProbe(system& sys, integer var, number threshold) :
	  TickTock(sys,sys.get_dt()),s(&sys), v(var),th(threshold) 
	{
	  reset();
	}
	
	void tick()
	{
	  if(!switched)
		if (s->get_current()[v] >th) 
		  {
		  switched=true;
		  switcht=get_time();
		}
	}

	const bool check(void)
	{
	  return switched;
	}
	
	const number switch_time(void)
	{
	  return switcht;
	}
	

	void reset()
	{
	  switched=false;
	}

	NO_COPY(SwitchProbe);	

  private:
	system* s;
	integer v;
	number th;
	
	bool switched;
	number switcht;
  };
  
}

#endif
