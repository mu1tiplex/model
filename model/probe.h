/***************************************************************************
                          probe.h  -  description
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

#ifndef PROBE_H
#define PROBE_H 

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "odesystem.h"
#include "numerictraits.h"

namespace MODEL 
{
  
  /** The base class for probes. Probe itself is already more
	  specialized and expects the routine to write a number of data
	  points on every call */

  class GenericProbe : public TickTock
  {
  public:
	/** Without ostream */
	GenericProbe(TimeFrame& T, const string& fn);

	/** With ostream */
	GenericProbe(TimeFrame & T, ostream& os);

	~GenericProbe(){}

	virtual void tick();

	virtual void probe(void)=0;

	/** Call this one in the destructor */
	void write_data(void);

	/** Implement this one to write the data. You can call it upon
		destruction of your derived object. Unfortunately, I have not
		yet figured out how to call a derived function inside a
		destructor (as a general rule, this is not possible: it will
		always call the local version. DUH. So it is implemented so
		people would be reminded */
	virtual void print(ostream& out)=0;

	/** Don't copy or assign */
	NO_COPY(GenericProbe);
	
  protected:
	string n;
	ostream* out;
  };
  
  /** Class to probe variables and parameters.
	  It writes out the collected data into a file upon destruction.

	  You just have to redefine probe() to write out the needed data
	  using add_data(). You can add as many datapoints as you like at
	  one time (and the number does not have to be constant).
  */

  class Probe : public TickTock
  {
  public:
	typedef vector<number> data_record;
	typedef vector< data_record > data_list;

	/** Constructor if we don't have an ostream handy.
		@param T the timeframe to which the probe should be added.
		@param fn the name of the file to which the data is written.
		Oh, a probe has the same resolution as the timeframe...
	*/
	Probe(TimeFrame & T, const string& fn);

	/** Constructor with ostream 
		@param T the timeframe to which the probe should be added.
		@param os the ostream (e.g. cout)
	*/
	Probe(TimeFrame & T, ostream& os);

	virtual ~Probe();	

	virtual void tick();

	virtual void probe(void)=0;
	
	void add_data(const number& d);
	
	/** Don't copy or assign */
	NO_COPY(Probe);
	
  private:
	string n;
	ostream* out;

	bool data_is_here;
	data_list* var_data;
	data_record* temp_data;
  };
  
  //------------------------------------------------------------

  /** A probe which records everything of note in a ODESystem */
template<integer dims, typename nelem=number, class NT =
NumericTraits<nelem,dims> >

class ODEProbe : public Probe 
  {
  public:
	typedef ODESystem<dims,nelem,NT> system;

	ODEProbe(TimeFrame& T, system& s, const string& n) 
	  : Probe(T,n), my_sys(&s){}

	ODEProbe(TimeFrame& T, system& s, ostream& o) 
	  : Probe(T,o), my_sys(&s){}

	virtual void probe(void)
	{
	  typename NT::vect current(my_sys->get_current());
	  for(counter i=0;i<dims;i++) add_data(current[i]);
	}
  
  private:
	system* my_sys;
};

  /** A probe which records a specific in a ODESystem,
	  averaging out ove a certain timescale */
template<integer dims, typename nelem=number, class NT =
NumericTraits<nelem,dims> >

class AvgProbe : public Probe 
  {
  public:
	typedef ODESystem<dims,nelem,NT> system;

	/** Constructor for an exponentially averaging probe.

		\param T The timeframe to which the probe is
		connected, at which specifies the rate at which the probe
		samples
		\param s The system where we will sample info from
		\param n The filename to write the output to
		\param no The variable to sample from
		\param avg The number of samples to average over (3dB point) 
		\param out The rate at which to output

		All the times and rates are expresses as multiples of the
		timeframe resolution.
	*/

	AvgProbe(TimeFrame& T, system& s, const string& n, integer no,
			 integer avg, integer out = 1) 
	  : Probe(T,n), my_sys(&s),whichone(no),n(avg),o(out),ocount(0),runningvalue(0){}

	AvgProbe(TimeFrame& T, system& s, ostream& o, integer no, integer
			 avg) 
	  : Probe(T,o), my_sys(&s), whichone(no), n(avg){}

	virtual void probe(void)
	{
	  // This was just plain silly
	  // const number factor=exp(-1./n);
	  // Should have been 
	  const number factor=1.-1./n;

	  runningvalue=factor*runningvalue+(1.-factor)*my_sys->get_current()[whichone];

	  if(++ocount==o){
		add_data(runningvalue);
		ocount=0;
	  }
	}
  
  private:
	system* my_sys;
	integer whichone;
	integer n;
	integer o;
	integer ocount;
	
	number runningvalue;
};
  
} // end namespace
#endif

/*********************************************************************
$Id: probe.h,v 1.2 2002-07-30 19:57:19 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.4  2002/02/06 12:49:39  mpeeters
Altered interface to make fudging easier. Will be reverted.

Revision 1.1.2.3  2001/10/04 12:57:20  mpeeters
Changed wrong calculation of 1st order digital filter for averaging.

Revision 1.1.2.2  2001/09/03 12:29:17  mpeeters
Bugfix: == should have been =. Bummer.

Revision 1.1.2.1  2001/09/03 10:06:07  mpeeters
Added parameter to allow for subsampling of output.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.4  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.3  2000/10/09 09:43:44  mpeeters
Added GenericProbe, because Probe was too specialized already.
Probe now only logs data when add_data is called, not at every timepoint.
This caused a huge memory leak, which has been shut now.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
