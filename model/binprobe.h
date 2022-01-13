/***************************************************************************
							 binprobe.h
                             -------------------
    begin                : Mon Oct 5 2000
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

#ifndef BINPROBE_H
#define BINPROBE_H 

#include "probe.h"
#include "bin.h"
#include "bin2D.h"

namespace MODEL
{

  /** This class collects data for stochastic processes. It takes a 
	  bin as an input (so you can add to it if you want to), and
	  writes out the bin at the end.

	  You just have to define ONE variable of which you would like a
	  histogram, and the ODESystem you are talking about.

	  Also, as always, you need a timeframe and a filename.
  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class BinProbe : public GenericProbe 
  {
  public:
	typedef ODESystem<dims,nelem,NT> system;
	
	BinProbe(TimeFrame& T,
			 const string& n, 
			 Bin& b, 
			 system& sys, 
			 integer v) : GenericProbe (T,n), my_sys(&sys), my_b(&b), var(v) 
	{
	}

	~BinProbe()
	{
	  // Otherwise, nothing gets written
	  write_data();
	}
	
	virtual void print(ostream& out)
	{
	  // Here, we write the data from the bin to the probe
	  vector<number> bi=my_b->get_bins(); 
	  vector<counter> hi=my_b->get_histo();
	
	  for(counter i=0;i<bi.size();i++)
		out << bi[i] << "\t" << hi[i] << endl;  
	}
	
	void probe(void)
	{
	  my_b->add_value(my_sys->get_current()[var]);
	}

	NO_COPY(BinProbe);
	
  private: 
	system* my_sys;
	Bin* my_b;
	integer var;
  };

  /** The same, but averages out during the integration to simulate a
	  detector with finite bandwidth */  

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class AvgBinProbe : public GenericProbe 
  {
  public:
	typedef ODESystem<dims,nelem,NT> system;
	

	/** no is the number of samples to average */

	AvgBinProbe(TimeFrame& T,
			 const string& n, 
			 Bin& b, 
			 system& sys, 
				integer v,integer no) : GenericProbe (T,n),
										my_sys(&sys), my_b(&b),
										var(v),n(no) 
	{
	  runningvalue=0.;
	}

	~AvgBinProbe()
	{
	  // Otherwise, nothing gets written
	  write_data();
	}
	
	virtual void print(ostream& out)
	{
	  // Here, we write the data from the bin to the probe
	  vector<number> bi=my_b->get_bins(); 
	  vector<counter> hi=my_b->get_histo();
	
	  for(counter i=0;i<bi.size();i++)
		out << bi[i] << "\t" << hi[i] << endl;  
	}
	
	void probe(void)
	{
	  // This was just plain silly
	  // const number factor=exp(-1./n);
	  // Should have been 
	  const number factor=1.-1./n;

	  runningvalue=factor*runningvalue+(1.-factor)*my_sys->get_current()[var];
	  my_b->add_value(runningvalue);
	}

	NO_COPY(AvgBinProbe);
	
  private: 
	system* my_sys;
	Bin* my_b;
	integer var;
	integer n;

	number runningvalue;
  };

  // ----------------------------------------------------------------------

  /** These are the equivalent 2D versions */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class BinProbe2D : public GenericProbe 
  {
  public:
	typedef ODESystem<dims,nelem,NT> system;
	
	BinProbe2D(TimeFrame& T,
			 const string& n, 
			 Bin2D& b, 
			 system& sys, 
			   integer vx, integer vy) : GenericProbe (T,n),
										 my_sys(&sys), my_b(&b),
										 varx(vx), vary(vy) 
	{
	}

	~BinProbe2D()
	{
	  // Otherwise, nothing gets written
	  write_data();
	}
	
	virtual void print(ostream& out)
	{
	  // Here, we write the data from the bin to the probe
	  vector<number> bi=my_b->get_bins(); 
	  Bin2D::hist2D hi=my_b->get_histo();
	
	  for(counter i=0;i<bi.size();i++)
		{
		  for(counter j=0;j<bi.size();j++)
			out << hi[i][j] << "\t";
		  out << endl;
		}
	}
	
	void probe(void)
	{
	  my_b->add_value(my_sys->get_current()[varx],my_sys->get_current()[vary]);
	}

	NO_COPY(BinProbe2D);
	
  private: 
	system* my_sys;
	Bin2D* my_b;
	integer varx,vary;
  };

} // end namespace
#endif

/*********************************************************************
$Id: binprobe.h,v 1.2 2002-07-30 19:57:18 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.3  2002/07/24 09:55:34  mpeeters
Better output format.

Revision 1.1.2.2  2002/07/01 13:34:42  mpeeters
Added 2D binning class and probe.

Revision 1.1.2.1  2001/10/04 12:57:20  mpeeters
Changed wrong calculation of 1st order digital filter for averaging.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.1  2000/10/13 11:20:22  mpeeters
Added two new probes (which were in the 0.9 tarball, but which I forgot
to add to the CVS). Removed history.

*********************************************************************/
