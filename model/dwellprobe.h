/***************************************************************************
							 dwellprobe.h
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

#ifndef DWELLPROBE_H
#define DWELLPROBE_H 

#include "probe.h"
#include "bin.h"

namespace MODEL
{

  /** This class collects dwell times. You provide ONE variable for
      which you want dwell times, and the ODESystem you are talking
      about.

      It then just compares the value to the previous value and the
      (provided) threshold level and determines if a switch has
      happened since last time. If so, it writes out a data point. 

      This one uses Probe, not GenericProbe.

      Also, as always, you need a timeframe and a filename.
  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
    class DwellProbe : public Probe 
    {
      public:
      typedef ODESystem<dims,nelem,NT> system;
	
      DwellProbe(TimeFrame& T,
		 const string& n, 
		 system& sys, 
		 integer v,
		 number threshold, number range) : Probe (T,n), my_sys(&sys),
      var(v),
      th(threshold), ra(range) 
      // range is the %of the threshold above or below to have a
      // definite switch - this is to avoid random paths thru the
      // middle, without a real switch taking place
      {
	// Set inital state
	more_or_less = (my_sys->get_current()[var] > th);
	last_t=T;
	avgn=0.;
      }

      DwellProbe(TimeFrame& T,
		 ostream& o, 
		 system& sys, 
		 integer v,
		 number threshold, number range) : Probe (T,o), my_sys(&sys),
      var(v),
      th(threshold), ra(range) 
      // range is the %of the threshold above or below to have a
      // definite switch - this is to avoid random paths thru the
      // middle, without a real switch taking place
      {
	// Set inital state
	more_or_less = (my_sys->get_current()[var] > th);
	last_t=T;
	avgn=0.;
      }

      const number& switches(void)
      {
	return avgn;
      }

      void probe(void)
      {
	bool now_more = (my_sys->get_current()[var] > (th*(1.+ra)));
	bool now_less = (my_sys->get_current()[var] < (th*(1.-ra)));
	  
	if(now_more!=now_less) // both false means in the middle
	if(now_more != more_or_less) 
	{
	  avgn++;
	  number dt=get_time()-last_t;
	  add_data(dt);
	  add_data(now_more?1.:-1.);
	  add_data(my_sys->get_current()[var]);

	  more_or_less=now_more;
	  last_t=get_time();
	}
      }

      NO_COPY(DwellProbe);
	
      private: 
      system* my_sys;
      integer var;
      number th;
      number ra;
	
      number avgn;

      bool more_or_less;
      number last_t;
	
    };

  /** This class collects bins dwell times (as described in the
      DwellProbe class). You provide ONE variable for
      which you want dwell times, and the ODESystem you are talking
      about.
  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
    class BinDwellProbe : public GenericProbe 
    {
      public:
      typedef ODESystem<dims,nelem,NT> system;
	
      BinDwellProbe(TimeFrame& T,
		    const string& n, Bin& b_up, Bin& b_down,
		    system& sys, 
		    integer v,
		    number threshold) : GenericProbe (T,n), my_sys(&sys),
      var(v), th(threshold),  
      my_b_up(&b_up), my_b_down(&b_down) 
      {
	// Set inital state
	more_or_less = (my_sys->get_current()[var] > th);
	last_t=T;

	pdfwrite=false;
      }

      /** call this with TRUE to make the output a pdf */
      void write_pdf(bool pdf=true){pdfwrite=pdf;}

      ~BinDwellProbe()
      {
	write_data();
      }
	
      virtual void print(ostream& out)
      {
	if (!pdfwrite){
	// Here, we write the data (i.e. the number of hits per bin) from the bin to the probe
	// First the up-data
	out << "# Dwell statistics for up state" << endl;
	  
	vector<number> bi=my_b_up->get_bins(); 
	vector<counter> hi=my_b_up->get_histo();
	
	for(counter i=0;i<bi.size();i++)
	out << bi[i] << "\t" << hi[i] << endl;  

	// Then the down-data; extra endl fo Gnuplot
	out << endl << endl << "# Dwell statistics for down state" << endl;

	bi=my_b_down->get_bins(); 
	hi=my_b_down->get_histo();
	
	for(counter i=0;i<bi.size();i++)
	out << bi[i] << "\t" << hi[i] << endl;
	}
	else // We have to write out the pdf
	{
	// Here, we write the data (i.e. the probability density in the middle of each bin ) from the bin to the probe
	// First the up-data
	out << "# Dwell statistics for up state" << endl;
	  
	vector<number> bi=my_b_up->get_bins(); 
	vector<number> hi=my_b_up->get_pdf();
	
	for(counter i=0;i<bi.size();i++)
	out << bi[i] << "\t" << hi[i] << endl;  

	// Then the down-data; extra endl fo Gnuplot
	out << endl << endl << "# Dwell statistics for down state" << endl;

	bi=my_b_down->get_bins(); 
	hi=my_b_down->get_pdf();
	
	for(counter i=0;i<bi.size();i++)
	out << bi[i] << "\t" << hi[i] << endl;  
	}  
      }

      void probe(void)
      {
	bool now_more = (my_sys->get_current()[var] > th);
	if(now_more != more_or_less) 
	{
	  number dt=get_time()-last_t;
	  avgn++;
		  
	  if(more_or_less) //bugfix 
	    my_b_up->add_value(dt);
	  else
	    my_b_down->add_value(dt);

	  more_or_less=now_more;
	  last_t=get_time();
	}
      }

      NO_COPY(BinDwellProbe);
	
      private: 
      system* my_sys;
      integer var;
      number th;
      number avgn;

      bool pdfwrite;
	
      Bin* my_b_up;
      Bin* my_b_down;
	
      bool more_or_less;
      number last_t;
	
    };

  /** AVGDwellProbe writes out the current running average of the dwell
      time; currently only for up switches (bottom state).  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
    class AVGDwellProbe : public Probe 
    {
      public:
      typedef ODESystem<dims,nelem,NT> system;
	
      /** By default, a 2x10% width (around the threshold) is assumed */
      AVGDwellProbe(TimeFrame& T,
		    const string& n, 
		    system& sys, 
		    integer v,
		    number threshold,
		    number width=0.1,bool below=true) : Probe (T,n), my_sys(&sys),
      var(v),
      th(threshold),
      w(width), lower(below)
      {
	// Set inital state
	/** \bug: assumed to be outside of the threshold */
	if(lower)
	  more_or_less = (my_sys->get_current()[var] > (th+w))?1:-1;
	else
	  more_or_less = (my_sys->get_current()[var] < (th-w))?-1:1;
	last_t=T;
	avg=0.;avgn=0.;
      }

      const number& switches(void)
      {
	return avgn;
      }

      const number& average(void)
      {
	return avg;
      }

      void probe(void)
      {
	bool more = (my_sys->get_current()[var] > (th+w));
	bool less = (my_sys->get_current()[var] < (th-w));
	integer now_more=more?1:(less?-1:0);		
	if(now_more != more_or_less) // if we have a change
	{
	  if(lower){			
	    /** \todo  We decouple the two kinds of probes (u & down). Agreed:
		this is wasteful and in the future I'll do the right
		thing) */
	    if (now_more==0) // inside DMZ
	      {
		// how did we enter
		from_top=(more_or_less==1);
	      }
	    if ((now_more==1) && !from_top) // and it is an up switch
	      {  
		number dt=get_time()-last_t; // get the time when we
		// swiched down
		avgn++;
		avg=((avgn-1.)*avg + dt)/avgn; // calculate average
		add_data(avg); // add the average
	      }
		  
	    if ((now_more==-1) && from_top) // down switch
	      {
		last_t=get_time(); // keep this point
	      }
	    more_or_less=now_more; // keep the current state 
	  }
	  else{
	    if (now_more==0) // inside DMZ
	      {
		// how did we enter
		from_bottom=(more_or_less==-1);
	      }
	    if ((now_more==-1) && !from_bottom) // and it is a down switch
	      {  
		number dt=get_time()-last_t; // get the time when we
		// swiched up
		avgn++;
		avg=((avgn-1.)*avg + dt)/avgn; // calculate average
		add_data(avg); // add the average
	      }
		  
	    if ((now_more==1) && from_bottom) // up switch
	      {
		last_t=get_time(); // keep this point
	      }
	    more_or_less=now_more; // keep the current state
	  }
	}
      }
      NO_COPY(AVGDwellProbe);
	
      private: 
      system* my_sys;
      integer var;
      number th;
      number w;

      number avg;
      number avgn;
	
      /** Tells you where the system is :
	  1) above threshold+w
	  -1) below threshold-w
	  0) inside DMZ	
      */
      integer more_or_less;
      bool from_top, from_bottom;
      bool lower;
      number last_t;
	
    };

  /** AVGAVGDwellProbe writes out the current running average of the dwell
      time for and averaged system; currently only for up switches (bottom state).  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
    class AVGAVGDwellProbe : public Probe 
    {
      public:
      typedef ODESystem<dims,nelem,NT> system;
	
      /** By default, a 2x10% width (around the threshold) is assumed */
      AVGAVGDwellProbe(TimeFrame& T,
		       const string& n, 
		       system& sys, 
		       integer v,
		       number threshold,
		       number width,
		       integer no,
		       bool below=true) : Probe (T,n), my_sys(&sys),
      var(v),
      th(threshold),
      w(width*threshold), lower(below),n(no)
      {
	// Set inital state
	runningvalue=my_sys->get_current()[var];

	/** \bug: assumed to be outside of the threshold */
	if(lower)
	  more_or_less = (my_sys->get_current()[var] > (th+w))?1:-1;
	else
	  more_or_less = (my_sys->get_current()[var] < (th-w))?-1:1;
	last_t=T;
	avg=0.;avgn=0.;
      }

      const number& switches(void)
      {
	return avgn;
      }

      const number& average(void)
      {
	return avg;
      }

      void probe(void)
      {
	const number factor=1.-1./n;

	runningvalue=factor*runningvalue+(1.-factor)*my_sys->get_current()[var];
	  
	bool more = (runningvalue > (th+w));
	bool less = (runningvalue < (th-w));
	integer now_more=more?1:(less?-1:0);		
	if(now_more != more_or_less) // if we have a change
	{
	  if(lower){			
	    /** \todo  We decouple the two kinds of probes (u & down). Agreed:
		this is wasteful and in the future I'll do the right
		thing) */
	    if (now_more==0) // inside DMZ
	      {
		// how did we enter
		from_top=(more_or_less==1);
	      }
	    if ((now_more==1) && !from_top) // and it is an up switch
	      {  
		number dt=get_time()-last_t; // get the time when we
		// swiched down
		avgn++;
		avg=((avgn-1.)*avg + dt)/avgn; // calculate average
		add_data(avg); // add the average
	      }
		  
	    if ((now_more==-1) && from_top) // down switch
	      {
		last_t=get_time(); // keep this point
	      }
	    more_or_less=now_more; // keep the current state 
	  }
	  else{
	    if (now_more==0) // inside DMZ
	      {
		// how did we enter
		from_bottom=(more_or_less==-1);
	      }
	    if ((now_more==-1) && !from_bottom) // and it is a down switch
	      {  
		number dt=get_time()-last_t; // get the time when we
		// swiched up
		avgn++;
		avg=((avgn-1.)*avg + dt)/avgn; // calculate average
		add_data(avg); // add the average
	      }
		  
	    if ((now_more==1) && from_bottom) // up switch
	      {
		last_t=get_time(); // keep this point
	      }
	    more_or_less=now_more; // keep the current state
	  }
	}
      }
      NO_COPY(AVGAVGDwellProbe);
	
      private: 
      system* my_sys;
      integer var;
      number th;
      number w;

      number avg;
      number avgn;
	
      /** Tells you where the system is :
	  1) above threshold+w
	  -1) below threshold-w
	  0) inside DMZ	
      */
      integer more_or_less;
      bool from_top, from_bottom;
      bool lower;
      number last_t;
	
      integer n;
      number runningvalue;
    };


}
#endif

/*********************************************************************
$Id: dwellprobe.h,v 1.2 2002-07-30 19:57:18 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.5  2002/07/30 08:39:25  mpeeters
Added pdf computation to binning stuff (bob).

Revision 1.1.2.3  2002/06/19 13:13:34  mpeeters
Added threshold functionality with DMZ.

Revision 1.1.2.2  2001/09/03 10:22:18  mpeeters
Added extra endl to allow for efficient gnuplot parsing.

Revision 1.1.2.1  2001/08/30 13:07:25  mpeeters
Fixed bug #456813, where the binning in BinDwellProbe was inverted.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.1  2000/10/13 11:20:22  mpeeters
Added two new probes (which were in the 0.9 tarball, but which I forgot
to add to the CVS). Removed history.

*********************************************************************/
