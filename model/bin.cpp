/***************************************************************************
                          bin.cpp  -  description
                             -------------------
    begin                : Tue Jul 25 2000
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
#include <stdexcept>
#include <iostream>
#include "bin.h"

namespace MODEL {

  Bin::Bin(number start, number end, counter nobins)
	: x0(start),x1(end),N(nobins)
  {
    if(x1>x0)
      {
	dx=(x1-x0)/number(N);
	histogram=new vector<counter>(0);
	histogram->resize(N+1,0);
      }
    else throw std::logic_error("Bin : start is larger or equal to first in  constructor");
  }

  /** copy constructor */

  Bin::Bin(const Bin& tocopy , counter firstbin, counter lastbin)
     {
       if(firstbin<=lastbin)
	 {
	   x0=tocopy.get_bins()[firstbin];
	   x1=tocopy.get_bins()[lastbin];
	   N=lastbin-firstbin+1;
	   dx=tocopy.get_width();
	   histogram=new vector<counter>(0);
	   histogram->resize(N+1,0);
	   (*histogram)[0]=0;
	   for(counter i=1;i<=N;i++)
	     { 
	       if(i+firstbin-1<=tocopy.get_numberofbins()) (*histogram)[i]=tocopy.get_histo()[firstbin+i-1];
	       else (*histogram)[i]=0;
	     };
	 }
       else throw std::logic_error("Bin : firstbin is larger dan lastbin in copy constructor");
     }


  Bin::~Bin()
  {
	delete histogram;
  }

  void
  Bin::add_value(number added)
  {
	counter bin;
        
	if(added<x0) {bin=0; cerr<<"warning : added number is smaller than the first bin"<<endl;}
	else if(added>=x1) {bin=N; cerr<<"warning : added number is larger than the last bin"<<endl;}
	else 
	  {
	  added-=x0;
	  added/=dx; 

	  bin=counter(added)+1;
	  }

	if(bin<1)bin=0;
	else if(bin>N)bin=N;
	(*histogram)[bin]++;
  }

  number  Bin::get_width() const
  {
    return dx;
  }

  counter Bin::totalbinned() const
  {
    counter total(0);
    for(counter i=1;i<=N;i++)
      {
	total+=(*histogram)[i];
      }
    return total;
  }


  vector<counter> 
  Bin::get_histo(void) const
  {
	return (*histogram);
  } 

  vector<number>
  Bin::get_bins(void) const
  {
	vector<number> bins(0);
	number x=x0-dx;
	for(counter i=0;i<N+1;i++)
	  {
				bins.push_back(x);
				x+=dx;
	  }
	return bins;
  }

  vector<number> Bin::get_pdf() const
  {
    vector<number> pdf(0);
    pdf.push_back(0);
    for(counter i=1;i<=N;i++)
    {  
      number prob = (*histogram)[i]/(this->totalbinned()*dx);
      pdf.push_back(prob);
    }
    return pdf;
  }
  
  counter Bin::get_numberofbins() const
  {
    return N;
  }

} // end namespace

/*********************************************************************
$Id: bin.cpp,v 1.2 2002-07-30 19:57:18 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.1  2002/07/30 08:39:24  mpeeters
Added pdf computation to binning stuff (bob).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
