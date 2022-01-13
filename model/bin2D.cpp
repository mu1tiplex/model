/***************************************************************************
                          bin2D.cpp  -  2D binninig
                             -------------------
    begin                : Tue Jul 1 2002
    copyright            : (C) 2002 by Michael Peeters
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

#include "bin2D.h"

namespace MODEL {

  Bin2D::Bin2D(number start, number end, counter nobins)
	: x0(start),x1(end),N(nobins)
  {
	const hist1D zero_element(N+1);

	dx=(x1-x0)/number(N);

	histogram=new hist2D(N+1,zero_element);
  }

  Bin2D::~Bin2D()
  {
	delete histogram;
  }

  void
  Bin2D::add_value(number addedx, number addedy)
  {
	counter binx,biny;

	// X
	if(addedx<x0) binx=0;
	else if(addedx>=x1) binx=N;
	else 
	  {
	  addedx-=x0;
	  addedx/=dx; 

	  binx=counter(addedx)+1;
	  }

	if(binx<1)binx=0;
	else if(binx>N)binx=N;

	// Y
	if(addedy<x0) biny=0;
	else if(addedy>=x1) biny=N;
	else 
	  {
	  addedy-=x0;
	  addedy/=dx; 

	  biny=counter(addedy)+1;
	  }

	if(biny<1)biny=0;
	else if(biny>N)biny=N;

	(*histogram)[binx][biny]++;
  }

  const Bin2D::hist2D&
  Bin2D::get_histo(void)
  {
	return (*histogram);
  } 

  Bin2D::start1D
  Bin2D::get_bins(void)
  {
	start1D bins(0);
	number x=x0-dx;
	for(counter i=0;i<N+1;i++)
	  {
				bins.push_back(x);
				x+=dx;
	  }
	return bins;
  }
} // end namespace


/*********************************************************************
$Id: bin2D.cpp,v 1.2 2002-07-30 19:57:18 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.1  2002/07/01 13:34:42  mpeeters
Added 2D binning class and probe.

*********************************************************************/
