/***************************************************************************
                          bin2D.h  -  2D binning
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

#ifndef BIN2D_H
#define BIN2D_H

#include "numerictypes.h"
#include "utility.h"
#include <vector>

/** The global namespace for this library.
	We should take care to use this consistently
*/

namespace MODEL {

  using namespace std;
  
  /**A class for binning: gathers large amounts 
	 of numeric 2D data into bins, for histogramming purposes.
	 */

	// 	define the binning:
	// 		0 : < start
	// 		1 : >= start, 			<start + inc
	// 		2 : >= start+inc,		<start + 2*inc
	// 		...
	// 		N : >= start+(N-1)*inc,	<end
	// 		N+1:>= end

  class Bin2D {
  public:
	/** typedefs for ease */
	typedef vector<counter> hist1D;
	typedef vector<hist1D> hist2D;

	typedef vector<number> start1D;
	typedef vector<start1D> start2D;

	/** Create a binning
		@param nobins the number of bins
	*/
	Bin2D(number start, number end, counter nobins=16);
	~Bin2D();

	/** add a value to the bin. */
	void	add_value(number addedx,number addedy);
	
	/** returns an array of numbers representing the histogram */
	const hist2D&	get_histo( );

	/** returns an array of bin starts */
	start1D get_bins(void);

	NO_COPY(Bin2D); 

  private:
	number x0;
	number x1;
	number dx;
	counter N;
	
	hist2D*	histogram;
  };

} // end MODEL

#endif


/*********************************************************************
$Id: bin2D.h,v 1.2 2002-07-30 19:57:18 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.1  2002/07/01 13:34:42  mpeeters
Added 2D binning class and probe.

*********************************************************************/
