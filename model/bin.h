/***************************************************************************
                          bin.h  -  description
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

#ifndef BIN_H
#define BIN_H

#include "numerictypes.h"
#include "utility.h"
#include <vector>

/** The global namespace for this library.
	We should take care to use this consistently
*/

namespace MODEL {

  using namespace std;
  
  /**A class for binning: gathers large amounts 
	 of numeric data into bins, for histogramming purposes.
	 */

	// 	define the binning:
	// 		0 : < start
	// 		1 : >= start, 			<start + inc
	// 		2 : >= start+inc,		<start + 2*inc
	// 		...
	// 		N : >= start+(N-1)*inc,	<end
	// 		N+1:>= end

  class Bin {
  public:
	
	/** Create a binning
		@param nobins the number of bins
	*/		
	Bin(number start, number end, counter nobins=16);
	~Bin();
        
	/** Copy constr1uctor, return a cropped (or enlarged) bin with the same bin width*/  
	Bin(const Bin& tocopy, counter firstbin, counter lastbin);

	/** add a value to the bin. */
	void	add_value(number added);
	
	/** returns the binwidth. */
	number get_width() const;

	/** returns the total number of binned numbers, those smaller than x0 not included*/
	counter totalbinned() const; 

	/** returns an array of numbers representing the histogram */
	vector<counter>	get_histo( ) const;
	
	/** returns an array of bin starts */
	vector<number> get_bins(void) const;

	/** returns an array of the pdf of the binned data (the probability density in the 
	    middle of the first bin is returned at place 1 in the vector), 
	    those smaller than x0 not included*/ 
        vector<number> get_pdf() const;

	/** returns the number of bins (N)*/
	counter get_numberofbins() const;

	NO_COPY(Bin); 
 
  private:
	number x0;
	number x1;
	number dx;
	counter N;
	
	vector<counter>*	histogram;
  };

} // end MODEL

#endif

/*********************************************************************
$Id: bin.h,v 1.2 2002-07-30 19:57:18 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.1  2002/07/30 08:39:25  mpeeters
Added pdf computation to binning stuff (bob).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.4  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.3  2000/10/09 09:46:56  mpeeters
Administrative changes (new files etc...)
Added NO_COPY to Bin

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
