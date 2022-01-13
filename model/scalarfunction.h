/***************************************************************************
                          scalarfunction.h  -  description
                             -------------------
    begin                : Tue May 2 2000
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

#ifndef SCALARFUNCTION_H
#define SCALARFUNCTION_H

#include "numerictypes.h"
#include "numerictraits.h"

namespace MODEL {

/**Base Class for scalar functions (of dims variables). Override the call() part to implement
  */

template <integer dims, class NT = NumericTraits<number,dims> >
class ScalarFunction
{
	public:
		typedef typename NT::number	numT;
        typedef typename NT::vect	vect;
		
		ScalarFunction() {}
		virtual ~ScalarFunction() {}
		
		/** Virtual Copy Constructor */
		virtual ScalarFunction* clone () const =0;
					
	/** Overloaded operator(), so the object can be presented as a function. */								
		numT	operator()(const vect& u)
			{ numT temp(0.); return function(temp,u); }
								
	/** Implement this function to create the return vector.
		fu is a reference to where the values should be stored.
		This is for efficiency in internal routines.
		Inherited classes should return fu too.*/
		virtual const numT& function(numT& fu,const vect& u) = 0;	
};

} //end namespace

#endif

/*********************************************************************
$Id: scalarfunction.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
