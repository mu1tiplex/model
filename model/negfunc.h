/***************************************************************************
                          cvr3dneg.h  -  description
                             -------------------
    begin                : Thu Apr 20 2000
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

#ifndef NEGFUNC_H
#define NEGFUNC_H

#include "vectorfunction.h"
#include "numvector.h"

namespace MODEL {

/**Returns the negative of the input for vectors
  */

// Prev: template <integer dims, class NT = NumericTraits<real,dims> >
template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
class NegFunc : public VectorFunction<dims,nelem,NT>{

	typedef	VectorFunction<dims,NT> 		base;
	typedef typename NT::vect				vect;
	
	public:				
		/** Virtual Copy Constructor */
		virtual NegFunc* clone () const
				{ return new NegFunc(*this); }

		virtual const vect& function(vect& fu, const vect& u)
				{
					fu=-u;
					return fu;
				}	
};

} // end namespace
#endif

/*********************************************************************
$Id: negfunc.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
