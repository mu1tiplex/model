/***************************************************************************
                          radpotfunction.h  -  description
                             -------------------
    begin                : Tue Apr 25 2000
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

#ifndef RADPOTFUNCTION_H
#define RADPOTFUNCTION_H

#include "vectorfunction.h"

namespace MODEL {


/**A testfunction for the all the others: returns a vector
that points away from the origin and scales with the
distance squared. Only for testing purposes.
  */

class RadPotFunction : public VectorFunction<3>  {
	public:	
    	typedef	VectorFunction<3>	base;
		
		RadPotFunction() {}
		virtual ~RadPotFunction() {}	

		RadPotFunction(const RadPotFunction& crvf) : base(crvf) {}
		const RadPotFunction&	operator=(const RadPotFunction& arvf)
								{if(this!=&arvf){
								base::operator=(arvf);} return *this;}
		
		/** Virtual Copy Constructor */
		virtual RadPotFunction* clone () const {return new RadPotFunction(*this);}
		
		/** Implement this function to create the return vector */
		virtual const vect& function( vect& fu, const vect& u)
						{	numT	r2=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
							fu=u; fu*=r2; return fu;}	
};

} // end namespace MODEL
#endif

/*********************************************************************
$Id: radpotfunction.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
