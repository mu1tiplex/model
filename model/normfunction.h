/***************************************************************************
                          normfunction.h  -  description
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

#ifndef NORMFUNCTION_H
#define NORMFUNCTION_H

#include "scalarfunction.h"

namespace MODEL{

/**Takes a vectorfunction and returns its norm instead
  */

template <integer dims, class NT = NumericTraits<number,dims> >
class NormFunction : public ScalarFunction<dims, NT>
{
	public:
		typedef	ScalarFunction<dims, NT>	base;
		typedef	typename NT::vf				vf;
		typedef typename base::numT numT;	
		typedef typename base::vect vect;

		/** Contructor which takes a copy
			@param own tells the const. if it is OK to take a copy
			@param scale scales the norm (handy in some cases)*/
		NormFunction(vf& vfu, bool own=false, numT scale=1.0)
						: owned(own), func((own)?(vfu.clone()):(&vfu)),values(0.),sc(scale) {}

		/** Copy constructor. Design decision: a copy off a Normfunction is
			NEVER owned (which is sort of logical, as we cannot access the
			original anymore anyway - no accessor for the moment) */
		NormFunction(const NormFunction& vfu)
			: base(vfu), owned(false),
				func(vfu.func), values(vfu.values), sc(vfu.sc) {}
			
		/** Assignment operator. Actually not sure if this is useful. Same
			strategy applies as for copy constructor */		
		const NormFunction& operator=(const NormFunction& vfu)
			{	if(this!=&vfu)
				{	
					if(owned)	delete func;
					func=vfu.func;
					owned=false;
					values=vfu.values;
					sc=vfu.sc;
				}
				return *this;}

		virtual ~NormFunction() { if(owned) delete func;}
	
		/** Virtual Copy Constructor */
		virtual NormFunction* clone (void) const
				{return new NormFunction(*this);}
		
					
	/** Implement this function to create the return vector.
		fu is a reference to where the values should be stored.
		This is for efficiency in internal routines.
		Inherited classes should return fu too.*/
	virtual const numT& function(numT& fu,const vect& u)
			{ values=(*func)(u); fu=sc*norm(values); return fu;}  	
		
			
	/** Return the value of the function during the last call */
	const vect&	function_value(void){return values;}
		
			
	private:
		bool	owned;
		vf*		func;
		/** Keeps the result of the last function call, so we
			do not have to call the function twice to get the values and the norm */
		vect	values;
		numT	sc;
};

} // end MODEL
#endif

/*********************************************************************
$Id: normfunction.h,v 1.2 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
