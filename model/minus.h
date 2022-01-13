/***************************************************************************
                          minus.h  -  description
                             -------------------
    begin                : Fri May 12 2000
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

#ifndef MINUS_H
#define MINUS_H

#include "vectorfunction.h"
#include "newtonroot.h"
#include <vector>

namespace MODEL {
  /**Makes a new vectorfunction, which is minus the input
	 The idea: use this to change stable states to unstable and vice-versa
	 */

  // Prev: template <integer dims, class NT = NumericTraits<number,dims> >
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class Minus : public VectorFunction<dims,nelem,NT>
  {
  public:
	/** Just create */
	Minus(vf& func,bool own=false)
	  : f((own)?func.clone():&func), owned(own) {}
	
	/** copy */
	Minus(const Minus& vfb) : f(vfb.f), owned(false){}
	/** assign */
	const Minus& operator=(const Minus& vfb)
	{ if (this!=&vfb){f=vfb.f;owned=false;} return *this;}
		
	/** Virtual Copy Constructor */
	virtual Minus* clone () const {return new Minus(*this);}
	
	/** Implement this function to create the return vector.
		fu is a reference to where the values should be stored.
		This is for efficiency in internal routines.
		Inherited classes should return fu too.*/
	virtual const vect& function(vect& fu,const vect& u)
	{			
	  f->function(fu,u);
	  fu*=-1.;
	  return fu;
	}	
		
  private:
	vf*		f;
	bool	owned;
  };

} // end namespace
#endif

/*********************************************************************
$Id: minus.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
