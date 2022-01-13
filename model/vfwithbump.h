/***************************************************************************
                          vfwithbump.h  -  description
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

#ifndef VFWITHBUMP_H
#define VFWITHBUMP_H

#include "vectorfunction.h"
#include "newtonroot.h"
#include <vector>

namespace MODEL {
  /**Makes a new vectorfunction, but with an extra hyperbump.
	 The idea: use this to locally remove zero's by dividing them out.
	 */

  // Prev: template <integer dims, class NT = NumericTraits<number,dims> >
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class VFwithBump : public VectorFunction<dims,nelem,NT>
  {
  public:
	  typedef VectorFunction<dims,nelem,NT> base;
	  typedef base vf;
	  typedef typename base::vect vect;
	  
  public:
	/** Just create a possibility to add bumps on top of an existing
		VectorFunction. */
	VFwithBump(vf& func,bool own=false)
	  : f((own)?func.clone():&func), owned(own) {}
	
	/** Create one, with a bump already specified */
	VFwithBump(vf& func,const vect& position,bool own=false)
	  : f((own)?func.clone():&func), owned(own) { AddBump(position);}
	~VFwithBump() {if(owned)delete f;}
	
	/** Adds a bump. */
	void	AddBump(const vect& position) {p.push_back(position);}	
		
	/** copy */
	VFwithBump(const VFwithBump& vfb) : f(vfb.f), owned(false), p(vfb.p){}
	/** assign */
	const VFwithBump& operator=(const VFwithBump& vfb)
	{ if (this!=&vfb){f=vfb.f;owned=false;p=vfb.p;} return *this;}
		
	/** Virtual Copy Constructor */
	virtual VFwithBump* clone () const {return new VFwithBump(*this);}
	
	/** Implement this function to create the return vector.
		fu is a reference to where the values should be stored.
		This is for efficiency in internal routines.
		Inherited classes should return fu too.*/
	virtual const vect& function(vect& fu,const vect& u)
	{			
	  // Here we have to iterate over all p's
	  typename vector<vect>::iterator	runner(p.begin());

	  // if we have no bumps, nothing should happen
	  if(runner==p.end()) {return f->function(fu,u); }
	  
	  // Otherwise
	  nelem	factor=1.;
	  while(runner!=p.end())
		{
		  // if we have a point,
		  
		  // Calculate delta u
		  vect bump(u);
		  bump-=(*runner);
		  
		  // Calculate how far we are from the point and add a
		  // hypersphere of zeroes
		  /** \todo It would be even better to make it a
			  hyperellipsoid instead of a sphere */
		  number distance=abs(length(bump)-1.5*NewtonRoot<dims,nelem,NT>::tolerancex)/length(bump);   
		  					
		  factor /= distance; // close to zero goes to zero
		  ++runner;
		}	
	  f->function(fu,u);
	  fu*=(1+factor);

	  return fu;
	}	
		
  private:
	vf*		f;
	bool	owned;
	vector<vect>	p; 	// not just one position
  };

} // end namespace
#endif

/*********************************************************************
$Id: vfwithbump.h,v 1.3 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2001/07/26 12:08:28  mpeeters
Removed a nasty error where the size of the hypersphere used to cut the zero out of existence would be absolute instead of relative. A better (relative) solution was found, but an even better one exists (see todo note).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
