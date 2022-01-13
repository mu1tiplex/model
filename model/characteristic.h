/***************************************************************************
                          characteristic.h  -  description
                             -------------------
    begin                : Wed May 24 2000
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

#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include "invariant.h"
#include "utility.h"

namespace MODEL{
  /**Calculates the characteristic function of a matrix
   */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class Characteristic
  {
  public:
	
	typedef typename NT::number	numT;
	typedef	typename NT::vect 	vect;
	typedef typename NT::matrix matrix;
	typedef typename NT::vf		vf;

	/** The Characteristic is fully calculated upon creation of the
		object. */
	Characteristic(const matrix& ma)
	{
	  Invariant<dims,NT>	inv(ma);
	  for(integer u=0;u<dims;u++) poly[u]=inv.calculate(u);
	}

	/** Calculates the characteristic polynomial at a certain point.
		\f[
		P(\lambda)=\lambda^n+a_0.\lambda^{n-1}+...+a_n
		\f]
		where the \f$a_x\f$ are the
		invariants of this particular matrix */
	template <typename n>	
	n	operator()(const n& l)
	{
	  n	res(0.);
	  n	power(1.);
	  numT	sign=((dims%2)==0)?1.:-1.;
	  for(integer u=dims-1;u>=0;u--)
		{
		  res+=sign*poly[u]*power;
		  sign*=-1.;
		  power*=l;
		}	
	  res+=power;
	  return res;
	}	
			
	/** No Copy or Assign */
	NO_COPY(Characteristic);
	
  private:
	vect				poly;
  };
} // end namespace
#endif

/*********************************************************************
$Id: characteristic.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
