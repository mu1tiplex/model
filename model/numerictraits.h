/***************************************************************************
                          numerictraits.h  -  description
                             -------------------
    begin                : Wed May 10 2000
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

#ifndef NUMERICTRAITS_H
#define NUMERICTRAITS_H

#include "numerictypes.h"

namespace MODEL
{

  // Forward definitions
  template <typename numT, integer dims>
  class NumericTraits;

  //Previously:  <integer dims, class NT>
  template<integer dims, typename nelem, class NT>
  class NumVector;

  // Previously template <integer dims, class NT>
  template<integer dims, typename nelem, class NT>
  class VectorFunction;


  template <integer dims, class NT>
  class ScalarFunction;

  /**Traits classes for numerics
   */

  template <typename numT, integer dims>
  class NumericTraits
  {
  public:
	/** Real Type */
	typedef	numT												number;
	/** Vector Type */
	typedef NumVector< dims , numT, NumericTraits >				vect;
	/** 2D Matrix Type */
	typedef	NumVector< dims , vect, NumericTraits<vect,dims> >	matrix;
	/** And a tensor ... */
	typedef	NumVector< dims , vect, NumericTraits<matrix,dims> > tensor;

		/** vect = f(vect)  */
	typedef	VectorFunction< dims , numT, NumericTraits >		vf;
	/** number = f(vect) */                             	
	typedef	ScalarFunction< dims , NumericTraits<number,dims> >	sf;		
		
	/** Dimensions */
	static const integer	N=dims;
  };


  /** Specialisation for complex type */

  /*
	template <typename numT, integer dims>
	class NumericTraits<std::complex<numT>,dims>
	{
	public:
	typedef	numT											number;
	typedef NumVector< dims , NumericTraits >				vect;
	typedef	NumVector< dims , NumericTraits<vect,dims> >	matrix;
	};
  */

}
#endif

/*********************************************************************
$Id: numerictraits.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
