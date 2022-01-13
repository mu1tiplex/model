/***************************************************************************
                   numvector.h  -  vector for numeric types
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

#ifndef NUMVECTOR_H
#define NUMVECTOR_H

#include <vector>
#include <stdexcept>
#include "numerictypes.h"
#include "numerictraits.h"

namespace MODEL {

	// prototyping friends - trying to fix error by defining template frined first

// Previously: template<integer dims, class NT>
template <integer dims, typename nelem, class NT >
	typename NT::number
	norm(const NumVector<dims,nelem,NT>& nv)
	{
	  return nv*nv;
	}

 // Previously: template<integer dims, class NT>
template <integer dims, typename nelem, class NT >
	typename NT::number
	length(const NumVector<dims,nelem,NT>& nv)
	{
	  return sqrt(norm(nv));
	}

// Previously: template<integer dims, class NT>
template <integer dims, typename nelem, class NT >
	inline
	NumVector<dims, nelem, NT>
	operator*(const typename NT::number& nvmul, const NumVector<dims,nelem,NT>& nv)
	{
	  NumVector<dims,nelem,NT> res(nv);
	  return res*=nvmul;
	}

// Previously: template<integer dims, class NT>
template <integer dims, typename nelem, class NT >
	inline
	NumVector<dims,nelem,NT>
	operator-(const NumVector<dims,nelem,NT>& m)
	{
	  NumVector<dims,nelem,NT>	v(m);
	  const typename NT::number	minus(-1.);
	  v*=minus;
	  return v;
	}



template <integer dims, typename nelem, class NT >
	inline
	NumVector<dims, nelem, NT>
	operator/(const NumVector<dims,nelem,NT>& nvmul, const typename NT::number&  nv)
	{
	  NumVector<dims,nelem,NT> res(nvmul);
	  return res/=nv;
	}

// Previously: template<integer dims, class NT>
template <integer dims, typename nelem, class NT >
	inline
	NumVector<dims,nelem,NT>
	operator+(const  NumVector<dims,nelem,NT>& nvmul, const NumVector<dims,nelem,NT>& nv)
	{
	  NumVector<dims,nelem,NT> res(nv);
	  return res+=nvmul;
	}

template <integer dims, typename nelem, class NT >
	inline
	NumVector<dims,nelem,NT>
	operator-(const  NumVector<dims,nelem,NT>& nvmul, const NumVector<dims,nelem,NT>& nv)
	{
	  NumVector<dims,nelem,NT> res(nvmul);
	  return res-=nv;
	}

/** Multiplier for real types
  */

 // Previously: template<integer dims, class NT>
 template <integer dims, typename nelem, class NT >
	inline
	typename NT::number
	operator*(const NumVector<dims,nelem,NT>& nv,const NumVector<dims,nelem,NT>& nvmul)
	{
	  typename NT::number r=0.0;
	  for(integer i=0;i<dims;i++) r += nv[i]*nvmul[i]; return r;
	}


  /**Class which adds functionality to Vector in case
	 the type of data is numeric.
	 */
  // Previously: template<integer dims, class NT = NumericTraits<number, dims> >
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
	class NumVector : public std::vector<typename NT::number>
  {
	public:
	typedef	typename	NT::number				number;
	typedef	typename	NT::number				numT;	// old code
	typedef	typename	NT::vect				vect;   // make life easier
	typedef	std::vector<number>	base;  	// base type

	/** Constuctor automatically resize for speed gain */
	NumVector(const numT& init=0.)
	{this->resize(dims,init);}
	virtual ~NumVector() {}

	/** Constructor to copy everything from an array
		resize has to be called, otherwise the objects[] do not exist !
	*/
	NumVector(const numT Tarray[])
	{ 	this->resize(dims);
		for(integer i=0;i<dims;i++) (*this)[i]=Tarray[i];}

	/** "virtual" copy constructor */
	NumVector*	clone(void){return new NumVector(*this);}

	/** Copy constructor calls to class*/
	NumVector(const NumVector& vcopy) : base(vcopy) {}

	/** Assignment operator calls top class */
	const NumVector&	operator=(const NumVector& vassign)
	{ if (this!=&vassign) base::operator=(vassign); return *this; }

	/** Overloaded accessor to catch errors (bound checking) */
	const numT&	operator[](integer i) const
	{
	  // removed bounds checking for speed
	  return std::vector<typename NT::number>::operator[](i);
	  // if (i>=0&&i<dims) return std::vector<typename NT::number>::operator[](i);
	  // else throw std::logic_error("array value out of bounds");
	  return (*this)[0]; // this should _never_ happen
	}

	/** Overloaded accessor to catch errors (bound checking) */
	numT&	operator[](integer i)
	{
	  if (i>=0&&i<dims) return std::vector<typename NT::number>::operator[](i);
	  else throw std::logic_error("array value out of bounds");
	  return (*this)[0]; // this should _never_ happen
	}

	// Arithmetic Operators
	// This will be taken over by PETE, undoubtedly

	NumVector& operator+=(const NumVector& nvplus)
	{ for(integer i=0;i<dims;i++) (*this)[i]+=nvplus[i]; return *this; }

	NumVector& operator-=(const NumVector& nvmin)
	{ for(integer i=0;i<dims;i++) (*this)[i]-=nvmin[i]; return *this; }

	NumVector& operator/=(const numT& nvdiv)
	{ for(integer i=0;i<dims;i++) (*this)[i]/=nvdiv; return *this; }

	NumVector& operator*=(const numT& nvmul)
	{ for(integer i=0;i<dims;i++) (*this)[i] *= nvmul; return *this; }

	friend NumVector operator*<dims, nelem, NT>(const numT& nvmul, const NumVector& nv);
	friend NumVector operator/<dims, nelem, NT>(const NumVector& nvmul, const numT& nv);


	friend NumVector operator+<dims, nelem, NT>(const NumVector& nvmul, const NumVector& nv);

	friend NumVector operator-<dims, nelem, NT>(const NumVector& nvmul, const NumVector& nv);

	/** Computes the scalar product. Wrong for complex types */
	/*
			friend numT operator*<dims,NT>(const NumVector& nv,
			const NumVector& nvmul);
	*/
	//		/** Compute the scalar product correct for complex types */
	//		friend real operator*(const NumVector< std::complex<numT> ,dims>& nv,
	//											const NumVector< std::complex<numT>,dims>& nvmul);

	/** Computes the norm of a vector.
	 */
	friend numT norm<dims,nelem,NT>(const NumVector& nv) ;

	/** Slightly more effort: the length. */
	friend numT length<dims,nelem,NT>(const NumVector& nv);

	friend NumVector operator-<dims,nelem,NT>(const NumVector& m);
  };



  ///** Multiplier for complex types
  //	*/
  //template<typename numT, integer dims>
  //inline
  //real operator*(const NumVector< std::complex<numT> ,dims>& nv,const NumVector< std::complex<numT>,dims>& nvmul)
  //{
  //	real r=0.0;
  //	for(integer i=0;i<dims;i++) r += nv[i].real()*nvmul[i].real()+nv[i].imag()*nvmul[i].imag(); return r;
  //}



} // end namespace MODEL
#endif

/*********************************************************************
$Id: numvector.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
