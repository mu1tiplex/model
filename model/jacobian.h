/***************************************************************************
                          jacobian.h  -  description
                             -------------------
    begin                : Fri Apr 21 2000
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

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "numerictypes.h"
#include "utility.h"
#include <stdexcept>

namespace MODEL {

  /**	Calculates the Jacobian of a VF using finite differences.
		i is always the component of the function in question,
		j is always the parameter we differentiate to.
  */

  template <integer dims, class NT = NumericTraits<number,dims> >
  class Jacobian {
  public:
	typedef typename NT::number	numT;
	typedef	typename NT::vect 	vect;
	typedef typename NT::matrix matrix;
	typedef typename NT::vf		vf;

	/** Constructor taking a copy of the function
		@param own tells te constructor if it should take a copy, or use the reference */
	Jacobian(vf& rvf, numT JEps=1E-4, bool own=false) : 
	  owned(own), 
	  f((own)?(rvf.clone()):(&rvf)),
	  epsilon(JEps) {}	
	~Jacobian() {if (owned) delete f;}
	/** copying is allowed */
	Jacobian(const Jacobian& cj) : owned(false), f(cj.f) {}
	
	/** Calculate one element of the Jacobian J(i,j) = df(i)/dj. Not implemented. */
	numT	calculate_didj(integer i, integer j, const vect& u);
	
	/** Calculate one row of the Jacobian. Not implemented. */
	vect	calculate_di(integer i, const vect& u);
	
 	/** Calculate full Jacobian using only 1 evaluation and a user supplied function value at u.
		Adapted from NRC */
	matrix	calculate(const vect& u, const vect& fu);

    /** Calculate the jacobian, using a more accurate scheme. This
		eats two evaluations, but is second order accurate, using
		central differences */
	matrix calculate_accurate(const vect& u, const vect& fu);
	
	/** Calculate full Jacobian using 2 evaluations. 
		j[i][j]=df[i]dx[j] */
	matrix	calculate(const vect& u)
	{ vect fu( (*f)(u) ); return calculate(u,fu);}

	/** Calculate full Jacobian using 3 evaluations. 
		j[i][j]=df[i]dx[j] */
	matrix	calculate_accurate(const vect& u)
	{ vect fu( (*f)(u) ); return calculate_accurate(u,fu);}

	/** What function are we working on? */
	vf& get_function(void)
	{
	  return *f;
	}
				
  private:
	/** Assignment operator is private and NOT implemented */
	PRIVATE_ASSIGN(Jacobian);
	
  private:
	bool	owned;
	vf*	f;	
	/** epsilon is the relative step size used in the finite difference jacobian. */
	number					epsilon;
  };

  // This returns:
  // # Jacobian
  // 0 1 0 = df[0]/dx[j]
  // 0 0 1 = df[1]/dx[j]
  // 1 0 0 = df[2]/dx[j]
  // so the matrix is organised along rows: matrix[i] gives all
  // the derivatives of one component. 
  // j[i][j]=df[i]dx[j]

  template <integer dims, class NT>
  typename NT::matrix
  Jacobian<dims,NT>::calculate(const vect& u, const vect& fu)
  {	
	// There might be a faster way to do this, but I haven found it yet
	matrix	jac(0.);
	
	for(integer j=0;j<dims;j++)
	  {
		vect	udu(u);
						
		numT	du = epsilon*abs(u[j]);		
		
		if(du==0.0) du = epsilon;		// avoid numerical error	
		udu[j] += du;			// "
		du = udu[j] - u[j];		// "
		

		vect	fudu( (*f)(udu) );
		
		for(integer i=0;i<dims;i++)
		  {
			jac[i][j]=(fudu[i]-fu[i])/du;	
		  }		
	  }
	return jac;
  }

  template <integer dims, class NT>
  typename NT::matrix
  Jacobian<dims,NT>::calculate_accurate(const vect& u, const vect& fu)
  {	
	// There might be a faster way to do this, but I haven found it yet
	matrix	jac(0.);
	
	for(integer j=0;j<dims;j++)
	  {
		vect	updu(u),umdu(u);
						
		numT	du = epsilon*abs(u[j])/2.,dup,dum;		
		
		if(du==0.0) du = epsilon/2.;		// avoid numerical error	

		updu[j] += du;			// "
		dup = updu[j] - u[j];		// "

		umdu[j] -= du;			// "
		dum = u[j] - umdu[j];		// "

		vect	fupdu( (*f)(updu) );
		vect	fumdu( (*f)(umdu) );
		
		for(integer i=0;i<dims;i++)
		  {
			jac[i][j]=(fupdu[i]-fumdu[i])/(dum+dup);	
		  }		
	  }
	return jac;
  }


} // end namespace
#endif

/*********************************************************************
$Id: jacobian.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
