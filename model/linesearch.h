/***************************************************************************
                          linesearch.h  -  description
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

#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "numerictypes.h"
#include "normfunction.h"
#include "utility.h"
#include <stdexcept>
// #include "debugmacro.h"

namespace MODEL{

  /**Implements a line search for real multidim functions
   */

  template <integer dims, class NT = NumericTraits<number,dims> >
  class LineSearch // Could we derive it off vectorfunction - that is what it is... ???
  {
  public:
	typedef typename NT::number			numT;
	typedef	typename NT::vect			vect;	
	typedef	NormFunction<dims,NT>		nf;
	typedef	typename NT::vf				vf;
		
	LineSearch(vf& fu, bool own=false)
	  :  func(own?fu.clone():&fu), fmin(nf(*func,false,0.5)), tooclose(false), owned(own) {};
	~LineSearch() {if(owned)delete func;}
		
	/** implements a line searching and backtracking algoritm for a real function
		of the <dims> variable (the norm, actually), as described in
		NRC.
		@param uold		the previous point
		@param fold		the function value
		@param grad		the gradient
		@param p		the newton direction (calculated from the Jacobian+LU decomposition)
		@param u		a new point, where the function has decreased (hopefully)
		@param f		the function value in the new point
		@param maxstep	the largest step allowed
	*/
	void operator()			(	const vect& uold, const numT fold,
								vect& grad, vect& p,
								vect& u, numT& f, numT maxstep);

	bool converged(void) {return tooclose;}
			
	/** Returns the internal NormFunction, so one can use this to calculate
		it without having to define an auxillary function */
	nf&		norm(void){return fmin;}
	
	LineSearch(const LineSearch& ls) :	func(ls.func), fmin(ls.fmin), tooclose(ls.tooclose),
										owned(false) {}
	// to be safe
	PRIVATE_ASSIGN(LineSearch);
	

  private:
	vf*		func;
	nf		fmin;
	bool	tooclose;
	bool	owned;
	
  public:
	/** Fraction of step to take at least */
	static const number	alfa=1e-4;	
	/** Tolerance to see if the new point is too close. Usually signals convergence */
	static const number	tolerance=1e-7;
  };

  template <integer dims, class NT>
  void LineSearch<dims,NT>::operator()			
	(	const vect& uold, const numT fold,
		vect& grad, vect& p,
		vect& u, numT& f, numT maxstep)
  {
	tooclose=false;				// in NRC called: check
	
	// rescale if too large
	number	pl(length(p)); // in NRC called: sum
	if(pl>maxstep) p*=maxstep/pl; 	

	// scalar poduct of the gradient and the direction	
	number	slope(grad*p);		
	
	// calculate minimal lambda using heuristic
	number 	lambdascale(0.0);	// in NRC called: test
	for (integer i=0;i<dims;i++)
	  {
		number temp;
		temp=abs(p[i])/max(abs(uold[i]),number(1.0));
		if (temp > lambdascale) lambdascale=temp;
	  }
	
	number	lambdamin(tolerance/lambdascale),lambda(1.);		// in NRC called alamin,alam
	
	number lambda2(0.),templambda(0.),f2(0.);
	while (true) // Search forever
	  {
			
		u=uold; u+=lambda*p;
		f=fmin(u);
			
		//		DEBUG_Microscopic << "linesearch:" << u << " | " << f << " | "
		//			  << lambda*p <<endl;
   			 			
		if (lambda < lambdamin) { tooclose=true; u=uold ;return;}	// wooha - we are on top of it
		else if (f <= fold+alfa*lambda*slope) return;			// step is OK (hopefully this happens)
		else
		  {
			// First step ? We use a quadratic model for lambda
			if (lambda == 1.0) templambda = -slope/(2.0*(f-fold-slope));
			// All other steps: use a cubic (a few more calculations)
			else
			  {	
				number	rhs1 = f-fold-lambda*slope;
				number	rhs2 = f2-fold-lambda2*slope;
   					
				number	a=(rhs1/(lambda*lambda)-rhs2/(lambda2*lambda2))/(lambda-lambda2);
				number	b=(-lambda2*rhs1/(lambda*lambda)
						   +lambda*rhs2/(lambda2*lambda2))/(lambda-lambda2);
   					
				if (a == 0.0) templambda = -slope/(2.0*b); // easy roots
				else
				  {	// use discriminant
					number	disc=b*b-3.0*a*slope;
					if (disc<0.0) templambda=0.5*lambda;
					else if (b <= 0.0) templambda=(-b+sqrt(disc))/(3.0*a);
					else templambda=-slope/(b+sqrt(disc));  												
				  }
				if (templambda>0.5*lambda) templambda=0.5*lambda; // don't overdo it.
			  }	
		  }
		lambda2=lambda;
		f2 = f;
		lambda=max(templambda,0.1*lambda);	// don't overdo it.
	  }
  }

} // end MODEL
#endif

		

/*********************************************************************
$Id: linesearch.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
