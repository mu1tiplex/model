/***************************************************************************
                          newtonroot.h  -  description
                             -------------------
    begin                : Mon May 8 2000
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

#ifndef NEWTONROOT_H
#define NEWTONROOT_H

#include "numerictypes.h"
#include "numerictraits.h"
#include <stdexcept>
#include "jacobian.h"
#include "linesearch.h"
#include "lusolve.h"
#include "utility.h"

// #include "debugmacro.h"

namespace MODEL {

  /**Finds a root of the vectorfunction.
   */

  // Previously: template <integer dims, class NT = NumericTraits<number,dims> >
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class NewtonRoot : public VectorFunction<dims,number,NT>
  {
  public:
	typedef typename NT::number	numT;
	typedef	typename NT::vect 	vect;
	typedef typename NT::matrix matrix;
	typedef typename NT::vf		vf;
	typedef VectorFunction<dims,number,NT>				base;
	
	NewtonRoot(vf& f, bool own=false)
	  : 	func( (own)? (f.clone()) : (&f) ), owned(own),
			wrongmin(false), ls(*func), fdjac(*func) {}	
				
	~NewtonRoot(){ if (owned) delete func;}

	/** Copy constructor needed for clone(). The function not copied ! */
	NewtonRoot(const NewtonRoot& nr)
	  : 	func(nr.func), owned(false),
			wrongmin(nr.wrongmin), ls(*func), fdjac(*func) {}	
	
	virtual NewtonRoot* clone () const {return new NewtonRoot(*this);}
	
	/** Implement this function to create the return vector.
		fu is a reference to where the values should be stored.
		This is for efficiency in internal routines.
		Inherited classes should return fu too.*/
	virtual const vect& function(vect& fu,const vect& startu);	
	
	/** Is this a correct minimum */
	bool	wrong_min(void){return wrongmin;}

	/** Have we run out of options */
	bool no_root(void)
	{
	  return noroot;
	}	
	// no assignment allowed
	PRIVATE_ASSIGN(NewtonRoot);
	
		
  private:
	vf*						func;
	bool					owned;
	bool					wrongmin;
	bool noroot;
	
	LineSearch<dims,NT>		ls;
	Jacobian<dims,NT>		fdjac;

  public:
	static const	integer		maxiterations=10000;      // NRC: MAXITS
	static const	numT	 	tolerancef=1.E-8;      	// NRC: TOLF
	static const	numT		tolerancemin=1.E-6;     // NRC: TOLMIN
	static const	numT		maxstep=100.;           // NRC: STPMX
	static const	numT		tolerancex=1.E-8;		// NRC: TOLX
  };

  // --- NON INLINED ---
//   template <integer dims, typename nelem, class NT >
//   NewtonRoot<dims,nelem,NT>::NewtonRoot(const NewtonRoot& nr) : ls(nr.ls),fdjac(nr.fdjac) {cout
// 	<< "NewtonRoot : BOOM" << endl; throw std::logic_error("Copying NewtonRoot is evil");}

  // Prev: template <integer dims, class NT >
  template <integer dims, typename nelem, class NT >
  const typename NewtonRoot<dims,nelem,NT>::vect&
  NewtonRoot<dims,nelem,NT>::function(vect& fu,const vect& startu)
  {	
	// We are looking
	noroot=false;

	// Calculate the initial function value and norm
	numT		f=ls.norm()(startu);
	vect		fvec=ls.norm().function_value();

	// alias for clarity
	vect&	u(fu);	// The return vector IS the final point where we end
	u=startu;   			// And we start fro, the initial point
	

	// Test if we are too close to a zero
	numT	testtol=norm(fvec);           // More efficient with square
	// for (integer i=0;i<dims;i++)
	// 	{ numT tm(0.); if ( (tm=abs(fvec[i]))> testtol) testtol=tm; }		
	
	if (testtol<0.01*tolerancef*tolerancef) return fu;      // close enough to zero already
	
	// Calculate maximal step
	numT	sum=norm(u);
	numT	stpmax=maxstep*max(numT(sqrt(sum)),numT(dims));
	
	for (integer its=0;its<maxiterations;its++)
	  {
		//	  DEBUG_Cellular << "newtonroot:" << u << endl;
	  
		// Calculate Jacobian
		matrix	j=fdjac.calculate(u,fvec);
	
		vect	gradient(0.);
		
		// Calculate Gradient
		for (integer i=0;i<dims;i++)
		  {
			numT	sum(0.);
			for (integer k=0;k<dims;k++) sum += j[k][i]*fvec[k];
			gradient[i]=sum;
		  }
			
		// Keep previous point
		vect		uold(u);
		numT		fold(f);
		
		// Descent direction
		vect		p(-fvec);
			
		// Solve system using LU decomposition
		// ludcmp(fjac,n,indx,&d);
		LUSolve<dims,NT>	lus(j);  // init
		
		lus.solve(p);	// solve by backsubstitution
		
		ls(uold,fold,gradient,p,u,f,stpmax); // Linesearch with this value
		
		fvec=ls.norm().function_value();		// Get value (without recalculating)
		
		wrongmin=ls.converged();				// Error coming from ls
				
		// Test if we are too close to a zero
		//testtol=0.;
		//for (integer i=0;i<dims;i++)
		//	{ numT tm(0.); if ( (tm=abs(fvec[i]) ) > testtol) testtol=tm; }		
		
		testtol=norm(fvec);
		if (testtol<tolerancef*tolerancef)
		  {
			//		  DEBUG_Grainy << "newtonroot: exited normally" << endl;
			wrongmin=false;
			return fu;
		  }      // close enough to zero already
		
		// Make sure we are not at a spurious convergence (grad f=0)
		if (wrongmin)
		  {  	
			// Besides, whenever we arrive here, something boogery has happened
			testtol=0.0;
			numT den=max(f,0.5*number(dims));
			numT tm(0.);
			for (integer i=0;i<dims;i++)
			  {
				tm=abs(gradient[i])*max(abs(u[i]),number(1.0))/den;
				if (tm > testtol) testtol=tm;
			  }
			wrongmin = (testtol < tolerancemin)||(!(norm(fvec)<tolerancef*tolerancef));
			// DEBUG_Grainy << "newtonroot: gradient?" << gradient <<
			// "zero:" << fvec <<endl;
			// DEBUG_Grainy<< "newtonroot: exited thru spurious
			// convergence"<<endl;
			return fu;
		  }
		
		testtol =0.0;	// check for converence on u
		numT tm(0.);
		for (integer i=0;i<dims;i++) {
		  tm=(abs(u[i]-uold[i]))/max(abs(u[i]),number(1.0));
		  if (tm > testtol) testtol=tm;
		}
		if (testtol < tolerancex)
		  {	
			//		  DEBUG_Grainy << "newtonroot: converged on u"<<endl;	
			return fu;
		  }
	  }
	noroot=true;
	throw std::logic_error("NewtonRoot: maxiteration exceeded!");
	return fu=vect(0.);     // BAD BAD BAD we should never come here
  }	



} // end MODEL
#endif

/*********************************************************************
$Id: newtonroot.h,v 1.3 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2002/07/30 19:57:18  mpeeters
Merge of 1.0.2 branch into trunk. From now on, the trunk will be a stable (meaning it compiles and has no quirks) branch and all special stuff will be done on separate branches waiting to be merged into the stable one. A tag STABLE will be made which moves with the trunk. Also, the main trunk will be considered to be 1.0.3.

Revision 1.1.2.2  2002/07/09 08:40:43  mpeeters
Added explicit type conversion for max(), some compilers (gcc 2.96) fail to convert  otherwise.

Revision 1.1.2.1  2002/02/06 12:48:04  mpeeters
Changes minor type error.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.4  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.3  2000/09/22 08:54:45  mpeeters
Change of default tolerances when looking for roots.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
