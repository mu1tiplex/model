/***************************************************************************
							 invariant.h
                             -----------

    begin                : Tue May 16 2000
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

#ifndef INVARIANT_H
#define INVARIANT_H

#include "lusolve.h"
#include "invariant.h"

namespace MODEL{
  /**Calculates the invariants of a matrix (sum over diagonal minors)
   */

  template <integer dims, class NT = NumericTraits<number,dims> >
  class Invariant
  {
  public:
	
	typedef typename NT::number	numT;
	typedef	typename NT::vect 	vect;
	typedef typename NT::matrix matrix;
	typedef typename NT::vf		vf;

	Invariant(const matrix& ma) : m(ma) {}
	
	/** Calculate the (n+1)-th order invariant. 1st order = trace, nth
		order = determinant. The others are the sums of the x-th order
		minor determinants. I agree: I must have been tripping when I
		wrote this code - that's recursive template programming for you. 
		@param done makes sure we don count elements
		twice */
	numT	calculate(integer i,integer done=0);	
	
	/** NO Copy or Assign */
	NO_COPY(Invariant);

  private:
	matrix				m;
  };

  template <integer dims, class NT>
  typename Invariant<dims,NT>::numT	
  Invariant<dims,NT>::calculate(integer i,integer done)
  {
	numT	result(0.);
	
	//cout << "M:" << m << endl;
	
    // Previously: typedef NumVector<dims-1,NumericTraits< NumVector<dims-1, NumericTraits<numT, dims-1> >,dims-1 > > 	smallermatrix;
	typedef NumVector<dims-1,numT,NT>	smallervector;
	typedef NumVector<dims-1,smallervector,NT>	smallermatrix;
	typedef Invariant<dims-1, NumericTraits< numT,dims-1> > 	smallerinv;
	
	if (i==(dims-1)) {LUSolve<dims,NT> l(m,true); return l.det();}
	else
	  {
        for(integer remove=done;remove<dims;remove++) // pick a row to be removed
		  {
   			smallermatrix	mr(0.);
        	integer s(0);	             			
        	for(integer j=0;j<dims;j++) // loop over rows
			  {
        		if (j!=remove)
				  {
        			integer t(0);
        			for (integer k=0;k<dims;k++) // loop over colums
					  {
        		    	if (k!=remove)
						  {	
        		    		mr[s][t] = m[j][k]; // make smaller matrix
        		    		t++;
						  }
					  }
        			s++;
				  }
			  }
        	//cout << "mr:" << mr << endl;
            smallerinv	inmo(mr);			// Take a copy
            result += inmo.calculate(i,remove);  // Calculate smaller one
		  }
	  }
 	return result;	
  }	

  /** Ending the recursion */
  template < class NT >
  class Invariant<1,NT>
  {
  public:
	
	typedef typename NT::number	numT;
	typedef	typename NT::vect 	vect;
	typedef typename NT::matrix matrix;
	typedef typename NT::vf		vf;

	Invariant(const matrix& ma) : m(ma) {}
	
	/** Calculate the n-th order invariant. 1st order = trace, nth order = determinant. */	
	numT	calculate(integer i,integer done=0);	
	
  private:
	/** No copy or assign */
	NO_COPY(Invariant);
		
  private:
	matrix				m;
  };

  template <class NT>
  typename Invariant<1,NT>::numT	
  Invariant<1,NT>::calculate(integer i,integer done)
  {
 	return m[0][0];
  }

} // end namespace

#endif

/*********************************************************************
CVS INFO: $Id: invariant.h,v 1.2 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.6  2000/09/22 08:51:45  mpeeters
Added NumericTraits NT template parameter to all declarations and
definitions. Not really needed, but cleaner this way.

Revision 1.5  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

Revision 1.4  2000/09/15 08:28:30  mpeeters
Still getting to know cvs. Testfile: invariant.h
	(which is why my version number is going nuts)

Revision 1.3  2000/09/14 13:54:51  mpeeters
Added code template file

Revision 1.2  2000/09/14 13:37:19  mpeeters
First step towards cleaner sources: CVS tags test

*********************************************************************/

