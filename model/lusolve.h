/***************************************************************************
                          lusolve.h  -  description
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

#ifndef LUSOLVE_H
#define LUSOLVE_H

#include "utility.h"
#include <stdexcept>
#include "numerictraits.h"
#include "numerictypes.h"

namespace MODEL {

  template <integer dims, class NT = NumericTraits<number,dims> >
  class LUSolve
  {
  public:
	typedef typename NT::number	numT;
	typedef	typename NT::vect	vect;	
	typedef typename NT::matrix	matrix;
		
	LUSolve(matrix& m, bool own=false) : M(own?m.clone():&m), owned(own) {decompose();};
	~LUSolve() {if(owned) delete M;}
		
	/** Solve it */
	vect	operator()(const vect& u) {vect res(u); solve(res); return
																  res;}

	/** Solve it, but with efficient return of vector */
	void	solve(vect& u);
	/** Determinant */
	numT	det(void);
 		
  private:
	// to avoid silly errors
	NO_COPY(LUSolve);
		
	// perform decomposition
	void	decompose(void);

  private:
	matrix*					M;		
	integer					pivotrows[dims];
	numT					pivotsign;
			
	bool	owned;

	static const numT	tiny=1.E-20;
  };

  template <integer dims, class NT >
  void LUSolve<dims,NT>::decompose(void)
  {
    matrix&		a(*M);

	pivotsign=1.0;	
	vect	rowscale(1.);
	// Test for singularity
	for (integer i=0;i<dims;i++)
	  {
		numT big=0.0;
		numT t(0.);
		for (integer j=0;j<dims;j++)
		  if ((t=fabs(a[i][j])) > big) big=t;
		if (big == 0.0) throw std::logic_error("Singular matrix in routine ludcmp");
		rowscale[i]/=big;
	  }
	
	numT sum(0.);
	for (integer j=0;j<dims;j++)
	  {
		for (integer i=0;i<j;i++)
		  {
			sum=a[i][j];
			for (integer k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		  }
		numT 	big(0.);
		integer imax(-1);
		numT    dum(0.);
		for (integer i=j;i<dims;i++)
		  {
			sum=a[i][j];
			for (integer k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=rowscale[i]*abs(sum) ) >= big) {
			  big=dum;
			  imax=i;
			}
		  }
		if (j != imax) {
		  for (integer k=0;k<dims;k++) {
			dum=a[imax][k];
			a[imax][k]=a[j][k];
			a[j][k]=dum;
		  }
		  pivotsign *= -1;
		  rowscale[imax]=rowscale[j];
		}
		pivotrows[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=tiny;
		if (j != (dims-1)) {
		  dum=1.0/(a[j][j]);
		  for (integer i=j+1;i<dims;i++) a[i][j] *= dum;
		}
	  }
  }

  template <integer dims, class NT>
  void LUSolve<dims,NT>::solve(vect& u)
  {
    matrix&		a(*M);
    numT		sum(0.);
	integer		ii(-1);
	bool	nonzero=false;
	
	for (integer i=0;i<dims;i++)
	  {		
		integer ip=pivotrows[i];
		sum=u[ip];
		u[ip]=u[i];
		
		if (nonzero) for (integer j=ii;j<i;j++) sum -= a[i][j]*u[j];
		else if (sum!=0.0) {nonzero=true; ii=i;}
		
		u[i]=sum;
	  }
	for (integer i(dims-1);i>=0;i--) {
	  sum = u[i];
	  for (integer j(i+1);j<dims;j++) sum -= a[i][j]*u[j];
	  u[i]=sum/a[i][i];
	}
  }

  template <integer dims, class NT>
  typename LUSolve<dims,NT>::numT	
  LUSolve<dims,NT>::det(void)
  {
	numT	d(pivotsign);
	for(integer j=0;j<dims;j++) d *= (*M)[j][j];
	return d;
  }


} // end MODEL;

#endif

/*********************************************************************
$Id: lusolve.h,v 1.3 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2002/07/30 19:57:18  mpeeters
Merge of 1.0.2 branch into trunk. From now on, the trunk will be a stable (meaning it compiles and has no quirks) branch and all special stuff will be done on separate branches waiting to be merged into the stable one. A tag STABLE will be made which moves with the trunk. Also, the main trunk will be considered to be 1.0.3.

Revision 1.1.2.1  2001/11/15 13:33:26  mpeeters
Last commit before moving everything to new gcc stdc++ headers (3.0).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.4  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.3  2000/09/22 08:54:08  mpeeters
Changed dependencies. Missed this one on the last commit

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
