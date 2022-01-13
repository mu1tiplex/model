/***************************************************************************
                          eigenvalues.h  -  description
                             -------------------
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

#ifndef EIGENVALUES_H
#define EIGENVALUES_H

#include "numerictypes.h"
#include "numerictraits.h"
#include "invariant.h"
#include "utility.h"
#include <stdexcept>

namespace MODEL
{

  /**Returns the eigenvalues of a nonsymmetrix matrix
   */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class Eigenvalues
  {
  public:
	
	typedef typename NT::number	numT;
	typedef	typename NT::vect 	vect;
	typedef typename NT::matrix matrix;
	typedef typename NT::vf		vf;

	typedef	typename NumericTraits<complex,dims>::vect	EVect;
	
	/** Calculate upon initialisation */
	Eigenvalues(const matrix& ma) : m(ma)
	{calculate();}
	
	/** Get the real part */
	const	vect&	real(void) {return wr;}	

	/** Get the imag part */
	const	vect&	imag(void) {return wi;}	
			
  private:
	/** Return the balanced matrix.

		Balances the matrix for more accurate calculation. From NRC */
	const matrix&	balance(void);

	/** Return het Hessenberg form.
		Converts the matrix to upper Hessenberg form thru Gaussian pivoting. From NRC */
	const matrix&	hessenberg(void);
		
   	/** Calculate the eigenvalues explicitly. From NRC [hqr] */
	void	calculate(void);

	/** No copy or assign */
	NO_COPY(Eigenvalues);
	
  private:
	matrix	m;
	
	vect	wr;			// real part
	vect	wi;			// imaginary part
	
	static const number radix=2.;
  };

  /** Calculate the eigenvalues. This is much too long: should be cut
	  into itti bitty bite size pieces */
  template <integer dims, typename nelem, class NT >
  void
  Eigenvalues<dims,nelem,NT>::calculate(void)
  {
	// Make sure we have balanced it and reduced to upper hessenberg
	balance();
	hessenberg();
	
	numT 	anorm(abs(m[0][0]));
	
	for (integer i=0;i<dims;i++)        	// compute matrix norm
	  for (integer j=max(i-1,0);j<dims;j++)   // Because we start from upper Hessenberg !
		anorm += abs(m[i][j]);
	
	integer nn(dims-1);
	numT 	t(0.);     		// Gets changed only by an exceptional shift.
	
	while (nn >= 0)        	// Search next eigenvalue
	  {
		integer	its(0);   	// iterations
		integer l(0);
		do
		  {
			for (l=nn;l>=1;l--)    	// look for single small subdiagonal
			  {
				numT s=abs(m[l-1][l-1])+abs(m[l][l]);
				if (s == 0.0) s=anorm;
				if ( ( abs(m[l][l-1]) + s ) == s) break;
			  }
			numT x(m[nn][nn]),y(0.),z(0.);
			if (l == nn)   			// One root found
			  {
				wr[nn]=x+t; 		// Compensate for shift
				wi[nn]=0.0;
				nn--;
			  }
			else
			  {
				y = m[nn-1][nn-1];
				numT w(m[nn][nn-1]*m[nn-1][nn]);
				if (l == (nn-1))
				  {  // Two roots found
					numT	p = 0.5*(y-x);
					numT	q=p*p+w;
					z=sqrt(abs(q));
					x += t;    		// Compensate for shift
					if (q >= 0.0)
					  {	// a real pair
						z= p + ((p>=0.)?abs(z):-abs(z));
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					  }
					else
					  { 	// a complex pair
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]= -(wi[nn]=z);
					  }
					nn -= 2;
				  }
				else
				  {	// no roots found, continue
					numT p(0.),q(0.),r(0.),s(0.);
			
					if (its == 30) throw std::logic_error("Too many iterations in eigenvalue loop");
					if (its == 10 || its == 20)
					  {   // do a exceptional shift
						t += x;
						for (integer i=0;i<=nn;i++) m[i][i] -= x;
						s = abs(m[nn][nn-1]) + abs(m[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					  }
					
					++its;
					
					integer mm(0);
					
					for (mm=(nn-2);mm>=l;mm--)
					  {	// do shift and look for 2 consequtive small elements
						z=m[mm][mm];
						r=x-z;
						s=y-z;
						
						p=(r*s-w)/m[mm+1][mm]+m[mm][mm+1];
						q=m[mm+1][mm+1]-z-r-s;
						
						r=m[mm+2][mm+1];
						s=abs(p)+abs(q)+abs(r);
						
						p /= s;
						q /= s;
						r /= s;
						
						if (mm == l) break;
						numT u= abs(m[mm][mm-1])*(abs(q)+abs(r));
						numT v= abs(p)*(abs(m[mm-1][mm-1])+abs(z)+abs(m[mm+1][mm+1]));
						if ((u+v) == v) break;
					  }
					
					for (integer i=mm+2;i<=nn;i++)
					  {
						m[i][i-2]=0.0;
						if (i != (mm+2)) m[i][i-3]=0.0;
					  }
					
					for (integer k=mm;k<=nn-1;k++)
					  {	// Double QR step
						if (k != mm)
						  {
							p=m[k][k-1];      	// Householder vector setup
							q=m[k+1][k-1];
							r=0.;
							if (k != (nn-1)) r=m[k+2][k-1];
							x=abs(p)+abs(q)+abs(r);
							if (x != 0.0)
							  {
								p /= x;    			// Scale to prevent overflow
								q /= x;
								r /= x;
							  }
						  }
						s=sqrt(p*p+q*q+r*r);
						s=((p>=0.)?s:-s);
						if (s != 0.0)
						  {
							if (k == mm)
							  {
								if (l != mm) m[k][k-1] = -m[k][k-1];
							  }
							else m[k][k-1] = -s*x;
							p += s;
							
							x= p/s;
							y= q/s;
							z= r/s;
							
							q /= p;
							r /= p;
							
							for (integer j=k;j<=nn;j++)
							  {
								p=m[k][j]+q*m[k+1][j];
								if (k != (nn-1))
								  {
									p += r*m[k+2][j];
									m[k+2][j] -= p*z;
								  }
								m[k+1][j] -= p*y;
								m[k][j] -= p*x;
							  }
							integer mmin = nn<k+3 ? nn : k+3;
							for (integer i=l;i<=mmin;i++)
							  {
								p=x*m[i][k]+y*m[i][k+1];
								if (k != (nn-1))
								  {
									p += z*m[i][k+2];
									m[i][k+2] -= p*r;
								  }
								m[i][k+1] -= p*q;
								m[i][k] -= p;
							  }
						  }
					  }
				  }
			  }
		  } while (l < nn);
	  }
  }

  template <integer dims, typename nelem, class NT >
  const typename Eigenvalues<dims,nelem,NT>::matrix&
  Eigenvalues<dims,nelem,NT>::balance(void)
  {
	integer last(0);	
	numT sqrdx(radix*radix);
	
	while (last == 0)
	  {
		last=1;
		for (integer i=0;i<dims;i++)
		  {
			// Calculate row and column norms.
			numT r(0.),c(0.),g(0.);
			for (integer j=0;j<dims;j++)
			  if (j != i)
				{
				  c += abs(m[j][i]);
				  r += abs(m[i][j]);
				}
				
			if ((c!=0.) && (r!=0.))
			  {
				// If both are nonzero,
				numT f(1.);
				numT s(c+r);
				
				g=r/radix;
				while (c<g)
				  {
					// nd the integer power of the machine radix that comes closest to balancing the matrix.
					f *= radix;
					c *= sqrdx;
				  }
				
				g=r*radix;
				while (c>g)
				  {
					f /= radix;
					c /= sqrdx;
				  }
				if (((c+r)/f) < (0.95*s))
				  {
					last=0;
					g=1.0/f;
					// numT a=m[0][0];
					for (integer k=0;k<dims;k++) m[i][k] *= g; // Apply similarity transforma- tion.
					// numT b=m[0][0];
					for (integer k=0;k<dims;k++) m[k][i] *= f;
				  }
			  }
		  }
	  }	
	return m;
  }

  /** just for me */
  template<typename T>
  void swap(T& a,T& b)
  {
	T	temp(a);
	a=b;
	b=temp;
  }

  template <integer dims, typename nelem, class NT >
  const typename Eigenvalues<dims,nelem,NT>::matrix&
  Eigenvalues<dims,nelem,NT>::hessenberg(void)
  {
	integer i(0);
	for (integer rp1=1;rp1<(dims-1);rp1++)
	  { 					//rp1 is r + 1in NRC.
		numT x(0.);
		i=rp1;
		
		for (integer j=rp1;j<dims;j++)
		  {  // Find the pivot.
			if (abs(m[j][rp1-1]) > abs(x))
			  {
				x=m[j][rp1-1];
				i=j;
			  }
		  }
		if (i != rp1)
		  { // Interchange rows and columns.
			for (integer j=rp1-1;j<dims;j++) swap(m[i][j],m[rp1][j]);
			for (integer j=0;j<dims;j++) swap(m[j][i],m[j][rp1]);
		  }
		if (x!=0.)
		  {
			//Carry out the elimination.
			for (integer k=rp1+1;k<dims;k++)
			  {
				numT	y(0.);
				if ((y=m[k][rp1-1]) != 0.0)
				  {
					y /= x;
					m[k][rp1-1]=y;
					for (integer j=rp1;j<dims;j++) m[k][j] -= y*m[rp1][j];
					for (integer j=0;j<dims;j++) m[j][rp1] += y*m[j][k];
				  }
			  }
		  }
	  }		
	return m;
  }

} // end namespace
#endif

/*********************************************************************
$Id: eigenvalues.h,v 1.2 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
