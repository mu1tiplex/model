/***************************************************************************
			reservoirfunction.h
			-----------
                             
    begin                : Wed Apr 25 00:47:40 CEST 2001
    author               : (C) 2001 by Michael Peeters
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

#ifndef RESERVOIR_H
#define RESERVOIR_H
#include "numerictypes.h"
#include "numvector.h"
#include "vectorfunction.h"
#include "scalarfunction.h"

namespace MODEL
{

  /** This will be a class for reservoir rate equations, on which we can
	  do noise simulations. There is always one dump reservoir, which
	  corresponds to the infinite universe. Use it (but don't fill it
	  up.)

	  \par What does it do?  For n dimensions, you have n+1 reservoirs
	  (include the universe which is indexed by 0), so we need to
	  define 0.5*n(n-1)+n exchange rates.  (The universe only eats
	  particles.)

	  \par You define these using the supplied macros for ease, but it
	  the future I might write a parser. The class then represents
	  these to the outside world as a normal VectorFunction, but of
	  course the extra consistency can be used to extract noise
	  response information.

	  \par How do I use it?  Derive your class from ReservoirFunction,
	  and implement all the required exchange rates (as friend
	  functions). If you miss one, the class will complain (I hope).

	  \code
	  class MySystem : ReservoirFunction<3>
	  {
	  ...
	  };

	  \endcode 
  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class ReservoirFunction : public VectorFunction<dims,nelem,NT>
  {
  public:
	typedef	typename NT::number					numT;
	typedef	typename NT::vect					vect;
	typedef typename NT::vf						vf;
	typedef typename NT::matrix					matrix;
	
	// typedef ScalarFunction<dims>* exchange_rate;
	typedef number exchange_rate_f(ReservoirFunction*,const vect& u);
	typedef exchange_rate_f* exchange_rate;

	// typedef number (ReservoirFunction::*exchange_rate)(const vect&);
	typedef exchange_rate*   exchange_rate_list;

	ReservoirFunction()
	{
	  for (integer i=0;i<dims+1;i++)
		{
		  R[i]=new exchange_rate[dims+1];
		  for(integer j=0;j<dims+1;j++)
			R[i][j]=0; // make sure all pointer are zero to start with  
		}
	}
	
	/** Copy constructor: copy parameters as well !
		\note: this should not work, as R is not copied, but strangely
		enough, it does. WHY? And it doesn't ...
		It used to work because the amputated version does not have a
		correct call to the base class, and just copied all  things
		binary. Of course, once you try to make a "real" copy, it fails.
		\todo: Clean out copy constructor dependencies and calling
		sequences in the VectorFunction Hierarchy 
		\bug: this leaks as much as an bucket without a bottom. We
		really should rewrite the destructor to mop up.
	*/

	ReservoirFunction(const
											   ReservoirFunction& rf)
	{
	  const map<string,number*>& subparlist=rf.get_parlist();
	  map<string,number*>::const_iterator run=subparlist.begin();
	  map<string,number*>::const_iterator endrun=subparlist.end();
	  
	  while(run!=endrun){
		define_parameter(run->first,*(run->second));
		++run;
	  }

	  // copy R functions...
	  for (integer i=0;i<dims+1;i++)
		{
		  R[i]=new exchange_rate[dims+1];
		  for(integer j=0;j<dims+1;j++)
			R[i][j]=rf.R[i][j]; // make sure all pointer are zero to start with  
		}
	  
	}

	virtual void add_exchange_rate(int from, int to,
								   exchange_rate er)
	{
	  R[from][to]=er;
	}

	/** This is the input vector. It might depend on the specific point */
	virtual const vect& input(vect& fu,const vect& u) = 0;

	/** This is the workhorse: it tries to take all the exchange rate
		functions and puts them into a coherent system which reacts
		like a vectorfunction. You have to
		make sure you defined them all.
	*/ 
	virtual const vect& function(vect& fu,const vect& u)	
	{
	  // set input (sources)
	  input (fu,u);
	  
	  for (integer i=0;i<dims;i++)
		{
		  // non-universe terms
		  // problem: you have to be exact/const when specifying the
		  // template parameters...
		  for (integer j=0;j<dims;j++)
			{
			  if(i!=j){
				
				if(R[i+1][j+1] && R[j+1][i+1]) // if this one is not specified
				  fu[i]+=-(R[i+1][j+1])(this,u)+(R[j+1][i+1])(this,u);
				else
				  cerr << "EXCHANGE RATE UNDEFINED: (" << i+1 << "," <<
					j+1 << ")" << endl;
			  }
			}
		  
		  //universe contributions
		  if(R[i+1][0]) // if this one is specified
			fu[i]-=(R[i+1][0])(this,u);
		  else
			cerr << "UNIVERSE EXCHANGE RATE UNDEFINED: (" << i+1 << ",0)" << endl;
		}

	  return fu;
	}

	number get_rate(integer from,integer to,const vect&
					u)
	{
	  if(R[from][to])
		return R[from][to](this,u);
	  else
		cerr << "EXCHANGE RATE UNDEFINED: (" << from << "," <<
		  to << ")" << endl;
	  return 0;
	  
	}
	

  private:
	exchange_rate_list R[dims+1];
  };
  
}

/** Utilities to make the definition of exchange rates easier */
#define EXCHANGE_RATE(a,b) static::MODEL::number R_a_b

#endif

/*********************************************************************
$Id: reservoirfunction.h,v 1.2 2001-07-26 12:01:22 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.1  2001/05/21 12:00:23  mpeeters
Reservoirfunction is a special kind of vectorfunction (at the moment, implemented like a heap of dung), vfamputation if a way to remove dependencies from certain vectorfuncions (to remove their dimensionality).

*********************************************************************/
