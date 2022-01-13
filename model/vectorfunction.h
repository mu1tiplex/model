/***************************************************************************
                          vectorfunction.h  -  description
                             -------------------
    begin                : Thu Apr 20 2000
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

#ifndef VECTORFUNCTION_H
#define VECTORFUNCTION_H

#include "numerictypes.h"
#include "numerictraits.h"
#include "numvector.h"
#include "parameter.h"
#include <stdexcept>
#include <map>
#include <string>
#include <algorithm>


namespace MODEL {

  using namespace std;

  /**	Base class for vectorfunctions: derived from this
		if you want to implement one. Just override/implement
		the function() part. The function has as input: a vector u.
		The time dependence of any parameter is left to the implementor (if needed).

		\todo Parameters should be copied in copy constructor
		-but this would break some existing stuff !
  */


  // Previously template <integer dims, class NT = NumericTraits<number,dims> >
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class VectorFunction
  {
  public:
	/** Really weird def. to be able to use NT outside */
	typedef NT NT_outside;
	typedef	typename NT::number					numT;
	typedef	typename NT::vect					vect;
	typedef typename NT::vf						vf;
	typedef typename NT::matrix					matrix;
       		
	VectorFunction() {}
	virtual ~VectorFunction() {};
	
	// Just to make sure no errors are made
	/** \bug gigantic BUG: parlist should never be copied - this is
		just a quick fix */
	VectorFunction(const VectorFunction& vf) {}
	//	const VectorFunction<dims,nelem> operator=(const VectorFunction<dims,nelem>& vf)
	// {if(this!=&vf){} return *this;}
			
	/** Virtual Copy Constructor */
	virtual VectorFunction* clone () const =0;
				
	/** Overloaded operator(), so the object can be presented as a
		function. */						
	vect	operator()(const vect& u)
	{ vect temp(0.); return function(temp,u); }
								
	/** Implement this function to create the return vector.
		fu is a reference to where the values should be stored.
		This is for efficiency in internal routines.
		Inherited classes should return fu too.*/
	virtual const vect& function(vect& fu,const vect& u) = 0;	

	void define_parameter(const string& name, number&
						  p){parlist[name]=&p;}

   	number& get_parameter(const string& name)
	{
	  map<string,number*>::iterator  end=parlist.end();
	  map<string,number*>::iterator  found=parlist.find(name);

	  if (found==end) throw std::logic_error("VectorFunction::Parameter not found");
	  
	  return *(found->second);
	}

   	number& get_parameter(const string& name) const
	{
	  map<string,number*>::iterator  end=parlist.end();
	  map<string,number*>::iterator  found=parlist.find(name);

	  if (found==end) throw
						std::logic_error("VectorFunction::Parameter not found");
	  
	  return *(found->second);
	}

  public:
 	map<string,number*>& get_parlist(void)
	{
	  return parlist;
	}

 	const map<string,number*>& get_parlist(void) const
	{
	  return parlist;
	}

  private:
	map<string,number*> parlist;	
  };

  // Used To make older code happy */
  // typedef VectorFunction<MODEL::real,3> RealVectorFunction;

} // end namespace MODEL
#endif

/*********************************************************************
$Id: vectorfunction.h,v 1.3 2002-07-30 19:57:19 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2.2.2  2002/07/30 08:59:15  mpeeters
Spring cleaning commit - administrative stuff.

Revision 1.2.2.1  2001/10/15 15:09:36  mpeeters
Finally removed multiline string literal.

Revision 1.2  2001/08/23 20:24:00  mpeeters
Modified docs to provide correct version info.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
