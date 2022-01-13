/***************************************************************************
			vfamputation.h
			-----------
                             
    begin                : Thu Apr 26 08:33:27 CEST 2001
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

#ifndef VFAMPUTATION_H
#define VFAMPUTATION_H

#include "vectorfunction.h"
#include "newtonroot.h"
#include "utility.h"
#include <vector>

namespace MODEL 
{

  /** A new vectorfunction comes out, and one of the dynamical
	  variables has been changed into a parameter */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class VFAmputation : public VectorFunction<dims,nelem,NT>
  {
  public:
	typedef typename NumericTraits<nelem,dims+1>::vf greatervf;
	typedef typename NumericTraits<nelem,dims+1>::vect greatervect;

  public:
	/** Which VF, what variable */
	VFAmputation(greatervf& func, const integer& var, bool own=false)
	  : f((own)?func.clone():&func), owned(own), dynpar(var) 
	{
	  define_parameter("amputated",paramvar);
	  // add subvariables
	  map<string,number*>& subparlist=func.get_parlist();
	  map<string,number*>::iterator run=subparlist.begin();
	  map<string,number*>::iterator endrun=subparlist.end();
	  
	  while(run!=endrun){
		// cerr << "ADDING PARAMETER : " << run->first << endl;
		define_parameter(run->first,*(run->second));
		++run;
	  }
	  
	}
	
	/** Create one, with a bump already specified */
	~VFAmputation() {if(owned)delete f;}
	
	/** copy 
		\bug parameters should be copied, too*/
	VFAmputation(const VFAmputation& vfb) : f(vfb.f), owned(false),
											dynpar(vfb.dynpar)  
	{
	  const map<string,number*>& subparlist=vfb.get_parlist();
	  map<string,number*>::const_iterator run=subparlist.begin();
	  map<string,number*>::const_iterator endrun=subparlist.end();
	  
	  while(run!=endrun){
		define_parameter(run->first,*(run->second));
		++run;
	  }
	}
		
	/** Virtual Copy Constructor */
	virtual VFAmputation* clone () const {return new VFAmputation(*this);}
	
	/** Implement this function to create the return vector.
		fu is a reference to where the values should be stored.
		This is for efficiency in internal routines.
		Inherited classes should return fu too.*/
	virtual const vect& function(vect& fu,const vect& u)
	{			
	  greatervect v=convert(u);
	  greatervect fv=convert(fu);
	  
	  fu=convert(f->function(fv,v));	  
	  return fu;	  
	}

	/** converts vectors n=>n+1. Moves parameter into variable*/
	greatervect convert(const vect& v)
	{
	  greatervect w;
	  for (integer i=0; i<dims+1; i++)
		if(i<dynpar) w[i]=v[i] ;
		else if (i==dynpar) w[i]=paramvar;
		else if (i>dynpar) w[i]=v[i-1] ;
	  
	  return w;
	}
	
	/** converts vectors n+1+>n. Updates parameter values ?*/
	vect        convert(const greatervect& v)
	{
	  vect w;
	  for (integer i=0; i<dims+1; i++)
		if(i<dynpar) w[i]=v[i] ;
		else if (i==dynpar) paramvar=v[i];
		else if (i>dynpar) w[i]=v[i+1] ;
	  
	  return w;
	}

  private:
	greatervf*		f;
	bool	owned;

	integer dynpar;
	number paramvar;

	PRIVATE_ASSIGN(VFAmputation);
	
  };


}// end namespace


#endif
/*********************************************************************
$Id: vfamputation.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************


#include "vectorfunction.h"
#include "newtonroot.h"
#include <vector>

$Log: not supported by cvs2svn $
Revision 1.1  2001/05/21 12:00:27  mpeeters
Reservoirfunction is a special kind of vectorfunction (at the moment, implemented like a heap of dung), vfamputation if a way to remove dependencies from certain vectorfuncions (to remove their dimensionality).

*********************************************************************/
