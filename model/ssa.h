/***************************************************************************
			ssa.h
			-----------
                             
    begin                : Tue Apr 24 17:56:51 CEST 2001
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

#ifndef SSA_H
#define SSA_H

#include "utility.h"
#include <stdexcept>
#include "numerictraits.h"
#include "numerictypes.h"
#include "jacobian.h"
#include "rootscan.h"

namespace MODEL {

  /** Feed me a Vectorfunction and I can give you a small signal
	  analysis for any parameter. Actually, I will calculate it for
	  all parameters */

  template <integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class SSA
  {
  public:
	typedef typename NT::number	numT;
	typedef	typename NT::vect	vect;	
	typedef typename NT::matrix	matrix;
	typedef typename Jacobian<2*dims>::matrix system;
	typedef typename Jacobian<2*dims>::vect input;	
	typedef typename NumericTraits<complex,dims>::vect response;
	typedef response cinput;
	typedef typename NT::vf						vf;
	
  public:
	/** For my_sys, for this input parameter */
	SSA(vf& my_sys, const string& parm);

	/** Set the stationary parameter value around which we will
		modulate. This will return the number of stable stationary points the
		system have been able to find. Use get_stat_points to get a
		list in ScanList format */
	counter set_stat_param(const number& value);
	
	/** Once you have some points, choose one */
	void set_stat_point(const counter& po=1);

	/** Or (DANGEROUS) set it explicitly. Nobody guarantees this will
		work (and give meaningfull results) if you do not have a valid point */
	void set_stat_point(const vect& right_there);

	/** get a list of the stationary points */
	const ScanList<dims,NT>& get_stat_points(void);
	
	/** get a single response point. For that detailed analysis */
	response calc_point(const number& omega);
	
	/** get a list of responses, equally spaced in a log scale */
	ScanList<dims,complex> calc_response(const number& from, const
										 number& to, const counter&
										 n);

	/** get a list of responses real,complex, equally spaced in a log scale */
	ScanList<dims> calc_responseRC(const number& from, const
										 number& to, const counter&
										 n);

	/** get a list of responses, but this this just the norm, no phase
		information, please */
	ScanList<dims> calc_response_norm(const number& from, const
									  number& to, const counter&
									  n);

	/** get a list of responses, but this this just the phase, no amplitude
		information, please */
	ScanList<dims> calc_response_phase(const number& from, const
									  number& to, const counter&
									  n);

	/** set the depency manually. This is needed if we want to use
		these routines for other purposes, like noise analysis, where
		the sources are not necessarily the inputs */
	void set_dep(const vect& input)
	{
	  linpar = 0.;
	  for (integer i=0;i<dims;i++)
		linpar[i]=input[i];
	}

	/** set the depency manually when it is complex. This is needed if we want to use
		these routines for other purposes, like noise analysis, where
		the sources are not necessarily the inputs */
	void set_dep(const cinput& input)
	{
	  linpar = 0.;
	  for (integer i=0;i<dims;i++)
		{
		  linpar[i]=input[i].real();
		  linpar[dims+i]=input[i].imag();
		}	  
	}
  
  private:
	/** Calculate the parameter dependency: partial FD */
	void calc_dep(void);

	/** fill omega depency in system matrix A */
	void put_omega(const number& omega);

  private:
	/** Tells all calculations that all is well and they can go ahead
	 */
	bool all_is_well;

	/** The function/jacobian in question */
	Jacobian<dims>* J; 

	/** The actual jacobian of the vectorfunction */
	matrix j; // small j: nice and confusing
	
	/** A matrix twice as large for the complex problem */
	system A;
	
	/** The rootscanner */
	RootScan<dims,nelem>* rs;

	/** The list of stationary points */
	ScanList<dims,nelem> statp;

	/** The stationary point */
	vect here;

	/** The parameter */
	ParameterP p;

	/** The linear approximation of the way the equations respond to
		the input. Also twice as large, but only the top half in used  */
	input linpar;

  private:
	// to avoid silly errors
	NO_COPY(SSA);
	
  };
  
  // ------------------------------------------------------------

  template <integer dims, typename nelem, class NT >
  SSA<dims,nelem,NT>::SSA(vf& my_sys, const string& parm)
  {
	all_is_well=false;
	
	// get the jacobian, copy the function as is
	J=new Jacobian<dims>(my_sys,1E-6,true);

	// get the parameter
	p=&(J->get_function().get_parameter(parm));

	// init the scanner
	rs=new RootScan<dims,nelem>(J->get_function());	
  }

  template <integer dims, typename nelem, class NT >
  counter SSA<dims,nelem,NT>::set_stat_param(const number& value)
  {
	*p=value;

	statp=rs->scan(0.,1.,1); // these are dummy values
	statp=statp.select(RootScan<dims,nelem>::stable);

	// statp.print_raw(cout);

	// cerr << "SSA: started number of solutions " << statp.get_data()[0.].size() << endl;

	return statp.get_data()[0.].size();
	// watch out here: you might find only one (wrong point), because
	// we are NOT scanning. This is a problem in RootScan.
  }
  
  template <integer dims, typename nelem, class NT >
  void SSA<dims,nelem,NT>::set_stat_point(const counter& po)
  {
	// go to correct point
	if(statp.get_data()[0.].size()>0)
	  {
		// take the first parameter value in the list
		typename ScanList<dims,nelem>::parpointlist::iterator
		  getit=statp.get_data()[0].begin();
	
		// go to the right value
		for(counter i=po;i>1;i--) getit++;
		
		// the first value is the point...
		here=(*getit)[0];
		
		// now that we know where we, calculate the approximations
		calc_dep();

		// now we can do something
		all_is_well=true;		
	  }
	
  }
  template <integer dims, typename nelem, class NT >
  void SSA<dims,nelem,NT>::set_stat_point(const vect& right_there)
  {
		here=right_there;
		
		// now that we know where we, calculate the approximations
		calc_dep();

		// now we can do something
		all_is_well=true;
  }
  
  template <integer dims, typename nelem, class NT >
  typename SSA<dims,nelem,NT>::response SSA<dims,nelem, NT>::calc_point(const number& omega)
  {
	if(all_is_well) {
	  
	  put_omega(omega);
	  LUSolve<2*dims> solA(A,true); // make a copy
	
	  // copy into correct form
	  input t1 (solA(linpar));
	  response t2;
	  for (integer i=0;i<dims;i++)
		{
		  t2[i]=t1[i]+I*t1[i+dims];
		}

	  return t2;
	}
	
  }

  template <integer dims, typename nelem, class NT >
  ScanList<dims,complex>  
  SSA<dims,nelem, NT>::calc_response(const number& from, const
									 number& to, const counter&
									 n)
  {
	if(all_is_well) {
	  
	  ScanList<dims,complex> res;

	  number scaler=pow(to/from,1./number(n));
	  for (number omega=from;omega<to;omega*=scaler)
		{
		  res.add_point(calc_point(omega));
		  res.add_param(omega);
		}
	  return res;
	}
	
  }

  template <integer dims, typename nelem, class NT >
  ScanList<dims>
  SSA<dims,nelem, NT>::calc_responseRC(const number& from, const
									 number& to, const counter&
									 n)
  {
	if(all_is_well) {

	  ScanList<dims> res;

	  number scaler=pow(to/from,1./number(n));
	  for (number omega=from;omega<to;omega*=scaler)
		{
			response full;
			full=calc_point(omega);

			vect re,im;

			for(int i=0;i<dims;i++)
			{
				re[i]=real(full[i]);
				im[i]=imag(full[i]);
			}

		  res.add_point(re);
		  res.add_point(im);
		  res.add_param(omega);
		}
	  return res;
	}

  }

template <integer dims, typename nelem, class NT >
  ScanList<dims>
  SSA<dims,nelem, NT>::calc_response_norm(const number& from, const
										  number& to, const counter&
										  n)
  {
	if(all_is_well) {
	  
	  ScanList<dims> res;

	  number scaler=pow(to/from,1./number(n));
	  for (number omega=from;omega<to;omega*=scaler)
		{
		  response t1=calc_point(omega);
		  vect t2;
		  for(integer i=0;i<dims;i++)
			t2[i]=norm(t1[i]);

		  res.add_point(t2);
		  res.add_param(omega);
		}
	  return res;
	}

  }

  template <integer dims, typename nelem, class NT >
  ScanList<dims>  
  SSA<dims,nelem, NT>::calc_response_phase(const number& from, const
										  number& to, const counter&
										  n)
  {
	if(all_is_well) {
	  
	  ScanList<dims> res;

	  number scaler=pow(to/from,1./number(n));
	  for (number omega=from;omega<to;omega*=scaler)
		{
		  response t1=calc_point(omega);
		  vect t2;
		  for(integer i=0;i<dims;i++)
			t2[i]=arg(t1[i]);

		  res.add_point(t2);
		  res.add_param(omega);
		}
	  return res;
	}
	
  }
  
  template <integer dims, typename nelem, class NT >
  void SSA<dims,nelem, NT>::calc_dep(void)
  {
	// generate the Jacobian accurately (slow)
	j=J->calculate_accurate(here);

	// And copy it to the larger system
	for (int i=0;i<dims;i++)
	  {
		for (int k=0;k<dims;k++)
		  {
			A[i][k] = -j[i][k];
			A[dims+i][dims+k] =-j[i][k];
		  }  
	  }
	
	// find the (linear) dependecies of the equations to the
	// (tiny/small/microscopic input). This is done in a way very
	// similar to the Jacobian.
	// it is store in linpar
	linpar=0.; // just to be safe

	// value now
	vf& f=J->get_function();
	vect    fu ( f(here) );

	number pdp=*p;

	/** \todo this should be modifiable */
	const number epsilon = 1;
	number dp=epsilon*abs(pdp);
		
	if(dp==0.0) dp = epsilon;		// avoid numerical error	
	pdp += dp;			// "
	dp = pdp - *p;		// "

	// get value at new point
	number oldp=*p;
	*p=pdp;

	Parameter herep = f.get_parameter("current");

	vect	fudp(f(here));
	
	*p=oldp; // restore parameter
		
	for(integer i=0;i<dims;i++)
	  {
		linpar[i]=(fudp[i]-fu[i])/dp;	
	  }  
  }
  
  template <integer dims, typename nelem, class NT >
  void SSA<dims,nelem, NT>::put_omega(const number& omega)
  {
	for (int i=0;i<dims;i++)
	  {
		A[i][dims+i] = -omega;
		A[dims+i][i] =omega;
	  }
  }
  
}
// end MODEL
#endif 

// CVS Test: hello Guy.

/*********************************************************************
$Id: ssa.h,v 1.4 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2002/11/18 14:10:30  mpeeters
Removed unnecessary ssa.cpp. Added DEBUG in rootscan catch{} routine. Added RC and phase info to the SSA routines.

Revision 1.2  2002/07/30 19:57:19  mpeeters
Merge of 1.0.2 branch into trunk. From now on, the trunk will be a stable (meaning it compiles and has no quirks) branch and all special stuff will be done on separate branches waiting to be merged into the stable one. A tag STABLE will be made which moves with the trunk. Also, the main trunk will be considered to be 1.0.3.

Revision 1.1.2.2  2001/08/30 14:48:27  mpeeters
Testing CVS Branching.

Revision 1.1.2.1  2001/08/30 14:01:22  mpeeters
Fixed stupid warning about extra token (due to #endif SSHA_H).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.1  2001/05/21 11:56:43  mpeeters
Commited lots of small changes (too many to tell) at the previous commit (accidentally). In this commit, we add the small signal analysis stuff. Whoohoo.

*********************************************************************/
