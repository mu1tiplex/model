/***************************************************************************
                          odesystem.h  -  description
                             -------------------
    begin                : Tue Jul 11 2000
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

#ifndef ODESYSTEM_H
#define ODESYSTEM_H

#include "rootscan.h"

#include "timeframe.h"
#include "random.h"
// #include "jacobian.h"
// #include "debugmacro.h"
#include "ticktock.h"
#include "integrator.h"
#include "cowner.h"
#include "cycler.h"
#include "numerictraits.h"
#include "numerictypes.h"

namespace MODEL {

  /**Integration of a system of  ordinary differential 
	 equations. The system is split into a deterministic and 
	 a stochastic part. Different integration methods take
	 these into account in different ways. You should move any
	 effect due to a nonzero mean of the noise into the deterministic part

	 \todo Things to be added (wish list):
	 Full solution of the stationary problem (using the NewtonRoot
	 solver), stability of a certain point (or the nearest stable one,
	 anyway).
  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class ODESystem : public TickTock
  {
  public:
	typedef	typename NT::number					numT;
	typedef	typename NT::vect					vect;
	typedef typename NT::vf						vf;
	typedef typename NT::matrix					matrix;
	
	
	/** dt is only a guiding value: Once linked to a timeframe, 
		will recalculate it so
		it ends up fitting a whole number of times inside the
		resolution of the timeframe */
	ODESystem( 	TimeFrame& T, vf& detpart, vf& stochpart, vect& initial,
				Integration it=Euler, time res=1E-3);

	/** The full constructor. You can pass your own integrator with parts, which
		will be deleted at the end of time */
	ODESystem( 	TimeFrame& T, vect& initial, Integrator<dims,nelem,NT>* i, time res=1e-3);

	/** The full constructor without starting point. You can pass your own integrator with parts, which
		will be deleted at the end of time */
	ODESystem( 	TimeFrame& T, Integrator<dims,nelem,NT>* i, time res=1e-3);

	/** A shortcut if you do not want to add a stochastic part */
	ODESystem( 	TimeFrame& T, vf& detpart, vect& initial,
				Integration it, time res);

	virtual ~ODESystem();

	/** Whatever should be done when the time changes. This is where
		the resolution of the ODESystem is passed to the Integrator*/
	virtual void tick()
	{
	  step(); // Take care of the TimeFrame aspect (the modulators...)
	  integ->step(current,dt); // Call the integrator
	}

	/** Get current values */
	const vect& get_current(void)
	{
	  return current;
	}

	void set_current(const vect& hereandnow)
	{
	  current = hereandnow;  
	}

	/** Overloaded for ease of use.
		\todo This should be removed as it promoted
		confusion. current() should be used instead.
	*/
	const vect&	operator()(void){return get_current();}

	/** A reset without specifying a new point */
	void reset(void) 
	{
	  reset(init);
	}
	
	
	/** Reset everything - find a stable state to start from */
	void	reset(const vect& fromhere)
	{  
	  current=fromhere; 
	  //Previously on MODEL:
	  //relax_independent();
	  
	  // But a much beter way to go is:
	  RootScan<dims,nelem,NT> find(*deter);
	  find.add_start(fromhere);
	  ScanList<dims,nelem,NT> roots=find.scan(0.,1.,1);
	  ScanList<dims,nelem,NT> goodone=roots.select(RootScan<dims>::stable);
	  
	  if (goodone.get_data()[0.].size()>0)
		{	  
		  current=goodone.get_data()[0.].front()[0];
		  cerr << "ODESystem: started at " << current << endl;
		}
	  else
		{
		  current=fromhere;
		  cerr << "ODEsystem - No stable starting point found @ " <<
			current << endl;
		}
	}

	/** Calculate statitionairy value for current parameters. Returns
		time taken to relax Note: You cannot relax to a level smaller than the
		noise level.  @param goal accurary required. */
	time relax(number goal=DEFAULT_RELAX);
	/** Use this function to relax to a stable state, without using
		the TimeFrame */
	time relax_independent(number goal=DEFAULT_RELAX);
	
	/** Scanning a certain parameter */

  private:
	// Variables
	static const number DEFAULT_RELAX=1E-7;

	/** Initial values */
	vect		init;

	/** Current values */
	vect current;

	/** Integrator */
	Integrator<dims,nelem,NT>*	integ;

	/** Deterministic Function */
	vf* deter;
	/** Stocastic Multiplier */
	vf* stoch;
  };

  // Member Functions
  //--------------------------------------------------------------------------------

  template<integer dims, typename nelem, class NT >
  ODESystem<dims,nelem,NT>::ODESystem( 	TimeFrame& T, vf& detpart, vf& stochpart, vect& initial,
										Integration it, time res) :
	TickTock(T,res), init(initial), integ(NULL), deter(&detpart), stoch(&stochpart)
  {
	// dt might have changed (because it has to fit integrally inside the top
	// TimeFrame)

	// This is not needed: dt is a member !!!!
	// dt=get_dt();
	// relax to stable state to start if possible

	current = init;
	reset();
	
	// Initialize integrator
	switch(it)
	  {
	  case Euler:	integ = new IEuler<dims,nelem,NT>(detpart,stochpart);
		break;
	  case RungeKutta: integ = new IRungeKutta<dims,nelem,NT>(detpart,stochpart);
		break;
	  case Milshtein: integ = new IMilshtein<dims,nelem,NT>(detpart,stochpart);
		break;

	  case Heun: integ = new IHeun<dims,nelem,NT>(detpart,stochpart);
		break;
	  
	  default:	throw std::logic_error("ODESystem::Integration method not implemented");
		break;
	  }
  }
  //-------------------------------------------------------------------------------

  /** with custom integrator, no starting point */

  template<integer dims, typename nelem, class NT >
	ODESystem<dims,nelem,NT>::ODESystem( 	TimeFrame& T, Integrator<dims,nelem,NT>* i, time res):
	TickTock(T,res), init(0.), integ(NULL)
  {
	// dt might have changed (because it has to fit integrally inside the top
	// TimeFrame)

	// This is not needed: dt is a member !!!!
	// dt=get_dt();
	// relax to stable state to start if possible

	integ = i;

	// get functions !!!
	deter = i->get_deterministic() ; stoch= i->get_stochastic();

	current = init;
  }

  //--------------------------------------------------------------------------------
  /** with custom integrator */

  template<integer dims, typename nelem, class NT >
	ODESystem<dims,nelem,NT>::ODESystem( 	TimeFrame& T, vect& initial, Integrator<dims,nelem,NT>* i, time res):
	TickTock(T,res), init(initial), integ(NULL)
  {
	// dt might have changed (because it has to fit integrally inside the top
	// TimeFrame)

	// This is not needed: dt is a member !!!!
	// dt=get_dt();
	// relax to stable state to start if possible

	integ = i;

	// get functions !!!
	deter = i->get_deterministic() ; stoch= i->get_stochastic();

	current = init;
	reset();	
  }

  //--------------------------------------------------------------------------------

  /** Without stoch */

  template<integer dims, typename nelem, class NT >
  ODESystem<dims,nelem,NT>::ODESystem( 	TimeFrame& T, vf& detpart, vect& initial,
										Integration it, time res) :
	TickTock(T,res), init(initial), integ(NULL), deter(&detpart)
  {
	// dt might have changed (because it has to fit integrally inside the top
	// TimeFrame)

	// This is not needed: dt is a member !!!!
	// dt=get_dt();
	// relax to stable state to start if possible

	current = init;
	reset();
	
	// Initialize integrator
	switch(it)
	  {
	  case Euler:	integ = new IEuler<dims,nelem,NT>(detpart,detpart);
		break;
	  case RungeKutta: integ = new IRungeKutta<dims,nelem,NT>(detpart,detpart);
		break;
	  
	  default:	throw std::logic_error("ODESystem::Integration method not implemented");
		break;
	  }
  }

  //--------------------------------------------------------------------------------

  template<integer dims, typename nelem, class NT >
  ODESystem<dims,nelem,NT>::~ODESystem()
  {
	delete integ;
  }

  //--------------------------------------------------------------------------------

  template<integer dims, typename nelem, class NT >
  time 
  ODESystem<dims,nelem,NT>::relax_independent(number goal)
  {
	// This will become a RungeKutta
	IRungeKutta<dims,nelem,NT> massage(*deter,*stoch); 
	ODESystem& system(*this); // alias for readability
	number change(0.);
	time chrono=0.;
	vect lastvalue(system());

	do 
	  {
		// One step
		chrono+=dt;
		vect newvalue=massage.step(current,dt);
		  		  
		// Have we converged?
		change=0.;
		for (counter i=0;i<dims;i++)
		  {
			number converge=abs((newvalue[i]-lastvalue[i])/lastvalue[i]);
			if(converge>change) change=converge; // find maximal change
		  }
		lastvalue=newvalue;
		// cout << change << "\t" << newvalue << endl;
	
	  } while (change>goal);

	return chrono; // return time spent
  }

  //--------------------------------------------------------------------------------
  template<integer dims, typename nelem, class NT >
  time 
  ODESystem<dims,nelem,NT>::relax(number goal)
  {
	ODESystem& system(*this); // alias for readability
	number change(0.);
	TimeFrame& t=get_timeframe();
	time init=t;
	vect converge(0.),lastvalue(system());

	do 
	  {
		// One step
		++t;
		vect newvalue=system();
		  
		// Have we converged?
		change=0.;
		for (counter i=0;i<dims;i++)
		  {
			number converge=abs((newvalue[i]-lastvalue[i])/lastvalue[i]);
			if(converge>change) change=converge; // find maximal change
		  }
		lastvalue=newvalue;

		
		//		DEBUGSTREAM << change <<endl;
	  } while (change>goal);

	return t-init; // return time spent
  }

} // end namespace
#endif

/*********************************************************************
$Id: odesystem.h,v 1.4 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2002/07/30 19:57:18  mpeeters
Merge of 1.0.2 branch into trunk. From now on, the trunk will be a stable (meaning it compiles and has no quirks) branch and all special stuff will be done on separate branches waiting to be merged into the stable one. A tag STABLE will be made which moves with the trunk. Also, the main trunk will be considered to be 1.0.3.

Revision 1.2.2.3  2002/07/24 09:59:34  mpeeters
Added possibility to start without predefined starting point.

Revision 1.2.2.2  2002/06/19 13:17:32  mpeeters
Possibility for custom integrator.

Revision 1.2.2.1  2001/10/15 15:08:24  mpeeters
Added staring from a different point to make restarting easier.

Revision 1.2  2001/07/26 11:53:28  mpeeters
Added constructor for non-stochastic case.
Made sure the initial stat point is found correctly.
Removed unneeded get_dt().

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.5  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.4  2000/09/29 08:55:07  mpeeters
Inserted check for stable state. Removed bug where system would choose wrong starting point.

Revision 1.3  2000/09/22 08:56:04  mpeeters
Bug fix: Rootscan was called for a 3D system always. Changed template
parameter to RootScan<dims>, as it should be.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
