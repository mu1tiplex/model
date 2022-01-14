/***************************************************************************
			singelmode.cpp
			-----------
                             
    begin                : 28/09/2000
    author               : (C) 2000 by Michael Peeters
	                                   Guy Van der Sande
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

/** This program demostrates the basic use of MODEL. Essentially, all
	you have to do is implement the vectorfunction describing the
	right-hand side of your system of ordinary differential equations.

	As system, we have chosen a single mode laser. These are the steps
	we will take:
	1) Define the system
	2) Use the routines to calculate the stationary solutions
	3) Use the routines to calculate the stability
	4) Dynamical response
	\todo Add Stochastics to this
*/
	
// 1) Define the system
// We need the headers 
#include "model/vectorfunction.h"
#include "model/odesystem.h"
#include "model/rootscan.h"
#include "model/timeframe.h"
#include "model/modulator.h"
#include "model/probe.h"
#include "model/lusolve.h"
#include "model/ssa.h"

#include <string>
#include <iostream>
#include <fstream>

/** This class describes a single mode laser. */
class SingleMode : public MODEL::VectorFunction<2> // A 2D system
{
public:
  /** Definition of base class for ease of use */
  typedef MODEL::VectorFunction<2> base;
  /** use the correct type number */
  typedef base::numT numT;

  /** Constructor - empty
	  We will hardcode all the parameter values, except the current */
  SingleMode() : j(0.) 
  {
	define_parameter("current", j);
  }
  /** Destructor - empty too - needed for good cleaning */
  ~SingleMode() {}

  /** Copy constructor - don't forget to call the base */
  SingleMode(const SingleMode& sm) : base(sm),j(sm.j) {
	define_parameter("current", j);}

  /** Virtual Copy Constructor */
  virtual base* clone () const 
  {
	return new SingleMode(*this);
  }

  /** Assignment operator */
  const SingleMode& operator=(const SingleMode& sm)
  {
	if(this != &sm) 
	  {
		j=sm.j;
	  }
	return *this;
  }
  
  /** Immplementing the function */
  virtual const vect& function(vect& fu,const vect& u)
  {
	// Making life easier
	const numT& p(u[0]); // or p=u[0] - photon density
	const numT& n(u[1]); // carrier density
	
	numT& dpdt(fu[0]);
	numT& dndt(fu[1]);
	
	// The parameters
	const numT rho = 1E-3; // ratio of timescales
	const numT g   = 1.1; // gain

	// The system
	dpdt = (1./rho)*(g*(n-1.)*p-p ) + 1E-9;
	dndt = j - n - n*p;
	
	return fu;
  }

private:
  /** The current */
  MODEL::number j;
};

/** A constant noise class */
  class ConstNoise : public MODEL::VectorFunction<2>
  {
  public:
	typedef	MODEL::VectorFunction<2>	base;
	
	ConstNoise(const vect& r=vect(0.) ) : cr(r) {}
	virtual ~ConstNoise() {};
		
	/** Copy constructor */
	ConstNoise(const ConstNoise& c) : base(c) {copy(c);}
		
	/** Assignment */
	const ConstNoise& operator=(const ConstNoise& c)
	{	if(this!=&c) {
	  base::operator=(c); copy(c);} return *this;}
		
	/** Virtual Copy Constructor */
	virtual ConstNoise* clone () const
	{ return new ConstNoise(*this); }

  private:
	void	copy(const ConstNoise& c) {cr=c.cr;}

	virtual const vect& function(vect& fu, const vect& u)	
	{
	  fu=cr;
	  return fu;
	}	

  private:
	vect	cr;
  };

// good practice - define a prototype
int main (void);

/** The main loop */
int main (void)
{
  using namespace MODEL;

  // 1) Define the system (again)
  cerr << "Defining system." << endl;
  
  SingleMode vcsel;

  // 2) Calculate the stationary solutions
  // These are roots of the system
  cerr << "Calculating stationary solutions." << endl;

  RootScan<2> rs(vcsel,"current");   // A root finder
  
  // the output is formatted in a scanlist
  ScanList<2> solutions=rs.scan(0.,6.,100); // from 0 to 2, 100 points
  
  // Save it all to a file (gnuplot readable)
  string filename="stationarystable.dat";
  ofstream data_out(filename.c_str());

  // raw and unformatted
  //solutions.print_raw(cout);

  // 3) Use the stability (which is already calculated !)
  cerr << "Selecting stable solutions." << endl;
  solutions.select(RootScan<2>::stable).print_list(data_out);

  // 4) A Dynamical simulation: a step from below to above threshold
  // a) Create a timeframe: specify the resolution
  cerr << "Doing a step (large signal dynamics)." << endl;

  TimeFrame t(2); // a resolution of 3 something (ns) 
  // This means t++ increase the time by 10 ns.
  // This also means any other activity linked to
  // this timeframe (measurements etc...) have this granularity   

  // b) create the ODESystem
  // An ODESystem has the following parameters:
  // i. a TimeFrame (t)
  // ii. a deterministic part (vcsel) 
  // iii. a multiplier for stochastics (noise)
  ConstNoise noise(0.);
  // iv. the initial values of the variables p and n (start)
  NumVector<2> start(0.); // superfluous
  // v. the method of integration (Euler=default, RungeKutta,
  //    Milshtein)
  // vi. the resolution of the integrator
  ODESystem<2> dynvcsel(t,vcsel,noise,start,RungeKutta,1E-3);
  //c) Define a modulation
  StepMod modulatecurrent(dynvcsel, vcsel, "current", 1.5,2.5, 10.);
  
  //d) Define measurement (all points are saved)
  ODEProbe<2> everything(dynvcsel,dynvcsel,"stepmodulation.dat");

  //e) Run the time. We write out the points for visual inspection,
  // but this is not necessary, actually

  while (t<30.)
	cerr << ++t << "\t" << dynvcsel() << endl;
  
  //5) Calculate small signal response, for different bias points
   cerr << "Calculating small signal response." << endl;

  SSA<2> small(vcsel,"current"); // what is the parameter we call bias
  ScanList<2> bode; // all results
  
  // one value
  for (number j=2.5;j<5.5;j+=0.5)
	{
	  if(small.set_stat_param(j)!=0)
		// set the param, this also find solutions
		small.set_stat_point(); // take the first point
	  else // find one in the list
		small.set_stat_point(*solutions.select(RootScan< 2
											   >::stable).get_solution(j).begin()->begin());

	  bode+=small.calc_response_norm(0.1,1e6,250); // add the response
	}

  string filenamessa="ssa.dat";
  ofstream ssa_out(filenamessa.c_str());
   
  bode.print_list(ssa_out); // write it out

  return 0;
}





/*********************************************************************
$Id: singlemode.cpp,v 1.4 2001-07-26 11:58:51 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2001/04/10 13:11:46  mpeeters
Added some plotting files. Changed simulation timescale.

Revision 1.1  2000/09/29 09:01:54  mpeeters
Small changes in cvr3dfunc.
Added a tutorial directory with an example.

*********************************************************************/
