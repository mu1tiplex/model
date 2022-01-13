/***************************************************************************
                          integrator.h  -  description
                             -------------------
    begin                : Tue Aug 12 2000
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

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "numerictypes.h"
#include "numerictraits.h"
#include "vectorfunction.h"
// #include "odesystem.h"
#include "random.h"
#include "jacobian.h"

namespace MODEL {

  // predef
  template<integer dims, typename nelem, class NT>
  class ODESystem;

  /** Integration types */
  enum Integration {Euler, RungeKutta, Milshtein,Heun};

  /** Integrator: take a step dt from here to the next point */
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class Integrator
  {
    friend class ODESystem<dims,nelem,NT>;

  public:
    typedef	typename NT::number					numT;
    typedef	typename NT::vect					vect;
    typedef typename NT::vf						vf;
    typedef typename NT::matrix					matrix;

    /** Set the integrator and tell it what to do. At the moment any
	special cases with for example correlated noise should be
	taken care of when we derive. I agree that storing a pointer
	to the function objects inside the base class might not be
	very kosher, but that's the way it is for now. Cleanup is NOT
	a priority.
	@param d is the deterministic vectorfunction.
	@param s is the multiplicative stochastic vectorfunction.
    */
    Integrator(vf& d, vf& s)
      :  deter(&d), stoch(&s){}

    virtual ~Integrator(){};

    /** Implement this one to take a step. takes the current point and
	modifies it.  Returns a ref to the modified point
	@param dt is the step to take. This makes it possible to generate a
	adaptive stepsize if needed.

	I wonder if the fact that an Integrator has a dt makes it a
	candidate for TickTock status? */
    virtual	vect& step(vect& current, time& dt)=0;

    /** Return current point */
    vect&	currentpoint(void) {return ODESystem<dims,nelem,NT>::current;}

  private:
    /** Get functions - not clean, but necessary */
    vf* get_deterministic(void)
    {
      return deter;
    }

    vf* get_stochastic(void)
    {
      return stoch;
    }

  protected:
    /** Storage for the functions */
    vf* deter;
    /** Storage for the functions */
    vf* stoch;

  };
  //------------------------------------------------------------

  /** The simplest integrator possible.
      If stochastic part is specified: Euler-Maryuama
      Used mainly for testing purposes.
      The function is evaluated with the old values */
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class IEuler : public Integrator<dims,nelem,NT>
  {
  public:
    typedef Integrator<dims,nelem,NT> base;
    typedef typename base::vf vf;
    typedef typename base::vect vect;

  public:
    /** Pass it on: nothing special needs to be done. */
    IEuler(vf& d, vf& s)
      : Integrator<dims,nelem,NT>(d,s) {}

    /** Agreed: this is very simple (and error prone) */
    virtual	vect& step(vect& current, time& dt)
    {
      vect noise;
      for(counter i=0;i<dims;i++) noise[i]=rnd();
      vect g=(*Integrator<dims,nelem,NT>::stoch)(current);
      for(counter i=0;i<dims;i++) g[i]*=noise[i];

      return current+=dt*(*Integrator<dims,nelem,NT>::deter)(current)+sqrt(dt)*g;
    }

    Normal& get_random()
    {
      return rnd;
    }

  private:
    /** Normal random */
    Normal rnd; // will be seeded with time()

  };


  //-------------------------------------------------------------------------
  /** A fourth order Runge-Kutta Integrator. Very robust for
      deterministic problems.  Even better would be a variable step size
      (but that defeats the problem...) */
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class IRungeKutta : public Integrator<dims,nelem,NT>
  {
  public:
    typedef Integrator<dims,nelem,NT> base;
    typedef typename base::vf vf;
    typedef typename base::vect vect;
  public:
    IRungeKutta(vf& d, vf& s)
      : Integrator<dims,nelem,NT>(d,s) {}

    virtual	vect& step(vect& current, time& dt)
    {
      // Runge-Kutta needs 4 evaluations:
      // 1) at the current point
      // 2) at a halfway point, using the values of 1)
      // 3) at a halfway point, using the values of 2)
      // 4) at the end point, using the values of 3)

      vect k1,k2,k3,k4;
      number dt2=dt/2.;

      // Step 1)

      k1=(*Integrator<dims,nelem,NT>::deter)(current);
      k1*=dt2;

      // Step 2)
      // halfway in t

      vect pos(k1);
      pos+=current; //halfway in u

      k2=(*Integrator<dims,nelem,NT>::deter)(pos);
      k2*=dt2;

      // Step 3)
      pos=k2;
      pos+=current;

      k3=(*Integrator<dims,nelem,NT>::deter)(pos);
      k3*=dt2;
      // Step 4)
      //another half step

      pos=k3; pos+=current; // BUG FIX: Forgot to add current point
      k4=(*Integrator<dims,nelem,NT>::deter)(pos);
      k4*=dt2;

      // Now calculate next point. Caveat
      // emptor, this is very slow (compared to what it could be).
      return current+=(k1+2.*k2+2.*k3+k4)/3.;
    }

  };

  //----------------------------------------------------------------------
  /** An integrator that can actually do stochastics - the counterpart
      of Euler
      @see Random */
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class IMilshtein : public Integrator<dims,nelem,NT>
  {
  public:
    typedef Integrator<dims,nelem,NT> base;
    typedef typename base::vf vf;
    typedef typename base::vect vect;
    typedef typename base::matrix matrix;

  public:
    IMilshtein(vf& d, vf& s)
      : Integrator<dims,nelem,NT>(d,s){}

    virtual	vect& step(vect& current, time& dt)
    {
      // The Milshtein method is the basic algo for noisy problems
      // => Needs a random number generator !

      vect oldcurrent(current);

      // Deterministix q
      vect q=(*Integrator<dims,nelem,NT>::deter)(oldcurrent);
      current+=dt*q;

      // Stochastix g
      // Generate numbers: noise is _NOT_ correlated between modes.
      vect noise;
      for(counter i=0;i<dims;i++) noise[i]=rnd();

      // Additive component (Euler)
      number sqdt=sqrt(dt);
      vect g=(*Integrator<dims,nelem,NT>::stoch)(oldcurrent);

      // Multiply each component with its noise
      vect ng(g);
      for(counter i=0;i<dims;i++) ng[i]*=noise[i];

      // Result
      current+=sqdt*ng;

      // Extra Millstein
      Jacobian<dims,NT> J(*Integrator<dims,nelem,NT>::stoch);
      matrix dgdy=J.calculate(oldcurrent,g);

      vect dgs(0.);
      for(counter i=0;i<dims;i++)
	for (counter j=0;j<dims;j++)
	  {
	    number gg=0.5*dgdy[i][j]*g[j];
	    dgs[i]-=gg; // Ito rules !
	    if(i==j)dgs[i]+=gg*(noise[i]*noise[j]);
	  }

      current+=dt*dgs;
      return current;
    }

    /** To give users a possibility to set the seed */
    Normal& get_random()
    {
      return rnd;
    }


  private:
    /** Normal random */
    Normal rnd; // will be seeded with time()
  };

  /** The IHeun integrator is a STRATONOVICH INTEGRATOR. However, you
      can force it to make the correction needed to integrate for the
      Ito solution by using IHeunIto
  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class IHeun : public Integrator<dims,nelem,NT>
  {
  public:
    typedef Integrator<dims,nelem,NT> base;
    typedef typename base::vf vf;
    typedef typename base::vect vect;
    typedef typename base::matrix matrix;
  public:
    IHeun(vf& d, vf& s)
      : Integrator<dims,nelem,NT>(d,s){}

    virtual	vect& step(vect& current, time& dt)
    {
      // The Heun method is the best easy algo for noisy problems
      // => Needs a random number generator !

      // u is N(0,1)

      // k=h.q(t,x)
      // l=sqrth.u.g(t,x)
      //
      // qavg = (1/2)*( q(t,x) + q(t+h,x+k+l) )
      // gavg = (1/2)*( g(t,x) + g(t+h,x+k+l) )
      //
      // Watch out: we might have to take a step for the modulators
      //
      // xnew = x + h*qavg + sqrth*u*gavg

      vect oldcurrent(current);

      // Some to be optimized out temps
      const number h=dt;
      const number sh=sqrt(dt);

      vect q = (*Integrator<dims,nelem,NT>::deter)(oldcurrent);
      vect g = (*Integrator<dims,nelem,NT>::stoch)(oldcurrent);

      // Generate numbers: noise is _NOT_ correlated between modes.
      vect u;
      for(counter i=0;i<dims;i++) u[i]=rnd();

      /** \todo Instead of generating noise multiplied by a funtion,
	  pass the noise to another function. This way, the function can
	  change correlations etc if needed. */

      // uneff: vect k=h*q;
      // uneff: vect l=sh*g;

      // generate next point
      for(counter i=0;i<dims;i++)
	{
	  oldcurrent[i]+=h*q[i]+sh*g[i]*u[i];
	  // The stochastic part should be sqrt(2D) when the autocorrelation
	  // is given by <F,F>=2D.delta; 2D is the variance (in fact)
	}

      // uneff:		  l[i]*=u[i];

      // uneff: oldcurrent+=k;
      // uneff: oldcurrent+=l;

      // vect qavg=0.5* h * ( q + (*deter)(oldcurrent) );
      // vect gavg=0.5* sh * ( g + (*stoch)(oldcurrent) );

      vect qn=(*Integrator<dims,nelem,NT>::deter)(oldcurrent);
      vect gn=(*Integrator<dims,nelem,NT>::stoch)(oldcurrent);

      for(counter i=0;i<dims;i++)
	{
	  current[i]+=0.5*( h*( q[i] + qn[i]) + sh*u[i]*(g[i]+gn[i]));

	  // uneff:  gavg[i]*=u[i];
	}

      // uneff: current += qavg;
      // uneff: current += gavg;

      return current;
    }

    /** To give users a possibility to set the seed */
    Normal& get_random()
    {
      return rnd;
    }

  private:
    /** Normal random */
    Normal rnd; // will be seeded with time()
  };

  /** The Heun integrator with the Ito correction (numerically, first
      order) added
      \todo: This is not yet finished
  */

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class IHeunIto : public Integrator<dims,nelem,NT>
  {
  public:
    typedef Integrator<dims,nelem,NT> base;
    typedef typename base::vf vf;
    typedef typename base::vect vect;
    typedef typename base::matrix matrix;
  public:

    IHeunIto(vf& d, vf& s)
      : Integrator<dims,nelem,NT>(d,s)
    {
      // Create the Jacobian for the stochastic part
      jacostoch = new Jacobian<dims,NT>(s);
    }

    virtual	vect& step(vect& current, time& dt)
    {
      // The Heun method is the best easy algo for noisy problems
      // => Needs a random number generator !

      // u is N(0,1)

      // k=h.q(t,x)
      // l=sqrth.u.g(t,x)
      //
      // qavg = (1/2)*( q(t,x) + q(t+h,x+k+l) )
      // gavg = (1/2)*( g(t,x) + g(t+h,x+k+l) )
      //
      // Watch out: we might have to take a step for the modulators
      //
      // xnew = x + h*qavg + sqrth*u*gavg


      vect oldcurrent(current);

      // Some to be optimized out temps
      const number h=dt;
      const number sh=sqrt(dt);

      vect q = (*Integrator<dims,nelem,NT>::deter)(oldcurrent);
      vect g = (*Integrator<dims,nelem,NT>::stoch)(oldcurrent);

      // first, we need to calculate the jacobian to be able to apply
      // the Ito correction (we also do this the second time)
      // We can use the code from the Millshtein algo to do this.

      matrix db;
      db=jacostoch.calculate(oldcurrent,g);

      // apply the correction
      // q-= ... ADD CODE HERE.


      // Generate numbers: noise is _NOT_ correlated between modes.
      vect u;
      for(counter i=0;i<dims;i++) u[i]=rnd();

      /** \todo Instead of generating noise multiplied by a funtion,
	  pass the noise to another function. This way, the function can
	  change correlations etc if needed. */

      // uneff: vect k=h*q;
      // uneff: vect l=sh*g;

      // generate next point
      for(counter i=0;i<dims;i++)
	{
	  oldcurrent[i]+=h*q[i]+sh*g[i]*u[i];
	  // The stochastic part should be sqrt(2D) when the autocorrelation
	  // is given by <F,F>=2D.delta; 2D is the variance (in fact)
	}

      // uneff:		  l[i]*=u[i];

      // uneff: oldcurrent+=k;
      // uneff: oldcurrent+=l;

      // vect qavg=0.5* h * ( q + (*deter)(oldcurrent) );
      // vect gavg=0.5* sh * ( g + (*stoch)(oldcurrent) );

      vect qn=(*Integrator<dims,nelem,NT>::deter)(oldcurrent);
      vect gn=(*Integrator<dims,nelem,NT>::stoch)(oldcurrent);

      for(counter i=0;i<dims;i++)
	{
	  current[i]+=0.5*( h*( q[i] + qn[i]) + sh*u[i]*(g[i]+gn[i]));

	  // uneff:  gavg[i]*=u[i];
	}

      // uneff: current += qavg;
      // uneff: current += gavg;

      return current;
    }

  private:
    /** Normal random */
    Normal rnd; // will be seeded with time()

    /** The jacobian */
    Jacobian<dims,NT>* jacostoch;
  };


  /** A Heun integrator with simple correlations : a vectorfunction*/

  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class IHeunSimpleCorr : public Integrator<dims,nelem,NT>
  {
  public:
    typedef Integrator<dims,nelem,NT> base;
    typedef typename base::vf vf;
    typedef typename base::vect vect;
    typedef typename base::matrix matrix;
  public:
    IHeunSimpleCorr(vf& d, vf& s, vf& transform)
      : Integrator<dims,nelem,NT>(d,s), trans(&transform){}

    virtual	vect& step(vect& current, time& dt)
    {
      // The Heun method is the best easy algo for noisy problems
      // => Needs a random number generator !

      // u is N(0,1)

      // k=h.q(t,x)
      // l=sqrth.u.g(t,x)
      //
      // qavg = (1/2)*( q(t,x) + q(t+h,x+k+l) )
      // gavg = (1/2)*( g(t,x) + g(t+h,x+k+l) )
      //
      // Watch out: we might have to take a step for the modulators
      //
      // xnew = x + h*qavg + sqrth*u*gavg

      vect oldcurrent(current);

      // Some to be optimized out temps
      const number h=dt;
      const number sh=sqrt(dt);

      vect q = (*Integrator<dims,nelem,NT>::deter)(oldcurrent);
      vect g = (*Integrator<dims,nelem,NT>::stoch)(oldcurrent);

      // Generate numbers: noise is _NOT_ correlated between modes.
      vect u,u0;
      for(counter i=0;i<dims;i++) u0[i]=rnd();

      // Do transform for eventual correlation
      trans->function(u,u0);

      // uneff: vect k=h*q;
      // uneff: vect l=sh*g;

      // generate next point
      for(counter i=0;i<dims;i++)
	{
	  oldcurrent[i]+=h*q[i]+sh*g[i]*u[i];
	}

      // uneff:		  l[i]*=u[i];

      // uneff: oldcurrent+=k;
      // uneff: oldcurrent+=l;

      // vect qavg=0.5* h * ( q + (*deter)(oldcurrent) );
      // vect gavg=0.5* sh * ( g + (*stoch)(oldcurrent) );

      vect qn=(*Integrator<dims,nelem,NT>::deter)(oldcurrent);
      vect gn=(*Integrator<dims,nelem,NT>::stoch)(oldcurrent);

      for(counter i=0;i<dims;i++)
	{
	  current[i]+=0.5*( h*( q[i] + qn[i]) + sh*u[i]*(g[i]+gn[i]));

	  // uneff:  gavg[i]*=u[i];
	}

      // uneff: current += qavg;
      // uneff: current += gavg;

      return current;
    }

  private:
    /** Normal random */
    Normal rnd; // will be seeded with time()

    /** storage for transformation */
    vf* trans;
  };

} // end namespace
#endif

/*********************************************************************
$Id: integrator.h,v 1.4 2002-11-21 10:08:51 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2002/11/21 09:47:54  mpeeters
Modified C++ code to compile under gcc 3.2 (implicit typenames removed, slist in different namespace).

Revision 1.2  2002/07/30 19:57:18  mpeeters
Merge of 1.0.2 branch into trunk. From now on, the trunk will be a stable (meaning it compiles and has no quirks) branch and all special stuff will be done on separate branches waiting to be merged into the stable one. A tag STABLE will be made which moves with the trunk. Also, the main trunk will be considered to be 1.0.3.

Revision 1.1.2.2  2002/07/08 14:49:02  mpeeters
Major bugfix: the milshtein algorithm (wrongly) used a stratonovich taylor expansion and thus converged to the stratonovich solution. This has been fixed.

Revision 1.1.2.1  2002/06/19 13:15:47  mpeeters
New integrator manipulators, Euler-Maryuama, transforms.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.4  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.3  2000/09/22 08:50:03  mpeeters
Changed dependencies causing trouble.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
