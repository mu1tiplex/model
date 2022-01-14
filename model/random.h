/***************************************************************************
                          random.h  -  description
                             -------------------
    begin                : Mon Jul 24 2000
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

#ifndef RANDOM_H
#define RANDOM_H

#include "numerictypes.h"
#include <iostream>
#include <stdexcept>
#include <time.h>
#include <unistd.h>

namespace MODEL
{
  /**Defines random generators which return numbers between 0 and 1
   */

  class Random {
  public:
	Random(counter initial_seed=0) {seed(initial_seed);}

	/** called to generate number */
	number operator()(void) {number temp; return generate(temp);}

	/** overload this one to define random number */
	virtual number& generate(number& storage) = 0;

	/** re-seed the generator */
	virtual void seed(integer seed)
	{
	  if(seed!=0)
		_seed=seed;
	  else
		{sleep(2);
		_seed=::time(NULL);
		/** \bug: A major problem here: when you call this more than
			once in the same second, a big problem occurs: the random
			generators have the same seed. Because our requirements are
			none to stringent, we wait a least two seconds before each
			call. This slows down the constructor.
			\todo: Use the kernel entropy pool to add a single sumber to
			the seed.
		*/
		}

	}
	/** grab the current seed status - to be able to continue if needed */
	counter get_seed(void) {return _seed;}

  protected:
	counter _seed;
  };

  /** Minimal uniform Random generator */
  class Uniform : public Random
  {
  public:
	Uniform(counter inital_seed=0);
	virtual number& generate(number& storage);
	const counter& seed_update(void);
	void seed(integer seed);

  private:
	// void Uniform::reshuffle(void);
	void reshuffle(void);

	static const counter IA=16807;
	static const counter IM=2147483647;
	static const number AM=1.0/IM;
	static const counter IQ=127773;
	static const counter IR=2836;
	static const counter tabelsize=32;
	static const counter NDIV = 1+(IM-1)/tabelsize;
	static const number RNMX;

	counter tabel[tabelsize];

	counter lastseed;
  };



  /** Gaussian (normal) distribution with stddev of 1.
	  This is not derived off Uniform, as it does not return numbers
	  between 0 and 1*/
  class Normal {
  public:
	Normal(counter seed=0);
	virtual ~Normal();
	virtual number& generate(number& storage);
	number operator()(void) {number temp; return generate(temp);}

	void seed(integer seed)
	{
	  rnd->seed(seed);
	}

  private:
	Uniform* rnd;
	bool regenerate2; //numbers are generated in pairs
	number saved;//the second one
  };

} // end namespace
#endif

/*********************************************************************
$Id: random.h,v 1.2 2002-07-30 19:57:19 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.2  2002/07/24 10:01:26  mpeeters
Changed routines to prevent duplication of seed (but not in an efficient way).

Revision 1.1.2.1  2002/07/05 09:04:15  mpeeters
Fixed problem with multiple (identical) random generators.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.3  2000/09/29 13:00:49  mpeeters
Added ChangeLog (using cvs2cl.pl script by Karl Fogel)
Changed makefile to take into account tutorial
Muddled around in random.h
Added tutorial section to the Introduction

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
