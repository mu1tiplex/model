/***************************************************************************
                          random.cpp  -  description
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

#include "random.h"

namespace MODEL {

  // Uniform
  //------------------------------------------------------------

	const number Uniform::RNMX=1.0-EPS;


  Uniform::Uniform(counter initial_seed) : Random(initial_seed)
  {
	reshuffle();
  }

  void Uniform::seed(integer seed)
  {
	Random::seed(seed);
	reshuffle();
  }

  void Uniform::reshuffle(void)
  {
	// Load the shuffle table after 8 warmups
	for (counter j=tabelsize+7;j>=0;j--)
	  if (j < tabelsize) tabel[j]=seed_update();
	lastseed=tabel[0];
  }
  const counter& Uniform::seed_update(void)
  {
	counter k=_seed/IQ;
	_seed=IA*(_seed-k*IQ)-IR*k;
	if(_seed<0) _seed+=IM;
	return _seed;
  }

  number& Uniform::generate(number& storage)
  {
	seed_update();
	counter shuffle=lastseed/NDIV;
	lastseed=tabel[shuffle];
	tabel[shuffle]=_seed;

	if ((storage=AM*lastseed)>RNMX) return storage=RNMX; // no end values
	else return storage;
  }

  // Normal
  //------------------------------------------------------------

  Normal::Normal(counter seed)
	: regenerate2(true)
  {
	rnd = new Uniform(seed);
  }

  Normal::~Normal()
  {
	delete rnd;
  }

  number& Normal::generate(number& storage)
  {
	number v1,v2,rsq,fac;
	if(regenerate2) // we need new numbers
	  {
		do { // until we find a pair inside the unit circle
		  v1=2.0*(*rnd)()-1.0;
		  v2=2.0*(*rnd)()-1.0;
		  rsq=v1*v1+v2*v2;
		} while (rsq >=1.0 || rsq==0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		saved=v1*fac; regenerate2=false;
		return storage=v2*fac;
	  }
	else
	  {
		regenerate2=true;
		return storage=saved;
	  }
  }

} // end namespace

/*********************************************************************
$Id: random.cpp,v 1.2 2002-07-30 19:57:19 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.1  2002/07/24 10:01:26  mpeeters
Changed routines to prevent duplication of seed (but not in an efficient way).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
