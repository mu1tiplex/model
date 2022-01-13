/***************************************************************************
			sinmod.h
			-----------
                             
    begin                : 20021119
    author               : (C) 2002 by Michael Peeters
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

#ifndef SINMOD_H
#define SINMOD_H

#include "model/modulator.h"

using namespace MODEL;

/** A modulator returning a sinewave */
class SinMod : public Modulator 
{
public:
  /** Constructor. 
	  \param T the timeframe
	  \param c the vectorfuction to which param belongs
	  \param param the name of the parameter
	  \param period the period of the modulation
	  \param amplitude the amplitude
	  \param bias an optional offset
  */
	template<integer dims, typename nelem, class NT  >
	SinMod( TimeFrame& T,
			VectorFunction<dims,nelem,NT>& v, 
			string param,
			MODEL::time period,
			number amplitude, number bias=0) : 
	  Modulator(T,v,param), 
	  a(amplitude), b(bias), p(period){}

  virtual number modulate(void);
  
private:
  number a;
  number b;
  MODEL::time p;
};

#endif
 
