/***************************************************************************
			sinmod.cpp
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

#include "sinmod.h"

number 
SinMod::modulate(void)
{
  MODEL::time t=get_time();

  MODEL::time phase=2*M_PI*t/p;
  return b+a*sin(phase);
}
