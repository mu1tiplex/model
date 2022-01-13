/***************************************************************************
                          numerictypes.h  -  description
                             -------------------
    begin                : Tue Apr 18 2000
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

#ifndef NUMTYPES_H
#define NUMTYPES_H

#include <complex>
#include <float.h>

namespace MODEL {

  // simple
  typedef long counter;
  typedef signed int integer;
  typedef long double number;
  typedef number	time;
  typedef std::complex<number> complex;

  const	number	EPS=1.1E-19; // actual value is 1.0842E-19
  const	complex	I=complex(0,1);

}
#endif // NUMTYPES_H

/*********************************************************************
$Id: numerictypes.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
