/***************************************************************************
                          parameter.h  -  description
                             -------------------
    begin                : Mon Aug 11 2000
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

#ifndef PARAMETER_H
#define PARAMETER_H

namespace MODEL{
  /** Parameter stub.
	  OK, OK, this is extremely dirty.
	  We should make a class out of it that encapsulated the pointer copying
	  For the moment however, I am more interested in working code.
	  There should be no change for the "end" user (me)

	  We are probably looking at something like
	  \code
	  class Parameter
	  {
	  public:
	     Parameter(number& changer) : n(&changer);
	     const Parameter& operator=(const number&);
	  private:
	     number* n;
	  };
	  \endcode
  */
  typedef number& Parameter;
  /** Even uglier... for use when they cannot be init.*/
  typedef number* ParameterP;
}
#endif

/*********************************************************************
$Id: parameter.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
