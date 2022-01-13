/***************************************************************************
                          utility.h  -  description
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

#ifndef UTILITY_H
#define UTILITY_H 

/** \file 
	Various macro's to ease programming
*/

/** Adds a private copy constructor.
	This way, we avoid idiotic default copiers
*/
#define PRIVATE_COPY(thisclass) private: thisclass (const thisclass & copy)

/** Adds a private assignment operator.
	This way, we avoid idiotic default assignments
*/
#define PRIVATE_ASSIGN(thisclass) private: const thisclass& operator=( const thisclass & assign)

/** Makes sure nothing can be copied. The functions are not defined and not public. 
 */
#define NO_COPY(thisclass) PRIVATE_COPY( thisclass ); PRIVATE_ASSIGN( thisclass)

#endif

/*********************************************************************
$Id: utility.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
