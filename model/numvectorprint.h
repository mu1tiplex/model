/***************************************************************************
                          numvectorprint.h  -  description
                             -------------------
    begin                : Fri May 12 2000
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

#ifndef NUMVECTORPRINT_H
#define NUMVECTORPRINT_H

#include "numerictraits.h"
#include <iostream>

using namespace MODEL;

/** for putting any kind of numvector in an ostream*/

// Previously: template<integer dims, class NT>
template <integer dims, typename nelem,class NT>
std::ostream& operator<<(std::ostream& out,const NumVector<dims,nelem,NT>& vprint)
{
	//out << " | ";
	for (integer i=0;i<(dims-1);i++)	out<< vprint[i] <<"\t";
										out<< vprint[dims-1];
	//out << ")";
	return out;
}

/** For Matrices **/
// Previously:template<integer dims, class NT>
template <integer dims, typename nelem,class NT>
std::ostream& operator<<( std::ostream& out,
		const NumVector<dims,NumVector<dims,nelem,NT>,NT >& vprint)
{
	//out << " | ";
	for (integer i=0;i<(dims-1);i++)	out<< vprint[i] <<std::endl;
										out<< vprint[dims-1];
	//out << ")";
	return out;
}

#endif

/*********************************************************************
$Id: numvectorprint.h,v 1.2 2002-07-30 19:57:18 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.1.2.1  2002/02/06 12:49:39  mpeeters
Altered interface to make fudging easier. Will be reverted.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.3  2000/09/22 08:51:45  mpeeters
Added NumericTraits NT template parameter to all declarations and
definitions. Not really needed, but cleaner this way.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
