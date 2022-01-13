/***************************************************************************
                          cowner.h  -  description
                             -------------------
    begin                : Tue May 9 2000
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

#ifndef COWNER_H
#define COWNER_H

#include <stdexcept>
#include <algorithm>
#include "cycler.h"

#ifdef __GNUC__
// gcc version < 3.1.0 ?
#if __GNUC__ < 3 || \
(__GNUC__ == 3 && __GNUC_MINOR__ < 1)
#include <slist.h>
using std::slist;
#else
#include <ext/slist>
using __gnu_cxx::slist;
#endif // GCC_VERSION
#else
using std::slist;
#endif // __GNUC__

namespace MODEL
{ 
  using namespace std; //for gcc 3.0
  
  /** This class is used to provide a group of Cycler. All the helping
	  hands are defined in here*/
  class COwner 
  {
	friend class Cycler;
  public:
	/** To make sure all the destructors are called */
	virtual ~COwner(){}
	/** Call all Cycler objects */
	void execute_list(void);
	/** What to do for each cycler - this is an extra layer of
		indirection for maximum flexibility. Normally, this will just
		call execute for each cycler */
	virtual void cycle(Cycler& c) {c.execute();}

  protected: 
	void add_to_list(Cycler& c);
	void remove(Cycler& c);

  private:
	slist<Cycler*> c_list;
  };
} // end namespace


#endif

/*********************************************************************
$Id: cowner.h,v 1.3 2002-11-21 10:08:51 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2002/11/21 09:47:54  mpeeters
Modified C++ code to compile under gcc 3.2 (implicit typenames removed, slist in different namespace).

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.3  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
