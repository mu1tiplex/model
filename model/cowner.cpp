/***************************************************************************
			cowner.cpp
			-----------
                             
    begin                : Fri Sep 15 12:22:12 CEST 2000
    author               : (C) 2000 by Michael Peeters
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

#include "cowner.h"

namespace MODEL
{
  void
  COwner::execute_list(void)
  {
	slist<Cycler*>::iterator runner(c_list.begin()),stop(c_list.end());  
	for(;runner!=stop;++runner) cycle(**runner);
  }

  void
  COwner::add_to_list(Cycler& c)
  {
	c_list.push_front(&c);
  }

  void 
  COwner::remove(Cycler& c)
  {
	slist<Cycler*>::iterator deadboy=find(c_list.begin(),c_list.end(),&c);
	if(deadboy!=c_list.end()) c_list.erase(deadboy);
	else throw(std::logic_error("Tried to remove a cycler that wasn't there in ther first place"));
  }
  
} // end namespace

/*********************************************************************
$Id: cowner.cpp,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
