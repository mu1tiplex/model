/***************************************************************************
                          cycler.cpp  -  description
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

#include "cycler.h"
#include "cowner.h"

namespace MODEL
{
  Cycler::Cycler(COwner& master)	
  {
	master.add_to_list(*this);
	theboss=&master;
  }

  Cycler::~Cycler()
  {
	theboss->remove(*this);
  }

  COwner* 
  Cycler::get_boss(void)
  {
	return theboss;
  }

}

/*********************************************************************
$Id: cycler.cpp,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
