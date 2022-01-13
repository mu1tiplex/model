/***************************************************************************
                          cycler.h  -  description
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

#ifndef CYCLER_H
#define CYCLER_H

namespace MODEL
{
  class COwner;

  /** This takes care of owner-listofelements type of
	  constructs. Cycler is the element, COwner is the master
  */
  class Cycler
  {
	friend class COwner;
	
  public:
	/** Create a Cycler, and tell it who it's COwner is*/
	Cycler(COwner& master);
	virtual ~Cycler();
	
	/** To be written for each one */
	virtual void execute(void)=0;
  protected:
	/** In case we need to remove ourself */
	COwner* get_boss(void);
	
  private:
	/** Not allowed */
	Cycler(const Cycler& c);

	COwner* theboss;
  };
  
  
}// end namespace


#endif

/*********************************************************************
$Id: cycler.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
