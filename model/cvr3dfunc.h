/***************************************************************************
                          cvr3dfunc.h  -  description
                             -------------------
    begin                : Thu Apr 20 2000
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

#ifndef CVR3DFUNC_H
#define CVR3DFUNC_H

#include "numerictypes.h"
#include "numvector.h"
#include "vectorfunction.h"

namespace MODEL
{


  /**Constant Function Class: returns a (constant real 3D vector),
	 whatever you do. Only for testing purposes.
	 */

  class CVR3DFunc : public VectorFunction<3>
  {
  public:
	typedef	VectorFunction<3>	base;
	
	CVR3DFunc(const vect& r=vect(0.) ) : cr(r) {}
	virtual ~CVR3DFunc() {};
		
	/** Copy constructor */
	CVR3DFunc(const CVR3DFunc& c) : base(c) {copy(c);}
		
	/** Assignment */
	const CVR3DFunc& operator=(const CVR3DFunc& c)
	{	if(this!=&c) {
	  base::operator=(c); copy(c);} return *this;}
		
	/** Virtual Copy Constructor */
	virtual CVR3DFunc* clone () const
	{ return new CVR3DFunc(*this); }

  private:
	void	copy(const CVR3DFunc& c) {cr=c.cr;}

	virtual const vect& function(vect& fu, const vect& u)	
	{
	  fu=cr;
	  return fu;
	}	

  private:
	vect	cr;
  };

} // end namespace
#endif

/*********************************************************************
$Id: cvr3dfunc.h,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
