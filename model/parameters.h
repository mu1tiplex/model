/***************************************************************************
			parameters.h
			-----------
                             
    begin                : 24/07/2001
    author               : (C) 2001 by Michael Peeters
    email                : Michael.Peeters@ieee.org
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PARM_H
#define PARM_H

#include "utility.h"
#include <string>
#include <iostream>
#include <stdexcept>

/** A class to hold parameters, Could
	possibly  be used in the future to replace the less than
	satisfactory  MODEL::Parameter class
*/

namespace MODEL 
{

using namespace std;

template <typename ParT>
class Parm 
{
public:
  typedef ParT& ParR;
  /** This defines a string for id */
  static const string id;

  /** Default constructor: no value is guanranteed */
  Parm() : n(id){}

  /** Constructor for parameter with internal storage */
  Parm(const string& name, const ParT& initvalue)
	: n(name), store(initvalue),
	  storage(&store),ext(false),extcopy(false) {}
  
  /** Constructor for parameter with external storage 
	  <copy> is not used, but will be needed to specify if subsequent
	  copies of this object should keep the parameter external. Also,
	  it is needed to distinguish between the two types*/
  Parm(const string& name, ParR p, bool copy) 
	: n(name), storage(&p), ext(true), extcopy(copy) 
  {}

  /** Copy */
  Parm(const Parm& p)
  {
	ext=p.extcopy;
	extcopy=p.extcopy;

	n=p.n;
	store=p.store;
	if(ext)
	  storage=p.storage;
	else
	  storage=&store;
  }

  /** Assignment operator */
  const Parm& operator=(const Parm& p)
  {
	ext=p.extcopy;
	extcopy=p.extcopy;

	n=p.n;
	store=p.store;
	if(ext)
	  storage=p.storage;
	else
	  storage=&store;
	
	return *this;
  }

  /**  Assignment to ParT */
  const Parm& operator=(const ParT& par)
  {
	*storage=par;
  }

  /** Conversion */
  operator ParT() {return *storage;}
  /** Conversion */
  operator ParT() const {return *storage;}

  /** This should not be done, but is anyway for efficiency */
  ParR get_parref() {return *storage;}
  const ParR get_parref() const {return *storage;}

  /** The nice way to set */
  void set(const ParT& p) 
  {
	*storage=p;
  }
  
  /** The nice way to get */
  ParT get(void)
  {
	return *storage;
  }
  
  /** If const, better like this */
  const ParT& get(void) const 
  {
	return *storage;
  }

  /** get name */
  const string& name(void) const
  {
	return n;
  }

  /** Output */
  friend std::ostream& operator<<<>(std::ostream& out, const
									Parm<ParT>& p);

private:
  string n;
  ParT store;
  ParT *storage;
  bool ext;
  bool extcopy;
};


/** Output */
template <typename ParT>
std::ostream& operator<<(std::ostream& out, const Parm<ParT>& p)
{
  out<< p.id << "\t" << p.n << "\t=\t" << *p.storage;
  return out;
}

/** Default ID */
template <typename ParT>
const string Parm<ParT>::id="id_dummy";

/** Specialization for string type */
const string Parm<string>::id="string";

/** Possibility for backward compatibility */
/** Specialization for long double type */
const string Parm<number>::id="number";

typedef Parm<number> Parameter;
typedef Parameter* ParameterP;
 
}
// end namespace MODEL

#endif 

/*********************************************************************
$Id: parameters.h,v 1.1 2001-08-23 21:52:15 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
*********************************************************************/
