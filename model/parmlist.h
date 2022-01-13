/***************************************************************************
			parmlist.h
			-----------
                             
    begin                : 27/07/2001
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

#ifndef PARMLIST_H
#define PARMLIST_H

#include "parameters.h"
#include "utility.h"
#include <map>
#include <iostream>
#include <fstream>
#include <stdexcept>


namespace MODEL 
{
  

/** A list of parameters. Can be written to a file and (hopefully) re-read */
template<typename ParT>
class ParmList
{
public:
  typedef Parm<ParT> par;
  typedef par::ParR ParR;
  typedef map<string,par> parlist;
  static const string file_header;

  /** Make an empty list */
  ParmList();

  /** Stream constructor. Allows for persistence */
  ParmList(std::istream& in);

  /** Load the list from file. It will add unknowns, and change known
	  parameters. This allows keeping externals (if the external is
	  known, of course */
  void load(std::istream& in);

  /** Write the list to file */
  void save(std::ostream& out);

  /** Check if parameter is defined */
  bool is_a_parameter(const string& p) 
  {
	return (the_pars.count(p)>0);
  }

  /** Add one to the list */
  void add(const Parm<ParT> p);

  /** Remove from list */
  void remove(const string& name);

  /** Get one from the list */
  par& get_param(const string& name);

  /** Get a naked one... */
  ParR operator[](const string& name)
  {
	return get_param(name).get_parref();
  }

private:
  parlist the_pars;

  NO_COPY(ParmList);
}; 

//
// ----------------------------------------------------------------------
// DECLARATIONS
// ----------------------------------------------------------------------

template<typename ParT>
const string ParmList<ParT>::file_header="# Parameter type: ";

template<typename ParT>
ParmList<ParT>::ParmList() 
{
}

template<typename ParT>
ParmList<ParT>::ParmList(std::istream& in)
{
  // Read header, check parameter type
  // First remove the file_header
  const unsigned int l=file_header.length();
  char dummy[l+1];
  in.get(dummy,l);
  
  // Read the type
  string rname;
  in >> rname;
  // cout << par::id << "=?=" << rname << endl;
  if (rname!=par::id) throw std::logic_error("ParmList::Wrong parameter type in file.");
  
  // The type is correct: read all entries
  string dumtype, newname, equals;
  ParT value;
  
  while (true) {
	in >> dumtype >>  newname >>  equals >> value;
	if (!in) break;
	// cout << "DUH: " << dumtype << newname << equals << value <<endl;
	par newone(newname,value);
	add(newone);
  }  
}

template<typename ParT>
void ParmList<ParT>::load(std::istream& in)
{
  // Read header, check parameter type
  // First remove the file_header
  const unsigned int l=file_header.length();
  char dummy[l+1];
  in.get(dummy,l);
  
  // Read the type
  string rname;
  in >> rname;
  // cout << par::id << "=?=" << rname << endl;
  if (rname!=par::id) throw std::logic_error("ParmList::Wrong parameter type in file.");
  
  // The type is correct: read all entries
  string dumtype, newname, equals;
  ParT value;
  
  while (true) {
	in >> dumtype >>  newname >>  equals >> value;
	if (!in) break;
	// cout << "DUH: " << dumtype << newname << equals << value <<endl;

	// This is the difference with the stream constructor
 	if(is_a_parameter(newname))
 	  get_param(newname)=value;
 	else
	  {
		par newone(newname,value);
		add(newone);
	  }
  }  
}


template<typename ParT>
void ParmList<ParT>::save(std::ostream& out)
{
  // First write ID string to avoid errors when reading
  out << file_header << par::id << endl;

  parlist::iterator run=the_pars.begin();
  parlist::iterator endrun=the_pars.end();
  
  while(run!=endrun){
	out << run->second << endl;
	++run;
  }
}

template<typename ParT>
void ParmList<ParT>::add(const Parm<ParT> p)
{
  the_pars[p.name()]=p;
}

template<typename ParT>
void ParmList<ParT>::remove(const string& name)
{
  throw std::logic_error("ParmList::Not Implemented");
}


template<typename ParT>
ParmList<ParT>::par& ParmList<ParT>::get_param(const string& name)
{
  	  parlist::iterator  end=the_pars.end();
	  parlist::iterator  found=the_pars.find(name);

	  if (found==end) 
		throw std::logic_error("ParmList::Parameter not found");
	  
	  return found->second;
}
 
}
// end namespace MODEL
#endif

/*********************************************************************
$Id: parmlist.h,v 1.1 2001-08-23 21:52:15 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
*********************************************************************/
