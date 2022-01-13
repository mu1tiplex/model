/***************************************************************************
			probe.cpp
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

#include "probe.h"

namespace MODEL
{
  using namespace std;

  // GenericProbe

  GenericProbe::GenericProbe(TimeFrame & T, const string& fn) 
	: TickTock(T,T.get_dt()), n(fn), out(NULL)
  {
  }

  GenericProbe::GenericProbe(TimeFrame & T, ostream& os) 
	: TickTock(T,T.get_dt()), out(&os)
  {
  }

  void 
  GenericProbe::tick()
  {
	probe();
  }

  void
  GenericProbe::write_data(void)
  {
	bool no_out=false;  
	if(no_out=(out==NULL)) out=new ofstream(n.c_str());

	// This does not work - doh !
	print(*out);

	if(no_out) {delete out; out=NULL;} // Make sure it is not used again
  }
  

  //Probe

  Probe::Probe(TimeFrame & T, const string& fn) 
	: TickTock(T,T.get_dt()), n(fn), out(NULL)
  {
	var_data=new data_list;
  }

  Probe::Probe(TimeFrame & T, ostream& os) 
	: TickTock(T,T.get_dt()), out(&os), data_is_here(false)
  {
	var_data=new data_list;
  }

  Probe::~Probe()
  { 
	bool no_out=false;  
	if(no_out=(out==NULL)) out=new ofstream(n.c_str());
	  
	(*out) << "# Probe Output" <<endl;
	  
	for( data_list::iterator printer=var_data->begin();
		 printer!=var_data->end();
		 printer++)
	  {
		if (printer->size()>1)
		  {
			for ( data_record::iterator clack=printer->begin();
				  clack!=(printer->end());
				  clack++)
			  (*out) << (*clack) << "\t";
			(*out) << endl;
		  }
		
	  }

	if(no_out) delete out;	
	delete var_data;
  }

  void 
  Probe::tick()
  {
	data_is_here=false;
	temp_data = new data_record;
	temp_data->push_back(get_time());

	probe();

	// the last element is the time
	if(data_is_here) 
	  {var_data->push_back(*temp_data); }
	delete temp_data; temp_data=NULL;

  }

  void 
  Probe::add_data(const number& d)
  {
	// As it is always push_back which is used, we just have to
	// look at the latest addition
	data_is_here=true;
	temp_data->push_back(d);
  }

} // end namespace

/*********************************************************************
$Id: probe.cpp,v 1.1 2001-05-22 10:54:55 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.4  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.3  2000/10/09 09:43:44  mpeeters
Added GenericProbe, because Probe was too specialized already.
Probe now only logs data when add_data is called, not at every timepoint.
This caused a huge memory leak, which has been shut now.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
