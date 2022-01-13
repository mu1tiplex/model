/***************************************************************************
			rootscan.h
			-----------
                             
    begin                : Fri Sep 15 12:22:11 CEST 2000
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

#ifndef ROOTSCAN_H
#define ROOTSCAN_H

#include "numvector.h"
#include "numvectorprint.h"
#include "numerictraits.h"
#include "vectorfunction.h"
#include <string>
#include "parameter.h"
#include <list>
#include <iostream>
#include <map>
#include "vfwithbump.h"
#include "eigenvalues.h"
#include "newtonroot.h"
#include "jacobian.h"

namespace MODEL 
{
  /** A class to store a parameter scan. Very ad hoc fro the moment
	  and not flexible at all 
  */
  template<integer dims, typename nelem=number, class NT =
  NumericTraits<nelem,dims> >
  class ScanList
  {
  public:
	typedef typename NT::vect vect;

	/** a parpoint is a list of vectors, any number of them.
		In practice, the first is the solutions, the second the
		real part of the eigenvalues and the third the imaginary part
	*/

	typedef vector< NumVector<dims,nelem,NT> > parpoint;
	
	/** a parpointlist is for example all the parpoints which belong
		to one specific parameter value. Hence it will be associated
		thru a map to a value */

	typedef list<parpoint> parpointlist;
	typedef map<number,parpointlist> datalist;

	/** a function which returns true when the point is "good" */
	typedef bool criteria(const parpoint&);

	/** Clear out everything
		\todo check if this doe not leak
	*/
	void clear() {dat.clear();}

	/** Find the closest solutions */
	parpointlist get_solution(number param)
	{
	  typename datalist::iterator here=dat.upper_bound(param);
	  return here->second;
	}

	/** Add a point to the newest parpoint*/
	void add_point(const vect& v)
	{
	  p.push_back(v);
	}

	/** Finalize: add the parpoint to the parpointlist for a certain
		parameter, and clear it*/
	void add_param(const number& param)
	{
	  dat[param].push_back(p);
	  p.clear();
	}
	
	/** select a part of the list, the rest is discarded. Returns the
		new list. Note: no use is made of remove_if from the STL,
		which would be the better way to go. Time... */
	ScanList select(criteria ifthis) const
	{
	  // make a copy
	  ScanList selected(*this);

	  // Run thru it, removing all bad elements
	  for(typename datalist::iterator p=selected.dat.begin();
		  p!=selected.dat.end();++p)
		{
		  typename parpointlist::iterator pp=p->second.begin();
		  while(pp!=p->second.end())
			if(!ifthis(*pp)) pp=p->second.erase(pp);
			else ++pp;
		}
	  
	  return selected;
	}

	/** write out a sorted list (due to the map), but raw for the
		rest. If you want continuous linked regions, look at
		print_list */
	void print_raw(ostream& out)
	{
	  out << "# Raw data -----------------------" << endl;
	  
	  for(typename datalist::iterator parr=dat.begin();
		  parr!=dat.end();
		  ++parr)
		{

		  for(typename parpointlist::iterator pointlist=parr->second.begin();
			  pointlist!=parr->second.end();
			  ++pointlist)
			{
			  out << parr->first; // Write parameter value
			  for(typename parpoint::iterator point=pointlist->begin();
				  point!=pointlist->end();
				  ++point)
				out << "\t" << *point;
			  out << endl;
			}
		}
	}

	/** Add two lists together. Good to generate plots for a certain
		parameter, varying another */
	ScanList operator+=(const ScanList<dims,nelem>& extra)
	{
	  const datalist& other=extra.dat;
	  
	  for(typename datalist::const_iterator parr=other.begin();
		  parr!=other.end();
		  ++parr)
		{

		  for(typename parpointlist::const_iterator pointlist=parr->second.begin();
			  pointlist!=parr->second.end();
			  ++pointlist)
			{
			  // out << parr->first; // Write parameter value
			  for(typename parpoint::const_iterator point=pointlist->begin();
				  point!=pointlist->end();
				  ++point)
				add_point(*point); // out << "\t" << *point;
			  add_param(parr->first); // out << endl;
			}
		}
	  return *this;
	}
	

	/** writes out a sorted list (gnuplot style - not clean, I know,
		and much too long to be included, but it works 
		@param accuracy defines the strictness of the sorter. It tells
		you by a factor of how much the roots and eigenvalues have to
		match
		\note FIXED: printlist eats all the data, and nothing is left 
		afterwards. So: we have to make a copy first 
		\bug We should look for the closest match !!!
	*/ 
	void print_list(ostream& out, const number& accuracy=1.)
	{
	  datalist toprint(dat);
	  
	  out << "# Scanlist output --------------------" << endl;
	  
	  counter pieces(0);
  
	  // 1) Find a starting parameter value
  
	  typename datalist::iterator parrunner(toprint.begin());
	  while(parrunner!=toprint.end()) { // while there are parameters left

		// 2) Find a starting point
		typename parpointlist::iterator varrunner(parrunner->second.begin());
		
		while (varrunner!=parrunner->second.end()) { // if there are any
		  typename datalist::iterator nextrunner(parrunner);
		  number param=nextrunner->first; // get first = KEY 

		  // write comment
		  out << "# -------------------------------" 
			  << " Start of list starting at " << param <<endl;

		  // We have a starting point, remove it from the list
		  parpoint previous(*varrunner);
		  number previousparam(param);

		  typename parpointlist::iterator diepiggy(varrunner);
		  
		  parrunner->second.erase(diepiggy); // first erase,
		  varrunner=parrunner->second.begin(); // then find another
	  
		  // Add it to our list:		  
		  out << param;
		  for(typename parpoint::iterator i=previous.begin();
			  i!=previous.end();++i)
			out << "\t" << *i;
		  out << endl;
	  
		  // Find a next point (a higher par value)
		  ++nextrunner;
		  while (nextrunner!=toprint.end()){  // loopy - we have points
			// find someone to love
			bool loved=false;
			number nextparam=nextrunner->first;

			typename parpointlist::iterator inrunner=nextrunner->second.begin();

			// calculate match requirements
			number maxchange=accuracy*nextparam/previousparam;

			while (inrunner!=nextrunner->second.end()) { // find match

			  // Ok, we have one
			  parpoint match(previous);
			  parpoint gotcha(*inrunner);
			  bool matched=true;
			  
			  // Is it any good: calculate all differences
			  for(typename parpoint::iterator 
					i=match.begin(), 
					j=gotcha.begin(),
					k=previous.begin();
				  i!=match.end();
				  ++i,++j,++k)
				{
				  *i -= *j;
				  *i /=length(*k); // normalized difference

				  // the abs is to make sure it works for complex
				  // types
				  /** \todo rewrite the length function for complex
					  types, so it returns a real number */
				  matched= matched && abs(length(*i)) <maxchange;
				}				
			 
			  /** \todo  actually, we should search the best match, but for
				  now, this should be sufficient */
			  if(matched) { // we have a match
				previous=gotcha;
			
				// Remove it from the list
				parrunner->second.erase(inrunner);
			
				// Add it to our list:
				out << nextparam;
				for(typename parpoint::iterator i=previous.begin();
					i!=previous.end();++i)
				  out << "\t" << *i;
				out << endl;

				loved=true;
				break; // We are happy, go to next parval		
			  }
			  else ++inrunner;
		  
			} // while find match

			if(loved) {
			  // try it with the next parval
			  ++nextrunner;
			}
			else {
			  // All exhausted, close our list
			  // cout << "#----------------------------" << endl;
			  ++nextrunner;
		  
			  // break;		  
			}

			previousparam=nextparam;
		
		  } // loopy
		  // close the list
		  out << "# -------------------------------" 
			  << " End of list starting at " << param <<endl<<endl<<endl;
		  ++pieces;
	  
		  // three endl's for gnuplot
	  
		} // if we had any
	
	
		// Go to next parameter value, as we are exhausted here
		++parrunner;
	
	  } // while parameters left  
	  out << "# Number of distict pieces: " << pieces << endl;
	}
	  
	datalist& get_data(void)
	{
	  return dat;
	}

  private:
	datalist dat;
	parpoint p;
  };
 
  //------------------------------------------------------------
 
  /** A class which will scan a vectorfunction over a parameter range
	  and return a ScanList containing the roots for each value of the
	  parameter, together with the eigenvalues

	  \todo For the rootfinder, it would be classy if we would be able
	  to specify the order of magnitude we expect the different vars
	  to have, and use this to better estimate starting points
  */
  template<integer dims, typename nelem=number, class NT = NumericTraits<nelem,dims> >
  class RootScan
  {
  public:
	typedef typename NT::vf vf;
	typedef typename NT::vect vect;

	/** Initialize */
	RootScan(vf& tobescanned, const string& param, vect starter=1.) : f(&tobescanned)
	{
	  Parameter p=f->get_parameter(param);
	  _p=&p;

	  // Generate some plausible starting points
	  def_starters.push_back(starter);		  
	}
  

  /** Added as a convenience: we want to find only the roots of the
	  system at the point it is right now: the parameter has been
	  set by other means */
  RootScan(vf& tobescanned) : f(&tobescanned)
  {
	_p=&dummy;
  }

  /** Add other possible starting points */
  void add_start(const vect start)
  {
	def_starters.push_back(start);
  }

  /** Scan */
  const ScanList<dims,nelem,NT>& scan(const number& from, 
									  const number& to,
									  const counter& n)
  {
	// Clear out previous scan
	roots.clear();

	number delta=(to-from)/number(n-1);

	// A starting value list
	list< NumVector<dims,nelem,NT> > starters(def_starters);

	for(*_p=from;*_p<=to;*_p+=delta)		{
#ifdef ROOTSCAN_DEBUG
	  cerr << "PARAMETER SCAN: " << *_p << endl;
#endif
	  // Add a bumpy layer, to find other zeroes
	  VFwithBump<dims,nelem,NT>	bumpy(*f);
	  // Find stationary solution (the zeroes)
	  NewtonRoot<dims,nelem,NT> stat(bumpy);
	  Jacobian<dims> J(*f,1E-6);

	  NumVector<dims,nelem,NT> solution;

	  // Get the list of starting points, defined from the
	  // previous parameter value. 
	  list< NumVector<dims,nelem,NT> > oldstarters(starters);

	  starters.clear();
	  starters=def_starters;

	  // Loop until no more points
		
		typename list< NumVector<dims,nelem,NT> >::iterator
		  start=oldstarters.begin();
		while(start!=oldstarters.end()){		  
		  while(true){		
			try {
			  solution=stat(*start); 
			}
			catch (logic_error le) 
			  {
#ifdef ROOTSCAN_DEBUG
				cerr<<le.what()<<endl;
#endif
				++start; 
				break;}
				
		  //debug : 
		  // cerr << *_p << ":\t" << solution << endl;
			  
		  if(stat.wrong_min()) {++start; break;}
		  if(stat.no_root()) {++start; break;}

		  // Oh goody, a root !
			  
		  // Find stability
		  NumVector<dims,nelem,NT> reals=Eigenvalues<dims,nelem,NT>(J.calculate(solution)).real();
		  NumVector<dims,nelem,NT> imags=Eigenvalues<dims,nelem,NT>(J.calculate(solution)).imag();
			
		  // Define the data
#ifdef ROOTSCAN_DEBUG
		  cerr << "STATPOINT:" << solution << endl;
#endif
		  roots.add_point(solution);
		  roots.add_point(reals);
		  roots.add_point(imags);
		  // Save the data
		  roots.add_param(*_p);

		  // add the solution to the bumplist
		  bumpy.AddBump(solution);

		  // add to next try for next point
		  starters.push_back( solution ); 

		  /** \todo this is a problem: add a boolean parameter to
			  select this. With a warning about extra roots.
		  */

		  // 			  // define new starting point
		  // 			  for(integer i=0;i<dims;i++)
		  //  				{
		  //  				  NumVector<dims,nelem,NT> nsplus(solution);
		  //  				  NumVector<dims,nelem,NT> nsmin(solution);
		  // 				  nsplus[i]+= (delta / *_p)*solution[i]; // Linear
		  // 				  // Interpol
		  // 				  nsmin[i]-= (delta / *_p)*solution[i]; // Linear Interpol 
		  // 				  starters.push_back(nsplus);
		  // 				  starters.push_back(nsmin);

		  // 				}
				  
		}
			
	  }
	}
	return roots;
  }
	
  // A few examples of criteria
  /** select all */
  static bool all(const typename ScanList<dims,nelem,NT>::parpoint& p) {return true;}
  /** select stable ones - this assumes the real parts of the
	  eigenvalues are the second element of the parpoint and that we
	  are working in a 3d system.*/
  static bool stable(const typename ScanList<dims,nelem,NT>::parpoint& p) 
  {
	bool stab=true;
	for(counter i=0;i<dims;++i) stab=stab&&(p[1][i]<0);
	  
	return stab;
  }
  /** select unstable ones - see note for stable */
  static bool unstable(const typename ScanList<dims,nelem,NT>::parpoint& p)
  {
	return !stable(p);
  }
  /** select positive ones - this assumes that the function is the
	  first element of the parpoint*/
  static bool positive(const  typename ScanList<dims,nelem,NT>::parpoint &p)
  {
	bool pos=true;
	for(counter i=0;i<dims;++i) pos=pos&&(p[0][i]>-1E-4);
	// Agree, this is still negative, but still... Just to remove
	// rounding errors which play hell on the routines
	  
	return pos;
  }

  static bool negative(const  typename ScanList<dims,nelem,NT>::parpoint &p)
  {
	return !positive(p);
  }

	
private:
  vf* f;
  number* _p;
  ScanList<dims,nelem,NT> roots;

  // if we don't have a parameter to scan
  number dummy;

public:
  /** A list of starting points to use */
  list< NumVector<dims,nelem,NT> > def_starters;
};
  




} // end namespace

#endif

/*********************************************************************
$Id: rootscan.h,v 1.6 2002-11-21 09:47:54 mpeeters Exp $
**********************************************************************

$Log: not supported by cvs2svn $
Revision 1.5  2002/11/18 14:10:19  mpeeters
Removed unnecessary ssa.cpp. Added DEBUG in rootscan catch{} routine. Added RC and phase info to the SSA routines.

Revision 1.4  2002/07/30 19:57:19  mpeeters
Merge of 1.0.2 branch into trunk. From now on, the trunk will be a stable (meaning it compiles and has no quirks) branch and all special stuff will be done on separate branches waiting to be merged into the stable one. A tag STABLE will be made which moves with the trunk. Also, the main trunk will be considered to be 1.0.3.

Revision 1.3.2.2  2002/06/19 13:19:07  mpeeters
Made the select member const, so one can make a subselection of a const list.

Revision 1.3.2.1  2001/08/30 14:46:04  mpeeters
Added #ifdef ROOTSCAN_DEBUG to surpress outputting ridiculous amounts of info to cerr.

Revision 1.3  2001/08/23 20:24:00  mpeeters
Modified docs to provide correct version info.

Revision 1.2  2001/07/26 12:01:22  mpeeters
Make start point selection easier (add_start). Removed some bugs.

Revision 1.1  2001/05/22 10:54:55  mpeeters
Moved sources and headers for libModel to model/ subdirectory, in an attempt to rationalize the source tree. This should make things "netter".

Revision 1.5  2001/05/21 11:53:16  mpeeters
Removed Makefile.in, which is automatically generated anyway.

Revision 1.4  2000/09/29 08:56:28  mpeeters
Inserted into print_raw: write the parameter on EVERY line, otherwise gnuplot flips.

Revision 1.3  2000/09/22 08:50:03  mpeeters
Changed dependencies causing trouble.

Revision 1.2  2000/09/15 10:26:31  mpeeters
Cleaned out header and added CVS tails to files, removed superfuous
@author comments, inserted dates

*********************************************************************/
