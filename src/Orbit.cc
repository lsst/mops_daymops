// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/* File: Orbit.cc
 * Author: Matthew Cleveland
 * Purpose: Implements Orbit object
 */

#include "lsst/mops/Orbit.h"


namespace lsst {
    namespace mops {


/**************************************************
 * Print Orbit data
 **************************************************/
void Orbit::print()
{
  std::cout << "perihelion: " << perihelion << "  eccentricity: "
	    << eccentricity << "  inclination: " << inclination 
	    << "  perihelion arg: " << perihelionArg << "  longitude: " 
	    << longitude << "  perihelion time: " << perihelionTime 
    	    << "  equinox: " << equinox 
	    << "  orbitID: " << orbitID << std::endl;
}




/**************************************************
 * Populate all data fields from formatted string
 **************************************************/
void Orbit::populateOrbitFromString(std::string values, int index)
{
  std::istringstream iss(values);

  //iss.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  
  //try{
  iss >> perihelion;
  iss >> eccentricity;
  iss >> inclination;
  iss >> perihelionArg;
  iss >> longitude;
  iss >> perihelionTime;
  iss >> equinox;
  
  //check for optional orbitID column, if not there, use
  //line number
  if(!iss.eof()){
    iss >> orbitID;
  }
  else{
    orbitID = index;
  }
  
  /*  }
      catch (...){
      throw LSST_EXCEPT(BadParameterException, 
      "Badly-formatted Orbit string\n");
      }*/
}



/**************************************************
 * SETTERS
 **************************************************/
void Orbit::setPerihelion(double _p)
{
  perihelion = _p;
}

void Orbit::setEccentricity(double _e)
{
  eccentricity = _e;
}

void Orbit::setInclination(double _i)
{
  inclination = _i;
}

void Orbit::setPerihelionArg(double _p)
{
  perihelionArg = _p;
}

void Orbit::setLongitude(double _l)
{
  longitude = _l;
}

void Orbit::setPerihelionTime(double _p)
{
  perihelionTime = _p;
}

void Orbit::setEquinox(double _e)
{
  equinox = _e;
}

void Orbit::setOrbitID(double _o)
{
  orbitID = _o;
}



/**************************************************
 * GETTERS
 **************************************************/
double Orbit::getPerihelion() const
{
  return perihelion;
}

double Orbit::getEccentricity() const
{
  return eccentricity;

}

double Orbit::getInclination() const
{
  return inclination;
}

double Orbit::getPerihelionArg() const
{
  return perihelionArg;
}

double Orbit::getLongitude() const
{
  return longitude;
}

double Orbit::getPerihelionTime() const
{
  return perihelionTime;
}

double Orbit::getEquinox() const
{
  return equinox;
}

double Orbit::getOrbitID() const
{
  return orbitID;
}





    }} //close lsst::mops
