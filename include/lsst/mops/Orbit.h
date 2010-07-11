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
 
/* File: Orbit.h
 * Author: Matthew Cleveland
 * Purpose: Represents and Orbit object
 */
#ifndef _ORBIT_H_
#define _ORBIT_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

//#include <stdlib.h>

//#include "Exceptions.h"
//#include <stdio.h>


namespace lsst{
namespace mops{

class Orbit{

 public:
  enum{PERIHELION, ECCENTRICITY, INCLINATION, PERIARG, LONGITUDE, PERITIME};

  void populateOrbitFromString(std::string, int);
  void print();


  /**************************************************
   * SETTERS
   **************************************************/
  void setPerihelion(double);
  void setEccentricity(double);
  void setInclination(double);
  void setPerihelionArg(double);
  void setLongitude(double);
  void setPerihelionTime(double);
  void setEquinox(double);
  void setOrbitID(double);

  /**************************************************
   * GETTERS
   **************************************************/
  double getPerihelion() const;
  double getEccentricity() const;
  double getInclination() const;
  double getPerihelionArg() const;
  double getLongitude() const;
  double getPerihelionTime() const;
  double getEquinox() const;
  double getOrbitID() const;

 private:
  
  double perihelion;
  double eccentricity;
  double inclination;
  double perihelionArg;
  double longitude;
  double perihelionTime;
  double equinox;
  double orbitID;
};

}} // close lsst::mops namespace

#endif
