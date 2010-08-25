// -*- LSST-C++ -*-
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
