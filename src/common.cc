// -*- LSST-C++ -*-

#include <algorithm> //for nth_element
#include <iostream>

#include "lsst/mops/common.h"
#include "lsst/mops/Exceptions.h"



namespace lsst {
namespace mops {

double timeElapsed(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}


    
double fastMedian(std::vector<double> fv) {
    unsigned int middleIndex = 0;
    unsigned int fvSize = fv.size();
    if (fvSize % 2 == 0) {
        middleIndex = (fvSize/2)-1;
    }
    else {
        middleIndex = fvSize/2;
    }
    std::nth_element(fv.begin(), fv.begin() + middleIndex, fv.end());
    return fv[middleIndex];
}





double minOfTwo(double a, double b)
{
    if (a <= b) {
        return a;
    }
    else {
        return b;
    }
}


double maxOfTwo(double a, double b)
{
    if (a > b) {
        return a;
    }
    else {
        return b;
    }
}
    

bool areEqual(double a, double b)
{
    Constants c;
    if (fabs(a - b) < c.epsilon()) {
        return true;
    }
    return false;
    
}


/* returns distance from p1 to p2, ASSUMES that p1 and p2 have
   size==dimensions*/
double euclideanDistance(std::vector<double> p1, std::vector<double> p2, 
                         unsigned int dimensions)
{
    double partialSum = 0.;
    for (unsigned int i = 0; i < dimensions; i++)
    {
        partialSum += pow(p1[i] - p2[i], 2); 
    }
    return sqrt(partialSum);
}   
    
    
    
    
double distance1D(double a, double b, GeometryType geo)
{
    if (geo == EUCLIDEAN) {
        return fabs(a - b);
    }
    else if (geo == CIRCULAR_DEGREES) {
        return circularShortestPathLen_Deg(a, b);
    }
    else if (geo == CIRCULAR_RADIANS) {
        return circularShortestPathLen_Rad(a, b);
    }
    else {
        throw LSST_EXCEPT(BadParameterException, "EE: distance1D: got unexpected geometry type.\n");
    }
}





double convertToStandardDegrees(double deg)
{            
    while (deg < 0) {
        deg += 360.;
    }
    while (deg >= 360) {
        deg -= 360.;
    }
    return deg;
}



/* return true iff regian ab encloses point c - all points
   in degrees */
bool circularIsBetween_Degrees(double a, double b, double c)
{
    a = convertToStandardDegrees(a);
    b = convertToStandardDegrees(b);
    c = convertToStandardDegrees(c);
    
    if (areEqual(a + 180., b)) {
        return true;
    }
    
    if (fabs(a - b) < 180)
    {
        if ((((a <= c) && (c <= b)) || (((a >= c) && (c >= b))))) {
            return true;
        }
        else {
            return false;
        }
    }
    else {
        /* shortest path from a to b crosses 0 point, so we
           really need to search two ranges */
        double lo = minOfTwo(a, b);
        double hi = maxOfTwo(a, b);
        if ((lo >= c) && (c >= 0)) {
            return true;
        }
        else if ((hi <= c) && (c < 360.)) {
            return true;
        }
        else {
            return false;
        }
    }
    
}


void printDoubleVec(std::vector<double> fv)
{
    std::vector<double>::iterator myIter;
    for (myIter = fv.begin();
         myIter != fv.end();
         myIter++)
    {
        std::cout << *myIter << " ";
    }
    std::cout << std::endl;
    
}



double circularShortestPathLen_Deg(double a, double b)
{
    double tmp;
    tmp = fabs(convertToStandardDegrees(a)
               - convertToStandardDegrees(b));
    if (tmp < 180.) {
        return tmp;
    }
    else {
        return fabs(360. - tmp);
    }
}
    

double circularShortestPathLen_Deg_signed(double a, double b)
{
    double tmp;
    tmp = convertToStandardDegrees(a)
        - convertToStandardDegrees(b);
    if ((tmp < 180.) && (tmp > -180.)) {
        return tmp;
    }
    else {
        if (tmp > 180.) {
            return tmp - 360.;
        }
        else {
            // tmp < -180
            return 360. + tmp;
        }
    }
}
    
double circularShortestPathLen_Rad(double a, double b)
{
    Constants c;
    return c.deg_to_rad()*circularShortestPathLen_Deg(c.rad_to_deg()*a, 
                                                      c.rad_to_deg()*b);
}




/* returns the great-circle distance between (RA0, Dec0) and (RA1, Dec1) - 
 * that is, the shortest path along the surface of a unit circle between
 * these two points. All units in degrees. */
double angularDistanceRADec_deg(double RA0, double Dec0, double RA1, double Dec1) 
{
    Constants c;
    double RADist = circularShortestPathLen_Deg(RA0, RA1);
    double DecDist = circularShortestPathLen_Deg(Dec0, Dec1);    
    //convert all factors to radians
    RADist = c.deg_to_rad()*convertToStandardDegrees(RADist);
    DecDist = c.deg_to_rad()*convertToStandardDegrees(DecDist);
    Dec0 = c.deg_to_rad()*convertToStandardDegrees(Dec0);
    Dec1 = c.deg_to_rad()*convertToStandardDegrees(Dec1);
    double r;
    r = 2*asin(sqrt( pow((sin(DecDist/2.)),2) + 
                     cos(Dec0)*cos(Dec1)*pow(sin(RADist/2),2)));
    //back to degrees
    return c.rad_to_deg()*r;
    
}

/* 
 * arcToRA: given a lattitude in declination (degrees) and an arc length
 * in degrees, return the number of degrees in RA which corresponding to
 * that arc length at the given declination.
 
 * E.G, if you are at RA,Dec = (100,0) and travel 1 degree to (101,0), you have
 * moved 1 degree in RA and 1 degree in 'actual' distance Ergo
 * arcToRA(0,1) == 1.  
 
 * however, very close to the north pole at RA,Dec = (0,89), you can walk halfway
 * 'around the earth' to (180,89) and only cross 2 degrees of arc. Ergo
 * arcToRA(89,2) == 180.
 */
double arcToRA(double dec, double arcDeg) 
{
    /* this is the haversine formula, with the dec taken as constant (so
     * delta dec == 0) and put in terms of change in RA. */
    
    Constants c;
    double decRad = c.deg_to_rad()*convertToStandardDegrees(dec);
    double arcRad = c.deg_to_rad()*convertToStandardDegrees(arcDeg);
    double inside = sqrt( pow(sin(arcRad/2),2) / pow(cos(decRad),2) );
    if ((inside >= 1.0) || (inside <= -1.0)) {
        return 180.0;
    }
    double x_Rad = 2*asin(inside);
    return c.rad_to_deg()*x_Rad;
}

        

bool regionsOverlap1D_unsafe(double aLo, double aHi, double bLo, double bHi) 
{
    
    if (aLo < bLo) {
        if (aHi >= bLo) {
            return true;
        }                            
    }
    else { //bmin <= amin
        if (bHi >= aLo) {
            return true;
        }
    }
    return false;
}
    

/*
 * when dealing with circles, this version ASSUMES that aLo < aHi and
 * fabs(aLo - aHi) < 180 - that is, that we are looking at a normal
 * range of angles, not something which crosses 0/360.
 * 
 * it is used by the tree classes, which can make this assumption safely
 * (because they are careful to maintain this property).
 */


bool regionsOverlap1D(double aLo, double aHi, double b1, double b2, 
                      GeometryType type)
{
    if (aLo > aHi) {
        throw LSST_EXCEPT(BadParameterException, "EE: common.cc:regionsOverlap1D:  we DEMAND that aLo < aHi.\n");
    }
    
    if (type == EUCLIDEAN)
    {
        double amin = aLo;
        double amax = aHi;
        double bmin = minOfTwo(b1, b2);
        double bmax = maxOfTwo(b1, b2);
        
        if (amin < bmin) {
            if (amax >= bmin) {
                return true;
            }                            
        }
        else { //bmin <= amin
            if (bmax >= amin) {
                return true;
            }
        }
    }
    
    else if (type == CIRCULAR_DEGREES) {
        if ((aLo > 360.) || (aLo < 0.) || (aHi > 360) || (aHi < 0))
        {
            throw LSST_EXCEPT(BadParameterException, "EE: unexpected condition in common.cc:regionsOverlap1D.  aLo or aHi are outside normal degree bounds (0-360).\n");
        }
        b1 = convertToStandardDegrees(b1);
        b2 = convertToStandardDegrees(b2);
        if (areEqual(fabs(b1 - b2), 180.)) {
            /* b1 b2 has two shortest paths - the two halves of the
             * semicircle */
            return true;
        }
        if (fabs(b1 - b2) < 180) {
            /* nothing funky here - just like euclidean */
            return regionsOverlap1D(aLo, aHi, b1, b2, EUCLIDEAN);
        }
        else {
            /* the shortest path from b1 to b2 crosses the 0 line */
            double bMin = minOfTwo(b1, b2);
            double bMax = maxOfTwo(b1, b2);
            return ((regionsOverlap1D(aLo, aHi, 0, bMin, EUCLIDEAN) ||
                     (regionsOverlap1D(aLo, aHi, bMax, 360, EUCLIDEAN))));
        }
    }
    else if (type == CIRCULAR_RADIANS)
    {
        /*just call recursively after converting to deg...*/
        Constants c;
        double RAD_TO_DEG = c.rad_to_deg();
        return regionsOverlap1D(aLo*RAD_TO_DEG, aHi*RAD_TO_DEG,
                                b1*RAD_TO_DEG, b2*RAD_TO_DEG,
                                CIRCULAR_DEGREES);
        
    }
    else {
        throw LSST_EXCEPT(BadParameterException, 
                          "EE: regionsOverlap1D: got unexpected value for geometry type, must be one of EUCLIDEAN, CIRCULAR_DEGREES or CIRCULAR_RADIANS\n");
    }
    // we should never reach this point.
    return false;
    
}
    
    
    
    
/*
 * this is a safe, general-purpose check for whether two ranges of
 * angles overlap.
 *
 * a region spanning 180 degrees is understood to overlap with ALL other
 * ranges.
 *
 */
    
    
bool angularRegionsOverlapSafe(double a1, double a2, double b1, double b2)
{
    
    double aLo = minOfTwo(a1, a2);
    double aHi = maxOfTwo(a1, a2);
    
    double bLo = minOfTwo(b1, b2);
    double bHi = maxOfTwo(b1, b2);
    
    if ((areEqual(fabs(aLo - aHi), 180.)) ||
        (areEqual(fabs(bLo - bHi), 180.))) {
        // take care of this weird degenerate case.                
        return true;
    }
    
    // make a and b each into ranges which are
    // can be treated as euclidean (do not cross the 0/360 line)
    while ( bHi - bLo > 180.) {
        bLo += 360.;
    }
    while ( aHi - aLo > 180.) {
        aLo += 360.;
    }
    
    double aLo2 = minOfTwo(aHi, aLo);
    double aHi2 = maxOfTwo(aHi, aLo);
    double bLo2 = minOfTwo(bHi, bLo);
    double bHi2 = maxOfTwo(bHi, bLo);
    
    // now get a and b on some contiguous 360-degree region.
    while ( bHi2 - aLo2 > 180) {
        aHi2 += 360;
        aLo2 += 360;
    }
    while (aHi2 - bLo2 > 180) {
        bHi2 += 360;
        bLo2 += 360;
    }
    
    //now we can treat them as euclidean.
    return  regionsOverlap1D_unsafe(aLo2, aHi2, bLo2, bHi2);
}
    





std::string boolToString(bool b) {
    if (b == true) {
        return std::string("true");
    }
    else {
        return std::string("false");
    }
}
    
    
// helper for guessBollFromStringOrGiveErr
std::string stringToLower(std::string strToConvert)
{
    for(unsigned int i=0; i<strToConvert.length(); i++)
    {
        strToConvert[i] = tolower(strToConvert[i]);
    }
    return strToConvert;//return the converted string
}
    
    
bool guessBoolFromStringOrGiveErr(std::string guessStr, std::string errStr) {
    if ((guessStr == "1") ||
        (stringToLower(guessStr) == "true") ||
        (stringToLower(guessStr) == "t")) {
        return true;
    }
    else if ((guessStr == "0") ||
             (stringToLower(guessStr) == "false") ||
             (stringToLower(guessStr) == "f")) {
        return false;
    }
    else {
        throw LSST_EXCEPT(CommandlineParseErrorException, errStr);
    }
    
}
        
        
void extendBounds(std::vector<double> &myBounds, 
                  const std::vector<double> &childBounds,
                  bool areUBounds) 
{
    if (myBounds.size() != childBounds.size()) {
        LSST_EXCEPT(BadParameterException, 
    "Expect argument vectors of same length.");
    }
    for (unsigned int i = 0; i < myBounds.size(); i++) {
        if (areUBounds) {
            if (childBounds.at(i) > myBounds.at(i)) {
                myBounds.at(i) = childBounds.at(i);
            }
        }
        else {
            if (childBounds.at(i) < myBounds.at(i)) {
                myBounds.at(i) = childBounds.at(i);
            }

        }
    }
}


void toCartesian_deg(double ra, double dec, double &x, double &y, double &z)
{
    Constants c;
    x = cos(c.deg_to_rad()*ra) * sin(c.deg_to_rad()*(90 - dec));
    y = sin(c.deg_to_rad()*ra) * sin(c.deg_to_rad()*(90 - dec));
    z = cos(c.deg_to_rad()*(90 - dec));        
}

void toRaDec_deg(double x, double y, double z, double &ra, double &dec)
{
    Constants c;
    ra  = c.rad_to_deg()*(atan2(y, x));
    dec = c.rad_to_deg()*(atan2(z, sqrt(x*x + y*y)));

}

        

}} // close namespace lsst::mops
