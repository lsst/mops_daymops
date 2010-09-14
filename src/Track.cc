#include <algorithm>

#include "lsst/mops/Track.h"
#include "lsst/mops/Exceptions.h"

namespace lsst { namespace mops {

Track::Track()
{
    epoch = 0;
}



void Track::addDetection(unsigned int detIndex, const std::vector<MopsDetection> & allDets)
{
    // TDB: might be nice to check that we haven't been given two dets from the same image time.
    // this will enforce some safety with quadratic fitting. not critical for now, 
    // since linkTracklets never does anything that dumb.
    componentDetectionIndices.insert(detIndex);
    componentDetectionDiaIds.insert(allDets.at(detIndex).getID());
}
	  
const std::set<unsigned int> Track::getComponentDetectionIndices() const
{
    const std::set<unsigned int>copy(componentDetectionIndices);
    return copy;
}

const std::set<unsigned int> Track::getComponentDetectionDiaIds() const
{
    const std::set<unsigned int>copy(componentDetectionDiaIds);
    return copy;
}


Track & Track::operator= (const Track &other) {
    componentTrackletIndices = other.componentTrackletIndices;
    componentDetectionIndices = other.componentDetectionIndices;
    componentDetectionDiaIds = other.componentDetectionDiaIds;
    epoch = other.epoch;
    raFunc = other.raFunc;
    decFunc = other.decFunc;
    return *this;
}








void Track::calculateBestFitQuadratic(const std::vector<MopsDetection> &allDets)
{
    std::vector<double> mjds;
    std::vector<double> ras;
    std::vector<double> decs;

    for (std::set<unsigned int>::const_iterator detIndIt = componentDetectionIndices.begin();
         detIndIt != componentDetectionIndices.end(); detIndIt++) {
        const MopsDetection* curDet = &allDets.at(*detIndIt);
        mjds.push_back(curDet->getEpochMJD());
        ras.push_back(curDet->getRA());
        decs.push_back(curDet->getDec());
    }

    epoch = *std::min_element(mjds.begin(), mjds.end());
    for (std::vector<double>::iterator mjdIt = mjds.begin();
         mjdIt != mjds.end(); mjdIt++) {
        *mjdIt = *mjdIt - epoch;
    }
    if (mjds.size() > 0) {
        
        bestFit1d(ras, mjds, raFunc);
        bestFit1d(decs, mjds, decFunc);
    }
}




void Track::bestFit1d(const std::vector<double> &X, 
                      const std::vector<double> &time, 
                      std::vector<double> & res) 
{
    /* jmyers sep 2010 
       
       copy-pasting Kubica's method for calculating quadratic fits,
       replacing dyv_ stuff with C++ style vectors.  use this compare
       performance with GSL.  It's probably quite a bit faster.
       
       see linkTracklets/obs.c line 651 for the original.
    */
    
    
    
    /* The Wise Dr. Kubica said:

       Computes the coefficients for the motion equation
       (WITHOUT too much computation... hopefully):
       
       X = M_0 + M_1 * t + M_2 * 0.5 * t^2
      
       will default to M_i = 0.0 if i > floor(log2(# points)) + 1
    */
    
    // jmyers - removed the 'if ... > NUM_FOR_QUAD. 
    // always use the quadratic fitting - we will never ask for a non-quadratic fit.
    
    double A, B, C, D, E, F, G;
    double a, b, c, t, x;
    double bot, w, sum;
    double dobN, tspread;
    res.resize(3,0);
    unsigned int i;
    int Nv;
    unsigned int N;

    /* Count the number of observations and the number of virtual
       observations (i.e. the number of distinct time steps) */
    N    = X.size();
    // Nv   = dyv_count_num_unique(time, 1e-6);
    // jmyers: we will trust that there are no redundant items per time
    
    dobN = (double)N;
    tspread = std::max_element( time.begin(), time.end() ) -
        std::min_element( time.begin(), time.end() );
     

    // jmyers: we will always do a real quadratic fit.
    // /* Very quickly filter out the 1 and 2 data points case */
    // if(Nv == 1) {
    //     sum = 0.0;
    //     bot = 0.0;
    //     for(i=0;i<N;i++) {
    //         w = 1.0;
    //         //sum += (dyv_ref(X,i) * w);
    //         sum += X.at(i) * w;
    //         bot += w;
    //     }
    //     //dyv_set(res,0,sum/bot);
    //     res.at(0) = sum/bot;
    // } else {
    A = 0.0; B = 0.0; C = 0.0;
    D = 0.0; E = 0.0; F = 0.0;
    G = 0.0;
    
    for(i=0;i<N;i++) {
        w = 1.0;	       
        t  = time.at(i); //jmyers: was dyv_ref(time,i);
        x  = X.at(i); //jmyers: was dyv_ref(X,i);
        A += ((t*t*t*t)*w);
        B += ((t*t*t)*w);
        C += ((t*t)*w);
        D += ((t)*w);
        E += ((x*t*t)*w);
        F += ((x*t)*w);
        G += ((x)*w);
    }
    
    // jmyers - always do a real quadratic. we will
    // never request a linear or any other model.
    
    // /* Default to linear for a few points or */
    // /* A very short arc (< 2 hours).         */
    // if ((Nv < NUM_FOR_QUAD)||(tspread < 0.1)) {
    //     bot = D*D - C*dobN;
    
    //     a = 0.0;    
    //     b = (G*D - dobN*F)/bot;
    //     c = (F*D - C*G)/bot;
    // } else {
    bot  = A*D*D - A*dobN*C + dobN*B*B;
    bot += C*C*C - 2.0*C*B*D;
    
    a  = C*C*G - G*B*D + D*D*E;
    a += B*dobN*F-dobN*E*C - D*F*C;
    a  = 2.0 * (a/bot);
    
    b  = D*A*G - D*E*C - B*G*C;
    b += B*dobN*E + F*C*C -dobN*A*F;
    b  = (b/bot);
    
    c  = E*C*C - E*B*D - A*G*C + A*D*F;
    c += B*B*G - C*B*F;
    c  = (c/bot);

    res.at(0) = c; //jmyers: was dyv_set(res,0,c);
    res.at(1) = b; //jmyers: was dyv_set(res,1,b);
    res.at(2) = a; //jmyers: was dyv_set(res,2,a);   
}


void Track::predictLocationAtTime(const double mjd, double &ra, double &dec) const
{
    double t = mjd - epoch;
    ra = raFunc.at(0) + raFunc.at(1) * t + .5 * raFunc.at(2) * t * t;
    dec = decFunc.at(0) + decFunc.at(1) * t + .5 * decFunc.at(2) * t * t;

}



void Track::getBestFitQuadratic(double &epoch, double &ra0, double &raV, double &raAcc,
                                double &dec0, double &decV, double &decAcc) const
{
    epoch = this->epoch;
    ra0 = raFunc.at(0);
    raV = raFunc.at(1);
    raAcc = raFunc.at(2);
    dec0 = decFunc.at(0);
    decV = decFunc.at(1);
    decAcc = decFunc.at(2);
}




}};
