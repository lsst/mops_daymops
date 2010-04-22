// -*- LSST-C++ -*-


/* jmyers 7/25/08
 * 
 * a class for mapping spatial points (presented as an arbitrary-length standard
 * vector) to 'values' of an arbitrary type.  This is to be used particularly in
 * spatial searching, where we would like to associate a point in space with
 * some value (e.g. RA 108.0, Dec -10.1 is the location of detection 1234).
 * 
 * TBD: It might really be wise to replace this entirely with a 
 *
 *  template <class T,K>
 *  pair<std::vector<T>, K> 
 * 
 * or possibly even be more flexible and only demand that the first element be
 * iterable and/or indexable.  I shudder to think how messy the KDTree code
 * will look then...!  We'll keep this for now.
 */

#ifndef LSST_POINT_AND_VALUE_H
#define LSST_POINT_AND_VALUE_H


#include <iostream>
#include <vector>

namespace lsst {
namespace mops {


    template <class T>
    class PointAndValue {

    public:
        void setPoint(std::vector <double> point) { myPoint = point; }
        void setValue(T value) { myValue = value; }
    
        std::vector <double> getPoint() const { return myPoint; }
        T getValue() const { return myValue; }

        void debugPrint() {
            std::cout << " my Point: [";            
            for (unsigned int i = 0; i < myPoint.size(); i++) {
                std::cout << myPoint.at(i);
            }
            std::cout << "\n my Value: [" << myValue << std::endl;
        }
    
    private:
        T myValue;
        std::vector <double> myPoint;
    };

}} // close namespace lsst::mops


#endif
