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
 
/* jonathan myers 
   8/10/09
*/

#ifndef LSST_TRACK_H
#define LSST_TRACK_H

#include <set>



class Track {
public:
    // the tracklets which were used to build this track, if any
    std::set<unsigned int> componentTrackletIndices;
    /* the set of detections associated with this track.  note that this is not
       merely the union of the detections associated with the tracklets
       associated with this track; if the track contains two tracklets which
       are in 'conflict' (contain multiple detections from the same image time)
       then only one of the possible detections is chosen.        
    */
    std::set<unsigned int> componentDetectionIndices;
    

    Track & operator= (const Track &other) {
        componentTrackletIndices = other.componentTrackletIndices;
        componentDetectionIndices = other.componentDetectionIndices;
        return *this;
    }

    bool operator==(const Track &other) const {
        bool toRet = componentDetectionIndices == other.componentDetectionIndices;
        //std::cout << " Track operator == called. returning " << (toRet ? "TRUE\n" : "FALSE\n");

        return toRet ;
    }

    bool operator!=(const Track &other) const {
        return ! (*this == other);
    }

    bool operator<(const Track &other) const {
        //std::cout << "Track operator< called." << std::endl;
        bool toRet = componentDetectionIndices < other.componentDetectionIndices;
        //std::cout << "Track operator< called. returning " << (toRet ? "TRUE" : "FALSE")  << std::endl;
        return toRet;
        
    }

};




#endif
