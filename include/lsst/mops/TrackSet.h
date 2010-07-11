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
   4/08/10
*/

#ifndef LSST_TRACKSET_H
#define LSST_TRACKSET_H

#include <set>
#include "Track.h"

class TrackSet {
public:
    std::set<Track> componentTracks;
    void insert(const Track &newTrack) {
        componentTracks.insert(newTrack);
    };
    unsigned int size() const {
        return componentTracks.size();
    };
    
    bool isSubsetOf(const TrackSet &other) const {
        if (other.size() < this->size()) {
            return false;
        }
        else {
            // create a copy TrackSet
            std::set<Track>::const_iterator tIter;
            for (tIter = componentTracks.begin();
                 tIter != componentTracks.end();
                 tIter++) {
                if (other.componentTracks.find(*tIter) 
                    == other.componentTracks.end()) {
                    return false;
                }
            }
            return true;
        }
    };

     bool operator==(const TrackSet other) const {
	  return (componentTracks == other.componentTracks);
     };

    bool operator!=(const TrackSet &other) const {
        return ! (*this == other);
    };

};





#endif
