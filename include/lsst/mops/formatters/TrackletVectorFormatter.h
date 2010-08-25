// -*- LSST-C++ -*-
/* jonathan myers 
   5/2/10
*/

#ifndef _MOPS_TRACKLET_VECTOR_FORMATTER_H
#define _MOPS_TRACKLET_VECTOR_FORMATTER_H

#include "lsst/pex/policy/Policy.h"

#include "lsst/daf/base.h"
#include "lsst/daf/persistence.h"

#include "lsst/daf/base/Persistable.h"


namespace lsst {
namespace mops {
namespace formatters {



class TrackletVectorFormatter : public lsst::daf::persistence::Formatter {
public:
    
    explicit TrackletVectorFormatter(lsst::pex::policy::Policy::Ptr const & p);

    virtual ~TrackletVectorFormatter();

    virtual void write(
        lsst::daf::base::Persistable const *,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );

    virtual lsst::daf::base::Persistable* read(
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );

    virtual void update(
        lsst::daf::base::Persistable*,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );
private:
    lsst::pex::policy::Policy::Ptr _policy;
};


}}} // close lsst::mops::formatters


#endif
