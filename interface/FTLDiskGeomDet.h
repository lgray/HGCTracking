#ifndef RecoTracker_FastTimeMatching_FTLDiskGeomDet_h
#define RecoTracker_FastTimeMatching_FTLDiskGeomDet_h

/// Class corresponding to one layer of FTL endcap

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/ForwardDetId/interface/FastTimeDetId.h"

class FTLDiskGeomDet : public GeomDet {
    public:
        FTLDiskGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float radlen, float xi) ;

        int subdet() const { return subdet_; }
        int zside() const { return zside_; }
        int layer() const { return layer_; }

        bool operator<(const FTLDiskGeomDet &other) const { 
            return (subdet_ == other.subdet_ ? layer_ < other.layer_ : subdet_ < other.subdet_);
        }

    protected:
        const int subdet_, zside_, layer_;
};

#endif
