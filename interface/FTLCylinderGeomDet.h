#ifndef RecoTracker_FastTimeMatching_FTLBarrelGeomDet_h
#define RecoTracker_FastTimeMatching_FTLBarrelGeomDet_h

/// Class corresponding to one layer of FTL barrel

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/ForwardDetId/interface/FastTimeDetId.h"

class FTLBarrelGeomDet : public GeomDet {
    public:
        FTLBarrelGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float zmin, float zmax, float radlen, float xi) ;

        int subdet() const { return subdet_; }
        int zside() const { return zside_; }
        int layer() const { return layer_; }

        bool operator<(const FTLBarrelGeomDet &other) const { 
            return (subdet_ == other.subdet_ ? layer_ < other.layer_ : subdet_ < other.subdet_);
        }

    protected:
        const int subdet_, zside_, layer_;
};

#endif
