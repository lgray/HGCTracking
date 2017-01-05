#ifndef RecoTracker_FastTimeMatching_FTLTracker_h
#define RecoTracker_FastTimeMatching_FTLTracker_h

/// Class holding all the layers of the FTL Tracker
/// With some more effort it could inherit from TrackingGeometry etc
// - holds (and owns) all the disks
// - implements navigation and lookup by id

#include "RecoTracker/FastTimeMatching/interface/FTLDiskGeomDet.h"
#include "Geometry/HGCalGeometry/interface/FastTimeGeometry.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

class FTLTracker {
    public:
        FTLTracker(const FastTimeGeometry *geom) ;
        ~FTLTracker() { 
            for (FTLDiskGeomDet * disk : disksPos_) delete disk;
            for (FTLDiskGeomDet * disk : disksNeg_) delete disk;
        }

        // disable copy constructor and assignment
        FTLTracker(const FTLTracker &other) = delete;
        FTLTracker & operator=(const FTLTracker &other) = delete;

        unsigned int numberOfDisks(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).size(); }
        const FTLDiskGeomDet * disk(int zside, int disk) const { return (zside > 0 ? disksPos_ : disksNeg_).at(disk); }
        const FTLDiskGeomDet * firstDisk(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).front(); }
        const FTLDiskGeomDet * lastDisk(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).back(); }
        const FTLDiskGeomDet * firstDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? firstDisk(zside) : lastDisk(zside); }
        const FTLDiskGeomDet * lastDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? lastDisk(zside) : firstDisk(zside); }

        // FIXME: interface to be extended for outside-in 
        const FTLDiskGeomDet * nextDisk(const FTLDiskGeomDet *from, PropagationDirection direction=alongMomentum) const ;

        const FTLDiskGeomDet * idToDet(DetId id) const ;

    protected:
        const FastTimeGeometry *geom_;
        std::vector<FTLDiskGeomDet *> disksPos_, disksNeg_;

        void makeDisks(int subdet, int disks) ;
        void addDisk(FTLDiskGeomDet *disk) { 
            (disk->zside() > 0 ? disksPos_ : disksNeg_).push_back(disk);
        }
};

#endif
