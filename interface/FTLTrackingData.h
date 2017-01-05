#ifndef RecoTracker_FastTimeMatching_FTLTrackingData_h
#define RecoTracker_FastTimeMatching_FTLTrackingData_h

#include "RecoTracker/FastTimeMatching/interface/FTLTrackingDiskData.h"
#include "RecoTracker/FastTimeMatching/interface/FTLTracker.h"
#include <map>

/// Class that holds the rechits of FTL Layer, and does the lookup
//  Almost like a MeasurementTrackerEvent in the tracker
//
//  note: the class may be expensive both to copy and to move


class FTLTrackingData {
    public:
        typedef FTLTrackingDiskData Disk;
        typedef Disk::TColl TColl;

        FTLTrackingData(const FTLTracker &tracker, const FTLTrackingBasicCPE *cpe) : tracker_(&tracker), cpe_(cpe) {}
        void addData(const edm::Handle<TColl> &data, int subdet) ;
        void addClusters(const edm::Handle<reco::CaloClusterCollection> &data) ;
        const Disk & diskData(const FTLDiskGeomDet *disk) const { 
            auto it = data_.find(disk);
            if (it == data_.end()) throwBadDisk_(disk);
            return it->second; 
        }

    private:
        const FTLTracker * tracker_;
        const FTLTrackingBasicCPE* cpe_;
        std::map<const FTLDiskGeomDet *, Disk> data_; // FIXME better type later

        void throwBadDisk_(const FTLDiskGeomDet *disk) const ;
};


#endif
