#ifndef RecoTracker_FastTimeMatching_FTLTrackingDiskData_h
#define RecoTracker_FastTimeMatching_FTLTrackingDiskData_h

/// Class that holds the rechits on one FTL endcap Layer, and does the lookup
//  Almost like a MeasurementDetWithData in the tracker
//
//  note 1: the class is typedeffed on the rechit collection, but is not a template
//  note 2: the class is cheap to move, but expensive to copy


#include "RecoTracker/FastTimeMatching/interface/FTLTrackingRecHit.h"
#include "RecoTracker/FastTimeMatching/interface/FTLTrackingClusteringRecHit.h"
#include "RecoTracker/FastTimeMatching/interface/FTLTracker.h"
#include "RecoTracker/FastTimeMatching/interface/FTLTrackingBasicCPE.h"

#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

#include "RecoTracker/FastTimeMatching/interface/ftldebug.h"


class FTLTrackingDiskData {
    public:
        typedef FTLRecHitCollection TColl;
        typedef TColl::value_type value_type;
        typedef TColl::const_iterator const_iterator;
        typedef edm::Ref<TColl> ref_type;
        typedef edm::RefVector<TColl> refvector_type;
        typedef FTLTrackingRecHit<ref_type> rechit_type;
        typedef FTLTrackingClusteringRecHit rechitcluster_type;

        FTLTrackingDiskData() {}
        FTLTrackingDiskData(const edm::Handle<TColl> &data, int subdet, int zside, int layer, const FTLTrackingBasicCPE *cpe) ;

        void addClusters(const edm::Handle<reco::CaloClusterCollection> &data, int subdet, int zside, int layer) ;

        unsigned int size() const { return index_.size(); }
        unsigned int nclusters() const { return clusterIndex_.size(); }

        std::vector<TrajectoryMeasurement> measurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest)  const ;
        std::vector<TrajectoryMeasurement> clusterizedMeasurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest, float rCut) const ;
        std::vector<TrajectoryMeasurement> clusterMeasurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest) const ;

        // for debugging
        void setTruth(const CaloTruthRevMap *truth) const { truthMap_ = truth; }

    private:
        const edm::Handle<TColl> *alldata_;
        std::vector<std::pair<FTLTrackingBasicCPE::PositionHint,const_iterator> > index_;

        const edm::Handle<reco::CaloClusterCollection> *clusterData_;
        std::vector<std::pair<FTLTrackingBasicCPE::PositionHint,reco::CaloClusterCollection::const_iterator> > clusterIndex_;

        const FTLTrackingBasicCPE *cpe_;

        // for debugging
        mutable const CaloTruthRevMap *truthMap_;

        void buildIndex_() {
            // could do some sorting here
        }
    
};

#endif
