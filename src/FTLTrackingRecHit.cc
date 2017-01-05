#include "RecoTracker/FastTimeMatching/interface/FTLTrackingRecHit.h"
#include "RecoTracker/FastTimeMatching/interface/FTLTrackingClusteringRecHit.h"

template<>
bool FTLTrackingRecHit<FTLRecHitRef>::sharesInput( const TrackingRecHit* other, SharedInputType what) const 
{
    if (typeid(*other) == typeid(FTLTrackingRecHit<FTLRecHitRef>)) {
        return objRef() == (static_cast<const FTLTrackingRecHit<FTLRecHitRef> *>(other))->objRef();
    } else if (typeid(*other) == typeid(FTLTrackingClusteringRecHit)) {
        if (what == TrackingRecHit::some) {
            const auto &refv = (static_cast<const FTLTrackingClusteringRecHit *>(other))->objRefs();
            return std::find(refv.begin(), refv.end(), objRef()) != refv.end();
        } else {
            return objRef() == (static_cast<const FTLTrackingClusteringRecHit *>(other))->objRef();
        }
    } else if (typeid(*other) == typeid(FTLTrackingRecHit<reco::CaloClusterPtr>)) {
        return other->sharesInput( this, what ); // see below
    } else {
        return false;
    }
}

template<>
bool FTLTrackingRecHit<reco::CaloClusterPtr>::sharesInput( const TrackingRecHit* other, SharedInputType what) const 
{
    if (typeid(*other) == typeid(FTLTrackingRecHit<reco::CaloClusterPtr>)) {
        return objRef() == (static_cast<const FTLTrackingRecHit<reco::CaloClusterPtr> *>(other))->objRef();
    } else if (typeid(*other) == typeid(FTLTrackingRecHit<FTLRecHitRef>)) {
        DetId otherId = (static_cast<const FTLTrackingRecHit<FTLRecHitRef> *>(other))->objRef()->id();
        for (const auto & pair : objref_->hitsAndFractions()) {
            if (pair.second == 0) continue;
            if (pair.first == otherId) return true;
        }
        return false;
    } else if (typeid(*other) == typeid(FTLTrackingClusteringRecHit)) {
        const auto &refv = (static_cast<const FTLTrackingClusteringRecHit *>(other))->objRefs();
        for (FTLRecHitRef ref : refv)  {
            DetId otherId = ref->id();
            for (const auto & pair : objref_->hitsAndFractions()) {
                if (pair.second == 0) continue;
                if (pair.first == otherId) return true;
            }
        }
        return false;
    } else {
        return false;
    }
}


