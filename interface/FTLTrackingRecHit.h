#ifndef RecoTracker_FastTimeMatching_FTLTrackingRecHit_h
#define RecoTracker_FastTimeMatching_FTLTrackingRecHit_h

/// Basic template class for a RecHit wrapping a Ref to an object

#include <cassert>
#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"

template<typename ObjRef>
class FTLTrackingRecHit : public RecHit2DLocalPos {
    public:
        typedef ObjRef ref_type;
        typedef typename ObjRef::value_type obj_type;

        FTLTrackingRecHit() {}
        FTLTrackingRecHit(DetId id, const ObjRef &objref, const LocalPoint &pos, const LocalError &err):
            RecHit2DLocalPos(id),
            objref_(objref),
            pos_(pos), err_(err) {}

        virtual FTLTrackingRecHit<ObjRef> * clone() const override { return new FTLTrackingRecHit<ObjRef>(*this); }

        virtual LocalPoint localPosition() const override { return pos_; }
        virtual LocalError localPositionError() const override { return err_; }

        const ObjRef & objRef() const { return objref_; }

        float energy() const { return objref_->energy(); }

        virtual bool sharesInput( const TrackingRecHit* other, SharedInputType what) const override { assert(false); }
    protected:
        ObjRef objref_;
        LocalPoint pos_;
        LocalError err_;
};


// Instantiations and specializations for FTLRecHitRef and reco::CaloClusterPtr
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
typedef FTLTrackingRecHit<FTLRecHitRef> FTLTrackingRecHitFromHit;
typedef FTLTrackingRecHit<reco::CaloClusterPtr> FTLTrackingRecHitFromCluster;

template<>
bool FTLTrackingRecHit<FTLRecHitRef>::sharesInput( const TrackingRecHit* other, SharedInputType what) const ;

template<>
bool FTLTrackingRecHit<reco::CaloClusterPtr>::sharesInput( const TrackingRecHit* other, SharedInputType what) const ;

#endif
