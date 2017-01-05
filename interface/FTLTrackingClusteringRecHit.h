#ifndef RecoTracker_FastTimeMatching_FTLTrackingClusteringRecHit_h
#define RecoTracker_FastTimeMatching_FTLTrackingClusteringRecHit_h

/// Class for a TrackingRecHit made from an on-the-fly clustering of FTLRecHits

#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

class FTLTrackingClusteringRecHit : public RecHit2DLocalPos {
    public:
        typedef FTLRecHitRef ObjRef;
        typedef edm::RefVector<ObjRef::product_type, ObjRef::value_type, ObjRef::finder_type> ObjRefVector;

        FTLTrackingClusteringRecHit() {}
        FTLTrackingClusteringRecHit(DetId id, const ObjRefVector &objrefs, float etot, const LocalPoint &pos, const LocalError &err):
            RecHit2DLocalPos(id),
            objrefs_(objrefs), etot_(etot),
            pos_(pos), err_(err) {}

        virtual FTLTrackingClusteringRecHit * clone() const override { return new FTLTrackingClusteringRecHit(*this); }

        virtual LocalPoint localPosition() const override { return pos_; }
        virtual LocalError localPositionError() const override { return err_; }

        ObjRef objRef() const { return objrefs_[0]; }
        const ObjRefVector & objRefs() const { return objrefs_; }

        float energy() const { return etot_; }

        virtual bool sharesInput( const TrackingRecHit* other, SharedInputType what) const override ;

    protected:
        ObjRefVector objrefs_;
        float etot_;
        LocalPoint pos_;
        LocalError err_;
};

#endif

