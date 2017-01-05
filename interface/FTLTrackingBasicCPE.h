#ifndef RecoTracker_FastTimeMatching_FTLTrackingBasicCPE_h
#define RecoTracker_FastTimeMatching_FTLTrackingBasicCPE_h

/// Class that estimates positions and uncertainties for a hit or cluster
//  Similar interface to the ClusterParameterEstimators of the tracking
#include "Geometry/HGCalGeometry/interface/FastTimeGeometry.h"

#include "RecoTracker/FastTimeMatching/interface/FTLTrackingRecHit.h"
#include "RecoTracker/FastTimeMatching/interface/FTLTrackingClusteringRecHit.h"

class FTLTrackingBasicCPE {
    public:
        struct PositionHint {
            PositionHint(float x_, float y_, float size_) : x(x_), y(y_), size(size_) {}
            const float x, y, size;
        };

        // For a circle, E<x^2> = 1/4 r^2 = A / (4 pi); HGcal sensors have 1cm area --> sqrt(E<x^2>) = 0.28
        FTLTrackingBasicCPE(const FastTimeGeometry* geom, float fixedError=0.28*1.2, float fixedErrorBH=3) ;

       std::pair<LocalPoint,LocalError> localParameters(const FTLRecHit &obj, const Surface &surf) const {
            float loce2 = (obj.id().det() == DetId::Forward ? fixedError2_ : fixedErrorBH2_);
            return std::make_pair(
                        surf.toLocal(getPosition(obj)),
                        LocalError(loce2,0,loce2)
                   );
        }
        std::pair<LocalPoint,LocalError> localParameters(const reco::CaloCluster & obj, const Surface &surf) const {
            float loce2 = (obj.hitsAndFractions().front().first.det() == DetId::Forward ? fixedError2_ : fixedErrorBH2_);
            return std::make_pair(
                        surf.toLocal(GlobalPoint(obj.position().X(), obj.position().Y(), obj.position().Z())),
                        LocalError(loce2 * obj.size(),0,loce2 * obj.size())
                   );
        }
        PositionHint hint(const FTLTrackingRecHitFromHit & hit) const { return hint(*hit.objRef()); }
        PositionHint hint(const FTLTrackingRecHitFromCluster & hit) const { return hint(*hit.objRef()); }
        PositionHint hint(const FTLRecHit & obj) const {
            float loce = (obj.id().det() == DetId::Forward ? fixedError_ : fixedErrorBH_);
            const GlobalPoint & gp = getPosition(obj);
            return PositionHint(gp.x(), gp.y(), loce);
        }
        PositionHint hint(const reco::CaloCluster & obj) const {
            float loce = (obj.hitsAndFractions().front().first.det() == DetId::Forward ? fixedError_ : fixedErrorBH_);
            const math::XYZPoint &gp = obj.position();
            return PositionHint(gp.X(), gp.Y(), loce * sqrt(obj.size()));
        }
        GlobalPoint getPosition(const FTLRecHit &obj) const { return getPosition(obj.id()); }
        GlobalPoint getPosition(const DetId id) const {
	  return geom_->getPosition(id);            
        }
    private:
        const FastTimeGeometry *geom_;
        const float fixedError_, fixedError2_;
        const float fixedErrorBH_, fixedErrorBH2_;
};

#endif
