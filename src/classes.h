#include "RecoTracker/FastTimeMatching/interface/FTLTrackingRecHit.h"
#include "RecoTracker/FastTimeMatching/interface/FTLTrackingClusteringRecHit.h"
#include "RecoTracker/FastTimeMatching/interface/TrajectorySeedFromTrack.h"

namespace RecoTracker_FastTimeMatching {
    struct dictionary {
        FTLTrackingRecHitFromHit dummy_hit;
        FTLTrackingRecHitFromCluster dummy_hit_cluster;
        FTLTrackingClusteringRecHit dummy_clustering_hit; 
    };
}
