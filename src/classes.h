#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/TrajectorySeedFromTrack.h"

namespace RecoParticleFlow_HGCTracking {
    struct dictionary {
        HGCTrackingRecHitFromHit dummy_hit;
        HGCTrackingRecHitFromCluster dummy_hit_cluster;
        HGCTrackingClusteringRecHit dummy_clustering_hit; 
    };
}
