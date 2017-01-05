import FWCore.ParameterSet.Config as cms

from RecoTracker.FastTimeMatching.ftlTrajectoryBuilder_cfi import *

ftlTracks = cms.EDProducer("FTLTracking",
    ftlTrajectoryBuilderPSet,

    ### Track collection to use for seeding
    srcTk = cms.InputTag("generalTracks"),
    cutTk = cms.string("0.0 < abs(eta) < 3.0 && quality('highPurity')"),

    ### Attempt a backwards refit of the HGC hits to get an unbiased track 
    doBackwardsRefit = cms.bool(False),

    ### MC Truth, for debugging
    srcTruth = cms.InputTag("mix","MergedCaloTruth"),

    # verbosity of printout (0 = none)
    debugLevel = cms.untracked.uint32(0)
)

