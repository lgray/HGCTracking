import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('TRY',eras.Phase2C2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D3Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('/store/cmst3/user/gpetrucc/l1phase2/081116/DR_PU0/2023D3_guns/ChargedPionGun_PU0_job1.FEVT.root'),
    fileNames = cms.untracked.vstring('file:/data/g/gpetrucc/ChargedPionGun_PU0_job1.FEVT.root'),
    secondaryFileNames = cms.untracked.vstring()
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# customisation of the process.

process.load('RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cfi')
process.load('RecoParticleFlow.HGCTracking.ftlTracks_cff')

process.ftlTracks.debugLevel = 3
process.ftlTracks.patternRecoAlgo = "hitsAndClusters"

process.run = cms.Path(process.ftlTracks)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("ftltracks_clusters.root"),
    outputCommands = cms.untracked.vstring("drop *", 
        "keep recoTracks_generalTracks_*_*",
        "keep *_HGCalRecHit_*_*", 
        "keep *_particleFlowClusterHGCal_*_*", 
        "keep *_hgcalLayerClusters_*_*",
        "keep *_ftlTracks_*_*",
    ),
)
process.e = cms.EndPath(process.out)

import sys
args = sys.argv[2:] if sys.argv[0] == "cmsRun" else sys.argv[1:]
if args:
    process.source.fileNames = [ args[0] ]
