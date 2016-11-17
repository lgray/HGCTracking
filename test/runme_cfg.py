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
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_8_1_0_pre16/RelValTenMuExtendedE_0_200/GEN-SIM-RECO/81X_upgrade2023_realistic_v3_2023D3-v1/10000/68CF8573-DBA6-E611-972F-0025905B855C.root'),
    secondaryFileNames = cms.untracked.vstring()
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# customisation of the process.

process.load('RecoLocalCalo.HGCalRecHitDump.imagingClusterHGCal_cfi')
process.load('RecoParticleFlow.HGCTracking.hgcTracks_cff')

process.hgcTracks.debugLevel = 3

process.run = cms.Path(process.imagingClusterHGCal + process.hgcTracks)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("hgctracks.root"),
    outputCommands = cms.untracked.vstring("drop *", 
        "keep recoTracks_generalTracks_*_*",
        "keep *_HGCalRecHit_*_*", 
        "keep *_particleFlowClusterHGCal_*_*", 
        "keep *_imagingClusterHGCal_*_*",
        "keep *_hgcTracks_*_*",
    ),
)
process.e = cms.EndPath(process.out)

import sys
args = sys.argv[2:] if sys.argv[0] == "cmsRun" else sys.argv[1:]
if args:
    process.source.fileNames = [ args[0] ]
