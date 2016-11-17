/** \class HGCTracking
 *   produce Candidates for each RecHit **/

#include <cassert>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cmath>

#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTracker.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingBasicCPE.h"
#include "RecoParticleFlow/HGCTracking/interface/TrajectorySeedFromTrack.h"
#include "RecoParticleFlow/HGCTracking/interface/TrajectoryCleanerBySharedEndpoints.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingData.h"
#include "RecoParticleFlow/HGCTracking/interface/hgcdebug.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/CkfPattern/interface/IntermediateTrajectoryCleaner.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryCleaning/interface/FastTrajectoryCleaner.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilterFactory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"


class HGCTracking : public edm::stream::EDProducer<> {
    public:
        explicit HGCTracking(const edm::ParameterSet& ps);
        ~HGCTracking() {}
        virtual void produce(edm::Event& evt, const edm::EventSetup& es) override;

    private:
        enum PatternRecoAlgo { SingleHitAlgo, ClusterizingAlgo, MixedAlgo };

        edm::EDGetTokenT<reco::TrackCollection> srcTk_;
        StringCutObjectSelector<reco::Track> cutTk_;
        edm::EDGetTokenT<HGCRecHitCollection> srcEE_, srcFH_, srcBH_;
        edm::EDGetTokenT<reco::CaloClusterCollection> srcClusters_;
        double lostHitRescaleErrorFactor_;
        PatternRecoAlgo algo_;
        std::string propName_, propNameOppo_;
        std::string estimatorName_;
        std::string updatorName_;
        std::string singleCleanerName_, trajectoryCleanerName_;
        bool fastCleaner_, endpointCleaner_;
        unsigned int theMaxCand;
        double theFoundHitBonus, theLostHitPenalty;
        unsigned int theMaxStartingEmptyLayers, theLayersBeforeCleaning;
        double deltaChiSquareForHits_, minChi2ForInvalidHit_;
        double clusterRadius_;
        bool hasMCTruth_;
        edm::EDGetTokenT<std::vector<CaloParticle>> srcTruth_;
        CaloTruthRevMap revTruthMap_;

        uint32_t geomCacheId_;
        std::unique_ptr<HGCTracker> disks_;
        edm::ESHandle<CaloGeometry> caloGeom_;
        edm::ESHandle<GlobalTrackingGeometry> trkGeom_;
        edm::ESHandle<Propagator> prop_, propOppo_;
        edm::ESHandle<MagneticField> bfield_;
        edm::ESHandle<Chi2MeasurementEstimatorBase> estimator_;
        edm::ESHandle<TrajectoryStateUpdator> updator_;
        edm::ESHandle<TrajectoryCleaner> trajectoryCleaner_;
        std::unique_ptr<TrajectoryFilter> trajFilter_;

        template<class Start>
        std::vector<TempTrajectory> advanceOneLayer(const Start &start, const TempTrajectory &traj, const HGCDiskGeomDet *disk, const HGCTrackingData &data) ;

        Trajectory bwrefit(const Trajectory &traj, float scaleErrors = 100.) const ;

        void declareOutput(const std::string &label); 
        void makeOutput(const std::vector<Trajectory> &trajs, edm::Event &out, const std::string &label); 
        void writeAllHitsByLayer(const HGCRecHitCollection &hits, edm::Event &out, const std::string &label, int layers);

        void trim(std::vector<TempTrajectory> &tempTrajectories) const {
            unsigned int nbad = 0;
            for (const auto &t : tempTrajectories) { if (!t.isValid()) nbad++; }
            if (hgctracking::g_debuglevel > 1) if (nbad) printf("\t\t\tpruning %d bad trajectories\n", nbad);
            tempTrajectories.erase(std::remove_if( tempTrajectories.begin(),tempTrajectories.end(), std::not1(std::mem_fun_ref(&TempTrajectory::isValid))), tempTrajectories.end());
        }
        void trim(std::vector<Trajectory> &trajectories) const {
            unsigned int nbad = 0;
            for (const auto &t : trajectories) { if (!t.isValid()) nbad++; }
            if (hgctracking::g_debuglevel > 1) if (nbad) printf("\t\t\tpruning %d bad trajectories\n", nbad);
            trajectories.erase(std::remove_if( trajectories.begin(),trajectories.end(), std::not1(std::mem_fun_ref(&Trajectory::isValid))), trajectories.end());
        }
        template<typename Traj> void printTraj(const Traj &t) const {
            TrajectoryStateOnSurface finalTSOS = t.lastMeasurement().updatedState();
            if (!finalTSOS.isValid()) finalTSOS = t.lastMeasurement().forwardPredictedState();
            GlobalPoint pos = finalTSOS.globalPosition(); GlobalVector mom = finalTSOS.globalMomentum();
            LocalTrajectoryError lte = finalTSOS.localError();
            //printf("found hits %3d, lost hits %3d, chi2 %8.1f ndf %3d outermost state pt = %6.1f +- %5.1f   eta %+5.2f phi %+5.2f x = %+7.2f +- %4.2f  y = %+7.2f +- %4.2f   z = %+7.2f",
            printf("found hits %3d, lost hits %3d, chi2 %8.1f ndf %3d outermost pt = %6.1f x = %+7.2f y = %+7.2f z = %+7.2f last-id %12d first-id %12d ",
                    t.foundHits(), t.lostHits(), t.chiSquared(), t.foundHits()*2-5, 
                    mom.perp(), //mom.perp() * mom.mag() * sqrt(lte.matrix()(0,0)),
                    //pos.eta(), float(pos.phi()), 
                    pos.x(), /*sqrt(lte.matrix()(3,3)),*/ pos.y(), /*sqrt(lte.matrix()(4,4)),*/
                    pos.z(),
                    t.lastMeasurement().recHit()->isValid() ? t.lastMeasurement().recHit()->geographicalId()() : 0,
                    t.firstMeasurement().recHit()->isValid() ? t.firstMeasurement().recHit()->geographicalId()() : 0
                    );
            // median energy over the last 4 hits
            if (t.foundHits() >= 3) {
                float etotal = 0;
                std::vector<float> energies;
                const auto &meas = t.measurements();
                for (int i = meas.size()-1; i >= 0; --i) {
                    if (meas[i].recHit()->isValid()) {
                        float energy = -1;
                        if (typeid(*meas[i].recHit()) == typeid(HGCTrackingRecHitFromHit)) {
                           energy = ((dynamic_cast<const HGCTrackingRecHitFromHit&>(*meas[i].recHit())).energy());
                        } else if (typeid(*meas[i].recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
                           energy = ((dynamic_cast<const HGCTrackingRecHitFromCluster&>(*meas[i].recHit())).energy());
                        } else if (typeid(*meas[i].recHit()) == typeid(HGCTrackingClusteringRecHit)) {
                           energy = ((dynamic_cast<const HGCTrackingClusteringRecHit&>(*meas[i].recHit())).energy());
                        } else {
                           throw cms::Exception("Invalid valid hit", typeid(*meas[i].recHit()).name());
                        }
                        if (energies.size() < 4) energies.push_back(energy);
                        etotal += energy;
                    }
                }
                if (energies.size()) {
                    std::sort(energies.begin(),energies.end());
                    float median = (energies.size() % 2 == 1) ? energies[energies.size()/2] : 0.5*(energies[energies.size()/2]+energies[energies.size()/2-1]);
                    //printf("   energy median %7.3f (%lu)", median, energies.size());
                    printf("   energy %6.3f et_tot %6.2f  ", median, etotal * pos.perp() / pos.mag());
                }
            }
            if (hasMCTruth_) {
                for (const auto & pair : truthMatch(t)) {
                    printf("  %.1f hits from %s pdgId %+d eid %d/%+d pt %.1f eta %+5.2f phi %+5.2f    ", pair.second, 
                        (pair.first->eventId().event()==0&&pair.first->eventId().bunchCrossing()==0 ? "SIGNAL" : "pileup"),
                        pair.first->pdgId(), pair.first->eventId().event(), pair.first->eventId().bunchCrossing(), pair.first->pt(), pair.first->eta(), pair.first->phi());
                }
            }
            printf("\n");
        }

        template<typename Traj>
        std::vector<std::pair<const CaloParticle *, float>> truthMatch(const Traj &t) const {
            std::vector<std::pair<const CaloParticle *, float>> ret;
            if (hasMCTruth_) { 
                std::map<const CaloParticle *,float> scores;       
                std::vector<const CaloParticle *> keys;
                std::vector<DetId> ids;
                for (const auto & tm : t.measurements()) {
                    if (!tm.recHit()->isValid()) continue;
                    ids.clear();
                    if (typeid(*tm.recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
                        for (auto & p : (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*tm.recHit())).objRef()->hitsAndFractions()) {
                            if (p.second > 0) ids.push_back(p.first);
                        }
                    } else {
                        ids.push_back(tm.recHit()->geographicalId());
                    }
                    for (DetId detid: ids) {
                        auto range = revTruthMap_.equal_range(detid.rawId());
                        for (; range.first != range.second; ++range.first) {
                            const auto &pair = *range.first;
                            if (std::find(keys.begin(), keys.end(), pair.second.first) == keys.end()) {
                                keys.push_back(pair.second.first);
                                scores[pair.second.first] += 1;
                            }
                            //scores[pair.second.first] += pair.second.second > 0.0 ? 1 : pair.second.second;
                        }   
                    }
                    keys.clear();
                }
                keys.clear();
                for (auto & pair : scores) keys.push_back(pair.first);
                std::sort(keys.begin(), keys.end(), [&scores](const CaloParticle *a, const CaloParticle *b) -> bool { return scores[b] < scores[a]; });
                for (auto & key : keys) {
                    ret.emplace_back(key, scores[key]);
                }
            }
            return ret;
        }
};


HGCTracking::HGCTracking(const edm::ParameterSet& ps) :
    srcTk_(consumes<reco::TrackCollection>(ps.getParameter<edm::InputTag>("srcTk"))),
    cutTk_(ps.getParameter<std::string>("cutTk")),
    srcEE_(consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcEE"))),
    srcFH_(consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcFH"))),
    srcBH_(consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcBH"))),
    srcClusters_(consumes<reco::CaloClusterCollection>(ps.getParameter<edm::InputTag>("srcClusters"))),
    lostHitRescaleErrorFactor_(ps.getParameter<double>("lostHitRescaleErrorFactor")),
    propName_(ps.getParameter<std::string>("propagator")),
    propNameOppo_(ps.getParameter<std::string>("propagatorOpposite")),
    estimatorName_(ps.getParameter<std::string>("estimator")),
    updatorName_(ps.getParameter<std::string>("updator")),
    trajectoryCleanerName_(ps.getParameter<std::string>("trajectoryCleaner")),
    fastCleaner_(ps.getParameter<bool>("fastCleaner")), 
    endpointCleaner_(ps.getParameter<bool>("endpointCleaner")), 
    theMaxCand(ps.getParameter<uint32_t>("maxCand")),
    theFoundHitBonus(ps.getParameter<double>("foundHitBonus")),
    theLostHitPenalty(ps.getParameter<double>("lostHitPenalty")),
    theMaxStartingEmptyLayers(ps.getParameter<uint32_t>("maxStartingEmptyLayers")),
    theLayersBeforeCleaning(ps.getParameter<uint32_t>("layersBeforeCleaning")),
    deltaChiSquareForHits_(ps.getParameter<double>("deltaChiSquareForHits")),
    minChi2ForInvalidHit_(ps.getParameter<double>("minChi2ForInvalidHit")),
    hasMCTruth_(ps.existsAs<edm::InputTag>("srcTruth")),
    srcTruth_(hasMCTruth_ ? consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("srcTruth")) : edm::EDGetTokenT<std::vector<CaloParticle>>()),
    geomCacheId_(0)
{
    hgctracking::g_debuglevel = ps.getUntrackedParameter<uint32_t>("debugLevel",0);

    std::string palgo = ps.getParameter<std::string>("patternRecoAlgo");
    if (palgo == "singleHit") {
        algo_ = SingleHitAlgo;
    } else if (palgo == "hitsAndClusters") {
        algo_ = MixedAlgo;
    } else if (palgo == "clusterizing") {
        algo_ = ClusterizingAlgo;
        clusterRadius_ = ps.getParameter<double>("clusterRadius");
    } else throw cms::Exception("Configuration") << "Unknown algo: knowns are 'singleHit', 'clusterizing'\n";
    auto pset = ps.getParameter<edm::ParameterSet>("trajectoryFilter");
    edm::ConsumesCollector && iC = consumesCollector();
    trajFilter_.reset( TrajectoryFilterFactory::get()->create(pset.getParameter<std::string>("ComponentType"), pset, iC) );

    declareOutput("");
    //declareOutput("bw");
    if (hgctracking::g_debuglevel > 0) {
        produces<std::vector<reco::CaloCluster>>("EE");
        produces<std::vector<reco::CaloCluster>>("FH");
        produces<std::vector<reco::CaloCluster>>("BH");
    }
}

void
HGCTracking::produce(edm::Event& evt, const edm::EventSetup& es)
{
    //Get Calo Geometry
    if (es.get<CaloGeometryRecord>().cacheIdentifier() != geomCacheId_) {
        es.get<CaloGeometryRecord>().get(caloGeom_);
        geomCacheId_ = es.get<CaloGeometryRecord>().cacheIdentifier();
        disks_.reset(new HGCTracker(caloGeom_.product()));
    } 

    es.get<GlobalTrackingGeometryRecord>().get(trkGeom_);
    es.get<IdealMagneticFieldRecord>().get(bfield_);
    es.get<TrackingComponentsRecord>().get(propName_, prop_);
    if (!propNameOppo_.empty()) es.get<TrackingComponentsRecord>().get(propNameOppo_, propOppo_);
    es.get<TrackingComponentsRecord>().get(estimatorName_,estimator_);
    es.get<TrackingComponentsRecord>().get(updatorName_,updator_);

    if (!trajectoryCleanerName_.empty()) es.get<TrajectoryCleaner::Record>().get(trajectoryCleanerName_, trajectoryCleaner_);
    TrajectoryCleanerBySharedEndpoints endpointCleaner(theFoundHitBonus, theLostHitPenalty);

    HGCTrackingBasicCPE cpe(&*caloGeom_); // FIXME better
    auto trajCandLess = [&](TempTrajectory const & a, TempTrajectory const & b) {
        return  (a.chiSquared() + a.lostHits()*theLostHitPenalty - a.foundHits()*theFoundHitBonus)  <
            (b.chiSquared() + b.lostHits()*theLostHitPenalty - b.foundHits()*theFoundHitBonus);
    };
    
    trajFilter_->setEvent(evt, es); 

    edm::Handle<std::vector<CaloParticle>> truth;
    if (hasMCTruth_) {
        evt.getByToken(srcTruth_, truth);
        for (const CaloParticle &p : *truth) {
            if (hgctracking::g_debuglevel > 0) {
                if (p.eventId().event() == 0 && p.eventId().bunchCrossing() == 0) {
                    printf("Considering truth particle of pdgId %+6d eid %d/%d pt %7.1f eta %+5.2f phi %+5.2f simclusters %4d \n", 
                            p.pdgId(), p.eventId().event(), p.eventId().bunchCrossing(), p.pt(), p.eta(), p.phi(), int(p.simClusters().size()));
                }
            }
            for (const auto & scr : p.simClusters()) {
                //printf("    simcluster pt %7.1f eta %+5.2f phi %+5.2f rechits %4d simhits %d \n",
                //    scr->pt(), scr->eta(), scr->phi(), scr->numberOfRecHits(), scr->numberOfSimHits());
                for (const auto & pair : scr->hits_and_fractions()) {
                    revTruthMap_.emplace(pair.first, std::make_pair(&p, pair.second));
                }
            }
        }
    }

    HGCTrackingData data(*disks_, &cpe);
    edm::Handle<HGCRecHitCollection> srcEE, srcFH, srcBH;
    evt.getByToken(srcEE_, srcEE); data.addData(srcEE, 3);
    evt.getByToken(srcFH_, srcFH); data.addData(srcFH, 4);
    evt.getByToken(srcBH_, srcBH); data.addData(srcBH, 5);
    if (hgctracking::g_debuglevel > 0) {
        writeAllHitsByLayer(*srcEE, evt, "EE", 28);
        writeAllHitsByLayer(*srcFH, evt, "FH", 12);
        writeAllHitsByLayer(*srcBH, evt, "BH", 12);
    }

    edm::Handle<reco::CaloClusterCollection> srcClusters;
    evt.getByToken(srcClusters_, srcClusters); data.addClusters(srcClusters);

    edm::Handle<reco::TrackCollection> tracks;
    evt.getByToken(srcTk_, tracks);
    std::vector<Trajectory> finaltrajectories;
    std::vector<TempTrajectory> myfinaltrajectories;
    unsigned int itrack = 0;
    for (const reco::Track &tk : *tracks) { ++itrack;
        if (!cutTk_(tk)) continue;
        if (hgctracking::g_debuglevel > 1) {
            printf("\n\nConsidering track pt %7.1f eta %+5.2f phi %+5.2f valid hits %d outer lost hits %d highPurity %1d\n", tk.pt(), tk.eta(), tk.phi(), tk.hitPattern().numberOfValidHits(), tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS), tk.quality(reco::Track::highPurity));
        }
        int zside = tk.eta() > 0 ? +1 : -1;
        const HGCDiskGeomDet *disk = disks_->firstDisk(zside);
        FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(tk, bfield_.product());
        fts.rescaleError(1.0 + lostHitRescaleErrorFactor_ * tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));
        std::vector<TempTrajectory> trajectories = advanceOneLayer(fts, TempTrajectory(alongMomentum, 0), disk, data);
        myfinaltrajectories.clear();
        unsigned int depth = 2;
        for (disk = disks_->nextDisk(disk); disk != nullptr; disk = disks_->nextDisk(disk), ++depth) {
            if (trajectories.empty()) continue;
            if (hgctracking::g_debuglevel > 1) {
                printf("   New destination: disk subdet %d, zside %+1d, layer %2d, z = %+8.2f\n", disk->subdet(), disk->zside(), disk->layer(), disk->toGlobal(LocalPoint(0,0,0)).z());
            }
            //printf("   Starting candidates: %lu\n", trajectories.size());
            std::vector<TempTrajectory> newCands;
            int icand = 0;
            for (TempTrajectory & cand : trajectories) {
                if (hgctracking::g_debuglevel > 1) {
                    printf("    Processing candidate %2d/%lu with ", ++icand, trajectories.size()); printTraj(cand); //found hits %3d, lost hits %3d, chi2 %8.1f\n", cand.foundHits(), cand.lostHits(), cand.chiSquared());
                }
                TrajectoryStateOnSurface start = cand.lastMeasurement().updatedState();
                if (!start.isValid()) start = cand.lastMeasurement().predictedState();
                std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, cand, disk, data);
                if (hisTrajs.empty() ) {
                    if (hgctracking::g_debuglevel > 1) printf("     --> stops here\n");
                    if (trajFilter_->qualityFilter(cand)) {
                        if (hgctracking::g_debuglevel > 1) printf("          --> passes filter, being retained for the moment \n");
                        myfinaltrajectories.push_back(cand);
                    }
                } else {
                    //printf("---- produced %lu trajectories \n", hisTrajs.size());
                    for ( TempTrajectory & t : hisTrajs ) { 
                        if (t.foundHits() == 0 && depth > theMaxStartingEmptyLayers) continue;
                        if (trajFilter_->toBeContinued(t)) {
                            newCands.push_back(t);
                        } else {
                            //printf("------    --> one is not to be continued.\n");
                            if (trajFilter_->qualityFilter(t)) {
                                //printf("------         --> but it is to be saved.\n");
                                myfinaltrajectories.push_back(std::move(t));
                            }
                        }
                    }
                }
            }
            //printf("   A total of %lu trajectories after this step\n", newCands.size());
            if (endpointCleaner_ && depth > theLayersBeforeCleaning) {
                unsigned int oldsize = newCands.size();
                endpointCleaner.clean(newCands); trim(newCands);
                if (hgctracking::g_debuglevel > 1) if (oldsize != newCands.size()) printf("    Reduced from %u to %lu trajectories after TrajectoryCleanerBySharedEndpoints\n", oldsize, newCands.size());
            }
            if (newCands.size() > theMaxCand && depth > theLayersBeforeCleaning) {
                std::sort(newCands.begin(), newCands.end(), trajCandLess);
                newCands.resize(theMaxCand);
                //printf("    Reduced to %lu trajectories after sorting and trimming\n", newCands.size());
            }
            trajectories.swap(newCands);
        }
        if (!trajectories.empty()) {
            if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories reached the end of the tracker from this track\n", trajectories.size());
            for (TempTrajectory & t : trajectories) {
                if (t.foundHits() > 0 && trajFilter_->qualityFilter(t)) {
                    myfinaltrajectories.push_back(std::move(t));
                }
            }
        }
        if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories found from this track\n", myfinaltrajectories.size());
        if (fastCleaner_) {
            FastTrajectoryCleaner cleaner(theFoundHitBonus/2,theLostHitPenalty); // the factor 1/2 is because the FastTrajectoryCleaner uses the ndof instead of the number of found hits
            cleaner.clean(myfinaltrajectories); trim(myfinaltrajectories);
            if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories after FastTrajectoryCleaner\n", myfinaltrajectories.size());
        }
        if (!myfinaltrajectories.empty()) {
            auto ostate = trajectoryStateTransform::outerStateOnSurface(tk, *trkGeom_, &*bfield_);
            auto pstate = trajectoryStateTransform::persistentState(ostate, DetId(tk.outerDetId()));
            boost::shared_ptr<TrajectorySeed> seed(new TrajectorySeedFromTrack(reco::TrackRef(tracks,itrack-1), pstate, alongMomentum));
            std::vector<Trajectory> toPromote;
            for (TempTrajectory &t : myfinaltrajectories) {
                while (!t.lastMeasurement().recHit()->isValid()) t.pop();
                if (hgctracking::g_debuglevel > 1) { printf("- TempTrajectory "); printTraj(t); }
                toPromote.push_back(t.toTrajectory());
            }
            if (toPromote.size() > 1 && !trajectoryCleanerName_.empty()) {
                trajectoryCleaner_->clean(toPromote); trim(toPromote);
                if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories after multiple cleaner\n", toPromote.size());
            }
            for (Trajectory &t : toPromote) { 
                if (hgctracking::g_debuglevel > 1) { printf("- Trajectory "); printTraj(t); }
                t.setSharedSeed(seed);
                finaltrajectories.push_back(std::move(t));
            }
            if (hgctracking::g_debuglevel > 1) printf("- from track pt %7.1f eta %+5.2f phi %+5.2f valid hits %2d outer lost hits %d highPurity %1d\n", tk.pt(), tk.eta(), tk.phi(), tk.hitPattern().numberOfValidHits(), tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS), tk.quality(reco::Track::highPurity));
        }
    }
    if (hgctracking::g_debuglevel > 0) printf("\n\nA total of %lu trajectories found in the event\n",finaltrajectories.size());
    if (!trajectoryCleanerName_.empty()) {
        trajectoryCleaner_->clean(finaltrajectories); trim(finaltrajectories);
        if (hgctracking::g_debuglevel > 0) printf("A total of %lu trajectories after multiple cleaner\n", finaltrajectories.size());
    }
    if (hgctracking::g_debuglevel > 0)  {
        std::set<const reco::Track *> done;
        for (const Trajectory &t : finaltrajectories) {
            printf("- Trajectory     "); printTraj(t);
            const reco::Track & tk  = * (dynamic_cast<const TrajectorySeedFromTrack &>(t.seed())).track();
            printf("\t from track pt %7.1f eta %+5.2f phi %+5.2f valid hits %2d outer lost hits %d highPurity %1d\n", tk.pt(), tk.eta(), tk.phi(), tk.hitPattern().numberOfValidHits(), tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS), tk.quality(reco::Track::highPurity));
            done.insert(&tk);
        }
        printf("Seed tracks not reconstructed:\n");
        for (const reco::Track &tk : *tracks) {
            if (!cutTk_(tk)) continue;
            if (done.count(&tk)) continue; // reconstructed
            printf("- Track pt %7.1f eta %+5.2f phi %+5.2f valid hits %2d outer lost hits %d highPurity %1d\n", tk.pt(), tk.eta(), tk.phi(), tk.hitPattern().numberOfValidHits(), tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS), tk.quality(reco::Track::highPurity));
        }
        printf("done.\n");
    }
    makeOutput(finaltrajectories, evt, "");

    if (!propNameOppo_.empty()) {
        //printf("Now going backwards\n");
        std::vector<Trajectory> bwfits;
        for (const Trajectory &t : finaltrajectories) {
            bwfits.push_back(bwrefit(t));
        }
        trim(bwfits);
        finaltrajectories.swap(bwfits);
        if (hgctracking::g_debuglevel > 0)  {
            for (const Trajectory &t : finaltrajectories) {
                printf("- TrajectoryBW "); printTraj(t);
            }
        }
    }

    revTruthMap_.clear();
}

template<class Start>
std::vector<TempTrajectory> 
HGCTracking::advanceOneLayer(const Start &start, const TempTrajectory &traj, const HGCDiskGeomDet *disk, const HGCTrackingData &data) {
    std::vector<TempTrajectory> ret;
    TrajectoryStateOnSurface tsos = prop_->propagate(start, disk->surface());
    if (!tsos.isValid()) { 
        if (hgctracking::g_debuglevel > 0)  {
            printf("        Destination disk subdet %d, zside %+1d, layer %2d, z = %+8.2f\n", disk->subdet(), disk->zside(), disk->layer(), disk->toGlobal(LocalPoint(0,0,0)).z());
            printf("         --> propagation failed.\n"); 
        }
        return ret; 
    }
    if (!disk->surface().bounds().inside(tsos.localPosition())) {
        if (hgctracking::g_debuglevel > 0)  {
            GlobalPoint gp = tsos.globalPosition();
            printf("        Prop point: global eta %+5.2f phi %+5.2f  x = %+8.2f y = %+8.2f z = %+8.2f rho = %8.2f\n", gp.eta(), float(gp.phi()), gp.x(), gp.y(), gp.z(), gp.perp());
            //LocalPoint lp = tsos.localPosition();
            //printf("            local                       x = %+8.2f y = %+8.2f z = %+8.2f rho = %8.2f\n", lp.x(), lp.y(), lp.z(), lp.perp());
            printf("         --> outside the bounds.\n"); 
        }
        return ret; 
    }
    const auto & diskData = data.diskData(disk);
    if (hasMCTruth_) diskData.setTruth(&revTruthMap_);
    //printf("        Looking for hits on disk subdet %d, zside %+1d, layer %2d: %u total hits\n", disk->subdet(), disk->zside(), disk->layer(), diskData.size());
    std::vector<TrajectoryMeasurement> meas;
    switch(algo_) {
        case SingleHitAlgo:    meas = diskData.measurements(tsos, *estimator_); break;
        case ClusterizingAlgo: meas = diskData.clusterizedMeasurements(tsos, *estimator_, clusterRadius_); break;
        case MixedAlgo:    
                meas = diskData.clusterMeasurements(tsos, *estimator_); 
                std::vector<TrajectoryMeasurement> meas2 = diskData.measurements(tsos, *estimator_); 
                if (meas.empty()) meas2.swap(meas);
            break;
    };
    std::sort(meas.begin(), meas.end(), TrajMeasLessEstim());
    if (hgctracking::g_debuglevel > 1)  printf("        Compatible hits: %lu\n", meas.size());
    for (const TrajectoryMeasurement &tm : meas) {
        if (deltaChiSquareForHits_ > 0) {
            if (meas.size() > 1  && !ret.empty() && tm.estimate() > meas.front().estimate() + deltaChiSquareForHits_) {
                //printf("        stop after the first %lu hits, since this chi2 of %.1f is too bad wrt the best one of %.1f\n", ret.size(), tm.estimate(), meas.front().estimate());
                break;
            }
        }
        TrajectoryStateOnSurface updated = updator_->update(tm.forwardPredictedState(), *tm.recHit());
        if (!updated.isValid()) { 
            if (hgctracking::g_debuglevel > 0)  {
                std::cout << "          Hit with chi2 = " << tm.estimate() << std::endl;
                std::cout << "              track state     " << tm.forwardPredictedState().localPosition() << std::endl;
                std::cout << "              rechit position " << tm.recHit()->localPosition() << std::endl;
                std::cout << "               --> failed update state" << std::endl;
            }
            continue;
        }
        //std::cout << "adding valid hit" << std::endl;
        ret.push_back(traj.foundHits() ? traj : TempTrajectory(traj.direction(),0)); // don't start with a lost hit
        ret.back().push(TrajectoryMeasurement(tm.forwardPredictedState(),
                                              updated,
                                              tm.recHit(),
                                              tm.estimate()), 
                        tm.estimate());
    }
    //std::cout << "adding invalid hit" << std::endl;
    if (minChi2ForInvalidHit_ > 0) {
        if (meas.size() > 0  && !ret.empty() && meas.front().estimate() < minChi2ForInvalidHit_) {
            //printf("        will not add the invalid hit after %lu valid hits, as the best valid hit has chi2 of %.1f\n", ret.size(), meas.front().estimate());
            return ret;
        }
    }
    ret.push_back(traj.foundHits() ? traj : TempTrajectory(traj.direction(),0)); // either just one lost hit, or a trajectory not starting on a lost hit
    ret.back().push(TrajectoryMeasurement(tsos, std::make_shared<InvalidTrackingRecHit>(*disk, TrackingRecHit::missing)));
    return ret;
}

Trajectory
HGCTracking::bwrefit(const Trajectory &traj, float scaleErrors) const {
    Trajectory ret(oppositeToMomentum);

    const Trajectory::DataContainer & tms = traj.measurements();

    ret.reserve(tms.size());

    TrajectoryStateOnSurface tsos = tms.back().updatedState();
    tsos.rescaleError(scaleErrors);

    for (int i = tms.size()-1; i >= 0; --i) {
        const HGCDiskGeomDet * det = disks_->idToDet(tms[i].recHit()->geographicalId());
        if (det == 0) {
            if (hgctracking::g_debuglevel > 0)  {
                printf(" ---> failure in finding det for step %d on det %d, subdet %d\n",i,tms[i].recHit()->geographicalId().det(),tms[i].recHit()->geographicalId().subdetId());
            }
            ret.invalidate();
            return ret;
        }
        TrajectoryStateOnSurface prop = propOppo_->propagate(tsos, det->surface());
        if (!prop.isValid()) {
            if (hgctracking::g_debuglevel > 0)  {
                printf(" ---> failure in propagation for step %d\n",i);
            }
            ret.invalidate();
            return ret;
        }
        if (tms[i].recHit()->isValid()) {
           auto pair = estimator_->estimate( prop, *tms[i].recHit() );
           if (hgctracking::g_debuglevel > 1)  {
               if (i % 7 == 0) {
               printf("   for step %2d pt %6.1f +- %5.1f   q/p %+7.4f +- %6.4f   x = %+7.2f +- %4.2f  y = %+7.2f +- %4.2f   dxdz = %+5.3f +- %4.3f dydz = %+5.3f +- %4.3f    adding hit %+7.2f %+7.2f  results in chi2 = %6.1f (old: %6.1f)\n", i, 
                    prop.globalMomentum().perp(), prop.globalMomentum().perp() * prop.globalMomentum().mag() * sqrt(prop.localError().matrix()(0,0)),
                    prop.globalParameters().signedInverseMomentum(), sqrt(prop.localError().matrix()(0,0)),
                    prop.localPosition().x(), sqrt(prop.localError().matrix()(3,3)), prop.localPosition().y(),  sqrt(prop.localError().matrix()(4,4)),
                    prop.localParameters().dxdz(), sqrt(prop.localError().matrix()(1,1)), prop.localParameters().dydz(), sqrt(prop.localError().matrix()(2,2)),
                    tms[i].recHit()->localPosition().x(), tms[i].recHit()->localPosition().y(),
                    pair.second, tms[i].estimate());
               }
           }
           TrajectoryStateOnSurface updated = updator_->update( prop, *tms[i].recHit() );
           if (!updated.isValid()) {
                if (hgctracking::g_debuglevel > 0) printf(" ---> fail in backwards update for step %d\n",i);
                ret.invalidate(); 
                return ret;
           } 
           //printf("       updated pt %6.1f +- %5.1f   q/p %+7.4f +- %6.4f   x = %+7.2f +- %4.2f  y = %+7.2f +- %4.2f   dxdz = %+5.3f +- %4.3f dydz = %+5.3f +- %4.3f\n",  
           //     updated.globalMomentum().perp(), updated.globalMomentum().perp() * updated.globalMomentum().mag() * sqrt(updated.localError().matrix()(0,0)),
           //     updated.globalParameters().signedInverseMomentum(), sqrt(updated.localError().matrix()(0,0)),
           //     updated.localPosition().x(), sqrt(updated.localError().matrix()(3,3)), updated.localPosition().y(),  sqrt(updated.localError().matrix()(4,4)),
           //     updated.localParameters().dxdz(), sqrt(updated.localError().matrix()(1,1)), updated.localParameters().dydz(), sqrt(updated.localError().matrix()(2,2)));
           ret.push( TrajectoryMeasurement(prop, updated, tms[i].recHit(), pair.second) ); 
           tsos = updated;
        } else {
           ret.push( TrajectoryMeasurement(prop, tms[i].recHit()) ); 
           tsos = prop;
        }
    }
    //printf("Reversed trajectory: "); printTraj(ret);
    return ret;
}

void HGCTracking::declareOutput(const std::string &label) 
{
    produces<std::vector<reco::TrackExtra>>(label);
    produces<std::vector<reco::Track>>(label);
    produces<std::vector<reco::Track>>(label+"Seed");
    produces<std::vector<reco::CaloCluster>>(label);
}

void HGCTracking::makeOutput(const std::vector<Trajectory> &trajs, edm::Event &out, const std::string &label) 
{
    std::unique_ptr<std::vector<reco::TrackExtra>> outTkEx(new std::vector<reco::TrackExtra>());
    std::unique_ptr<std::vector<reco::Track>> outTk(new std::vector<reco::Track>());
    std::unique_ptr<std::vector<reco::CaloCluster>> outCl(new std::vector<reco::CaloCluster>());
    std::unique_ptr<std::vector<reco::Track>> outTkSeed(new std::vector<reco::Track>());
    auto rTrackExtras = out.getRefBeforePut<reco::TrackExtraCollection>(label);

    for (const Trajectory &t : trajs) {
        reco::CaloCluster hits;
        bool first = true;
        float energy = 0;
        GlobalPoint pos; GlobalVector mom; int q = 1;
        GlobalPoint outpos; GlobalVector outmom; 
        for (const auto & tm : t.measurements()) {
            if (!tm.recHit()->isValid()) continue;
            if (first) {
                pos = tm.updatedState().globalPosition();
                mom = tm.updatedState().globalMomentum();
                q = tm.updatedState().charge();
                hits.setPosition(math::XYZPoint(pos.x(), pos.y(), pos.z()));
                first = false;
            } else if (tm.updatedState().isValid()) {
                outpos = tm.updatedState().globalPosition();
                outmom = tm.updatedState().globalMomentum();
            }
            if (typeid(*tm.recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
                for (auto & p : (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*tm.recHit())).objRef()->hitsAndFractions()) {
                    if (p.second > 0) hits.addHitAndFraction(p.first, p.second);
                    energy += p.first * p.second;
                }
            } else if (typeid(*tm.recHit()) == typeid(HGCTrackingClusteringRecHit)) {
                for (auto & hr : (dynamic_cast<const HGCTrackingClusteringRecHit&>(*tm.recHit())).objRefs()) {
                    hits.addHitAndFraction(hr->id(), 1.0);
                }
                energy += ((dynamic_cast<const HGCTrackingClusteringRecHit&>(*tm.recHit())).energy());
            } else {
                energy += ((dynamic_cast<const HGCTrackingRecHitFromHit&>(*tm.recHit())).energy());
                hits.addHitAndFraction(tm.recHit()->geographicalId(), 1.0);
            }
        }
        hits.setEnergy(energy);
        outCl->push_back(hits);

        reco::TrackExtra extra(reco::Track::Point(outpos.x(), outpos.y(), outpos.z()), reco::Track::Vector(outmom.x(), outmom.y(), outmom.z()), true,  
                               reco::Track::Point(pos.x(), pos.y(), pos.z()), reco::Track::Vector(mom.x(), mom.y(), mom.z()), true,
                               reco::Track::CovarianceMatrix(), 0,  reco::Track::CovarianceMatrix(), 0, t.direction());
        outTkEx->push_back(extra);

        reco::Track trk(t.chiSquared(), t.foundHits()*2, hits.position(),
                        reco::Track::Vector(mom.x(), mom.y(), mom.z()), q, reco::Track::CovarianceMatrix());
        for (const auto & tm : t.measurements()) {
            int subdet = 5, layer = 0;
            if (tm.recHit()->geographicalId().det() == DetId::Forward) {
                subdet = tm.recHit()->geographicalId().subdetId();
                layer = HGCalDetId(tm.recHit()->geographicalId()).layer();
            } else {
                layer = HcalDetId(tm.recHit()->geographicalId()).depth();
            }
            // convert subdet to something representable in a hitpattern
            switch (subdet) {
                case 3: subdet = PixelSubdetector::PixelEndcap; break;
                case 4: subdet = StripSubdetector::TID; break;
                case 5: subdet = StripSubdetector::TEC; break;
            }
            trk.appendTrackerHitPattern(subdet, layer/2, layer%2, tm.recHit()->type());
        }
        if (hasMCTruth_) {
            auto mresult = truthMatch(t);
            trk.setNLoops(mresult.empty() ? 0 : ceil(mresult.front().second));
            if (!mresult.empty() && mresult.front().first->eventId().event()==0 && mresult.front().first->eventId().bunchCrossing()==0) {
                trk.setQuality(reco::Track::confirmed);
            }
        }
        outTk->push_back(trk);

        const reco::Track & seedtk  = * (dynamic_cast<const TrajectorySeedFromTrack &>(t.seed())).track();
        outTkSeed->push_back(seedtk);
    }
    for (unsigned int i = 0, n = outTk->size(); i < n; ++i) {
        (*outTk)[i].setExtra(reco::TrackExtraRef(rTrackExtras, i));
    }
    out.put(std::move(outTkEx), label);
    out.put(std::move(outTk), label);
    out.put(std::move(outCl), label);
    out.put(std::move(outTkSeed), label+"Seed");
}

void HGCTracking::writeAllHitsByLayer(const HGCRecHitCollection &hits, edm::Event &out, const std::string &label, int layers) 
{
    std::unique_ptr<std::vector<reco::CaloCluster>> outCl(new std::vector<reco::CaloCluster>(2*layers));

    for (const HGCRecHit &hit : hits) {
        int zside, layer;
        if (hit.id().det() == DetId::Forward) {
            HGCalDetId parsed(hit.id());
            zside = parsed.zside(); layer = parsed.layer();
        } else {
            HcalDetId parsed(hit.id());
            zside = parsed.zside(); layer = parsed.depth();
        }
        (*outCl)[(zside>0)+2*(layer-1)].addHitAndFraction(hit.id(), 1.0);
    }
    for (int ilayer = 0; ilayer < layers; ++ilayer) {
       (*outCl)[2*ilayer  ].setEnergy(ilayer+1);
       (*outCl)[2*ilayer+1].setEnergy(ilayer+1.5);
       math::RhoEtaPhiVectorF pointP(1, +0.1*(ilayer+1), 0);
       math::RhoEtaPhiVectorF pointN(1, -0.1*(ilayer+1), 0);
       (*outCl)[2*ilayer  ].setPosition(math::XYZPoint(pointN.x(),pointN.y(),pointN.z()));
       (*outCl)[2*ilayer+1].setPosition(math::XYZPoint(pointP.x(),pointP.y(),pointP.z()));
    }

    out.put(std::move(outCl), label);
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( HGCTracking );
