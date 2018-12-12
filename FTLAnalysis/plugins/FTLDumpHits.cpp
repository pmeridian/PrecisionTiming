#ifndef _FTL_DUMP_HITS_
#define _FTL_DUMP_HITS_

#include "TMath.h"

#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDRingForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetRing.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "PrecisionTiming/FTLAnalysis/interface/FTLHitsTree.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "RecoMTD/TransientTrackingRecHit/interface/MTDTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderWithPropagator.h"
#include "RecoTracker/TransientTrackingRecHit/interface/Traj2TrackHits.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"

using namespace std;

float mmToCm=1./10.;

auto cmp = [](const unsigned one, const unsigned two) -> bool { return one < two; };

class FTLDumpHits : public edm::EDAnalyzer
{
public:
  explicit FTLDumpHits(const edm::ParameterSet& pSet);
  ~FTLDumpHits() {};
  
  //---utils
  
  //---methods
  virtual void beginJob() override {};
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
  virtual void endJob() override {};
  
  std::string PrintPosition(const GlobalPoint& gp);
  std::string PrintPosition(const LocalPoint& lp);
  
private:
  const MTDGeometry* mtdGeometry_;
  
  //---inputs
  edm::Handle<reco::GenParticleCollection> genParticlesHandle_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<std::vector<PSimHit> > simHitsHandle_;
  edm::EDGetTokenT<std::vector<PSimHit> > simHitsToken_;
  edm::Handle<FTLRecHitCollection> recHitsHandle_;
  edm::EDGetTokenT<FTLRecHitCollection> recHitsToken_;    
  edm::Handle<FTLClusterCollection> clustersHandle_;
  edm::EDGetTokenT<FTLClusterCollection> clustersToken_;    
  edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_;
  edm::Handle<edm::View<reco::Track> > tracksHandle_;
  edm::Handle<MTDTrackingDetSetVector> tHitsHandle_;
  edm::EDGetTokenT<MTDTrackingDetSetVector> tHitsToken_;    
  edm::EDGetTokenT<vector<SimVertex> >                 genVtxToken_;
  edm::Handle<vector<SimVertex> >                      genVtxHandle_;    
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  edm::Handle<reco::BeamSpot > bsHandle_; 
  edm::ESHandle<GlobalTrackingGeometry> gtg;
  std::unique_ptr<TrackTransformer> theTransformer;
  std::unique_ptr<MeasurementEstimator> theEstimator;
  edm::ESHandle<TransientTrackBuilder> builder;
  edm::ESHandle<TransientTrackingRecHitBuilder> hitbuilder;

  //---options
  BTLDetId::CrysLayout crysLayout_;
  double track_hit_DRMax_;
  double track_hit_distMax_;
  bool verbosity_;
  
  //---outputs
  FTLHitsTree outTree_;
  edm::Service<TFileService> fs_;  
  
};



FTLDumpHits::FTLDumpHits(const edm::ParameterSet& pSet):
  genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
  simHitsToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsTag"))),
  recHitsToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("recHitsTag"))),
  clustersToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersTag"))),
  tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  tHitsToken_(consumes<MTDTrackingDetSetVector>(pSet.getUntrackedParameter<edm::InputTag>("tHitsTag"))),
  genVtxToken_(consumes<vector<SimVertex> >(pSet.getUntrackedParameter<edm::InputTag>("genVtxTag"))),    
  bsToken_(consumes<reco::BeamSpot>(pSet.getUntrackedParameter<edm::InputTag>("bsTag"))),    
  crysLayout_((BTLDetId::CrysLayout)(pSet.getUntrackedParameter<int>("crysLayout"))),
  track_hit_DRMax_(pSet.getParameter<double>("track_hit_DRMax")),
  track_hit_distMax_(pSet.getParameter<double>("track_hit_distMax")),
  verbosity_(pSet.getParameter<bool>("verbosity"))
{
  outTree_ = FTLHitsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "FTLHits tree for FTL studies");
  theTransformer = std::make_unique<TrackTransformer>(pSet.getParameterSet("TrackTransformer"));
  float theMaxChi2=25.;
  float theNSigma=3.;
  theEstimator = std::make_unique<Chi2MeasurementEstimator>(theMaxChi2,theNSigma);
}



void FTLDumpHits::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
  outTree_.Reset();
  
  setup.get<GlobalTrackingGeometryRecord>().get(gtg);
  theTransformer->setServices(setup);

  //---get the MTD geometry
  edm::ESHandle<MTDGeometry> geoHandle;
  setup.get<MTDDigiGeometryRecord>().get(geoHandle);
  mtdGeometry_ = geoHandle.product();
  
  edm::ESHandle<MTDDetLayerGeometry> layerGeo;
  setup.get<MTDRecoGeometryRecord>().get(layerGeo);
  
  //--- get the B field
  edm::ESHandle<MagneticField> theField;
  setup.get<IdealMagneticFieldRecord>().get(theField);


  PropagationDirection dir(alongMomentum);
  SteppingHelixPropagator* propagator = new SteppingHelixPropagator(theField.product(),dir);
  propagator -> setMaterialMode(false);
  propagator -> setNoErrorPropagation(false);

  PropagatorWithMaterial* propagatorWMa= new PropagatorWithMaterial(anyDirection,0.13957018,theField.product(),1.6,false,0.1,true);
  
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  setup.get<TransientRecHitRecord>().get("MTDRecHitBuilder",hitbuilder);

  //--- load gen particles
  event.getByToken(genParticlesToken_, genParticlesHandle_);
  auto genParticles = *genParticlesHandle_.product();
  
  //---load sim hits
  event.getByToken(simHitsToken_, simHitsHandle_);
  auto simHits = *simHitsHandle_.product();
  
  //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
  event.getByToken(recHitsToken_, recHitsHandle_);
  auto recHits = FTLRecHitCollection();
  if(recHitsHandle_.isValid())
    recHits = *recHitsHandle_.product();

  event.getByToken(clustersToken_, clustersHandle_);
  // auto clusters = FTLClusterCollection();
  // if(clustersHandle_.isValid())
  auto clusters = *clustersHandle_.product();
  
  event.getByToken(tHitsToken_, tHitsHandle_);
  auto mtdTHits = *tHitsHandle_.product();

  event.getByToken(genVtxToken_, genVtxHandle_);    
  const SimVertex* genPV = NULL;
  if(genVtxHandle_.isValid())
    genPV = &(genVtxHandle_.product()->at(0));

  event.getByToken(bsToken_,bsHandle_);
  const auto& bs = *bsHandle_.product();

  //---load tracks
  event.getByToken(tracksToken_,tracksHandle_);
  auto tracks = *tracksHandle_.product();
  
  int nMods = (crysLayout_ == BTLDetId::CrysLayout::barzflat) ? 14 : 18;

  //---fill the tree - simHits
  outTree_.simHits_n = 0;
  for(auto simHit : simHits)
  {
    BTLDetId id = simHit.detUnitId();
    DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+nMods*(id.modType()-1),0,1);
    const auto& det = mtdGeometry_ -> idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
    
    double energy = simHit.energyLoss()*1000.;
    double time   = simHit.tof();
    
    int RR = id.mtdRR();
    int module = id.module();
    int modType = id.modType();
    int crystal = id.crystal();
    int ieta = id.ieta(crysLayout_);
    int iphi = id.iphi(crysLayout_);
    
    LocalPoint lp_entry(simHit.entryPoint().x()*mmToCm, simHit.entryPoint().y()*mmToCm, simHit.entryPoint().z()*mmToCm);
    LocalPoint lp_exit ( simHit.exitPoint().x()*mmToCm,  simHit.exitPoint().y()*mmToCm,  simHit.exitPoint().z()*mmToCm);
    GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,id.row(topo.nrows()),id.column(topo.nrows())));
    GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit,id.row(topo.nrows()),id.column(topo.nrows())));
    
    outTree_.simHits_n += 1;
    
    outTree_.simHits_energy->push_back(energy);
    outTree_.simHits_time->push_back(time);
    outTree_.simHits_rr->push_back(RR);
    outTree_.simHits_module->push_back(module);
    outTree_.simHits_modType->push_back(modType);
    outTree_.simHits_crystal->push_back(crystal);
    outTree_.simHits_ieta->push_back(ieta);
    outTree_.simHits_iphi->push_back(iphi);
    outTree_.simHits_entry_local_x->push_back(lp_entry.x());
    outTree_.simHits_entry_local_y->push_back(lp_entry.y());
    outTree_.simHits_entry_local_z->push_back(lp_entry.z());
    outTree_.simHits_entry_global_R->push_back(sqrt(gp_entry.perp2()));
    outTree_.simHits_exit_local_x->push_back(lp_exit.x());
    outTree_.simHits_exit_local_y->push_back(lp_exit.y());
    outTree_.simHits_exit_local_z->push_back(lp_exit.z());
    outTree_.simHits_exit_global_R->push_back(sqrt(gp_exit.perp2()));
  }
  
  
  
  //---fill the tree - recHits
  outTree_.recHits_n = 0;
  for(auto recHit : recHits)
  {
    BTLDetId id = recHit.id();
    DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+nMods*(id.modType()-1),0,1);
    const auto& det = mtdGeometry_ -> idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
    
    double energy = recHit.energy();
    double time   = recHit.time();
    
    MeasurementPoint mp(recHit.row()+0.5,recHit.column()+0.5);
    LocalPoint lp = topo.localPosition(mp);
    GlobalPoint gp = det->toGlobal(lp);
    
    int RR = id.mtdRR();
    int module = id.module();
    int modType = id.modType();
    int crystal = id.crystal();
    int ieta = id.ieta(crysLayout_);
    int iphi = id.iphi(crysLayout_);
    
    outTree_.recHits_n += 1;
    
    outTree_.recHits_energy->push_back(energy);
    outTree_.recHits_time->push_back(time);
    outTree_.recHits_rr->push_back(RR);
    outTree_.recHits_module->push_back(module);
    outTree_.recHits_modType->push_back(modType);
    outTree_.recHits_crystal->push_back(crystal);
    outTree_.recHits_ieta->push_back(ieta);
    outTree_.recHits_iphi->push_back(iphi);
    outTree_.recHits_local_x->push_back(lp.x());
    outTree_.recHits_local_y->push_back(lp.y());
    outTree_.recHits_local_z->push_back(lp.z());
    outTree_.recHits_global_R->push_back(sqrt(gp.perp2()));
  }


  //---fill the tree - recHits
  outTree_.clusters_n = 0;
  
  for (auto clusIt : clusters)
    {    
      DetId id = clusIt.detId();
      //    DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+nMods*(id.modType()-1),0,1);
      const auto& det = mtdGeometry_ -> idToDet(id);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      for ( auto cluster : clusIt)
	{
	  float energy = cluster.energy();
	  float time   = cluster.time();
	  float x=cluster.x();
	  float y=cluster.y();
	  int size=cluster.size();
	  int sizeX=cluster.sizeX();
	  int sizeY=cluster.sizeY();
	  float seed_energy=cluster.seed().energy;
	  float seed_time=cluster.seed().time;
	  int seed_x=cluster.seed().x;
	  int seed_y=cluster.seed().y;
	  
	  MeasurementPoint mp(cluster.x(),cluster.y());
	  LocalPoint lp = topo.localPosition(mp);
	  GlobalPoint gp = det->toGlobal(lp);
	  
	  MTDDetId mtdId(id);
	  int RR = 0;
	  int module = 0;
	  int modType = 0;
	  int crystal = 0;
	  int ieta = 0;
	  int iphi = 0;
	  
	  if ( mtdId.mtdSubDetector() == MTDDetId::BTL )
	    {
	      BTLDetId btlId(id);
	      RR = btlId.mtdRR();
	      module = btlId.module();
	      modType = btlId.modType();
	      crystal = btlId.crystal();
	      ieta = btlId.ieta(crysLayout_);
	      iphi = btlId.iphi(crysLayout_);
	    }
	  
	  outTree_.clusters_n += 1;    
	  outTree_.clusters_size->push_back(size);
	  outTree_.clusters_size_x->push_back(sizeX);
	  outTree_.clusters_size_y->push_back(sizeY);
	  outTree_.clusters_energy->push_back(energy);
	  outTree_.clusters_time->push_back(time);
	  outTree_.clusters_rr->push_back(RR);
	  outTree_.clusters_module->push_back(module);
	  outTree_.clusters_modType->push_back(modType);
	  outTree_.clusters_crystal->push_back(crystal);
	  outTree_.clusters_ieta->push_back(ieta);
	  outTree_.clusters_iphi->push_back(iphi);
	  outTree_.clusters_x->push_back(x);
	  outTree_.clusters_y->push_back(y);
	  outTree_.clusters_seed_energy->push_back(seed_energy);
	  outTree_.clusters_seed_time->push_back(seed_time);
	  outTree_.clusters_seed_x->push_back(seed_x);
	  outTree_.clusters_seed_y->push_back(seed_y);
	  outTree_.clusters_local_x->push_back(lp.x());
	  outTree_.clusters_local_y->push_back(lp.y());
	  outTree_.clusters_local_z->push_back(lp.z());
	  outTree_.clusters_global_R->push_back(sqrt(gp.perp2()));
	}
    }

  //--- fill the tree - tracks
  int idx=0;
  for(unsigned iTrack = 0; iTrack < tracks.size(); ++iTrack)
  {
    auto track = tracks.at(iTrack);
    // skip neutrals
    if( track.charge() == 0 ) continue;
    if( fabs(track.eta()) > 1.5 ) continue;
    if( track.pt() < 0.7 ) continue;
    
    // match with gen particles
    float DRMin = 999999;
    int genPdgId = 0;
    float genEta = -999.;
    float genPhi = -999.;
    float genPt = -999.;
    float genTime = -999.;
    for(auto& genPart : genParticles)
    {
      if( genPart.status() != 1 ) continue;
      if( genPart.charge() == 0 ) continue;
      
      float Deta = track.eta()-genPart.eta();
      float Dphi = deltaPhi(track.phi(),genPart.phi());
      float DR   = sqrt(Deta*Deta+Dphi*Dphi);
      
      if( DR < DRMin )
      {
        DRMin = DR;
        
        genPdgId = genPart.pdgId();
        genEta   = genPart.eta();
        genPhi   = genPart.phi();
        genPt    = genPart.pt();
      }
    }
    
    if (genPV)
      {
	genTime = genPV->position().t();
        outTree_.track_mcMatch_genX -> push_back(genPV->position().x());
        outTree_.track_mcMatch_genY -> push_back(genPV->position().y());
        outTree_.track_mcMatch_genZ -> push_back(genPV->position().z());
      }
    else
      {
	outTree_.track_mcMatch_genX -> push_back(-999.);
        outTree_.track_mcMatch_genY -> push_back(-999.);
        outTree_.track_mcMatch_genZ -> push_back(-999.);
      }
    outTree_.track_idx -> push_back(idx);
    outTree_.track_pt -> push_back(track.pt());
    outTree_.track_p -> push_back(track.p());
    outTree_.track_eta -> push_back(track.eta());
    outTree_.track_phi -> push_back(track.phi());
    outTree_.track_x -> push_back(track.vx());
    outTree_.track_y -> push_back(track.vy());
    outTree_.track_z -> push_back(track.vz());

    outTree_.track_energy -> push_back(sqrt(track.momentum().mag2()));
    outTree_.track_normalizedChi2 -> push_back(track.normalizedChi2());
    outTree_.track_numberOfValidHits -> push_back(track.numberOfValidHits());
    outTree_.track_numberOfLostHits -> push_back(track.numberOfLostHits());
    outTree_.track_isHighPurity -> push_back(track.quality(reco::TrackBase::TrackQuality::highPurity));
    outTree_.track_mcMatch_genPdgId -> push_back(genPdgId);
    outTree_.track_mcMatch_genPt -> push_back(genPt);
    outTree_.track_mcMatch_genEta -> push_back(genEta);
    outTree_.track_mcMatch_genPhi -> push_back(genPhi);
    outTree_.track_mcMatch_genTime -> push_back(genTime);
    outTree_.track_mcMatch_DR -> push_back(DRMin);
        
    if( verbosity_ ) std::cout << "*** track " << iTrack << " / " << tracks.size() << "   pt: " << track.pt() << "   eta: " << track.eta() << "   phi: " << track.phi() << std::endl;
    if( verbosity_ ) std::cout << "*** match with gen particle   DR: " << DRMin << "   gen pdgId: " << genPdgId << "   gen eta: " << genEta << "   gen phi: " << genPhi << "   genPt: " << genPt << std::endl;
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    outTree_.matchedSimHits_n->resize(idx+1);
    outTree_.matchedRecHits_n->resize(idx+1);
    outTree_.matchedClusters_n->resize(idx+1);
    
    outTree_.matchedSimHits_idx->resize(idx+1);
    outTree_.matchedSimHits_energy->resize(idx+1);
    outTree_.matchedSimHits_energyCorr->resize(idx+1);
    outTree_.matchedSimHits_time->resize(idx+1);
    outTree_.matchedSimHits_rr->resize(idx+1);
    outTree_.matchedSimHits_module->resize(idx+1);
    outTree_.matchedSimHits_modType->resize(idx+1);
    outTree_.matchedSimHits_crystal->resize(idx+1);
    outTree_.matchedSimHits_ieta->resize(idx+1);
    outTree_.matchedSimHits_iphi->resize(idx+1);
    outTree_.matchedSimHits_row->resize(idx+1);
    outTree_.matchedSimHits_col->resize(idx+1);
    outTree_.matchedSimHits_entry_local_x->resize(idx+1);
    outTree_.matchedSimHits_entry_local_y->resize(idx+1);
    outTree_.matchedSimHits_entry_local_z->resize(idx+1);
    outTree_.matchedSimHits_entry_global_R->resize(idx+1);
    outTree_.matchedSimHits_exit_local_x->resize(idx+1);
    outTree_.matchedSimHits_exit_local_y->resize(idx+1);
    outTree_.matchedSimHits_exit_local_z->resize(idx+1);
    outTree_.matchedSimHits_exit_global_R->resize(idx+1);
    outTree_.matchedSimHits_track_Deta->resize(idx+1);
    outTree_.matchedSimHits_track_Dphi->resize(idx+1);
    outTree_.matchedSimHits_track_DR->resize(idx+1);
    outTree_.matchedSimHits_track_Dz->resize(idx+1);
    outTree_.matchedSimHits_track_RDphi->resize(idx+1);
    outTree_.matchedSimHits_track_dist->resize(idx+1);
    
    outTree_.matchedRecHits_idx->resize(idx+1);
    outTree_.matchedRecHits_energy->resize(idx+1);
    outTree_.matchedRecHits_energyCorr->resize(idx+1);
    outTree_.matchedRecHits_time->resize(idx+1);
    outTree_.matchedRecHits_rr->resize(idx+1);
    outTree_.matchedRecHits_module->resize(idx+1);
    outTree_.matchedRecHits_modType->resize(idx+1);
    outTree_.matchedRecHits_crystal->resize(idx+1);
    outTree_.matchedRecHits_ieta->resize(idx+1);
    outTree_.matchedRecHits_iphi->resize(idx+1);
    outTree_.matchedRecHits_local_x->resize(idx+1);
    outTree_.matchedRecHits_local_y->resize(idx+1);
    outTree_.matchedRecHits_local_z->resize(idx+1);
    outTree_.matchedRecHits_global_R->resize(idx+1);
    outTree_.matchedRecHits_track_Deta->resize(idx+1);
    outTree_.matchedRecHits_track_Dphi->resize(idx+1);
    outTree_.matchedRecHits_track_DR->resize(idx+1);
    outTree_.matchedRecHits_track_Dz->resize(idx+1);
    outTree_.matchedRecHits_track_RDphi->resize(idx+1);
    outTree_.matchedRecHits_track_dist->resize(idx+1);
    outTree_.matchedRecHits_sietaieta->resize(idx+1);
    outTree_.matchedRecHits_siphiiphi->resize(idx+1);

    outTree_.matchedClusters_idx->resize(idx+1);
    outTree_.matchedClusters_energy->resize(idx+1);
    outTree_.matchedClusters_energyCorr->resize(idx+1);
    outTree_.matchedClusters_time->resize(idx+1);
    outTree_.matchedClusters_rr->resize(idx+1);
    outTree_.matchedClusters_module->resize(idx+1);
    outTree_.matchedClusters_modType->resize(idx+1);
    outTree_.matchedClusters_crystal->resize(idx+1);
    outTree_.matchedClusters_ieta->resize(idx+1);
    outTree_.matchedClusters_iphi->resize(idx+1);
    outTree_.matchedClusters_size->resize(idx+1);
    outTree_.matchedClusters_size_x->resize(idx+1);
    outTree_.matchedClusters_size_y->resize(idx+1);
    outTree_.matchedClusters_local_x->resize(idx+1);
    outTree_.matchedClusters_local_y->resize(idx+1);
    outTree_.matchedClusters_local_z->resize(idx+1);
    outTree_.matchedClusters_global_R->resize(idx+1);
    outTree_.matchedClusters_track_Deta->resize(idx+1);
    outTree_.matchedClusters_track_Dphi->resize(idx+1);
    outTree_.matchedClusters_track_DR->resize(idx+1);
    outTree_.matchedClusters_track_Dz->resize(idx+1);
    outTree_.matchedClusters_track_RDphi->resize(idx+1);
    outTree_.matchedClusters_track_dist->resize(idx+1);
        
    // track extrapolation
    const Surface::RotationType dummyRot;
    std::vector<float> cyl_R;
    if( (crysLayout_ == BTLDetId::CrysLayout::tile) )
      {
	cyl_R.push_back(117.450);
	cyl_R.push_back(117.456);
	cyl_R.push_back(117.473);
	cyl_R.push_back(117.501);
	cyl_R.push_back(117.540);
	cyl_R.push_back(117.590);
	cyl_R.push_back(117.653);
	cyl_R.push_back(117.723);
	cyl_R.push_back(117.810);
      }
    if( crysLayout_ == BTLDetId::CrysLayout::bar )
      {
	cyl_R.push_back(117.450);
	cyl_R.push_back(117.541);
	cyl_R.push_back(117.812);
      }
    if( crysLayout_ == BTLDetId::CrysLayout::barzflat )
      {
	cyl_R.push_back(117.450);
	cyl_R.push_back(117.450); cyl_R.push_back(117.451); cyl_R.push_back(117.453); cyl_R.push_back(117.456); cyl_R.push_back(117.459); cyl_R.push_back(117.462); cyl_R.push_back(117.467); cyl_R.push_back(117.473);
	cyl_R.push_back(117.479); cyl_R.push_back(117.485); cyl_R.push_back(117.493); cyl_R.push_back(117.500); cyl_R.push_back(117.510); cyl_R.push_back(117.519); cyl_R.push_back(117.530); cyl_R.push_back(117.540);
	cyl_R.push_back(117.552); cyl_R.push_back(117.565); cyl_R.push_back(117.578); cyl_R.push_back(117.591); cyl_R.push_back(117.606); cyl_R.push_back(117.621); cyl_R.push_back(117.637); cyl_R.push_back(117.654);
	cyl_R.push_back(117.671); cyl_R.push_back(117.689); cyl_R.push_back(117.708); cyl_R.push_back(117.727); cyl_R.push_back(117.747); cyl_R.push_back(117.768); cyl_R.push_back(117.789); cyl_R.push_back(117.812);
      }
    
    //study propagator and track refitting
    reco::TransientTrack ttrack(track,theField.product(),gtg);
    auto trajs = theTransformer->transform(track);
    auto thits = theTransformer->getTransientRecHits(ttrack);

    for( const auto& trj : trajs ) {
      std::cout << "original track chi2: " << trj.chiSquared() 
		<< " ndof: " << trj.ndof() << std::endl;
    }

    int cylIt = 0;
    std::map<int,GlobalPoint> gp_ext;

    for(auto val : cyl_R)
      {
	Cylinder::ConstCylinderPointer theTargetCylinder = Cylinder::build(val,Surface::PositionType(0.,0.,0.),dummyRot);
	
	std::pair<TrajectoryStateOnSurface,double> aTsosPath_outer;
	if (trajs.front().direction() == alongMomentum )
	  aTsosPath_outer=propagator->propagateWithPath(trajs.front().lastMeasurement().updatedState(),(*theTargetCylinder));
	else
	  aTsosPath_outer=propagator->propagateWithPath(trajs.front().firstMeasurement().updatedState(),(*theTargetCylinder));
	
	gp_ext[cylIt] = GlobalPoint(0.,0.,0.);
	if( aTsosPath_outer.first.isValid() )
	  {
	    GlobalPoint temp(aTsosPath_outer.first.globalPosition().x(),aTsosPath_outer.first.globalPosition().y(),aTsosPath_outer.first.globalPosition().z());
	    gp_ext[cylIt] = GlobalPoint(temp);
	    if( verbosity_ ) std::cout << "*** track extrapolation: " << PrintPosition(temp) << std::endl;
	    
	    if( cylIt == 0 )
	      {
		outTree_.track_eta_atBTL -> push_back(gp_ext[cylIt].eta());
		outTree_.track_phi_atBTL -> push_back(gp_ext[cylIt].phi());
	      }
	  }
	else
	  {
	    if( cylIt == 0 )
	      {
		outTree_.track_eta_atBTL -> push_back(-999.);
		outTree_.track_phi_atBTL -> push_back(-999.);
	      }
	  }
      
	++cylIt;
      }
    
    bool hasMTD = false;
    TransientTrackingRecHit::ConstRecHitContainer mtdhits;
    auto tTrack = builder->build(track);
    
    // try propagation to BTL layers and find compatible hits
    TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
    const vector<const DetLayer*>& layers = layerGeo->allBTLLayers();
    
    for (auto ilay = layers.begin(); ilay!=layers.end(); ++ilay) {
      const MTDTrayBarrelLayer* layer = (const MTDTrayBarrelLayer*) (*ilay);
      pair<bool, TrajectoryStateOnSurface> comp = layer->compatible(tsos,*propagatorWMa,*theEstimator);
      if ( !comp.first )
	continue;
      hasMTD = true;
      vector<DetLayer::DetWithState> compDets = layer->compatibleDets(tsos,*propagatorWMa,*theEstimator);
      for( const auto& detWithState : compDets ) 
	{
	  auto range = mtdTHits.equal_range(detWithState.first->geographicalId(),cmp);
	  for( auto detitr = range.first; detitr != range.second; ++detitr ) {
	    auto best = detitr->end();
	    double best_chi2 = std::numeric_limits<double>::max();
	    for( auto itr = detitr->begin(); itr != detitr->end(); ++itr ) 
	      {
		auto est =  theEstimator->estimate(detWithState.second,*itr);
		if( est.first && est.second < best_chi2 ) { // just take the best chi2
		  best = itr;
		  best_chi2 = est.second;
		}
	      }
	    if( best != detitr->end() ) {
	      mtdhits.push_back(hitbuilder->build(&*best));
	    }
	  }    
	}     
    } //end loop btl layers
    
    hasMTD = mtdhits.size();
    if (hasMTD)
      std::cout << "Found MTD matching hits " << mtdhits.size() << std::endl;
        
    //bool validpropagation = true;
    //    double pathlength = 0.;
    double timeAtBTL=-999.;
    double pathlength1 = 0.;
    double pathlength2 = 0.;
    double bsToFirst1 = 0.;
    double bsToFirst2 = 0.;
    double lastToMTD=0.;

    reco::Track refitTrack;
    math::XYZPoint  refit_pos( 0, 0. ,0. );
    math::XYZVector  refit_mom( 0, 0. ,0. );

    if (hasMTD)
      {
	bool outsideIn = false;

	if (!thits.empty()){
	  GlobalPoint first = gtg->idToDet(thits.front()->geographicalId())->position();
	  GlobalPoint last = gtg->idToDet(thits.back()->geographicalId())->position();
	  
	  // maybe perp2?
	  auto rFirst = first.mag2();
	  auto rLast  = last.mag2();
	  if(rFirst > rLast) outsideIn = true;
	}

	if (!outsideIn)
	  for( auto& ahit : mtdhits ) thits.push_back(ahit);
	else
	  {
	    std::reverse(mtdhits.begin(),mtdhits.end());
	    for( auto& ahit : thits ) mtdhits.push_back(ahit);
	    thits.swap(mtdhits);
	  }
	auto trajwithmtd = theTransformer->transform(ttrack,thits);
	std::cout << "refitting resulted in " << trajwithmtd.size() << " trajectories!" << std::endl;

	TSCBLBuilderWithPropagator tscblBuilder(*propagatorWMa);	  

	TrajectoryStateOnSurface stateForProjectionToBeamLineOnSurface = 
	  trajwithmtd.front().closestMeasurement(GlobalPoint(bs.x0(),bs.y0(),bs.z0())).updatedState();
	TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(*stateForProjectionToBeamLineOnSurface.freeState(),bs);

	TrajectoryStateOnSurface stateForProjectionToBeamLineOnSurface_ori = 
	  trajs.front().closestMeasurement(GlobalPoint(bs.x0(),bs.y0(),bs.z0())).updatedState();
	TrajectoryStateClosestToBeamLine tscbl_ori = tscblBuilder(*stateForProjectionToBeamLineOnSurface_ori.freeState(),bs);
	
	if (trajwithmtd.front().direction() == alongMomentum)
	  {
	    for (auto it=trajwithmtd.front().measurements().begin(); it!=trajwithmtd.front().measurements().end()-1; ++it) 
	      {
		const auto &propresult = propagatorWMa->propagateWithPath(it->updatedState(), (it+1)->updatedState().surface());
		double layerpathlength = std::abs(propresult.second);
		pathlength1 += layerpathlength;
	      }
	    bsToFirst1 = std::abs(propagatorWMa->propagateWithPath(tscbl.trackStateAtPCA(), trajwithmtd.front().firstMeasurement().updatedState().surface()).second);
	  }
	else
	  {
	    for (auto it=trajwithmtd.front().measurements().rbegin(); it!=trajwithmtd.front().measurements().rend()-1; ++it) 
	      {
		const auto &propresult = propagatorWMa->propagateWithPath(it->updatedState(), (it+1)->updatedState().surface());
		double layerpathlength = std::abs(propresult.second);
		pathlength1 += layerpathlength;
	      }
	    bsToFirst1= std::abs(propagatorWMa->propagateWithPath(tscbl.trackStateAtPCA(), trajwithmtd.front().lastMeasurement().updatedState().surface()).second);
	  }

	if (trajs.front().direction() == alongMomentum)
	  {
	    for (auto it=trajs.front().measurements().begin(); it!=trajs.front().measurements().end()-1; ++it)
	      {
		const auto &propresult = propagatorWMa->propagateWithPath(it->updatedState(), (it+1)->updatedState().surface());
		double layerpathlength = std::abs(propresult.second);
		pathlength2 += layerpathlength;
	      }
	    bsToFirst2 = std::abs(propagatorWMa->propagateWithPath(tscbl_ori.trackStateAtPCA(), trajs.front().firstMeasurement().updatedState().surface()).second);
	  }
	else
	  {
	    for (auto it=trajs.front().measurements().rbegin(); it!=trajs.front().measurements().rend()-1; ++it)
	      {
		const auto &propresult = propagatorWMa->propagateWithPath(it->updatedState(), (it+1)->updatedState().surface());
		double layerpathlength = std::abs(propresult.second);
		pathlength2 += layerpathlength;
	      }
	    bsToFirst2 = std::abs(propagatorWMa->propagateWithPath(tscbl_ori.trackStateAtPCA(), trajs.front().lastMeasurement().updatedState().surface()).second);
	  }
	
	for (auto it=trajwithmtd.front().measurements().begin(); it!=trajwithmtd.front().measurements().end(); ++it) {
	  if (it->recHit()->geographicalId().det() == DetId::Forward && ForwardSubdetector(it->recHit()->geographicalId().subdetId()) == FastTime)
	    {
	      const MTDTrackingRecHit *mtdhit = static_cast<const MTDTrackingRecHit*>(it->recHit()->hit());
	      // if (mtdhit->time()<timeAtBTL)
	      timeAtBTL = mtdhit->time();
	      std::cout << "timeAtBTL " << timeAtBTL << std::endl; 
	      if (trajs.front().direction() == alongMomentum )
		lastToMTD=std::abs(propagatorWMa->propagateWithPath(trajs.front().lastMeasurement().updatedState(), (it)->updatedState().surface()).second);
	      else
		lastToMTD=std::abs(propagatorWMa->propagateWithPath(trajs.front().firstMeasurement().updatedState(), (it)->updatedState().surface()).second);
	    }
	}
	
	std::cout << "bsToFirst1 " << bsToFirst1 << std::endl;
	std::cout << "bsToFirst2 " << bsToFirst2 << std::endl;
	std::cout << "lastToMTD "  << lastToMTD << std::endl;
	pathlength1 += bsToFirst1 ;
	pathlength2 += bsToFirst2 ;
	pathlength2 += lastToMTD ;
	std::cout << "pathLength1 " << pathlength1 << std::endl;
	std::cout << "pathLength2 " << pathlength2 << std::endl;

	GlobalPoint v = tscbl.trackStateAtPCA().position();
	refit_pos=math::XYZPoint( v.x(), v.y(), v.z() );
	GlobalVector p = tscbl.trackStateAtPCA().momentum();
	refit_mom=math::XYZVector( p.x(), p.y(), p.z() );
	int ndof = trajs.front().ndof();

	refitTrack = reco::Track(trajs.front().chiSquared(),
		     int(ndof),//FIXME fix weight() in TrackingRecHit
		     refit_pos, refit_mom, tscbl.trackStateAtPCA().charge(), 
		     tscbl.trackStateAtPCA().curvilinearError(),
		     track.algo(),reco::TrackBase::undefQuality,0,0,1.0,1.0);

      }

    outTree_.track_pathLength->push_back(pathlength1);
    outTree_.track_pathLength2->push_back(pathlength2);
    outTree_.track_time_atBTL-> push_back(timeAtBTL);

    if (hasMTD)
      {
	outTree_.track_mtdrefit_pt -> push_back(refitTrack.pt());
	outTree_.track_mtdrefit_p -> push_back(refitTrack.p());
	outTree_.track_mtdrefit_eta -> push_back(refitTrack.eta());
	outTree_.track_mtdrefit_phi -> push_back(refitTrack.phi());
	outTree_.track_mtdrefit_normalizedChi2 -> push_back(refitTrack.normalizedChi2());
	outTree_.track_mtdrefit_x -> push_back(refitTrack.vx());
	outTree_.track_mtdrefit_y -> push_back(refitTrack.vy());
	outTree_.track_mtdrefit_z -> push_back(refitTrack.vz());
      }
    else
      {
	outTree_.track_mtdrefit_pt -> push_back(0.);
	outTree_.track_mtdrefit_p -> push_back(0.);
	outTree_.track_mtdrefit_eta -> push_back(0.);
	outTree_.track_mtdrefit_phi -> push_back(0.);
	outTree_.track_mtdrefit_normalizedChi2 -> push_back(0.);
	outTree_.track_mtdrefit_x -> push_back(0.);
	outTree_.track_mtdrefit_y -> push_back(0.);
	outTree_.track_mtdrefit_z -> push_back(0.);
      }
    /*
    //---get compatible layers
    const vector<const DetLayer*>& layers = layerGeo -> allBTLLayers();
    
    GlobalTrajectoryParameters gtp(vtx_inner,vec_inner,track.charge(),theField.product());
    SteppingHelixPropagator prop(theField.product(),alongMomentum);
    
    float theMaxChi2 = 25.;
    float theNSigma = 3.;
    theEstimator = new Chi2MeasurementEstimator(theMaxChi2,theNSigma);
    
    if( verbosity_ ) std::cout << "*** get compatible layers" << std::endl;
    int it = 0;
    for(auto ilay = layers.begin(); ilay!=layers.end(); ++ilay)
    {
      const MTDTrayBarrelLayer* layer = (const MTDTrayBarrelLayer*)(*ilay);
      const BoundCylinder& cyl = layer -> specificSurface();
      TrajectoryStateOnSurface tsos(gtp,cyl);
      
      std::pair<bool,TrajectoryStateOnSurface> comp = layer -> compatible(tsos,prop,*theEstimator);
      
      if( verbosity_ ) std::cout << ">>> layer " << it << "   cylinder R: " << cyl.radius() << " cm   cylinder half length: " << cyl.bounds().length()/2. << " is compatible: " << comp.first << std::endl;
      
      std::vector<DetLayer::DetWithState> compDets = layer->compatibleDets(tsos,prop,*theEstimator);
      
      if( verbosity_ )
      {
        std::cout << ">>>>>> number of compatibleDets: " << compDets.size() << std::endl;
        
        int it2 = 0;
        for(auto compIt : compDets)
        {
          std::cout << ">>>>>>>>> compatibleDet " << it2 << ":   final state pos: " << PrintPosition(compIt.second.globalPosition())          << std::endl;
          std::cout << ">>>>>>>>> compatibleDet " << it2 << ":           det pos: " << PrintPosition(compIt.first->position())                << std::endl;
          std::cout << ">>>>>>>>> compatibleDet " << it2 << ":          distance: " << (tsos.globalPosition()-compIt.first->position()).mag() << std::endl;
        }
      }
      
      ++it;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;    
    */
    
    
    
    //---get associated simHits
    if( verbosity_ ) std::cout << "*** simHits - n tot: " << simHits.size() << std::endl;
    int simHitIt = 0;
    for(auto simHit : simHits)
    {
      BTLDetId id = simHit.detUnitId();
      DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+nMods*(id.modType()-1),0,1);
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      
      double energy = simHit.energyLoss()*1000.;
      double time   = simHit.tof();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = id.crystal();
      int ieta = id.ieta(crysLayout_);
      int iphi = id.iphi(crysLayout_);
      int row=id.row(topo.nrows());
      int col=id.column(topo.nrows());

      LocalPoint lp_entry(   simHit.entryPoint().x()*mmToCm,   simHit.entryPoint().y()*mmToCm,   simHit.entryPoint().z()*mmToCm);
      LocalPoint lp_mid  (simHit.localPosition().x()*mmToCm,simHit.localPosition().y()*mmToCm,simHit.localPosition().z()*mmToCm);
      LocalPoint lp_exit (    simHit.exitPoint().x()*mmToCm,    simHit.exitPoint().y()*mmToCm,    simHit.exitPoint().z()*mmToCm);
      GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,id.row(topo.nrows()),id.column(topo.nrows())));
      GlobalPoint gp_mid   = det->toGlobal(topo.pixelToModuleLocalPoint(lp_mid,id.row(topo.nrows()),id.column(topo.nrows())));
      GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit,id.row(topo.nrows()),id.column(topo.nrows())));
      
      float eta = gp_mid.eta();
      float phi = gp_mid.phi();
      
      GlobalPoint gp_track = gp_ext[abs(id.row(topo.nrows())-int(topo.nrows()/2))];
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp_mid.z()      : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp_mid-gp_track).mag()      : -999.;
      
      if( DR < 0.05 && DR > 0. )
      {
        if( verbosity_ )  std::cout << ">>> topology:   nRows: " << topo.nrows() << "   nColumns: " << topo.ncolumns() << "   pitchx: " << topo.pitch().first << "   pitchy: " << topo.pitch().second << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
                                   << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << id.row(topo.nrows()) << "   column: " << id.column(topo.nrows()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   entryPoint:   local: " << PrintPosition(lp_entry)        << "   global: " << PrintPosition(gp_entry) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":     midPoint:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid)   << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":    exitPoint:   local: " << PrintPosition(lp_exit)         << "   global: " << PrintPosition(gp_exit)  << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ": hit position:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid) << "   DR: " << DR << "   dist: " << (gp_track-gp_mid).mag() << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.matchedSimHits_n->at(idx) += 1;
      
      outTree_.matchedSimHits_idx->at(idx).push_back(idx);
      outTree_.matchedSimHits_energy->at(idx).push_back(energy);
      outTree_.matchedSimHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.matchedSimHits_time->at(idx).push_back(time);
      outTree_.matchedSimHits_rr->at(idx).push_back(RR);
      outTree_.matchedSimHits_module->at(idx).push_back(module);
      outTree_.matchedSimHits_modType->at(idx).push_back(modType);
      outTree_.matchedSimHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedSimHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedSimHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedSimHits_row->at(idx).push_back(row);
      outTree_.matchedSimHits_col->at(idx).push_back(col);
      outTree_.matchedSimHits_entry_local_x->at(idx).push_back(lp_entry.x());
      outTree_.matchedSimHits_entry_local_y->at(idx).push_back(lp_entry.y());
      outTree_.matchedSimHits_entry_local_z->at(idx).push_back(lp_entry.z());
      outTree_.matchedSimHits_entry_global_R->at(idx).push_back(sqrt(gp_entry.perp2()));
      outTree_.matchedSimHits_exit_local_x->at(idx).push_back(lp_exit.x());
      outTree_.matchedSimHits_exit_local_y->at(idx).push_back(lp_exit.y());
      outTree_.matchedSimHits_exit_local_z->at(idx).push_back(lp_exit.z());
      outTree_.matchedSimHits_exit_global_R->at(idx).push_back(sqrt(gp_exit.perp2()));
      outTree_.matchedSimHits_track_Deta->at(idx).push_back(Deta);
      outTree_.matchedSimHits_track_Dphi->at(idx).push_back(Dphi);
      outTree_.matchedSimHits_track_DR->at(idx).push_back(DR);
      outTree_.matchedSimHits_track_Dz->at(idx).push_back(Dz);
      outTree_.matchedSimHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedSimHits_track_dist->at(idx).push_back(dist);
      
      ++simHitIt;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    
    //---find associated recHits
    float sieie=0, sipip=0;
    float ss_hit_count=0;
    int recHitIt = 0;
    if( verbosity_ ) std::cout << "*** recHits - tot: " << recHits.size() << std::endl;
    for(auto recHit : recHits)
    {
      BTLDetId id = recHit.id();
      DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+nMods*(id.modType()-1),0,1);
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      
      double energy = recHit.energy();
      double time   = recHit.time();
      
      MeasurementPoint mp(id.row(topo.nrows())+0.5,id.column(topo.nrows())+0.5);
      LocalPoint lp = topo.localPosition(mp);
      GlobalPoint gp = det->toGlobal(lp);

      float eta = gp.eta();
      float phi = gp.phi();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = id.crystal();
      int ieta = id.ieta(crysLayout_);
      int iphi = id.iphi(crysLayout_);
      
      GlobalPoint gp_track = gp_ext[abs(id.row(topo.nrows())-int(topo.nrows()/2))];
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
      
      // if( DR < 0.2 && DR > 0. )
      {
        if( verbosity_ )  std::cout << ">>> topology:   nRows: " << topo.nrows() << "   nColumns: " << topo.ncolumns() << "   pitchx: " << topo.pitch().first << "   pitchy: " << topo.pitch().second << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
                                   << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << recHit.row() << "- " << id.row(topo.nrows()) << "   column: " << recHit.column() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << "   DR: " << DR << "   dist: " << (gp_track-det->position()).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": hit position:   local: " << PrintPosition(lp) << "   global: " << PrintPosition(gp) << "DR: " << DR << "   dist: " << (gp_track-gp).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": extrapolated track " << abs(recHit.row()-int(topo.nrows()/2)) << " : " << PrintPosition(gp_track) << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.matchedRecHits_n->at(idx) += 1;
      
      outTree_.matchedRecHits_idx->at(idx).push_back(idx);
      outTree_.matchedRecHits_energy->at(idx).push_back(energy);
      outTree_.matchedRecHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.matchedRecHits_time->at(idx).push_back(time);
      outTree_.matchedRecHits_rr->at(idx).push_back(RR);
      outTree_.matchedRecHits_module->at(idx).push_back(module);
      outTree_.matchedRecHits_modType->at(idx).push_back(modType);
      outTree_.matchedRecHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedRecHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedRecHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedRecHits_local_x->at(idx).push_back(lp.x());
      outTree_.matchedRecHits_local_y->at(idx).push_back(lp.y());
      outTree_.matchedRecHits_local_z->at(idx).push_back(lp.z());
      outTree_.matchedRecHits_global_R->at(idx).push_back(sqrt(gp.perp2()));
      outTree_.matchedRecHits_track_Deta->at(idx).push_back(Deta);
      outTree_.matchedRecHits_track_Dphi->at(idx).push_back(Dphi);
      outTree_.matchedRecHits_track_Dz->at(idx).push_back(Dz);
      outTree_.matchedRecHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedRecHits_track_DR->at(idx).push_back(DR);
      outTree_.matchedRecHits_track_dist->at(idx).push_back(dist);
      
      if(recHit.energy() > 0.5)
      {
        sieie += energy*pow(eta-gp_track.eta(),2);
        sipip += energy*pow(phi-gp_track.phi(),2);
        ss_hit_count += energy;
      }
      
      ++recHitIt;
    }
  
    //---find associated recHits
    int clusterIt = 0;

    for (auto clusIt : clusters)
      {    
	DetId id = clusIt.detId();
	//    DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+nMods*(id.modType()-1),0,1);
	const auto& det = mtdGeometry_ -> idToDet(id);
	const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
	const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
	for ( auto cluster : clusIt)
	  {
	    
	    MTDDetId mtdId(id);
	    int RR = 0;
	    int module = 0;
	    int modType = 0;
	    int crystal = 0;
	    int ieta = 0;
	    int iphi = 0;
	    
	    if ( mtdId.mtdSubDetector() == MTDDetId::BTL )
	      {
		BTLDetId btlId(id);
		RR = btlId.mtdRR();
		module = btlId.module();
		modType = btlId.modType();
		crystal = btlId.crystal();
		ieta = btlId.ieta(crysLayout_);
		iphi = btlId.iphi(crysLayout_);
	      }
	    
	    double energy = cluster.energy();
	    double time   = cluster.time();
	    int size=cluster.size();
	    int sizeX=cluster.sizeX();
	    int sizeY=cluster.sizeY();
	    
	    MeasurementPoint mp(cluster.x(),cluster.y());
	    LocalPoint lp = topo.localPosition(mp);
	    GlobalPoint gp = det->toGlobal(lp);
	    
	    float eta = gp.eta();
	    float phi = gp.phi();
	    
	    GlobalPoint gp_track = gp_ext[abs(cluster.seed().x-int(topo.nrows()/2))];
	    
	    float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
	    float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
	    float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
	    float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
	    float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
	    float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
	    
	    // if( DR < 0.2 && DR > 0. )
	    //{
	    // if( verbosity_ )  std::cout << ">>> topology:   nRows: " << topo.nrows() << "   nColumns: " << topo.ncolumns() << "   pitchx: " << topo.pitch().first << "   pitchy: " << topo.pitch().second << std::endl;
	    // if( verbosity_ ) std::cout << ">>> " << clusterIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
	      // 				 << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
	      // 				 << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << cluster.row() << "- " << id.row(topo.nrows()) << "   column: " << cluster.column() << std::endl;
	      // if( verbosity_ ) std::cout << ">>> " << clusterIt << ":  detPosition:  global: " << PrintPosition(det->position()) << "   DR: " << DR << "   dist: " << (gp_track-det->position()).mag() << std::endl;
	      // if( verbosity_ ) std::cout << ">>> " << clusterIt << ": hit position:   local: " << PrintPosition(lp) << "   global: " << PrintPosition(gp) << "DR: " << DR << "   dist: " << (gp_track-gp).mag() << std::endl;
	      // if( verbosity_ ) std::cout << ">>> " << clusterIt << ": extrapolated track " << abs(cluster.row()-int(topo.nrows()/2)) << " : " << PrintPosition(gp_track) << std::endl;
	    //}
	    
	    
	    if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
	    if( gp_track.mag() <= 0. ) continue;
	    
	    outTree_.matchedClusters_n->at(idx) += 1;
	    
	    outTree_.matchedClusters_idx->at(idx).push_back(idx);
	    outTree_.matchedClusters_energy->at(idx).push_back(energy);
	    outTree_.matchedClusters_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
	    outTree_.matchedClusters_time->at(idx).push_back(time);
	    outTree_.matchedClusters_rr->at(idx).push_back(RR);
	    outTree_.matchedClusters_module->at(idx).push_back(module);
	    outTree_.matchedClusters_modType->at(idx).push_back(modType);
	    outTree_.matchedClusters_crystal->at(idx).push_back(crystal);
	    outTree_.matchedClusters_ieta->at(idx).push_back(ieta);
	    outTree_.matchedClusters_iphi->at(idx).push_back(iphi);
	    outTree_.matchedClusters_size->at(idx).push_back(size);
	    outTree_.matchedClusters_size_x->at(idx).push_back(sizeX);
	    outTree_.matchedClusters_size_y->at(idx).push_back(sizeY);
	    outTree_.matchedClusters_local_x->at(idx).push_back(lp.x());
	    outTree_.matchedClusters_local_y->at(idx).push_back(lp.y());
	    outTree_.matchedClusters_local_z->at(idx).push_back(lp.z());
	    outTree_.matchedClusters_global_R->at(idx).push_back(sqrt(gp.perp2()));
	    outTree_.matchedClusters_track_Deta->at(idx).push_back(Deta);
	    outTree_.matchedClusters_track_Dphi->at(idx).push_back(Dphi);
	    outTree_.matchedClusters_track_Dz->at(idx).push_back(Dz);
	    outTree_.matchedClusters_track_RDphi->at(idx).push_back(RDphi);
	    outTree_.matchedClusters_track_DR->at(idx).push_back(DR);
	    outTree_.matchedClusters_track_dist->at(idx).push_back(dist);
	    
	    ++clusterIt;
	  }
      }
  
  
    if( verbosity_ ) std::cout << "---\n\n\n" << std::endl;
    
    outTree_.matchedRecHits_sietaieta->at(idx).push_back( ss_hit_count>0 ? sqrt(sieie)/ss_hit_count : -999. );
    outTree_.matchedRecHits_siphiiphi->at(idx).push_back( ss_hit_count>0 ? sqrt(sipip)/ss_hit_count : -999. );
  
    ++idx;
  }
  
  outTree_.GetTTreePtr()->Fill();
}


std::string FTLDumpHits::PrintPosition(const GlobalPoint& gp)
{
  std::stringstream output;
  
  output << "(";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.x() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.y() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.z();
  output << ") cm";
  
  output << "   R: " << std::fixed << std::setprecision(3) << std::setw(7) << gp.perp();
  output << " cm";
  
  output << "   eta: " << std::setprecision(3) << std::setw(6) << gp.eta(); 
  output << "   phi: " << std::setprecision(3) << std::setw(6) << gp.phi();
  
  return output.str();
}

std::string FTLDumpHits::PrintPosition(const LocalPoint& lp)
{
  std::stringstream output;
  
  output << "(";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.x() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.y() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.z();
  output << ") cm";
  
  return output.str();
}
DEFINE_FWK_MODULE(FTLDumpHits);

#endif
