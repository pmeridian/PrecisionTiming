import FWCore.ParameterSet.Config as cms

FTLDumpHits = cms.EDAnalyzer('FTLDumpHits',
                             genParticlesTag = cms.untracked.InputTag("genParticles"),
                             simHitsTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsBarrel"),
                             recHitsTag = cms.untracked.InputTag("mtdRecHits:FTLBarrel"),
                             clustersTag = cms.untracked.InputTag("mtdClusters:FTLBarrel"),
                             tracksTag = cms.untracked.InputTag("generalTracks"),
                             tHitsTag = cms.untracked.InputTag("mtdTrackingRecHits"),
                             genVtxTag = cms.untracked.InputTag("g4SimHits"),
                             bsTag = cms.untracked.InputTag("offlineBeamSpot"),
                             crysLayout = cms.untracked.int32(0),
                             track_hit_DRMax = cms.double(0.05),
                             track_hit_distMax = cms.double(99999.),
                             treeName = cms.untracked.string("DumpHits"),
                                 TrackTransformer = cms.PSet(
        DoPredictionsOnly = cms.bool(False),
        Fitter = cms.string('KFFitterForRefitInsideOut'),
        #TrackerRecHitBuilder = cms.string('WithTrackAngleAndTemplate'),
        TrackerRecHitBuilder = cms.string('WithTrackAngle'),
        Smoother = cms.string('KFSmootherForRefitInsideOut'),
        MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
        MTDRecHitBuilder = cms.string('MTDRecHitBuilder'),
        RefitDirection = cms.string('alongMomentum'),
        RefitRPCHits = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyRKOpposite')
        ),
                             verbosity = cms.bool(True)
                             )
