import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('datasets',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Input dataset(s)")
options.register('eosdirs',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "files location(s) on EOS")
options.register('eosusdirs',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "files location(s) on EOS-US")
options.register('pattern',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "pattern of file names to be processed")
options.register('antipattern',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "pattern of file names not to be processed")
options.register('output',
                 'test.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "output file name")
options.register('useparent',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Load data from parent datasets")
options.register('debug',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.register('crysLayout',
                 'barzflat',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "crystal layout (tile, barphi, barzflat)")
options.register('nThreads',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "# threads")
options.maxEvents = -1
options.parseArguments()

#from Configuration.StandardSequences.Eras import eras
#if 'tile' in options.crysLayout:
#    myera=eras.Phase2_timing_layer_tile
#if 'barphi' in options.crysLayout:
#    myera=eras.Phase2_timing_layer_bar
#if 'barzflat' in options.crysLayout:
#    myera=eras.Phase2C4_timing_layer_bar
from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('MTD4DVertexingStudies', Phase2C9)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(options.nThreads),
    numberOfStreams=cms.untracked.uint32(0),
    wantSummary = cms.untracked.bool(True)
    )

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# import of standard configurations
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Geometry
# if 'tile' in options.crysLayout:
#     process.load('Configuration.Geometry.GeometryExtended2023D24Reco_cff')
#     process.load('Configuration.Geometry.GeometryExtended2023D24_cff')
# if 'barphi' in options.crysLayout:
#     process.load('Configuration.Geometry.GeometryExtended2023D25Reco_cff')
#     process.load('Configuration.Geometry.GeometryExtended2023D25_cff')
# if 'barzflat' in options.crysLayout:
#     process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
#     process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
    
# process.load('Configuration.StandardSequences.SimIdeal_cff')

# process.load("Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cfi")
# process.load("Geometry.MTDNumberingBuilder.mtdTopology_cfi")
# process.load("Geometry.MTDGeometryBuilder.mtdGeometry_cfi")
# process.load("Geometry.MTDGeometryBuilder.mtdParameters_cfi")
# process.mtdGeometry.applyAlignment = cms.bool(False)

files = []
secondary_files = []
for dataset in options.datasets:
    print('>> Creating list of files from: \n'+dataset)
    for instance in ['global', 'phys03']:
        query = "-query='file dataset="+dataset+" instance=prod/"+instance+"'"
        if options.debug:
            print(query)
        lsCmd = subprocess.Popen(['dasgoclient '+query+' -limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        str_files, err = lsCmd.communicate()
        files.extend(['root://cms-xrd-global.cern.ch/'+ifile for ifile in str_files.split("\n")])
        files = [k for k in files if options.pattern in k]
        files.pop()
        if options.useparent:
            print('>> Creating list of secondary files from: \n'+dataset)
            for file in files:
                query = "-query='parent file="+file[len('root://cms-xrd-global.cern.ch/'):]+" instance=prod/"+instance+"'"
                if options.debug:
                    print(query)
                lsCmd = subprocess.Popen(['dasgoclient '+query+' -limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                str_files, err = lsCmd.communicate()
                secondary_files.extend(['root://cms-xrd-global.cern.ch/'+ifile for ifile in str_files.split("\n")])        
                secondary_files.pop()

for eosdir in options.eosdirs:
    if eosdir[-1] != '/':
        eosdir += '/'
    print('>> Creating list of files from: \n'+eosdir)
    
    command = '/bin/find '+eosdir+' -type f | grep root | grep -v failed'
    str_files = subprocess.check_output(command,shell=True).splitlines()
    print str_files
    files.extend(['file:'+ifile for ifile in str_files])

for eosusdir in options.eosusdirs:
    if eosusdir[-1] != '/':
        eosusdir += '/'
    print('>> Creating list of files from: \n'+eosusdir)
    
    command = 'eos root://cmseos.fnal.gov ls '+eosusdir+' | grep root | grep -v failed | grep '+options.pattern+' | grep -v '+options.antipattern
    str_files = subprocess.check_output(command,shell=True).splitlines()
    print str_files
    files.extend(['root://cms-xrd-global.cern.ch/'+eosusdir+'/'+ifile for ifile in str_files])
        
if options.debug:
    for ifile in files:
        print(ifile)        
    for ifile in secondary_files:
        print(ifile)        

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(files),
    secondaryFileNames = cms.untracked.vstring(secondary_files)
    )
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# Recreate trackingParticles
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.MixingModule.trackingTruthProducerSelection_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import *
process.quickTrackAssociatorByHits = quickTrackAssociatorByHits
process.trackingParticleRecoTrackAsssociation = cms.EDProducer("TrackAssociatorEDProducer",
                                                               associator = cms.InputTag('quickTrackAssociatorByHits'),
                                                               label_tp = cms.InputTag("mix","MergedTrackTruth"),
                                                               label_tr = cms.InputTag("generalTracks"),
                                                               ignoremissingtrackcollection=cms.untracked.bool(False)
)
process.mix.playback = cms.untracked.bool(True)
process.mix.digitizers = cms.PSet(
         mergedtruth = cms.PSet(process.trackingParticles)
    )
for a in process.aliases: delattr(process, a)

# Analyzer
process.load('PrecisionTiming.FTLAnalysis.MTD4DVertexingAnalyzer_cfi')
MTDDumper = process.MTD4DVertexingAnalyzer
MTDDumper.trackingParticlesTag = cms.untracked.InputTag("mix", "MergedTrackTruth", "MTD4DVertexingStudies")
MTDDumper.trackAndTrackingParticlesAssociatorMapTag = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation", "", "MTD4DVertexingStudies")

# Output TFile
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.output)
    )

process.runseq = cms.Sequence(process.mix*process.tpClusterProducer*process.quickTrackAssociatorByHits*process.trackingParticleRecoTrackAsssociation+MTDDumper)
process.path = cms.Path(process.runseq)
process.schedule = cms.Schedule(process.path)
