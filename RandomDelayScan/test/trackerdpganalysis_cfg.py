import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser                                                                                                                                        
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('delayStep',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,
                  'integer numer to identify the delay xml for each partition');

options.register ('eventToSkip',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,
                  'interget number to indicate how many events to skip')

options.register ('secondaryFiles',[],VarParsing.multiplicity.list,VarParsing.varType.string,
                  'list of secondary files')

options.register ('ouputFileName',"trackerDPG.root",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  'name of the outtput root file')

options.register ('jsonFile',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  'json file to apply in case one wants to ....')

options.register ('isRawFile',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                  'when running as input a raw file instead of a standard FEVT')

options.register ('isDatFile',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                  'when running as input a dat file file instead of a standard FEVT')

options.register ('dropAnalyzerDumpEDM',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                  'dump all the collections produced by the config')

options.register ('inputDirectory',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  'directory where the xml files with the delays are stored')


options.parseArguments()

### start cmssw job
import FWCore.PythonUtilities.LumiList as LumiList
from Configuration.StandardSequences.Eras import eras

process = cms.Process("clusterAnalysis",eras.Run2_2018)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
readFiles = cms.untracked.vstring(options.inputFiles)
secFiles  = cms.untracked.vstring(options.secondaryFiles) 

if not options.isDatFile:
    process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
    process.source.skipEvents = cms.untracked.uint32(options.eventToSkip)
else:
    process.source = cms.Source("NewEventStreamFileReader",
                                fileNames = readFiles)
    process.source.skipEvents = cms.untracked.uint32(options.eventToSkip)


if options.jsonFile != "": ### to be checked / created by hand when the runs are take
    process.source.lumisToProcess = LumiList.LumiList(filename = options.jsonFile).getVLuminosityBlockRange()

# Conditions (Global Tag is used here):
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Express_v7', '')

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(4),
    numberOfStreams = cms.untracked.uint32(4))


#Geometry and magnetic field to be loaded
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometryDB_cfi")
process.load("TrackingTools.RecoGeometry.RecoGeometries_cff")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

### to select a specific set of triggers --> if one runs on the streamer files at T0 is no needed
process.skimming = cms.EDFilter("PhysDecl",
  applyfilter = cms.untracked.bool(False),
  debugOn = cms.untracked.bool(False),
  HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
)

process.hltfiter = cms.EDFilter("HLTHighLevel",
  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
  HLTPaths = cms.vstring("HLT_ZeroBias*"),
  throw = cms.bool(True),
  andOr = cms.bool(True) ## logical OR between trigger bits                               
)


### output file definition
if not options.dropAnalyzerDumpEDM:
    process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(options.ouputFileName)
                                       )
else:
    process.out = cms.OutputModule("PoolOutputModule",    
                                   fileName = cms.untracked.string(options.ouputFileName),
                                   outputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_*_*_*clusterAnalysis*',            
            ))
    
    process.output = cms.EndPath(process.out)


### when one runs on a RECO file and wants to re-fit tracks from already existing ones
if not options.isRawFile and not options.isDatFile: 

    process.load('RecoTracker.DeDx.dedxEstimators_cff')
    process.load("RecoTracker.TrackProducer.TrackRefitter_cfi")
    process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
    process.TrackRefitter.NavigationSchool = ''
    process.TrackRefitter.Fitter = 'FlexibleKFFittingSmoother'
    process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True
    process.generalTracksFromRefit = process.TrackRefitter.clone(
        src = cms.InputTag("generalTracks")
        )
    
    process.doAlldEdXEstimators += process.dedxMedian
    process.dedxTruncated40.tracks = "generalTracksFromRefit"
    process.dedxHarmonic2.tracks   = "generalTracksFromRefit"
    process.dedxMedian.tracks      = "generalTracksFromRefit"
    process.dedxHitInfo.tracks     = "generalTracksFromRefit"
    process.refit = cms.Sequence(process.MeasurementTrackerEvent*process.generalTracksFromRefit*process.doAlldEdXEstimators)

### run the tracking step on top of raw files
if options.isRawFile or options.isDatFile:

    process.load('Configuration.StandardSequences.RawToDigi_cff')    
    process.load('Configuration.StandardSequences.Reconstruction_cff')
    process.load('RecoTracker.DeDx.dedxEstimators_cff')

    ## save trajectories everywhere since they are used by the TrackerDpgAnalyzer
    process.generalTracks.copyTrajectories = cms.untracked.bool(True);
    process.mergedDuplicateTracks.TrajectoryInEvent = cms.bool(True);
    process.preDuplicateMergingGeneralTracks.copyTrajectories = cms.untracked.bool(True);
    process.earlyGeneralTracks.copyTrajectories = cms.untracked.bool(True);
    process.initialStepTracks.TrajectoryInEvent = cms.bool(True);
    process.initialStepTracksPreSplitting.TrajectoryInEvent = cms.bool(True);
    process.jetCoreRegionalStepTracks.TrajectoryInEvent = cms.bool(True);
    process.lowPtTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.highPtTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.lowPtQuadStepTracks.TrajectoryInEvent = cms.bool(True);
    process.pixelPairStepTracks.TrajectoryInEvent = cms.bool(True);
    process.detachedTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.mixedTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.pixelLessStepTracks.TrajectoryInEvent = cms.bool(True);
    process.tobTecStepTracks.TrajectoryInEvent = cms.bool(True);
    process.muonSeededTracksInOut.TrajectoryInEvent = cms.bool(True);
    process.muonSeededTracksOutIn.TrajectoryInEvent = cms.bool(True);
    process.detachedQuadStepTracks.TrajectoryInEvent = cms.bool(True);

### Main Analyzer
process.analysis = cms.EDAnalyzer('TrackerDpgAnalysis',
                                  ClustersLabel = cms.InputTag("siStripClusters"),
                                  PixelClustersLabel = cms.InputTag("siPixelClusters"),
                                  TracksLabel   = cms.VInputTag( cms.InputTag("generalTracksFromRefit")),
                                  vertexLabel   = cms.InputTag('offlinePrimaryVertices'),
                                  pixelVertexLabel = cms.InputTag('pixelVertices'),
                                  beamSpotLabel = cms.InputTag('offlineBeamSpot'),
                                  DeDx1Label    = cms.InputTag('dedxHarmonic2'),
                                  DeDx2Label    = cms.InputTag('dedxTruncated40'),
                                  DeDx3Label    = cms.InputTag('dedxMedian'),
                                  L1Label       = cms.InputTag('gtDigis'),
                                  HLTLabel      = cms.InputTag("TriggerResults::HLT"),
                                  InitalCounter = cms.uint32(1),
                                  keepOntrackClusters  = cms.untracked.bool(True),
                                  keepOfftrackClusters = cms.untracked.bool(False),
                                  keepPixelClusters    = cms.untracked.bool(False),
                                  keepPixelVertices    = cms.untracked.bool(False),
                                  keepMissingHits      = cms.untracked.bool(False),
                                  keepTracks           = cms.untracked.bool(True),
                                  keepVertices         = cms.untracked.bool(True),
                                  keepEvents           = cms.untracked.bool(True),
                                  DelayFileNames = cms.untracked.vstring(
        "TI_27-JAN-2010_2_delayStep_"+str(options.delayStep)+".xml",
        "TO_30-JUN-2009_1_delayStep_"+str(options.delayStep)+".xml",
        "TP_09-JUN-2009_1_delayStep_"+str(options.delayStep)+".xml",
        "TM_09-JUN-2009_1_delayStep_"+str(options.delayStep)+".xml",
        ),
                                  PSUFileName =  cms.untracked.string("PSUmapping.csv")
                                  )

if  options.isRawFile or options.isDatFile: 
    process.analysis.TracksLabel = cms.VInputTag(cms.InputTag("generalTracks"))


if options.inputDirectory != "":
    for string in range(len(process.analysis.DelayFileNames)) :
        process.analysis.DelayFileNames[string] = options.inputDirectory+"/"+process.analysis.DelayFileNames[string];
    process.analysis.PSUFileName = options.inputDirectory+"/PSUmapping.csv";

### Define the Path to be executed
if options.isRawFile or options.isDatFile:

    process.doAlldEdXEstimators += process.dedxMedian

    if not options.dropAnalyzerDumpEDM:
        process.p = cms.Path(
                        process.RawToDigi*
                        process.skimming*
                        process.hltfiter*  ## HLT skim                                                                                                                                            
                        process.reconstruction_trackingOnly* ## local and gloabl reco                                                                                                              
                        process.doAlldEdXEstimators*
                        process.analysis
                        )

    else: ## drop the analyzer
        process.p = cms.Path(
            process.RawToDigi* 
            process.skimming*
            process.hltfiter*  ## HLT skim
            process.reconstruction_trackingOnly* ## local and gloabl reco
            process.doAlldEdXEstimators*
            process.dedxMedian)

else: ## in this case one just need to re-fit the tracks
    if not options.dropAnalyzerDumpEDM:
        process.p = cms.Path(process.skimming*
                             process.hltfiter*  ## HLT skim
                             process.refit*
                             process.analysis)
    else:
        process.p = cms.Path(process.skimming*
                             process.hltfiter*  ## HLT skim
                             process.refit)


process.schedule = cms.Schedule(process.p);

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
