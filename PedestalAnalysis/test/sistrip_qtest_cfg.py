import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
#prepare options

options = VarParsing.VarParsing("analysis")

options.register ('globalTag',
                  "DONOTEXIST",
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "GlobalTag")

options.register ('runNumber',
                  1,
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Run Number")

options.parseArguments()

process = cms.Process("SiStripQualityStatJob")

### to store output text file from MessageLogger
process.MessageLogger = cms.Service("MessageLogger",
                                    destinations = cms.untracked.vstring('cout','cerr','SiStripQualityStatSummary'), #Reader, cout
                                    categories = cms.untracked.vstring('SiStripQualityStatistics'),
                                    cerr = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')
                                                              ),
                                    cout = cms.untracked.PSet(threshold = cms.untracked.string('WARNING'),
                                                              default = cms.untracked.PSet(limit=cms.untracked.int32(0))
                                                              ),
                                    SiStripQualityStatSummary = cms.untracked.PSet(threshold = cms.untracked.string('INFO'),
                                                                                   default = cms.untracked.PSet(limit=cms.untracked.int32(0)),
                                                                                   SiStripQualityStatistics = cms.untracked.PSet(limit=cms.untracked.int32(100000))
                                                                                   )
                                    
                                    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

### source file
process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    # The RunInfo for this run is NOT in the globalTag
    firstValue = cms.uint64(options.runNumber),
    lastValue = cms.uint64(options.runNumber),
    interval = cms.uint64(1)
)

### Geometry info
process.load("Configuration.Geometry.GeometryRecoDB_cff")


## connect to condDB
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')

process.siStripQualityESProducer.ListOfRecordToMerge=cms.VPSet(
    cms.PSet(record=cms.string('SiStripDetCablingRcd'),tag=cms.string(''))
    , cms.PSet(record=cms.string('SiStripBadChannelRcd'),tag=cms.string(''))
    , cms.PSet(record=cms.string('SiStripBadModuleRcd' ),tag=cms.string(''))
    , cms.PSet(record=cms.string('SiStripBadFiberRcd'),tag=cms.string(''))
    , cms.PSet(record=cms.string('SiStripBadStripRcd' ),tag=cms.string(''))
    , cms.PSet(record=cms.string('RunInfoRcd'),tag=cms.string(''))
)

process.siStripQualityESProducer.ReduceGranularity = cms.bool(False)
process.siStripQualityESProducer.PrintDebugOutput = cms.bool(True)
process.siStripQualityESProducer.UseEmptyRunInfo = cms.bool(False)

#### Add these lines to produce a tracker map
process.load("DQMServices.Core.DQMStore_cfg")
process.TkDetMap = cms.Service("TkDetMap")
process.SiStripDetInfoFileReader = cms.Service("SiStripDetInfoFileReader")
####

process.stat = cms.EDAnalyzer("SiStripQualityStatistics",
                              TkMapFileName = cms.untracked.string('TkMapBadComponents_offline.png'),
                              SaveTkHistoMap = cms.untracked.bool(True),                              
                              dataLabel = cms.untracked.string('test'),
                              PrintStripLevelInfo = cms.untracked.bool(True)
                              )

process.out = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.stat)
process.ep = cms.EndPath(process.out)

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
