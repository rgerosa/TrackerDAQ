import FWCore.ParameterSet.Config as cms
### CMSSW command line parameter parser                                                                                                                                        
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

process = cms.Process("Counter")

process.source = cms.Source("NewEventStreamFileReader",
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )

process.count = cms.EDFilter("countEvents");
process.path = cms.Path(process.count);

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10
