import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import glob

options = VarParsing.VarParsing ('analysis')

options.register ('inputPath',
	'.',
	VarParsing.VarParsing.multiplicity.singleton,
	VarParsing.VarParsing.varType.string,
	"Directory containing raw files"
)

options.parseArguments()

process = cms.Process("SIDQM")

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring()
)
process.source.fileNames.extend(cms.untracked.vstring())
for f in options.inputFiles :
	process.source.fileNames.extend(cms.untracked.vstring('file:'+options.inputPath+"/"+f))


process.load("DQMServices.Components.EDMtoMEConverter_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.CommissioningSourceHarvester = cms.EDAnalyzer("SiStripCommissioningSourceHarvester",
	RootFileName = cms.untracked.string("SiStripCommissioningSource.root"),
)

process.p = cms.Path(process.EDMtoMEConverter*process.CommissioningSourceHarvester)


