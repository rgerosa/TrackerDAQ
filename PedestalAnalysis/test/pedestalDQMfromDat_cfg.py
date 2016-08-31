import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import glob
import os

options = VarParsing.VarParsing ('analysis')

options.register ('runNumber',
	276956,
	VarParsing.VarParsing.multiplicity.singleton,
	VarParsing.VarParsing.varType.int,
	"Run Number"
)

options.register ('partition',
		  #'TI_27-JAN-2010_2',
		  #'TM_09-JUN-2009_1',
		  #'TO_30-JUN-2009_1',
		  'TP_09-JUN-2009_1',
	VarParsing.VarParsing.multiplicity.singleton,
	VarParsing.VarParsing.varType.string,
	"Partition name"
)

options.register ('isOnEOS',False,
	VarParsing.VarParsing.multiplicity.singleton,
	VarParsing.VarParsing.varType.bool,
	"Tell if the location of the files is on cern eos"
)

options.register ('inputPath',
	'/opt/cmssw/Data',
	VarParsing.VarParsing.multiplicity.singleton,
	VarParsing.VarParsing.varType.string,
	"Directory containing raw files"
)

options.register ('doFEDErr',
    False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "Make the map of bad channels"
)

options.parseArguments()

process = cms.Process("SISRC")

process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.load("DQM.SiStripCommon.DaqMonitorROOTBackEnd_cfi")

process.load("OnlineDB.SiStripConfigDb.SiStripConfigDb_cfi")
process.SiStripConfigDb.UsingDb = True                   
process.SiStripConfigDb.ConfDb = 'user/password@account'
process.SiStripConfigDb.Partitions.PrimaryPartition.PartitionName = options.partition
process.SiStripConfigDb.Partitions.PrimaryPartition.RunNumber = options.runNumber
process.SiStripConfigDb.TNS_ADMIN = cms.untracked.string('/etc')

process.SiStripCondObjBuilderFromDb = cms.Service("SiStripCondObjBuilderFromDb")
process.FedCablingFromConfigDb = cms.ESSource("SiStripFedCablingBuilderFromDb",
    CablingSource = cms.untracked.string('UNDEFINED')
)
process.SiStripDetInfoFileReader = cms.Service("SiStripDetInfoFileReader")
process.PedestalsFromConfigDb = cms.ESSource("SiStripPedestalsBuilderFromDb")
process.NoiseFromConfigDb = cms.ESSource("SiStripNoiseBuilderFromDb")
process.sistripconn = cms.ESProducer("SiStripConnectivity")

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring()
)

if not options.isOnEOS:
	fnames = glob.glob(options.inputPath+"/"+str(options.runNumber)+"/RAW*.root")
	for f in fnames :
		process.source.fileNames.extend(cms.untracked.vstring('file:'+f))
else:	
	os.system('/afs/cern.ch/project/eos/installation/cms/bin/eos.select find '+options.inputPath+" -name \"*.root\" > file.tmp");
	fnames = open('file.tmp','r');
	for f in fnames :
		f = f.replace('\n','');
		if f == "" or not "root" in f: continue;
		process.source.fileNames.extend(cms.untracked.vstring('root://eoscms.cern.ch//'+f))
	os.system("rm file.tmp");
	 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("EventFilter.SiStripRawToDigi.FedChannelDigis_cfi")
process.FedChannelDigis.UnpackBadChannels = cms.bool(True)
process.FedChannelDigis.DoAPVEmulatorCheck = cms.bool(True)
process.FedChannelDigis.ProductLabel = cms.InputTag('rawDataCollector')
process.FedChannelDigis.LegacyUnpacker = cms.bool(True)

process.load("DQM.SiStripCommissioningSources.CommissioningHistos_cfi")
process.CommissioningHistos.CommissioningTask = 'PEDS_FULL_NOISE'
process.CommissioningHistos.PedsFullNoiseParameters.NrEvToSkipAtStart = 20
process.CommissioningHistos.PedsFullNoiseParameters.NrEvForPeds = 2000
process.CommissioningHistos.PedsFullNoiseParameters.FillNoiseProfile = True

process.p = cms.Path(process.FedChannelDigis*process.CommissioningHistos)

process.load("DQMServices.Components.MEtoEDMConverter_cff")
process.out = cms.OutputModule("PoolOutputModule",
	outputCommands = cms.untracked.vstring(
		'drop *',
		'keep *_MEtoEDMConverter_*_*'
		),
	fileName = cms.untracked.string('SiStripCommissioningSourceDQM_'+options.partition+'_Run'+str(options.runNumber)+'.root')
)
process.e = cms.EndPath(process.MEtoEDMConverter*process.out)

if options.doFEDErr : 
    process.load("Configuration.Geometry.GeometryIdeal_cff")
    process.load("DQM.SiStripMonitorHardware.siStripFEDMonitor_cfi")
    process.siStripFEDMonitor.RawDataTag = cms.untracked.InputTag("rawDataCollector")
    process.siStripFEDMonitor.FedIdVsApvIdHistogramConfig.Enabled = cms.untracked.bool(True)
    process.siStripFEDMonitor.FedEventSizeHistogramConfig.Enabled = cms.untracked.bool(False)
    process.siStripFEDMonitor.DataPresentHistogramConfig.Enabled = False
    process.siStripFEDMonitor.AnyFEDErrorsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.AnyDAQProblemsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.AnyFEProblemsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.CorruptBuffersHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadChannelStatusBitsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadActiveChannelStatusBitsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FEOverflowsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FEMissingHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadMajorityAddressesHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadMajorityInPartitionHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FeMajFracTIBHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FeMajFracTOBHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FeMajFracTECBHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FeMajFracTECFHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FETimeDiffTIBHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FETimeDiffTOBHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FETimeDiffTECBHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FETimeDiffTECFHistogramConfig.Enabled = False
    process.siStripFEDMonitor.ApveAddressHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FeMajAddressHistogramConfig.Enabled = False
    process.siStripFEDMonitor.MedianAPV0HistogramConfig.Enabled = False
    process.siStripFEDMonitor.MedianAPV1HistogramConfig.Enabled = False
    process.siStripFEDMonitor.DataMissingHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadIDsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadDAQPacketHistogramConfig.Enabled = False
    process.siStripFEDMonitor.InvalidBuffersHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadDAQCRCsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadFEDCRCsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.TkHistoMapHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FEOverflowsDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FEMissingDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadMajorityAddressesDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.BadAPVStatusBitsDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.APVErrorBitsDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.APVAddressErrorBitsDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.UnlockedBitsDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.OOSBitsDetailedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDErrorsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDDAQProblemsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDsWithFEProblemsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDCorruptBuffersHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nBadChannelStatusBitsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nBadActiveChannelStatusBitsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDsWithFEOverflowsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDsWithMissingFEsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDsWithFEBadMajorityAddressesHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nUnconnectedChannelsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nAPVStatusBitHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nAPVErrorHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nAPVAddressErrorHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nUnlockedHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nOutOfSyncHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nTotalBadChannelsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nTotalBadActiveChannelsHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nTotalBadChannelsvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nTotalBadActiveChannelsvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDErrorsvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDCorruptBuffersvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nFEDsWithFEProblemsvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nAPVStatusBitvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nAPVErrorvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nAPVAddressErrorvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nUnlockedvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.nOutOfSyncvsTimeHistogramConfig.Enabled = False
    process.siStripFEDMonitor.FedMaxEventSizevsTimeHistogramConfig.Enabled = False

    process.f = cms.Path(process.siStripFEDMonitor)

