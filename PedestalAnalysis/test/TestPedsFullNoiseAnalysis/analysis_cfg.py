import FWCore.ParameterSet.Config as cms
process = cms.Process("SiStripCommissioningOfflineDbClient")
process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.load("DQM.SiStripCommon.DaqMonitorROOTBackEnd_cfi")

process.load("OnlineDB.SiStripConfigDb.SiStripConfigDb_cfi")
process.SiStripConfigDb.UsingDb = False ### cause we don't have access to the db 
process.SiStripConfigDb.ConfDb  = 'overwritten/by@confdb'  
#process.SiStripConfigDb.Partitions.PrimaryPartition.PartitionName = 'TP_09-JUN-2009_1'
process.SiStripConfigDb.Partitions.PrimaryPartition.PartitionName = 'TM_09-JUN-2009_1'
#process.SiStripConfigDb.Partitions.PrimaryPartition.PartitionName = 'TO_30-JUN-2009_1'
#process.SiStripConfigDb.Partitions.PrimaryPartition.PartitionName = 'TI_27-JAN-2010_2'
#process.SiStripConfigDb.Partitions.PrimaryPartition.RunNumber     = 290534
process.SiStripConfigDb.Partitions.PrimaryPartition.RunNumber     = 290532
#process.SiStripConfigDb.Partitions.PrimaryPartition.RunNumber     = 290536
#process.SiStripConfigDb.Partitions.PrimaryPartition.RunNumber     = 290538
process.SiStripConfigDb.TNS_ADMIN = '/etc'

process.source = cms.Source("EmptySource") 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) ) 

### trick to get the FedCabling object
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data_GRun', '')


process.load("DQM.SiStripCommissioningDbClients.OfflineDbClient_cff")
#process.db_client.FilePath         = cms.untracked.string('/home/rgerosa/TrackerDAQ/PedestalAnalysis_BadChannels/WorkflowValidation/SourceFiles/TECP/')
process.db_client.FilePath         = cms.untracked.string('/home/rgerosa/TrackerDAQ/PedestalAnalysis_BadChannels/WorkflowValidation/SourceFiles/TECM/')
#process.db_client.FilePath         = cms.untracked.string('/home/rgerosa/TrackerDAQ/PedestalAnalysis_BadChannels/WorkflowValidation/SourceFiles/TOB/')
#process.db_client.FilePath         = cms.untracked.string('/home/rgerosa/TrackerDAQ/PedestalAnalysis_BadChannels/WorkflowValidation/SourceFiles/TIB/')
#process.db_client.RunNumber        = cms.untracked.uint32(290534)
process.db_client.RunNumber        = cms.untracked.uint32(290532)
#process.db_client.RunNumber        = cms.untracked.uint32(290536)
#process.db_client.RunNumber        = cms.untracked.uint32(290538)
process.db_client.UseClientFile    = cms.untracked.bool(False)
process.db_client.UploadHwConfig   = cms.untracked.bool(False)
process.db_client.UploadAnalyses   = cms.untracked.bool(False)
process.db_client.DisableDevices   = cms.untracked.bool(False)
process.db_client.DisableBadStrips = cms.untracked.bool(False)
process.db_client.SaveClientFile   = cms.untracked.bool(True)
#process.db_client.OutputRootFile   = cms.untracked.string("SiStripCommissioningClient_TECP")
process.db_client.OutputRootFile   = cms.untracked.string("SiStripCommissioningClient_TECM")
#process.db_client.OutputRootFile   = cms.untracked.string("SiStripCommissioningClient_TOB")
#process.db_client.OutputRootFile   = cms.untracked.string("SiStripCommissioningClient_TIB")

process.p = cms.Path(process.db_client)

