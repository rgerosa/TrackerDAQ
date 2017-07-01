Pedestal Analysis

CMSSW Setup:

      cmsrel CMSSW_9_2_4 ;
      cd CMSSW_9_2_4/src ;
      cmsenv;
      git-cms-init;
      git-cms-addpkg DataFormats/SiStripCommon/;
      git-cms-addpkg DPGAnalysis/SiStripTools/;
      git-cms-addpkg CommonTools/TrackerMap;
      git-cms-addpkg CondFormats/SiStripObjects	
      git cms-addpkg DQM/SiStripCommissioningAnalysis;
      git cms-addpkg DQM/SiStripCommissioningClients;
      git cms-addpkg DQM/SiStripCommissioningDbClients;
      git cms-addpkg DQM/SiStripCommissioningSources;	 
      git cms-addpkg DQM/SiStripCommissioningSummary;
      git clone git@github.com:rgerosa/TrackerDAQAnalysis.git -b Pedestal_Analysis_2017;
      scramv1 b -j 8;

Producing DQM files from Pedestal Runs:

      cd TrackerDAQAnalysis/PedestalAnalysis/test;
      cmsRun pedestalDQMfromDat_cfg.py partition=<partition> inputPath=<input path for the dat files> doFEDErr=<make map of bad channels> runNumber=<Run number to pick files in the /opt/cmssw directory>


Producing Source file from:
	  cd TrackerDAQAnalysis/PedestalAnalysis/test;
	  cmsRun pedestalSourcefromDQM_cfg.py inputFiles=<list of files> inputPath=<input directory> 	  