##################################
### Random Delay Scan Analysis ###
##################################

### CMSSW Setup:

    cmsrel CMSSW_10_1_1_patch1 ;
    cd CMSSW_10_1_1_patch1/src ;
    cmsenv;		      
    git-cms-init;
    git-cms-addpkg DPGAnalysis/SiStripTools/;
    git clone git@github.com:rgerosa/TrackerDAQAnalysis.git -b delayScan_10X;
    cp TrackerDAQAnalysis/RandomDelayScan/plugins/TrackerDpgAnalysis.cc DPGAnalysis/SiStripTools/plugins;
    cp TrackerDAQAnalysis/RandomDelayScan/test/trackerdpganalysis_cfg.py DPGAnalysis/SiStripTools/test;
    cp TrackerDAQAnalysis/RandomDelayScan/test/crab/*xml DPGAnalysis/SiStripTools/test;
    scramv1 b -j 4;					 

### To run the TrackerDPGAnalyzer locally:

    cd  TrackerDAQAnalysis/RandomDelayScan/test/
    cmsRun trackerdpganalysis_cfg.py delayStep=<int corresponding to the delay xml chosen>

Please: modify by hand the cfg to give an example file to run on
    
### To run crab jobs: 

    TrackerDAQAnalysis/RandomDelayScan/test/crab/json*txt = example of json file to select a particular run
    TrackerDAQAnalysis/RandomDelayScan/test/crab/crab_*py = example of json crab config to analyze the delay run

Typically after that one andom delay run is taken, it appears as a part of Express stream dataset on DAS. Otherwise, you can wait untill it appears inside the prompt reco dataset, typically Minimum Bias prompt reco during the commissioning time.


### To run lxbatch jobs: 

    TrackerDAQAnalysis/RandomDelayScan/test/crab/submitBatchJobs.py 
   
    
### Event Skim:

Run Locally:

    cd TrackerDAQAnalysis/RandomDelayScan/macros/;
    root -l;
    .L skimTrees.C;
    skimTrees(<input file>, <outputfile>, <isBon = rule the selection string written inside the code>);

After running jobs, the output files belonging to a single run must be skimmed applying the desired analysis selection and dropping tracks/vertxes and event information. Once these selections are applied, we just need the delayMap and the clusters tree for the offiline analysis.

LXBATCH submission (assumes files are stored on cern EOS):
    
    cd TrackerDAQAnalysis/RandomDelayScan ;
    pythn scripts/submitTreeSkim.py  --inputDIR <directory with all the files for a given run, produced by crab is ok> --outputDIR <output location on Cern EOS> --outputBaseName <base name for the output root file> --isBOn (in case you want to apply bOn selections) --batchMode --jobDIR <JOBDIR> --queque <QUEQUE> --submit
    

### Copy from EOS to a local machine:

    cd TrackerDAQAnalysis/RandomDelayScan ;
    python scripts/copyFilesEOS.py --inputDIR <directory where skimmed trees are located> --outputDIR <local directory to be copied>

Once skimmed trees are ready, to copy them to a local machine you could use:
    

### Merge trees belonging to a given run:

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L TreeMerge.C;
    TreeMerge(<destination file>, < one reference to take the readoutMap>, <directory where all the single files are located>, <if you want to cancel single inputs after merging)


### Fit Cluster Charge or S/N for each detId:

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L fitChargeDistribution.C;
    fitChargeDistribution(<merged file>,<output dir>,<observable name (branch name)>, delayMin, delayMax, <apply or not path lenght correction>, <store outputs>)
    
The code produces a root file with some fits, just to see if they are making sense. Then, a text file is produced to be displayed on the tracker map: <detId,peak of landau conv gaussian shape>. To display it: http://test-stripdbmonitor.web.cern.ch/test-stripdbmonitor/PrintTrackerMap/print_TrackerMap.php

 
### Delay analysis per layer/ring/partition:

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayValidation.C;
    delayValidation(<merged file>,<no correction file stored in ../data/nocorrection.root>,<observable name (branch name)>, <plotParitions: to analyze the delay per partions>, <plotLayers: to find the delay per layer>, <plotSlices: to find the best delay per slices>, <outputDIR: name and path of the output directory>)

The code produces a set of root files with canvases for the different plots: profile fits and best delya vs partition, rings or layers (TEC divide by thin and thick sensors).

### Delay analysis per module

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayValidationPerModule.C;
    delayValidationPerModule(<input directory where all the merged files for different runs are located>,<no correction file stored in ../data/nocorrection.root>,<postfix: substring to be find to be sure to run on the merged files>, <observable name (branch name)>, <outputDIR: name and path of the output directory>,<saveMeanCanvas: store some gaussian fits of mean charge vs delay>, <saveMPVCanvas: save MPV fit canvases vs delay >, <saveCorrectionTree: to save the delay per channel in a TTree format. Can be analyzed then through the tkCommissioner>

Run the delay analysis over a set of different runs with different random delay configuration (fine time calibration per module).

### Delay correction per module

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayCorrectionPerModule.C;
    delayCorrectionPerModule(<input file from delayValidationPerModule.C>,<output directory>,<outout root file name>)

### Compare different runs:
    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L compareClustersVsRun.C;
    compareClustersVsRun(<inputDIR : directory with all the dpg trees for all the runs> <runList: list of run numbers to be considered> <outputDIR: where plots and root files are created> <postfix: string to grep to access to input root files>);

The code produces a lighter output file with a tree called again delayCorrection. It contains three branches: Detid, fedCh and te correction in units of 25/24ns, to be propagated to the DB.

#######################################################
### Check delays with the Offile new tkCommissioner ###
#######################################################

Run the script inside the installation directory called INSTALL_QTROOT.sh to properly install QTROOT (please make sure that no cmsenv command run before, use the libraries installed with yum on your machine)

    mkdir <working directory with ~3GB space>
    cd <working directory with ~3GB space>
    cp <>/TrackerDAQAnalysis/RandomDelayScan/installation/INSTALL_QTROOT.sh ./
    sh INSTALL_QTROOT.sh

It will take ~2 hour to complie QT package, than ROOT 5.34 and finally QTRoot

Then copy the offline tkCommissioner folder and set the environment:

     scp -r /afs/cern.ch/user/r/rgerosa/work/public/OfflineCommissioner ./
     cp <>/TrackerDAQAnalysis/RandomDelayScan/installation/set_environment.sh ./
     source set_environment.sh

To run the final analysis:

    newCommissioner < root file > < tree name >


