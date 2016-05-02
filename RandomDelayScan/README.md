##################################
### Random Delay Scan Analysis ###
##################################

### CMSSW Setup:

    cmsrel CMSSW_8_0_4 ;
    cd CMSSW_8_0_4/src ;
    cmsenv;
    git-cms-init;
    git-cms-addpkg DPGAnalysis/SiStripTools/;
    git clone git@github.com:rgerosa/TrackerDAQAnalysis.git;
    cp TrackerDAQAnalysis/RandomDelayScan/plugins/TrackerDpgAnalysis.cc DPGAnalysis/SiStripTools/plugins;
    cp TrackerDAQAnalysis/RandomDelayScan/test/trackerdpganalysis_cfg.py DPGAnalysis/SiStripTools/test;
    cp TrackerDAQAnalysis/RandomDelayScan/test/crab/*xml DPGAnalysis/SiStripTools/test;
    scramv1 b -j 8;

### To run the DPGAnalyzer locally (modify by hand the cfg to give an example file to run on):

    cd  DPGAnalysis/SiStripTools/test/
    cmsRun trackerdpganalysis_cfg.py delayStep=<int corresponding to the delay xml chosen>
    
### To run crab jobs: typically one a random delya run is taken, it appears as a part of a dataset on DAS (Express physics stram, Miniu Bias prompt reco).

    TrackerDAQAnalysis/RandomDelayScan/test/crab/json*txt = example of json file to select a particular run
    TrackerDAQAnalysis/RandomDelayScan/test/crab/crab_*py = example of json crab config to analyze the delay run
    
### After running jobs, the output files belonging to a single run must be skimmed applying the desired analysis selection and dropping tracks/vertxes and event information. Once the selections are applied, we just need the delayMap and the clusters tree:

    Run Locally:

    cd TrackerDAQAnalysis/RandomDelayScan/macros/;
    root -l;
    .L skimTrees.C;
    skimTrees(<input file>, <outputfile>, <isBon = rule the selection string written inside the code>);

    
    LXBATCH submission (assumes files are stored on cern EOS):
    
    cd TrackerDAQAnalysis/RandomDelayScan ;
    pythn scripts/submitTreeSkim.py  --inputDIR <directory with all the files for a given run, produced by crab is ok> --outputDIR <output location on Cern EOS> --outputBaseName <base name for the output root file> --isBOn (in case you want to apply bOn selections) --batchMode --jobDIR <JOBDIR> --queque <QUEQUE> --submit
    

### Once skimmed trees are ready, to copy them to a local machine you could use:

    cd TrackerDAQAnalysis/RandomDelayScan ;
    python scripts/copyFilesEOS.py --inputDIR <directory where skimmed trees are located> --outputDIR <local directory to be copied>
    

### Merge trees belonging to a given run:

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L TreeMerge.C;
    TreeMerge(<destination file>, < one reference to take the readoutMap>, <directory where all the single files are located>, <if you want to cancel single inputs after merging)


### Run the fitting script to analyze charge or S/N for each detId:

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L fitChargeDistribution.C;
    fitChargeDistribution(<merged file>,<output dir>,<observable name (branch name)>, delayMin, delayMax, <apply or not path lenght correction>, <store outputs>)
    
The code produces a root file with some fits, just to see if they are making sense. Then, a text file is produced to be displayed on the tracker map: <detId,peak of landau conv gaussian shape>. 
To display it: http://test-stripdbmonitor.web.cern.ch/test-stripdbmonitor/PrintTrackerMap/print_TrackerMap.php

 
### Run the delay analysis on a single run to extract the delay per layer/ring/partition:

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayValidation.C;
    delayValidation(<merged file>,<no correction file stored in ../data/nocorrection.root>,<observable name (branch name)>, <plotParitions: to analyze the delay per partions>, <plotLayers: to find the delay per layer>, <plotSlices: to find the best delay per slices>, <outputDIR: name and path of the output directory>)

The code produces a set of root files with canvases for the different plots: profile fits and best delya vs partition, rings or layers (TEC divide by thin and thick sensors).

### Runt the delay analysis over a set of different runs with different random delay configuration (fine time calibration per module):

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayValidationPerModule.C;
    delayValidationPerModule(<input directory where all the merged files for different runs are located>,<no correction file stored in ../data/nocorrection.root>,<postfix: substring to be find to be sure to run on the merged files>, <observable name (branch name)>, <outputDIR: name and path of the output directory>,<saveCanvas: store some gaussian fits of chraged TProfile vs delay>, <saveCorrectionTree: to save the delay per channel in a TTree format. Can be analyzed then through the tkCommissioner>

