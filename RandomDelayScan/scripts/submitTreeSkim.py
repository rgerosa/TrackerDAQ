### A priori assumption: run on cern batch system with files stored on EOS
### To run: python scripts/submitTreeSkim.py --inputDIR /store/user/rgerosa/TrackerDAQ/DELAYSCAN/ExpressPhysics/crab_delaySCAN_Express_Run271332/ --outputDIR /store/user/rgerosa/TrackerDAQ/DELAYSCAN/ExpressPhysics/SkimmedTrees/Run271332 --outputBaseName tree --batchMode --jobDIR JOB_Run_271332 --queque 1nh --submit
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT

from optparse import OptionParser
from subprocess import Popen

############################################                                                                                                                                   
#            Job steering                  #                                                                                                                                   
############################################                                                                                                                                    
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                 
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",     default="",   help="input directory where files are contained")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",    default="",   help="output DIR")
parser.add_option('--outputBaseName', action="store", type="string", dest="outputBaseName", default="",   help="output base name")
parser.add_option('--isBOn',        action="store_true", dest="isBOn", help="rule the track/event/vertex/clusters selection in the analysis")
##  for submitting jobs in lxbatch                                                                                                                                              
parser.add_option('--batchMode',    action="store_true",           dest="batchMode",                  help="batchMode")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="",        help="queque for LSF")
parser.add_option('--submit',       action="store_true",           dest="submit",                     help="submit")

(options, args) = parser.parse_args()

if __name__ == '__main__':


   print "################################";
   print "##### Start job submission #####";
   print "################################";
   
   currentDIR = os.getcwd();
   ## generate binary file                                                                                                                                                      
   ROOT.gROOT.ProcessLine(".L macros/skimTrees.C");

   os.system("rm -r "+options.jobDIR);
   
   ## make the file list ... typically all the files of a given run
   fileList = [];
   os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+options.inputDIR+" -name \"*.root\" > file_temp.txt");
   fs = open("file_temp.txt","r");
   for line in fs:
      if line == "": continue;
      if ".root" in line:
         line = line.replace('\n','');
         fileList.append(line);
   os.system("rm file_temp.txt");

   isBOn = 0;
   if options.isBOn:
       isBOn = 1;

   ## loop on the file list to create and submit jobs
   iFile = 1;
   for ifile in fileList:
       os.system("mkdir -p "+options.jobDIR)
       os.system("mkdir -p "+options.jobDIR+"/"+"JOB_"+str(iFile));

       nameList = ifile.split("/");
       while True:
          if nameList[len(nameList)-1] != '':
             break
          else:
             nameList.pop();
       fileName    = nameList[len(nameList)-1];
     

       ## write job sh file                                                                                                                                             
       jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_"+str(iFile)),'w')
       jobmacro.write("{\n");
       jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/skimTrees.C\");\n");
       jobmacro.write("gROOT->ProcessLine(\""+"skimTrees(\\\"%s\\\",\\\"%s\\\",%i)\");\n"%(fileName,options.outputBaseName+"_"+str(iFile)+".root",isBOn));
       jobmacro.write("}\n");
       jobmacro.close();

       jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_"+str(iFile)),'w')
       jobscript.write('cd %s \n'%currentDIR)
       jobscript.write('eval ` scramv1 runtime -sh ` \n')
       jobscript.write('cd - \n')
       jobscript.write("xrdcp -f root://eoscms.cern.ch//"+ifile+" ./\n")
       jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_"+str(iFile)))
       jobscript.write('root -l -b -q job.C\n');
       jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR+"/\n");
       jobscript.write("xrdcp -f "+options.outputBaseName+"_"+str(iFile)+".root root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/");       
       os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_"+str(iFile)));
       
       if options.submit:
           os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_"+str(iFile),currentDIR,options.jobDIR,"JOB_"+str(iFile),currentDIR,options.jobDIR,"JOB_"+str(iFile)));

       iFile += 1;
