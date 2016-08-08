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

def foo_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()

############################################                                                                                                                                 
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="",             help="to be used to correctly set the working area")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="",             help="output dir where workspaces are copied")
parser.add_option('--jsonFile',     action="store", type="string", dest="jsonFile",      default=""         ,    help="json file to be applied")
parser.add_option('--eventsPerJob', action="store", type=int,      dest="eventsPerJob",  default=500        ,    help="number of events for each job")
parser.add_option('--delayStep',    action="store", type=int,      dest="delayStep",     default=0          ,    help="used to pickup the right delay file")
parser.add_option('--isRawFile',    action="store_true", dest="isRawFile", help="isRawFile")
parser.add_option('--isDatFile',    action="store_true", dest="isDatFile", help="isDatFile")

############################################                                                                                                                                  
parser.add_option('--batchMode',    action="store_true", dest="batchMode", help="batchMode")
parser.add_option('--submit',       action="store_true", dest="submit",    help="submit")
parser.add_option('--jobDIR',       action="store",      type="string", dest="jobDIR", default="",  help="directory for job")
parser.add_option('--queque',       action="store",      type="string", dest="queque", default="",  help="queque for LSF")

(options, args) = parser.parse_args()
############################################                                                                                                                                  
if __name__ == '__main__':

  currentDIR = os.getcwd();

  if options.isRawFile and options.isDatFile:
    sys.exit('isRawFile and isDatFile cannnot be set to true at the same time');

  isRawFile = False;
  if options.isRawFile:
    isRawFile = True;

  isDatFile = False;
  if options.isDatFile:
    isDatFile = True;
  

  listOfFiles = [];
  os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+options.inputDIR+" -name *.root > file.temp ")
  file = open("file.temp","r")
  for line in file:
    if line == "" or line == "\n": continue;
    if not ".root" in line: continue;
    listOfFiles.append(line.replace("\n",""));
        
  os.system("rm file.temp");  
  os.system("mkdir -p "+options.jobDIR);
  os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select  mkdir -p "+options.outputDIR);

  ## open each file and split according to eventsPerJob
  ifile = 0;
  njob  = 0;
  for file in listOfFiles:
    print file
    tfile = ROOT.TFile.Open("root://eoscms.cern.ch//"+file);
    tree  = tfile.Get("Events");
    ## loop on events and save the starting one
    starEvent = [];
    for event in range(tree.GetEntries()):
      if event % options.eventsPerJob == 0:
        starEvent.append(event);
    print starEvent
    nJobs = len(starEvent);
    for ijob in range(nJobs):
      job = open('%s/job_file_%d_sub_%d.sh'%(currentDIR+"/"+options.jobDIR,ifile,ijob),'w')
      job.write('cd '+currentDIR+"\n");
      job.write('eval `scramv1 runtime -sh` \n');
      job.write('cd - \n');
      job.write('cp '+currentDIR+'/trackerdpganalysis_cfg.py ./ \n');
      job.write('cp '+currentDIR+'/*delaystep*'+str(options.delayStep)+'*xml ./ \n');
      job.write('cp '+currentDIR+'/'+options.jsonFile+' ./ \n');
      job.write('cmsRun trackerdpganalysis_cfg.py inputFiles=root://eoscms.cern.ch//'+file+' delayStep='+str(options.delayStep)+' eventToSkip='+str(starEvent[ijob])+' maxEvents='+str(options.eventsPerJob)+' ouputFileName=trackerDPG_'+str(njob)+".root jsonFile="+options.jsonFile+" isRawFile="str(isRawFile)+" isDatFile="+str(isDatFile)+"\n");
      job.write("xrdcp -f trackerDPG_"+str(njob)+".root root://eoscms.cern.ch//eos/cms/"+options.outputDIR+"\n");
      os.system('chmod a+x %s/job_file_%d_sub_%d.sh'%(options.jobDIR,ifile,ijob))
      if options.submit:
        os.system('bsub -q %s -o %s/job_file_%d_sub_%d.log -e %s/job_file_%d_sub_%d.err %s/job_file_%d_sub_%d.sh'%(options.queque,currentDIR+"/"+options.jobDIR,ifile,ijob,currentDIR+"/"+options.jobDIR,ifile,ijob,currentDIR+"/"+options.jobDIR,ifile,ijob));

      njob = njob+1;
    ifile = ifile + 1;
      

