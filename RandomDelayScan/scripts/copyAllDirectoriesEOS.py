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

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="",   help="input directory where files are located on eos")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="",   help="director to be copied")
parser.add_option('--toEOS',        action="store_true", dest = "toEOS", help="copy to eos system")
parser.add_option('--grepName', action="callback", type="string", dest="grepName", default="", callback=foo_callback, help="grep a set of names in the directory") 
parser.add_option('--skipName', action="callback", type="string", dest="skipName", default="", callback=foo_callback, help="drop a set of names in the directory")
parser.add_option('--nStreams',     action="store", type=int,      dest="nStreams",      default=1,    help="number of parallel streams")
parser.add_option('--nParallel',    action="store", type=int,      dest="nParallel",     default=2,    help="number of parallel streams")

(options, args) = parser.parse_args()

if __name__ == '__main__':
    

    if options.toEOS:
         os.system('/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select  mkdir -p '+options.outputDIR);
    else:
        os.system('mkdir -p '+options.outputDIR);

    if 'eos/cms/' in options.inputDIR or '/eos/cms/' in options.inputDIR:        
        options.inputDIR = options.inputDIR.replace('/eos/cms/','');
        options.inputDIR = options.inputDIR.replace('eos/cms/','');

    if 'eos/cms/' in options.outputDIR or '/eos/cms/' in options.outputDIR:
        options.outputDIR = options.outputDIR.replace('/eos/cms/','');
        options.outputDIR = options.outputDIR.replace('eos/cms/','');

    #### read all the subdirectories in the inputdit
    dirList = [];
    if not options.toEOS:   
        command = "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls "+options.inputDIR+" | grep -v txt | grep -v root | grep -v failed ";
        for name in options.grepName:
            command += " | grep "+name;
        for name in options.skipName:
            command += " | grep -v "+name;
        command += " > dir_list.txt";
        print command
        os.system(command);
        fs = open("dir_list.txt","r");
        for line in fs:
            line = line.replace('\n','');
            dirList.append(line);
        os.system("rm dir_list.txt");
    else:
        command = "ls "+options.inputDIR+" | grep -v txt | grep -v root | grep -v failed";
        for name in options.grepName:
            command += " | grep "+name;
        for name in options.skipName:
            command += " | grep -v "+name;
        command += " > dir_list.txt";
        print command
        os.system(command);
        fs = open("dir_list.txt","r");
        for line in fs:
            line = line.replace('\n','');
            dirList.append(line);
        os.system("rm dir_list.txt");    


    #### loop on sub directories
    for dir in dirList:
        if options.toEOS:
          os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR+"/"+dir);
          os.system('xrdcp -r -f -S '+str(options.nStreams)+' --parallel '+str(options.nParallel)+' '+options.inputDIR+"/"+dir+' root://eoscms.cern.ch//eos/cms'+options.outputDIR+"/"+dir);        
        else:
          os.system("mkdir -p "+options.outputDIR+"/"+dir);
          os.system('xrdcp -r -f -S '+str(options.nStreams)+' --parallel '+str(options.nParallel)+' root://eoscms.cern.ch//eos/cms'+options.inputDIR+"/"+dir+' '+options.outputDIR+"/"+dir);
