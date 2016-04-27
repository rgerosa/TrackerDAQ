#python ../../scripts/copyFilesEOS.py --inputDIR /store/user/rgerosa/MONOJET_ANALYSIS/Production-14-11-2016/PhotonJets/Jets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 --outputDIR /home/rgerosa/MONOJET_ANALYSIS/Production-14-1-2016/PhotonJets/

# file rename: ls | grep root | awk '{print "mv "$1"  "$1}' | sed 's/ScalarWH/MonoW_Scalar/2' | /bin/sh
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
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="",   help="input directory where files are located on eos")
parser.add_option('--nStreams',     action="store", type=int,      dest="nStreams",      default=3,    help="number of parallel streams")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="",   help="director to be copied")
parser.add_option('--toEOS',        action="store_true", dest = "toEOS", help="copy to eos system")
(options, args) = parser.parse_args()

if __name__ == '__main__':

    if options.toEOS:
         os.system('/afs/cern.ch/project/eos/installation/cms/bin/eos.select  mkdir -p '+options.outputDIR);
    else:
        os.system('mkdir -p '+options.outputDIR);

    if 'eos/cms/' in options.inputDIR or '/eos/cms/' in options.inputDIR:
        options.inputDIR = options.inputDIR.replace('/eos/cms/','');
        options.inputDIR = options.inputDIR.replace('eos/cms/','');

    if 'eos/cms/' in options.outputDIR or '/eos/cms/' in options.outputDIR:
        options.outputDIR = options.outputDIR.replace('/eos/cms/','');
        options.outputDIR = options.outputDIR.replace('eos/cms/','');

    if options.toEOS:
        os.system('xrdcp -r -f -S '+str(options.nStreams)+' '+options.inputDIR+' root://eoscms.cern.ch//eos/cms'+options.outputDIR);
        
    else:
        os.system('xrdcp -r -f -S '+str(options.nStreams)+' root://eoscms.cern.ch//eos/cms'+options.inputDIR+' '+options.outputDIR);
