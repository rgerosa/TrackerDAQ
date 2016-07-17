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
parser.add_option('--inputDelayFile',  action="store", type="string", dest="inputDelayFile",   default="",  help="file from delay analysis")
parser.add_option('--inputTimingFile',  action="store", type="string", dest="inputTimingFile",   default="",  help="file from delay analysis")
(options, args) = parser.parse_args()

if __name__ == '__main__':

    delays = {};
    delayFile = open(options.inputDelayFile,'r');
    for line in delayFile:
        line = line.replace('\n','')
        line = line.replace('\t',' ')
        lines = line.split(" ");
        lines = filter(None,lines) 
        if len(lines) != 2 : 
            sys.exit("problem with the delay file structure -- return");
        delays[lines[0]] = float(lines[1]);
        
    timing = {};
    timingFile = open(options.inputTimingFile,'r');
    for line in timingFile:
        line = line.replace('\n','')
        line = line.replace('\t',' ')
        lines = line.split(" ");
        lines = filter(None,lines)
        if lines[0] == "Device" or lines[0] == "Detid"  : continue;
        if len(lines) < 4: 
             sys.exit("problem with the timing file structure -- return");                         
        timing[lines[3]] = float(lines[1])

    outputfile = open("difference.txt",'w');
    for keys,values in delays.items():
        if keys in timing.keys():
            outputfile.write(str(keys)+" "+str(timing[keys]-values)+"\n");
    outputfile.close();
