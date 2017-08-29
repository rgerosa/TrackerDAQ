#python scripts/mergeFile.py --inputDIR /home/rgerosa/MONOJET_ANALYSIS/Production-14-1-2016/QCD/ --outputName tree
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
parser.add_option('--inputDIR',    action="store", type="string", dest="inputDIR",      default="",   help="input directory where files are located")
parser.add_option('--outputName',  action="store", type="string", dest="outputName",    default="",   help="input directory where files are located")
(options, args) = parser.parse_args()

if __name__ == '__main__':

    originalDIR = os.getcwd();

    os.chdir(options.inputDIR);
    currentDIR = os.getcwd();
    print currentDIR

    
    os.system("ls "+options.inputDIR+" | grep -v txt | grep -v root  > file_temp.txt");

    fs = open("file_temp.txt","r");
    for line in fs:
        line = line.replace('\n','');
        print "go into directory --> ",line;
        print ;
        os.chdir(currentDIR+"/"+line);
        os.system("find . -name \"*.root\" | grep -v failed | grep .root  > file_2_temp.txt");        
        fs_2 = open("file_2_temp.txt","r");        
        print "create output file with name --> "+options.outputName+"_"+line+".root";
        print ;
        command = "hadd -f "+options.outputName+"_"+line+".root";
            
        for line_2 in fs_2:
            line_2 = line_2.replace('\n','');
            command += " "+line_2;
        print command;        
        print ;

        os.system(command);
        print "mv "+options.outputName+"_"+line+".root "+currentDIR
        print ;
        os.system("mv "+options.outputName+"_"+line+".root "+currentDIR);
        os.system("rm file_2_temp.txt");
        os.chdir(currentDIR);

    os.system("rm file_temp.txt")
     
