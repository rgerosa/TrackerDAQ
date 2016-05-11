#!/bin/bash
export QTDIR=/home/rgerosa/CommissionerSoftware/QtRoot/Qt-4.8.6/
export QMAKESPEC=
export ROOTSYS=/home/rgerosa/CommissionerSoftware/QtRoot/root
export QTROOTSYSDIR=/home/rgerosa/CommissionerSoftware/QtRoot/root
export QXT=/home/rgerosa/CommissionerSoftware/QtRoot/libqxt
#export LD_LIBRARY_PATH=/opt/xdaq/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$QTDIR/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$QTDIR/bin:$PATH
export LD_LIBRARY_PATH=$QTROOTSYSDIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$QXT/lib

