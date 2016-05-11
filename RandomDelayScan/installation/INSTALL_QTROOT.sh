#!/bin/bash 
#
#######################################################
##
# $Id: INSTALL_QTROOT.sh 3547 2010-05-29 21:36:38Z fine@BNL.GOV $
#
# This macro tries to install a complete ROOT + QT + QTROOT(cvs) + COIN3D
# build environment.
#
# Original instructions here: http://root.bnl.gov/QtRoot/How2Install4Unix.html
#
# the **COMPLETE** stuff is installed in the current working directory.
# so you can also test different versions if you feel adventurous...
#
# you will need:
#  -  curl and/or wget to download the Qt and ROOT files.
#  -  svn to check out the Coin3D, QtRoot. root packages
# you should read also ALL the various licenses.
#
# about 100' on a AMD 3800+, ~3 Gb disk space.
##
#######################################################
# Author:                    L.Bardelli    [ bardelli  at fi  dot infn.it   ]
# MacOS patch (10.12.2009):  SIZUN Patrick [patrick dot sizun at cea dot fr ]
#######################################################

# the versions we will use:
QT_V=4.8
QT_VERSION=4.8.6
QT_WHERE=everywhere 
ROOT_VERSION_MAJOR=5.34
ROOT_VERSION_MINOR=00
ROOT_VERSION_MAJOR_GIT=5-34

# ROOT_VERSION_PATCHES=patches
ROOT_VERSION=$ROOT_VERSION_MAJOR.$ROOT_VERSION_MINOR
TAR_EXTENSION=tar.gz
MAKE_TOOL=make
NCPUS=0
copycommand="ln -s "
qtRootDir=qtRoot
COIN3DINSTALLDIR=coin3d
build_type=debug

#######################################################
# stop at the first error.
trap  "echo ; echo SOME FATAL ERROR DURING EXECUTION, SORRY... ; echo ; exit;" ERR


## http://en.wikipedia.org/wiki/Uname
PLATFORM=$(uname -s )
IV_PLATFORM=$(echo $PLATFORM)

uname -s | grep CYGWIN && PLATFORM=Win



if [ "$PLATFORM" != "Linux" ]; then
   if [ "$PLATFORM" != "Darwin" ]; then
      if [ "$PLATFORM" != "Win" ]; then
         echo "This macro is not tested outside Linux, Windows or MacOSX! Good luck..."
         sleep 5;
      fi
   fi
fi

##=================================================================
QT_CONF_OPTS=" -opensource -no-exceptions  -confirm-license $QT_CONF"
echo  -----  QT_CONF_OPT=$QT_CONF------  
echo  -----  QT_CONF_OPT=${QT_CONF_OPTS}------  
##=================================================================
if [ "$PLATFORM" == "Linux" ]; then
#--    QT_CONF_OPTS=-platform linux-g++-32 $QT_CONF_OPTS
    QT_CONF_OPTS=" -no-rpath -opengl $QT_CONF_OPTS"
    QT_PLATFORM=x11
    NCPUS=$(grep -e 'cpu[0-9]' /proc/stat | grep -c .)
fi
##=================================================================
if [ "$PLATFORM" == "Darwin" ]; then
    QT_PLATFORM=mac
    QT_CONF_OPTS=" -no-framework $QT_CONF_OPTS "
    ROOT_PLATFORM=macosx
    NCPUS=$(/usr/sbin/system_profiler SPHardwareDataType | grep Cores: | sed s/.*Cores://g )
    QMAKESPEC=macx-g++
    MACOSX_MAJOR_MINOR=$(/usr/bin/sw_vers -productVersion | awk -F. '{print $1 "." $2}')
fi
##=================================================================
if [ "$PLATFORM" == "Win" ]; then
    NCPUS=$(grep -e 'cpu[0-9]' /proc/stat | grep -c .)
    QT_PLATFORM=win
    TAR_EXTENSION=zip
    MAKE_TOOL=nmake
    copycommand="cp -rf "
#--
#--  Make sure the link from cygwin does not hide the link.exe from Visual C++
     linkexebase=$(dirname "$(which link)")
     mtexebase=$(dirname "$(which mt)")
     cygwinbase=$(dirname $(which cp))
     if  [ "$linkexebase" == "$cygwinbase" ]; then
         echo 
         echo FATAL EROOR:
         echo The cygwin command \"link\" from "$linkexebase" 
         echo HIDES the link.exe from Microsoft. We need the later to proceed.
         echo
         echo Please, DELETE !!! \"`which link`\"  and restart the script
         exit 1
     fi
     if  [ "$mtexebase" == "$cygwinbase" ]; then
         echo 
         echo FATAL EROOR:
         echo The cygwin command \"mt\" from "$mtexebase" 
         echo HIDES the mt.exe from Microsoft. We need the later to proceed.
         echo
         echo Please, DELETE !!! \"`which mt`\"  and restart the script
         exit 1
     fi
fi
##=================================================================
if  [ "$STAR" != "" ]; then
   IV_PLATFORM=$(echo $STAR_HOST_SYS)
fi

if [ "$NCPUS" -le "0" ]; then
    NCPUS=1;
fi
#if [ "$NCPUS" -ge "4" ]; then
# --    NCPUS=4;
#fi

##=================================================================
## this is the working directory where EVERYTHING will go.
MYWD=$(pwd)

##=================================================================
## If QTDIR was pre-installed ask the user whether he /she wants to re-use that installation
QTPREINSTALLED=$QTDIR
if [ "$QT_PLATFORM" == "win" ]; then
   if [ "$QTDIR" != "" ]; then 
      QTPREINSTALLED=`cygpath "$QTDIR"`
   fi
fi

if [ -f "$QTPREINSTALLED/include/Qt/qglobal.h" ]; then
  echo "--------------------$QTDIR -------------------"
  ls $QTPREINSTALLED
  echo "------------------------------------------------------"
  echo ""
  echo "The existing version of Qt package has been found under QTDIR=$QTDIR directory"
  echo ""
  read -p "Do you want to use it? (yes/no) " ANS
  if [ "$ANS" == "yes" ]; then
      touch done.qt
      touch download.qt
      export QTDIR=$QTPREINSTALLED
  else
      unset QTDIR
      export QTDIR=$MYWD/Qt-$(echo $QT_VERSION)
      echo "The Qt package will be re-built from the scratch and installed under $QTDIR"
  fi
else
   export QTDIR=$MYWD/Qt-$(echo $QT_VERSION)
fi

export ROOTSYS=$MYWD/root
export QTROOTSYSDIR=$ROOTSYS

##=================================================================

if [ -f done.root ]; then
   ANS=yes
   if [ -f done.coin3d ]; then
      USE_coin=yes
      echo "It seems to me the QT+ROOT(+COIN3D)  had been installed" 
   else
      echo "It seems to me the QT+ROOT with no COIN3D  had been installed"
   fi
   echo "If you want to start over then invoke this script from the empty directory"
else
   echo "========================================="
   echo "(L.B.)"
   echo "This will compile and install QT+ROOT(+COIN3D) in the current directory!!"
   echo ""
   echo "Current config: Qt v.$QT_VERSION opensource, ROOT v.$ROOT_VERSION"
   echo "                dir=$MYWD"
   echo "                PLATFORM=$PLATFORM, make will use $NCPUS cpu(s)"
   echo ""
   echo "It will require ~1-2 hours and ~3 Gb on disk. Current dir disk space is:"
   df -h $MYWD 
   echo "========================================="
   echo ""
   read -p "Do you want to proceed? (yes/no) " ANS

   if [ "$ANS" != "yes" ]; then
      exit 0;
   fi

   read -p "Do you want to install COIN3D also? (yes/no) " USE_coin
fi

if [ "$USE_coin" == "yes" ]; then
        ## no ending "/" here!!!
   COIN_VERSION=3     
   COIN3D_NEW_DIR=$MYWD/Coin3D/$COIN3DINSTALLDIR/$IV_PLATFORM/coin3d-${COIN_VERSION}
   if [ -f "$IVROOT/include/Inventor/Qt/SoQt.h" ]; then
      echo "--------------------$IVROOT -------------------"
      ls $IVROOT
      echo "------------------------------------------------------"
      echo ""
      echo "The existing version of Coin3D package has been found under IVROOT=$IVROOT directory"
      echo ""
      read -p "Do you want to use it? (yes/no) " ANS
      if [ "$ANS" == "yes" ]; then
         touch $MYWD/done.coin3d
         touch $MYWD/download.coin3d
         if [ "$QT_PLATFORM" == "win" ]; then
            export IVROOT=`cygpath -m $IVROOT`
         fi    
      else
         unset IVROOT
         export IVROOT=$(echo $COIN3D_NEW_DIR)
         if [ "$QT_PLATFORM" == "win" ]; then
            export IVROOT=`cygpath -m $(echo $COIN3D_NEW_DIR)`
         fi
         echo "The Coin3D package will be re-built from the scratch and installed under $IVROOT"
      fi
   else
      export IVROOT=$(echo $COIN3D_NEW_DIR)
      if [ "$QT_PLATFORM" == "win" ]; then
         export IVROOT=`cygpath -m $(echo $COIN3D_NEW_DIR)`
      fi
   fi
 
   if  [ "$STAR" == "" ]; then
      export LD_LIBRARY_PATH=$IVROOT/lib:$LD_LIBRARY_PATH
   fi
   if  [ "$STAR" != "" ]; then
      export LD_LIBRARY_PATH=$IVROOT/lib:`dropit coin -p $LD_LIBRARY_PATH`
   fi
fi


if  [ "$STAR" == "" ]; then
  #  remove the possible version of qt and root on PATH and LD_LIBRARY_PATH
  if [ "$USE_coin" == "yes" ]; then
      export LD_LIBRARY_PATH=$IVROOT/lib:$LD_LIBRARY_PATH
  fi
  export LD_LIBRARY_PATH=$ROOTSYS/lib:$QTDIR/lib:$LD_LIBRARY_PATH
  export PATH=$ROOTSYS/bin:$QTDIR/bin:$PATH
else
  if [ "$USE_coin" == "yes" ]; then
      export LD_LIBRARY_PATH=`dropit coin -p $LD_LIBRARY_PATH`
      export LD_LIBRARY_PATH=$IVROOT/lib:`dropit Coin -p $LD_LIBRARY_PATH`
  fi
  export LD_LIBRARY_PATH=`dropit qt -p $LD_LIBRARY_PATH`
  export LD_LIBRARY_PATH=$ROOTSYS/lib:$QTDIR/lib:`dropit root -p $LD_LIBRARY_PATH`
  export PATH=`dropit qt $PATH`
  export PATH=$ROOTSYS/bin:$QTDIR/bin:`dropit root $PATH`
fi

          echo "Make sure your ROOT version is known bug free ;)"
          echo "-----------------------------------------------"
          echo "Bug:     https://savannah.cern.ch/bugs/\?59014"
          echo "Patch: https://savannah.cern.ch/bugs/download.php\?file_id=12905"
          echo ""
          echo "Bug:     https://savannah.cern.ch/bugs/\?65402"
          echo "Patch:   https://savannah.cern.ch/bugs/download.php\?file_id=13273"
          echo "-----------------------------------------------"

echo ""
echo ">>>>> COMPILATION STARTS. Be patient..."
sleep 2 ## last chance for CTRL+C




## wget or curl? I prefer wget
GET="wget -c"
which wget || GET="curl -C - -O "


## ============= DOWNLOADS ==========================================
## with ls I check than the "true" file exists, and not a useless
## redirect from the server (like maintenance of similar...)
if [ "x$QT_WHERE" == "x" ]; then 
   QT_WHERE=${QT_PLATFORM}
fi
QTPKG=qt-$(echo $QT_WHERE)-opensource-src-$(echo $QT_VERSION)
QTPKGFILE=$(echo $QTPKG).$(echo $TAR_EXTENSION)

if [ ! -f download.qt ]; then
  echo ">>>>> DOWNLOADING  ........ Qt $QTPKGFILE.  Be patient..."
  $GET https://download.qt.io/archive/qt/$QT_V/$QT_VERSION/$QTPKGFILE
#  ftp://ftp.trolltech.com/qt/source/$QTPKGFILE
  ls -lh $QTPKGFILE
  touch download.qt
fi

if [ ! -f $MYWD/download.root ]; then
   if [ -d root ]; then
      echo ATTENTION !!! Removing the existent ROOT source tree . . . 
      rm -rf root
      if [ -f $MYWD/done.unpackroot ]; then
         rm $MYWD/done.unpackroot
      fi
   fi

   echo ">>>>> DOWNLOADING  ........ ROOT $ROOT_VERSION.$ROOT_VERSION_PATCHES  Be patient..."
   if [ "x$ROOT_VERSION" == "xtrunk" ]; then 
      svn -q co https://root.cern.ch/svn/root/trunk root
   else 
       git clone git@github.com:root-mirror/root.git -b v${ROOT_VERSION_MAJOR_GIT}-${ROOT_VERSION_MINOR} root
 #    $GET ftp://root.cern.ch/root/root_v$(echo $ROOT_VERSION).source.tar.gz
 #    ls -lh root_v$(echo $ROOT_VERSION).source.tar.gz
   fi
   if [ -d root ]; then
      touch $MYWD/done.unpackroot   
   fi
   touch $MYWD/download.root
fi

##OLDTAR $GET http://root.bnl.gov/QtRoot/downloads/qtFullRoot.tar.gz
##OLDTAR ls -lh qtFullRoot.tar.gz

if [ ! -f download.qtroot ]; then
  echo ">>>>> CHECKING OUT   ........ QtRoot.  Be patient..."
  echo "ATTENTION.  QtRoot CVS repository was replaced with SVN one"
#  svn -q co https://svn.bnl.gov/root/trunk qtRoot"
#  svn co https://qtroot.svn.sourceforge.net/svnroot/qtroot/trunk qtRoot
  svn co https://svn.code.sf.net/p/qtroot/code/trunk qtRoot
#  cvs -d :pserver:cvsuser:cvsuser@cvs.bnl.gov:/data01/CVS login
#  cvs -Q -d :pserver:cvsuser@cvs.bnl.gov:/data01/CVS co -Pd $qtRootDir root
  #- $GET http://root.bnl.gov/QtRoot/downloads/qtFullRoot.tar.gz
  #- ls -lh qtFullRoot.tar.gz
  #- tar -xzf qtFullRoot.tar.gz
  touch $MYWD/download.qtroot
fi



## ============= ENVIRONMENT.sh =====================================
rm -rf set_environment.*
echo "#!/bin/bash" > set_environment.sh
echo export QTDIR=$QTDIR >> set_environment.sh
echo export QMAKESPEC=$QMAKESPEC >> set_environment.sh
## echo export QT_CONF_OPTS=$QT_CONF_OPTS >> set_environment.sh
echo export ROOTSYS=$ROOTSYS >> set_environment.sh
echo export QTROOTSYSDIR=$QTROOTSYSDIR >> set_environment.sh
if [ "$USE_coin" == "yes" ]; then
    echo export IVROOT=$IVROOT >>  set_environment.sh
    echo export LD_LIBRARY_PATH=\$ROOTSYS/lib:\$IVROOT/lib:\$QTDIR/lib:\$LD_LIBRARY_PATH>> set_environment.sh
    echo export PATH=\$ROOTSYS/bin:\$QTDIR/bin:\$IVROOT/bin:\$PATH>> set_environment.sh
else
    echo export LD_LIBRARY_PATH=\$ROOTSYS/lib:\$QTDIR/lib:\$LD_LIBRARY_PATH>> set_environment.sh
    echo export PATH=\$ROOTSYS/bin:\$QTDIR/bin:\$PATH>> set_environment.sh
fi
if [ "$QTROOTSYSDIR" != "$ROOTSYSDIR" ]; then
    echo export LD_LIBRARY_PATH=\$QTROOTSYSDIR/lib:\$LD_LIBRARY_PATH>> set_environment.sh
fi
if [ "x$QT_PLATFORM" == "xmac" ]; then
   echo export  DYLD_LIBRARY_PATH=\$LD_LIBRARY_PATH>> set_environment.sh
fi
chmod -x  set_environment.sh ## user must "source" and not exec!


## ============= ENVIRONMENT.csh =====================================
cat set_environment.sh | sed s/=/" "/g | sed s/^export/set/g | sed s_"/bin/bash"_"/bin/csh"_g > set_environment.csh
chmod -x set_environment.csh

## ============= ENVIRONMENT.tcsh =====================================
cat set_environment.sh | sed s/=/" "/g | sed s/^export/setenv/g | sed s_"/bin/bash"_"/bin/tcsh"_g > set_environment.tcsh
chmod -x set_environment.tcsh

if [ "$PLATFORM" == "Win" ]; then
## ============= ENVIRONMENT.cmd =====================================
#cat set_environment.sh | sed s/^export/SET/g | sed s_"#!/bin/bash"_"rem Windows command file"_g | sed s_\$_\%_g > set_environment.cmd
vcbinpath="`which vcvars32.bat`"
vcbinpath=`cygpath -w $vcbinpath`

echo echo Set the VC environment                 > set_environment.cmd
echo CALL \"$vcbinpath\"                        >> set_environment.cmd
echo SET QTDIR=`cygpath -w $QTDIR`              >> set_environment.cmd
echo SET QMAKESPEC=$QMAKESPEC                   >> set_environment.cmd
echo SET QT_CONF_OPTS=\"$QT_CONF_OPTS\"         >> set_environment.cmd
echo SET ROOTSYS=`cygpath -w $ROOTSYS`          >> set_environment.cmd
echo SET QTROOTSYSDIR=`cygpath -w $QTROOTSYSDIR`>> set_environment.cmd
  if [ "$USE_coin" == "yes" ]; then
     echo SET IVROOT=`cygpath -w $IVROOT`                             >> set_environment.cmd
     echo "SET LIB=%ROOTSYS%\\lib;%IVROOT%\\lib;%QTDIR%\\lib;%LIB%"   >> set_environment.cmd
     echo "SET PATH=%QTDIR%\\bin;%ROOTSYS%\\bin;%IVROOT%\\bin;%PATH%" >> set_environment.cmd
  fi
  if [ "$USE_coin" != "yes" ]; then
     echo "SET LIB=%ROOTSYS%\\lib;%QTDIR%\\lib;%LIB%"                >> set_environment.cmd
     echo "SET PATH=%ROOTSYS%\\bin;%QTDIR%\\bin;%PATH%"             >> set_environment.cmd
  fi
  chmod -x set_environment.cmd
fi


## ============= for MAC-OSX: environment.plist ======================
echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?> 
<!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\"> 
<plist version=\"1.0\"> 
<dict> 
        <key>DYLD_LIBRARY_PATH</key> 
        <string>$ROOTSYS/lib</string> 
        <key>QTROOTSYSDIR</key> 
        <string>$ROOTSYS</string> 
        <key>ROOTSYS</key> 
        <string>$ROOTSYS</string> 
</dict> 
</plist> " > environment.plist

##===================== QT ==========================
cd $MYWD
if [ ! -f done.qt ]; then
    rm -rf $QTDIR
    if [ "$QT_PLATFORM" == "win" ]; then
	    echo 'Answer "no" to the question: "replace $QTPKG/configure?"'
       unzip -q $QTPKGFILE
       ls -lh $QTPKG
	   # There is no "configure -prefix" option for Windows 
       mv $QTPKG $(echo $QTDIR)
       cd $(echo $QTDIR)
    fi
    if [ "$QT_PLATFORM" != "win" ]; then
       tar xzf $QTPKGFILE
# -- use prefix     ln -s  $QTPKG $(echo $QTDIR)
       cd $(echo $QTPKG)
    fi
# -- use prefix      cd $(echo $QTDIR)
    if [ "$QT_PLATFORM" == "win" ]; then
	     chmod +x configure.exe
       ./configure $QT_CONF_OPTS
       $MAKE_TOOL
    fi
    if [ "$QT_PLATFORM" != "win" ]; then
        ./configure --prefix=$QTDIR $QT_CONF_OPTS  
        $MAKE_TOOL -j $NCPUS
    fi
    $MAKE_TOOL install
    echo We need no INSTALL to $QTDIR pwd $MAKE_TOOL install
    cd $MYWD
    touch done.qt 
fi

##=====================coin3d=============================
if [ "$USE_coin" == "yes" ]; then

    cd $MYWD
    if [ ! -f done.coin3d ]; then
        rm -f  done.qtroot
        if [ ! -d Coin3D/srcdir ]; then
            mkdir -p Coin3D/srcdir
        fi
        cd Coin3D/srcdir
        if [ ! -f $MYWD/download.coin ]; then
## server may have outdated certificate. (t)emporarily accept it.
    echo "t
t
t
t
t
t
t
t
t
t
t
t
t " | ../../qtRoot/qtgl/qtcoin/InstallCoin3D/download.sh
           touch $MYWD/download.coin
        fi
        rm -rf $IVROOT
        ../../qtRoot/qtgl/qtcoin/InstallCoin3D/installCoin3D${COIN_VERSION}.sh $COIN3DINSTALLDIR $build_type
        cd $MYWD
        touch $MYWD/done.coin3d
   fi
fi


##===================== ROOT ==========================
cd $MYWD
if [ ! -f $MYWD/done.unpackroot ]; then
    if [  -f $MYWD/done.root ]; then
       rm $MYWD/done.root
       rm $MYWD/done.qtroot
    fi
    rm -rf root
    tar xzf root_v$(echo $ROOT_VERSION).source.tar.gz
    touch $MYWD/done.unpackroot
fi    

## note1: DON\'T use --prefix!
## note2: mac seems to need "./configure macosx" to work properly...
if [ ! -f $MYWD/done.winpatch ]; then
    if [ "$QT_PLATFORM" == "win" ]; then
         if [  -f $MYWD/done.root ]; then
            rm $MYWD/done.root
            rm $MYWD/done.qtroot
         fi
         cd root
# --------- On Win32 we have to patch the core ROOT first  --------
#
# -- root
#
          $copycommand ../$qtRootDir/MyModules.mk .
          $copycommand ../$qtRootDir/MyMakefile.depend  .
          $copycommand ../$qtRootDir/root.diff .
          $copycommand ../$qtRootDir/plugins etc/
#
#  -- fix configure
#
        if [ -f root.diff/configure.$ROOT_VERSION_MAJOR ]; then 
           echo "-- fix configure"
           mv configure configure.root   
           $copycommand root.diff/configure.$ROOT_VERSION_MAJOR configure    
        fi
#
#  -- add qt to root
#
          rootdir=qt
          if [ -d graf2d ]; then 
            rootdir=graf2d/$rootdir
          fi
          mv $rootdir $rootdir.root
          $copycommand ../$qtRootDir/qt       $rootdir/

          rootdir=qtroot
          if [ -d gui/gui ]; then 
            rootdir=gui/$rootdir
          fi
          mv  $rootdir $rootdir.root
          $copycommand ../$qtRootDir/qtroot   $rootdir/

          rootdir=.
          $copycommand ../$qtRootDir/qtgui    $rootdir/
          $copycommand ../$qtRootDir/qt4ged   $rootdir/
          $copycommand ../$qtRootDir/qtged    $rootdir/
          $copycommand ../$qtRootDir/qtthread $rootdir/
          $copycommand ../$qtRootDir/qtgl     $rootdir/
          $copycommand ../$qtRootDir/qtimage  $rootdir/
#
# -- base
#
          echo " ---- Patching  build/win . . ."
          $copycommand ../$qtRootDir/root.diff/build/win/compiledata.sh build/win/
#
# -- base
#
          echo " ---- Patching  base/src . . ."
          rootdir=base
          if [ -d core ]; then 
            rootdir=core/$rootdir
          fi

          $copycommand ../$qtRootDir/root.diff/base/inc/TPadEditorHelper.h   $rootdir/inc
          $copycommand ../$qtRootDir/root.diff/base/src/TPadEditorHelper.cxx $rootdir/src
### Patch to activate the graphics lib under Windows in advance. 
          if [ -f root.diff/base/src/TApplication.$ROOT_VERSION_MAJOR.cxx ]; then 
             $copycommand ../$qtRootDir/root.diff/base/src/TApplication.$ROOT_VERSION_MAJOR.cxx $rootdir/src/TApplication.cxx
          fi
          echo "Patching core/base/src/TAttText.cxx"
          echo "Bug:     https://savannah.cern.ch/bugs/?59014"
          echo "Patch:   https://savannah.cern.ch/bugs/download.php?file_id=13273"
          if [ -f root.diff/base/src/TAttText.$ROOT_VERSION_MAJOR.cxx ]; then 
             $copycommand ../$qtRootDir/root.diff/base/src/TAttText.$ROOT_VERSION_MAJOR.cxx $rootdir/src/TAttText.cxx
          fi
# -- gui
#
          rootdir=gui
          if [ -d gui/gui ]; then 
            rootdir=gui/$rootdir
          fi
##
##    int maxPendingCounter = 30;
##    while (ProcessOneEvent() && maxPendingCounter--)          
##    There was a workaround agaisnt of infinite loop under Windows. it seems we do not need this protection anymore (5.22)
##  5.22         $copycommand ../$qtRootDir/root.diff/gui/src/TGClient.cxx       $rootdir/src
#
# -- clib
#
          echo " ---- Patching clib/src"
          rootdir=clib
          if [ -d core ]; then 
            rootdir=core/$rootdir
          fi
## There is some protection against of the infinite loop during the keyborad reading
          $copycommand ../$qtRootDir/root.diff/clib/src/Getline.c         $rootdir/src
#
# -- gpad -- fix the ROOT 5.14 bug (see TPad.diff for the futher details)
#
          rootdir=graf
          if [ -d graf2d ]; then 
            rootdir=graf2d/$rootdir
          fi
          echo "Patching    graf2d/graf/src/TPaveLabel.cxx"
          echo "Bug:    https://savannah.cern.ch/bugs/?59014"
          echo "Patch:  https://savannah.cern.ch/bugs/download.php?file_id=12905"
          if [ -f root.diff/graf2d/src/TPaveLabel.$ROOT_VERSION_MAJOR.cxx ]; then 
              $copycommand ../$qtRootDir/root.diff/graf2d/src/TPaveLabel.$ROOT_VERSION_MAJOR.cxx   $rootdir/src/TPaveLabel.cxx
          fi
#--          $copycommand ../$qtRootDir/root.diff/gpad/src/TPad.cxx      $rootdir/src

#
# -- thread -- fix the Win32 CPP flag
#
          if [ "$ROOT_VERSION" == "5.18.00" ]; then
             rootdir=thread
            if [ -d core ]; then 
              rootdir=core/$rootdir
            fi
## There was a compilation issue. It was fixed and the patch is not needed anymore
###            $copycommand ../$qtRootDir/root.diff/thread1/Module.mk           $rootdir/
###            $copycommand ../$qtRootDir/root.diff/thread1/src/TWin32Mutex.cxx $rootdir/src
          fi
# -- test
          echo " ---- Install extra tests"
          cd test
          $copycommand ../../$qtRootDir/test/qt .
          $copycommand ../../$qtRootDir/test/qtRootShower .
          cd ..
# -- config
        if [ -f root.diff/config/rootrc.in.$ROOT_VERSION_MAJOR ]; then 

           echo " ---- Install config"

           cd config
           mv rootrc.in rootrc.in.hold
           $copycommand ../root.diff/config/rootrc.in.$ROOT_VERSION_MAJOR  rootrc.in
           cd ..
        fi
#
# -- winnt
#
          echo " ---- Install winnt"
          rootdir=winnt
          if [ -d core ]; then 
            rootdir=core/$rootdir
          fi
## The patch to create the single threaded appliction.
## The recent ROOT version is single threaded and the patch is not needed .          
# ---   5.22                  $copycommand  ../$qtRootDir/root.diff/winnt/src/TWinNTSystem.$ROOT_VERSION_MAJOR.cxx $rootdir/src
        if [ -f root.diff/winnt/src/TWinNTSystem.$ROOT_VERSION_MAJOR.cxx ]; then 
           echo "-- patching TWinNTSystem"
           mv $rootdir/src/TWinNTSystem.cxx $rootdir/src/TWinNTSystem.cxx.root
           $copycommand root.diff/winnt/src/TWinNTSystem.$ROOT_VERSION_MAJOR.cxx   $rootdir/src/TWinNTSystem.cxx
        fi

#
# -- Tutorials
#          
          echo " ---- Install tutorials "
          cd tutorials
          $copycommand  ../../$qtRootDir/root.diff/tutorials/rootlogon.C .
          cd ..
          echo " ---- Install QT examples  "
          $copycommand ../$qtRootDir/qtExamples .
#
# -- There is no separate step to build QtRoot for Windows yet
#
         touch  $MYWD/done.qtroot
   fi
   touch $MYWD/done.winpatch
fi

#
# -- Configure and build ROOT
#

cd $MYWD
if [ ! -f done.root ]; then 
   cd root
   echo  "Replacing \"native\" type of GUI plugin  with \"qt\" and \"qtgui\""
   mv config/rootrc.in config/rootrc.in.native
   cat config/rootrc.in.native | sed   s/^Gui\.*Backend:\.*native\$/Gui\.Backend\:\ \ \ qt/   | sed   s/^Gui\.*Factory:\.*native\$/Gui\.Factory:\ \ \ \ \ qtgui/ >config/rootrc.in
   if [ "$QT_PLATFORM" = "win" ]; then
      if [ ! -f root.${ROOT_VERSION_MAJOR}.tar.gz ]; then
      cd ..
         echo Create a backup copy root.${ROOT_VERSION_MAJOR}.tar.gz to rebuild the project
         tar -czf root.${ROOT_VERSION_MAJOR}.tar.gz root  
         cd root
      fi
      ./configure $ROOT_PLATFORM --build=$build_type  --enable-qt --enable-table --disable-xrootd 
      make -j $NCPUS
#      make all-build
      make
      if [ ! -f $MYWD/done.qtroot ]; then 
		   touch $MYWD/done.qtroot
		fi
   fi
   if [ "$QT_PLATFORM" = "mac" ]; then
       sed -i.bak 's/\-ge \$decref98/-ge $decref97/g' net/xrootd/src/xrootd/configure.classic
   fi
   if [ "$QT_PLATFORM" != "win" ]; then
      ./configure  $ROOT_PLATFORM --build=$build_type  --enable-table --disable-xrootd --disable-alien
       make -j $NCPUS
       echo no prefix - no make install       make  install
   fi
   if [ "$QT_PLATFORM" = "mac" ]; then
	find lib -name "*.so" -exec bash -c 'ln -s `basename {}` `dirname {}`/`basename -s .so {}`.dylib' \;
   fi
   touch $MYWD/done.root
fi

##======================  QT-ROOT (2/2) ==========================
#cd $MYWD
#if [ ! -f $MYWD/done.qtroot ]; then
#    QGLVIEWER_DIR=qtRoot/qtgl/qglviewer/QGLViewer
##-- Build QGLViewer first
#    cd $QGLVIEWER_DIR
#if [ "$PLATFORM" != "Darwin" ]; then
#    qmake CONFIG+=$build_type PREFIX=$QTROOTSYSDIR
#else
#    qmake CONFIG+=$build_type PREFIX=$QTROOTSYSDIR QMAKE_MACOSX_DEPLOYMENT_TARGET=${MACOSX_MAJOR_MINOR}
#fi
#    make 
#    make install

# -- Build the QtRoot
if [ ! -f $MYWD/done.qtroot ]; then
   cd $MYWD
   cd qtRoot
   if [ "$PLATFORM" != "Darwin" ]; then
        qmake CONFIG+=$build_type
   else
        qmake CONFIG+=$build_type QMAKE_MACOSX_DEPLOYMENT_TARGET=${MACOSX_MAJOR_MINOR}
   fi
   make
   make install 
   cd $MYWD
#    ln -f -s $MYWD/qtRoot/qtgl/qglviewer/QGLViewer/libQGLViewer* $ROOTSYS/lib/
   if [ "$PLATFORM" == "Darwin" ]; then
      find root/lib -name "*.dylib" -exec bash -c 'ln -s `basename {}` `dirname {}`/`basename -s .dylib {}`.so 2> /dev/null' \;
      sed -i .bak "s/QT_CONF_OPTS=\(.*\)$/QT_CONF_OPTS=\"\1\"/g" set_environment.sh
      sed -i .bak "s/^features=\"\(.*\)\"/features=\"qt \1\"/g" root/bin/root-config
   fi
   touch  $MYWD/done.qtroot
fi

cd $MYWD
# ROOTRC=root/etc/system.rootrc
# --- ROOTRC=.rootrc
# --- if [ -f  $ROOTRC ]; then
# ---    rm $ROOTRC
# --- fi
# --- echo "To Activate Qt layer you need the custom rootrc file to be present either in yiour HOME or current directory"
# --- echo "" >  $ROOTRC
# --- echo "## added by INSTALL_QTROOT.sh ##" >>  $ROOTRC
# --- echo "" >>  $ROOTRC
# --- cat qtRoot/qtExamples/QtGBrowser/rootrcqtgui  >>   $ROOTRC


##===================== example test ======================
cd $MYWD
if [ "$QT_PLATFORM" != "win" ]; then
  cd qtRoot/qtExamples/HelloCanvas
  qmake
  make
fi

##====================== the end ==========================
cd $MYWD
cp -rf set_environment.* root/bin/
echo ""
echo "====================== DONE! ========================="
echo ""
echo "The local $build_type ROOT copy is Qt-enabled (patched rootrc)"
echo ""
if [ "USE_coin" == "yes" ]; then
    echo "you can use qtRoot/qtExamples/macros/rootgeom_coin.C to test COIN3D"
    echo ""
fi
echo "environment variables saved into $MYWD/set_environment.sh"
echo "                             and  $ROOTSYS/bin/set_environment.sh"
echo "            use it with:  source $MYWD/set_environment.sh"
echo "            or append it to your .bash_profile"
echo ""
echo "if your shell is csh use $MYWD/set_environment.csh"
echo ""
if [ "$PLATFORM" == "Win" ]; then
  echo ""
  echo "if your shell is Windows Command Prompt use `cygpath -w $MYWD/set_environment.cmd`"
  echo "                                         or `cygpath -w $ROOTSYS/bin/set_environment.cmd`"
  echo ""
fi
if [ "$PLATFORM" == "Darwin" ]; then
## mac-specific warnings
    echo "  MAC-OSX: you should now create a folder in your home named .MacOSX"
    echo "           and copy there the file environment.plist"
    echo ""
fi
echo "====================== BYE BYE! ======================"
echo ""
