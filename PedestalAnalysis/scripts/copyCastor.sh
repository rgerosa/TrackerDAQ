eos mkdir -p /store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157935/
eos mkdir -p /store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157936/
eos mkdir -p /store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157932/
eos mkdir -p /store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157927/

#nsls /castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/935/ | awk '{print "xrdcp root://castorcms.cern.ch//castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/935/"$1" root://eoscms.cern.ch//store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157935/"}' | /bin/sh

#nsls /castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/936/ | awk '{print "xrdcp root://castorcms.cern.ch//castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/936/"$1" root://eoscms.cern.ch//store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157936/"}' | /bin/sh

#nsls /castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/932/ | awk '{print "xrdcp root://castorcms.cern.ch//castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/932/"$1" root://eoscms.cern.ch//store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157932/"}' | /bin/sh

nsls /castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/927/ | awk '{print "xrdcp root://castorcms.cern.ch//castor/cern.ch/cms/store/streamer/SiStripCommissioning11-source/000/157/927/"$1" root://eoscms.cern.ch//store/caf/user/rgerosa/TrackerDAQ/PEDFULLNOISE/Source157927/"}' | /bin/sh
