#ifndef BadTrackerStrip_H
#define BadTrackerStrip_H

class BadTrackerStrip {

 public:
  
  BadTrackerStrip(){};
  ~BadTrackerStrip(){};
 BadTrackerStrip(uint32_t detid, uint16_t lldCh, uint16_t apvid, uint32_t stripid):
  detid_(detid),
    lldCh_(lldCh),
    apvid_(apvid),
    stripid_(stripid){};

  bool operator == (const BadTrackerStrip & a) const {
    if(detid_    == a.detid_   and
       lldCh_    == a.lldCh_   and
       apvid_    == a.apvid_   and
       stripid_  == a.stripid_ )
      return true;
    else
      return false;
  }

  uint32_t detid_;
  uint16_t lldCh_;
  uint16_t apvid_;
  uint32_t stripid_;

};

#endif
