#ifndef TrackerStrip_H
#define TrackerStrip_H

/// easy data-format to define a tracker strip

class TrackerStrip {

 public:
  
  TrackerStrip(){};
  ~TrackerStrip(){};
 TrackerStrip(uint16_t fecCrate, uint16_t fecSlot, uint16_t fecRing, uint16_t ccuAdd, uint16_t ccuCh, uint32_t fedKey, uint16_t lldCh, uint32_t detid, uint16_t apvid, uint32_t stripid):
  fecCrate_(fecCrate),
    fecSlot_(fecSlot),
    fecRing_(fecRing),
    ccuAdd_(ccuAdd),
    ccuCh_(ccuCh),
    fedKey_(fedKey),
    lldCh_(lldCh),
    detid_(detid),
    apvid_(apvid),
    stripid_(stripid){};

  bool operator == (const TrackerStrip & a) const {
    if(fecCrate_ == a.fecCrate_ and
       fecSlot_  == a.fecSlot_ and
       fecRing_  == a.fecRing_ and
       ccuAdd_   == a.ccuAdd_  and
       ccuCh_    == a.ccuCh_   and
       fedKey_   == a.fedKey_  and
       lldCh_    == a.lldCh_   and
       detid_    == a.detid_   and
       apvid_    == a.apvid_   and
       stripid_  == a.stripid_ )
      return true;
    else
      return false;
  }
  uint16_t fecCrate_;
  uint16_t fecSlot_;
  uint16_t fecRing_;
  uint16_t ccuAdd_;
  uint16_t ccuCh_;
  uint32_t fedKey_;
  uint16_t lldCh_;
  uint32_t detid_;
  uint16_t apvid_;
  uint32_t stripid_;

};

#endif
