#ifndef __StPileupUtil_h__
#define __StPileupUtil_h__

#include <vector>
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include <string>
#include <TH1.h>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class TH1F;

//______________________________________________________________________________
class StPileupUtil {
  public:
    StPileupUtil ();
    virtual ~StPileupUtil(); /// Default destructor

    int init();
    int initEvent(const StPicoDst* pico);

    int get_centrality16() const;
    int get_centrality9() const;
    float get_centralityWeight() const;

    bool isPileupEPD(int option=0) const; //0 = default, 1 = 50% pileup rejection, 2 = 80%, 3 = 90%

    int   get_refMultPrim() const; 
    float get_nmipsum() const; 
    int   get_nch() const; 

  private:
    // Functions
    int  read() ; 
    int get_nchBin_a() const;
    int get_nchBin_b() const;

    bool  isGoodEvent           (const StPicoEvent* event) const;
    bool  isGoodTrackRefMultPrim(const StPicoTrack *ptrk, const StPicoDst* pico) const;
    bool  isGoodTrackNch        (const StPicoTrack *ptrk, const StPicoDst* pico) const;
    float getTofBeta            (StPicoTrack const* const trk, const StPicoDst* pico) const;

    int m_refMultPrim;
    int m_nch;
    float m_nmipsum;
    std::vector<int> m_goodRuns;
    std::vector<int> m_centCuts;
    std::vector<int> m_nchaBins;
    std::vector<int> m_nchbBins;

    TH1F* m_hpucuts_default;
    TH1F* m_hpucuts_fifty;
    TH1F* m_hpucuts_eighty;
    TH1F* m_hpucuts_ninty;
    TH1F* m_hreweight;

    bool kBadEventCentrality;
    bool kBadEventPileup;

    StEpdGeom *mEpdGeom;

    ClassDef(StPileupUtil, 0)
};

#endif

