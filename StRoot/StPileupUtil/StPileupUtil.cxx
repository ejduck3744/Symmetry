#include <algorithm>
#include <string>
#include "StPileupUtil.h"

#include <TFile.h>
#include <TH1F.h>
#include <StThreeVectorF.hh>
#include <StHelix.hh>
#include <TLorentzVector.h>

#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StEpdUtil/StEpdGeom.h"
#include "phys_constants.h"

ClassImp(StPileupUtil)

using namespace std ;

//______________________________________________________________________________
// Default constructor
StPileupUtil::StPileupUtil() {

   m_refMultPrim=0; m_nch = 0; m_nmipsum = 0;
}

//______________________________________________________________________________
// Default destructor
StPileupUtil::~StPileupUtil() {

   delete m_hpucuts_default;
   delete m_hpucuts_fifty;
   delete m_hpucuts_eighty;
   delete m_hpucuts_ninty;
   delete m_hreweight;
}

//______________________________________________________________________________
int StPileupUtil::init() {

   int runsa[] = {19151031, 19151034, 19151036, 19151039, 19151041, 19151043, 19151044, 19151045,19151046, 19151047, 19151048, 19151049, 19151050, 19151052, 19151053, 19151054,19151055, 19151056, 19151066, 19151067, 19151068, 19151069, 19151070, 19151071,19151072, 19151082, 19151083, 19151084, 19152002, 19152003, 19152008, 19152009,19152010, 19152014, 19152016, 19152021, 19152023, 19152024, 19152025, 19152027,19152028, 19152029, 19152030, 19152031, 19152032, 19152033, 19152034, 19152035,19152036, 19152037, 19152038, 19152039, 19152040, 19152041, 19152042, 19152043,19152044, 19152045, 19152046, 19152048, 19152051, 19152052, 19152053, 19152054,19152055, 19152071, 19152073, 19152074, 19152075, 19152076, 19152081, 19153001,19153002, 19153003, 19153004, 19153007, 19153009, 19153010, 19153011, 19153012,19153013, 19153014, 19153015, 19153016, 19153017, 19153018, 19153019, 19153020,19153021, 19153022, 19153024, 19153025, 19153027, 19153028,19153029, 19153031, 19153033, 19153034, 19153035, 19153036, 19153037, 19153042,19153043, 19153044, 19153050, 19153051, 19153052, 19153053, 19153054, 19153055,19153056, 19153057, 19153058, 19153059, 19153061, 19153062, 19153063, 19153064,19153066, 19154001, 19154002, 19154005, 19154007, 19154027, 19154028, 19154029, 19154030, 19154031, 19154032, 19154036, 19154037, 19154038, 19154039, 19154040, 19154041, 19154044, 19154045, 19154046, 19154047, 19154048, 19154049, 19154052,19154053, 19154054, 19154055, 19154056, 19154057, 19154058, 19154061, 19154063,19154064, 19154065, 19154066, 19154067, 19155001, 19155003, 19155004, 19155005,19155006, 19155008, 19155009, 19155010, 19155011, 19155016, 19155017, 19155018, 19155019, 19155020, 19155021, 19155022};
                  
   int NRUNSA = 170; 
   for (int ir=0;ir<NRUNSA;ir++) m_goodRuns.push_back(runsa[ir]);

   int NCENT = 17;
   int centcuts16[] = {3,4,6,9,12,16,21,27,34,42,51,63,76,92,110,133,210};
   for (int ic=0;ic<NCENT;ic++) m_centCuts.push_back(centcuts16[ic]);

   int NCHA = 20, NCHB = 11;
   float nchb_bins[] = {0,8,10,12,14,18,22,28,34,40,48,300};
   float ncha_bins[] = {0,8,10,12,14,18,22,28,34,40,48,52,56,60,64,68,72,76,80,300};
   for (int in=0;in<NCHA;in++) m_nchaBins.push_back(ncha_bins[in]);
   for (int in=0;in<NCHB;in++) m_nchbBins.push_back(nchb_bins[in]);

   int isStatusOk = read();

   mEpdGeom = new StEpdGeom();

   return isStatusOk;;
}

//______________________________________________________________________________
int StPileupUtil::read() {

   TFile* fpileup = new TFile("StRoot/StPileupUtil/EPDPileUpCutsNch.root");
   m_hpucuts_default = (TH1F*)((TH1F*)fpileup->Get("hpucuts_default"))->Clone("m_hpucuts_default");  m_hpucuts_default->SetDirectory(0);
   m_hpucuts_fifty   = (TH1F*)((TH1F*)fpileup->Get("hpucuts_fiftyp"))->Clone("m_hpucuts_fifty");     m_hpucuts_fifty->SetDirectory(0);
   m_hpucuts_eighty  = (TH1F*)((TH1F*)fpileup->Get("hpucuts_eightyp"))->Clone("m_hpucuts_eighty");   m_hpucuts_eighty->SetDirectory(0);
   m_hpucuts_ninty   = (TH1F*)((TH1F*)fpileup->Get("hpucuts_nintyp"))->Clone("m_hpucuts_ninty");     m_hpucuts_ninty->SetDirectory(0);
   fpileup->Close();

   TFile* fcent = new TFile("StRoot/StPileupUtil/RefMultPrimGlauberFit.root");
   m_hreweight = (TH1F*)((TH1F*)fcent->Get("Ratioglaubertodata"))->Clone("m_hreweight"); m_hreweight->SetDirectory(0);
   fcent->Close();

   return 0;
}

int StPileupUtil::initEvent(const StPicoDst* pico) {

   StPicoEvent* picoEvent = (StPicoEvent*)pico->event();
   if( !picoEvent ){
      std::cout << " No PicoEvent! Skip! " << std::endl;
      return 1; kBadEventCentrality = 1; kBadEventPileup = 1;
   }

   if (!isGoodEvent(picoEvent)) {
      return 1; kBadEventCentrality = 1; kBadEventPileup = 1;
   }

   m_refMultPrim = 0;
   m_nch = 0;
   m_nmipsum = 0;

   const Int_t nTrack = pico->numberOfTracks();

   for (Int_t itr=0;itr<nTrack;itr++) {
      const StPicoTrack *ptrk = (StPicoTrack*)pico->track(itr);
      if(!ptrk)  continue;

      if (isGoodTrackRefMultPrim(ptrk,pico)) m_refMultPrim++;
      if (isGoodTrackNch(ptrk,pico))         m_nch++;
   }

   Int_t nepdHits = pico->numberOfEpdHits();
   if(nepdHits < 75) {return 1; kBadEventPileup = 1;}
   
   for(Int_t iHit=0; iHit<nepdHits; iHit++){
      StPicoEpdHit *epdHit = pico->epdHit(iHit);
      if( !epdHit) continue;
      float mip = epdHit->nMIP();

      TVector3 StraightLine_center;
      StraightLine_center = mEpdGeom->TileCenter(epdHit->id()) - picoEvent->primaryVertex();  
      double eta_epd_center = StraightLine_center.Eta();

      if ( (eta_epd_center>-4.0 && eta_epd_center<-2.5) && (mip>0.3 && mip<=4.0)) m_nmipsum  += mip;
   }

   return 0;
}

int StPileupUtil::get_centrality16() const {

   int centb = -1;
   if (kBadEventCentrality) return centb;
   
   for (unsigned int ic=m_centCuts.size()-1;ic>0;ic--) {
      if (m_refMultPrim >= m_centCuts[ic-1] && m_refMultPrim < m_centCuts[ic]) {centb = ic-1; break;}
   }

   return centb;
}

int StPileupUtil::get_centrality9() const {

   int centb = -1;
   if (kBadEventCentrality) return centb;

   int centb_16 = get_centrality16();
   if (centb_16 == 0  || centb_16 == 1 ) centb = 0;
   if (centb_16 == 2  || centb_16 == 3 ) centb = 1;
   if (centb_16 == 4  || centb_16 == 5 ) centb = 2;
   if (centb_16 == 6  || centb_16 == 7 ) centb = 3;
   if (centb_16 == 8  || centb_16 == 9 ) centb = 4;
   if (centb_16 == 10 || centb_16 == 11) centb = 5;
   if (centb_16 == 12 || centb_16 == 13) centb = 6;
   if (centb_16 == 14)                   centb = 7;
   if (centb_16 == 15)                   centb = 8;

   return centb;
}

float StPileupUtil::get_centralityWeight() const {

   if (kBadEventCentrality) return -1;

   return m_hreweight->GetBinContent(m_refMultPrim);
}

bool StPileupUtil::isPileupEPD(int option) const {

   if (kBadEventPileup) return 1;

   TH1F* hTemp;
   int nchbin;

   if (option==0)      hTemp = m_hpucuts_default;
   else if (option==1) hTemp = m_hpucuts_fifty;
   else if (option==2) hTemp = m_hpucuts_eighty;
   else if (option==3) hTemp = m_hpucuts_ninty;
   else {std::cout << "StPileupUtil::isPileup Unkown option" << std::endl; return 1;}

   if (option==0) nchbin = get_nchBin_a();
   else           nchbin = get_nchBin_b();

   if (nchbin<0) nchbin = 0;

   return m_nmipsum >= hTemp->GetBinContent(nchbin+1);
}

int StPileupUtil::get_refMultPrim() const {
   
   return m_refMultPrim;
}

float StPileupUtil::get_nmipsum() const {

   return m_nmipsum;
}

int StPileupUtil::get_nch() const {

   return m_nch;
}

bool StPileupUtil::isGoodEvent(const StPicoEvent* event) const {

   bool goodflag = 1;

   if (std::find(m_goodRuns.begin(),m_goodRuns.end(),event->runId()) == m_goodRuns.end()) goodflag = 0;

   if(!(event->isTrigger(620052))) goodflag = 0;

   Float_t vx=event->primaryVertex().X();
   Float_t vy=event->primaryVertex().Y();
   Float_t vz=event->primaryVertex().Z();
   if(vz < 198 || vz > 202)                goodflag = 0;
   if(( (vx*vx) + (vy+2.0)*(vy+2.0) ) > 4) goodflag = 0;

   return goodflag;
}

int StPileupUtil::get_nchBin_a() const {

   int centb = -1;

   for (unsigned int ic=0;ic<m_nchaBins.size()-1;ic++) {
      if (m_nch >= m_nchaBins[ic] && m_nch < m_nchaBins[ic+1]) {centb = ic; break;}
   }

   return centb;
}

int StPileupUtil::get_nchBin_b() const {

   int centb = -1;

   for (unsigned int ic=0;ic<m_nchbBins.size()-1;ic++) {
      if (m_nch >= m_nchbBins[ic] && m_nch < m_nchbBins[ic+1]) {centb = ic; break;}
   }

   return centb;
}

bool StPileupUtil::isGoodTrackRefMultPrim(const StPicoTrack *ptrk, const StPicoDst* pico) const {

   const Float_t dca       = ptrk->gDCA( pico->event()->primaryVertex() ).Mag();
   const Int_t   nHitsFit  = ptrk->nHitsFit();
   const Int_t   nHitsPoss = ptrk->nHitsMax();
   const Float_t quality   = (Float_t)nHitsFit/(Float_t)nHitsPoss;

   if (fabs(dca)>3.0 )     return false;
   if (quality < 0.52)     return false;
   if (nHitsFit < 10 )     return false;
   if (!ptrk->isPrimary()) return false;

   return true;
}

bool StPileupUtil::isGoodTrackNch(const StPicoTrack *ptrk, const StPicoDst* pico) const {

   const Float_t eta       = ptrk->pMom().PseudoRapidity();
   const Float_t dca       = ptrk->gDCA( pico->event()->primaryVertex() ).Mag();
   const Int_t   nHitsFit  = ptrk->nHitsFit();
   const Int_t   nHitsPoss = ptrk->nHitsMax();
   const Float_t quality   = (Float_t)nHitsFit/(Float_t)nHitsPoss;

   if (eta < -1.5 || -0.5 < eta) return false;
   if (fabs(dca)>3.0 )           return false;
   if (nHitsFit < 15 )           return false;
   if (quality < 0.52 )          return false;

   bool TofMatch = kFALSE;
   int tofIndex = ptrk->bTofPidTraitsIndex();
   StPicoBTofPidTraits* tofPidTraits;
   if (tofIndex >= 0)   tofPidTraits = pico->btofPidTraits(tofIndex);
   if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
   if (!TofMatch)                return false;

   float beta = getTofBeta(ptrk, pico);
   bool tofAvailable = !isnan(beta) && beta > 0;
   if (!tofAvailable)            return false;

   return true;                     
}

float StPileupUtil::getTofBeta(StPicoTrack const* const trk, const StPicoDst* pico) const {
   int index2tof = trk->bTofPidTraitsIndex();

   TVector3 vtxTV3 = pico->event()->primaryVertex();
   float beta = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0) {
      StPicoBTofPidTraits const* const tofPid = pico->btofPidTraits(index2tof);

      if (tofPid) {
         beta = tofPid->btofBeta();
         if (beta < 1e-4) {
            TVector3 const btofHitPosTV3 = tofPid->btofHitPos();
            StThreeVectorF const vtx(vtxTV3.X(),vtxTV3.Y(),vtxTV3.Z());
            StThreeVectorF btofHitPos(btofHitPosTV3.X(),btofHitPosTV3.Y(),btofHitPosTV3.Z());

            StPicoPhysicalHelix helix = trk->helix(pico->event()->bField());
            float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }
   return beta;
}
