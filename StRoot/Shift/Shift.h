#ifndef Shift_hh
#define Shift_hh
//
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StMaker.h"
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom3.h>
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
//
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StEpdEpFinder;
class StPicoTrack;
class StPileupUtil;
class TFile;
class TTree;
class TH1;
class TH2;
class TProfile;
class TProfile2D;
class TProfile3D;

// const Int_t CONST_VZ_BINS = 4; //epd vz bins {-145,-70,0,70,145}
const Int_t CONST_VZ_BINS = 4; //epd vz bins {positive, negative}
#ifndef ST_NO_NAMESPACES
using std::string;
#endif
//
//  The class declaration. It innherits from StMaker.
class Shift : public StMaker {

public:
  Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid);   // constructor
  virtual ~Shift(){};                                 // destructor

  virtual void Clear(Option_t *option=""); // called after every event to cleanup
  virtual Int_t  Init();                   // called once at the beginning of your job
  virtual Int_t  Make();                   // invoked for every event
  virtual Int_t  Finish();                 // called once at the end


  // My functions
  Int_t GetRunIndex( const Int_t run );
  Double_t Pileup_rejection(Double_t gRefmult, Int_t gNtofMatch, Double_t gvz,const StPicoEvent *event);
  bool   isGoodTrack(const StPicoTrack *ptrk);
  //bool   isGoodTrackSooraj(const StPicoTrack *ptrk);
  bool   isGoodEvent(const StPicoEvent *event);
  bool   isGoodTrigger(const StPicoEvent *event);
  bool	 isProton(const StPicoTrack *ptrk);
  bool   isPion(const StPicoTrack *ptrk);
  int    EP_group(int run_number, int cent);
  //double phi_correction(double phi, int which_one);//which one being which sub event plane.
  
private:
	//because I don't know how to use phys_constants:
	double protonMass;
	double y_cm;
	int do_eta_weighting;
	int fdEdXMode;
	
	TRandom3 *r1;

	//primary vertex
	TH2F *hist_VyVx;
	
	TH1D *hist_Vz;
	TH1D *hist_vpd_vz;
	TH2F *hist_vz_vpd_tpc;
	
	//QA profiles
	TProfile* run_vs_refmult;
	TProfile* run_vs_tofmult;
	TProfile* run_vs_vz_diff;
	TProfile* run_vs_avg_eta;
	TProfile* run_vs_sdca_xy;
	TH2F*     run_vs_sdca_xy_hist;
	TProfile* run_vs_avg_phi;
	TProfile* run_vs_avg_dedx;
	TProfile* run_vs_avg_pt;
	TProfile* run_vs_avg_nhits;
	TProfile* run_vs_avg_zdcx;
	TProfile* run_vs_avg_bbcx;

	// TH1D *hist_cent_weighted;
	//dedx:
	TH2F *hist_dEdx;
	TH2F *hist_dEdx_w_tof;
	

	TH2F *hist_tof_pid;

	TH2F *hist_phi_vs_eta;
	TH2F *hist_phi_vs_eta_cluster;
	TH2F *hist_phi_vs_eta_TOF;
	TH2F *hist_eta_vs_pt;

	// 0 - dEdX, 1 - dEdX positive tracks, 2 - dEdX negative tracks, 3 - dEdX tracks with ToF, 4 - ToF PID, 
	//5 - PV errors vs N tracks, 6 - PV errors vs N PV tracks, 7 dedx_sans_p_pi_k , 8 mass2, 9 mass2 sans_deuterons,
	//  10 hist_p_v_gb, 11 hist_p_v_gb_prim,12 hist_p_v_gb_dca
	TH2F* mass2_pt_pi_broad;
	TH2F* mass2_pt_proton_corrected_eta[10];
	TH2F* mass2_pt_proton_uncorrected_eta[10];
	TH2F* mass2_pt_pion_uncorrected_eta[10];
	TH2F* mass2_pt_kaon_uncorrected_eta[10];
	TH2F* mass2_inv_pt_proton_uncorrected_eta[10];
	
	TH2F* mass2_pt_proton_corrected_vz[14];
	TH2F* mass2_pt_proton_uncorrected_vz[14];
	TH2F* mass2_pt_pion_uncorrected_vz[14];
	TH2F* mass2_pt_kaon_uncorrected_vz[14];
	
	TH2F* mass2_pt_proton_corrected_mean_vz[40];
	TH2F* mass2_pt_proton_uncorrected_mean_vz[40];
	TH2F* mass2_pt_pion_uncorrected_mean_vz[40];
	TH2F* mass2_pt_kaon_uncorrected_mean_vz[40];
	TH2F* mass2_pt_proton_uncorrected_mean_vz_no_membrane[40];
	
	TH2F* mass2_phi_proton[12];
	TH2F* mass2_phi_proton_eta[10];
	TH2F* mass2_phi_pion[12];
	TH2F* mass2_phi_kaon[12];
	TH2F* mass2_phi_proton_eta_pos[12];
	TH2F* mass2_phi_proton_eta_neg[12];
	//TProfile* p_v_gb_full;
	TF1* delta_m2_over_m2;
	TTree* PID_tree;
	Int_t event_number;
	Int_t nHitsDedx;
	Int_t nHitsFit;
	float vx_tpc;
	float vy_tpc;
	float vz_tpc;
	float vz_tpc_error;
	float vz_vpd;
	float pt;
	float px;
	float py;
	float pz;
	float dca_x;
	float dca_y;
	float dca_z;
	Int_t q;
	float m2tof;
	float tof;
	float gb2_prim;
	float gb2_dca;
	float btofYLocal;
	float btofZLocal;
	float nsigma_He4;
	float nsigma_He3;
	float nsigma_tri;
	float nsigma_He4_old;
	float nsigma_He3_old;
	float nsigma_tri_old;
	float magn;
	float rescale;
	float dedx;
	float dedx_error;


	//PID histograms
	static const int NTrackHistoFolders = 18;
	TH2F* fHistodEdXTracks[NTrackHistoFolders];
	TH2F* fHistodEdXwithToFTracks[NTrackHistoFolders];
	TH2F* fHistoTofPIDTracks[NTrackHistoFolders];
	TH1D* fHistoMomentumTracks[NTrackHistoFolders];
	TH2F* fHistodEdXPull[NTrackHistoFolders];
	TH2F* fHistodEdXZ[NTrackHistoFolders];
	TH2F* fHistodEdX_of_tofpid[NTrackHistoFolders];
	TH2F* fHistotof_of_dEdXpid[NTrackHistoFolders];
	TH2F* fHisto_accuracy[NTrackHistoFolders];
	TH2F* fHisto_mass2[NTrackHistoFolders];
	TH1D* fh1d_mass2[NTrackHistoFolders];
	TH2F* fHisto_particle_acceptance[NTrackHistoFolders];
	// TH2F* hist_p_v_gb[NTrackHistoFolders];
	// TProfile* p_v_gb_particle[NTrackHistoFolders];
	// TH2F* hist_pt_m2[NTrackHistoFolders];
	// TProfile* profile_pt_m2[NTrackHistoFolders];
	// TH2F* hist_y_m2[NTrackHistoFolders];
	// TProfile* profile_y_m2[NTrackHistoFolders];
	
	
	//keeping these because they're definitely useful
	StPicoDstMaker *mPicoDstMaker;
	StPicoDst   *mPicoDst;
	StPicoEvent *picoEvent;
	StEpdGeom *mEpdGeom;
	StEpdEpFinder* mEpFinder;   
	StEpdEpInfo *mEpdEpInfo;
	StPileupUtil* mPileupTool;
	StRefMultCorr *mRefMultCorrUtil;

	TFile *File;
	
	TString mout_shift;
	
	
	ClassDef(Shift,0);
};
#endif
