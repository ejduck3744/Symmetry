#include <TFile.h>
#include <TTree.h>
#include <StMessMgr.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <StThreeVectorF.hh>
#include <StHelix.hh>
#include <TLorentzVector.h>

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
// #include "StPicoEvent/StPicoETofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StPileupUtil/StPileupUtil.h"
#include "StBichsel/Bichsel.h"
#include "StBichsel/StdEdxModel.h"
#include "phys_constants.h"

#include "../run/run.h"
#include "Shift.h"



//////if you have any issues with refmultcorr, see zuowen's code at:
///// /direct/gpfs01/star/pwg/liuzw/19_6GeV/piKp_v1/shift/StRoot/Shift
ClassImp(Shift)

    //__________________________________________________________________________________
    Shift::Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid) : StMaker(name) {
        mPicoDstMaker = picoMaker;
        mPicoDst = 0;
        mout_shift=Form("%s.root", jobid);
    }

//__________________________________________________________________________________
Int_t Shift::Init() {
    cout << "Init" << endl;
	do_eta_weighting = 1;
	fdEdXMode = 1;
    const float pi = acos(-1.0);                     
	mRefMultCorrUtil = new StRefMultCorr("refmult");
    mPileupTool = new StPileupUtil();
	mEpdGeom = new StEpdGeom();  
    mPileupTool->init();
    cout << "Define Histograms" << endl;
	r1 = new TRandom3();
	r1->GetSeed();
    //===============================                          
    //  Define Histograms                           
    //===============================                                            

    File=new TFile(mout_shift.Data(),"RECREATE");
	//because I don't know how to use phys_constants
	protonMass = 0.938;
	cout<<dataset<<endl;
	
	//primary vertex
	hist_VyVx = new TH2F("hist_VyVx","V_{Y} [cm] vs. V_{X} [cm];V_{X} [cm];V_{Y} [cm]",500,-5.0,5.0,500,-5.0,5.0);
	
	if(dataset == "4p5_gev_FXT_2020"){hist_Vz = new TH1D("hist_vz","V_{Z} (post-cut);V_{Z} (cm)",200,198,202);}
	else{hist_Vz = new TH1D("hist_vz","V_{Z} (post-cut);V_{Z} (cm)",600,-150,150);}
	
	hist_vpd_vz = new TH1D("hist_vpd_vz","V_{Z,vpd} (post-cut);V_{Z,vpd} (cm)",600,-150,150);
	hist_vz_vpd_tpc = new TH2F("hist_vz_vpd_tpc","V_{Z,vpd} vs. V_{Z,tpc};V_{Z,vpd} (cm);V_{Z,tpc} (cm)",600,-150,150,600,-150,150);
	
	TString p_name;
	TString p_title;
	
	//dedx:
	hist_dEdx = new TH2F("hist_dEdx","dE/dx vs |p|/q;|p|/q (GeV/c);dE/dx (keV/cm)",500,-10.0,10.0,1000,0.0,50.0);
	hist_dEdx_w_tof = new TH2F("hist_dEdx_w_tof","dE/dx vs |p|/q for ToF Tracks;|p|/q (GeV/c);dE/dx (keV/cm)",500,-10.0,10.0,1000,0.0,50.0);
	
	//TOF:
	hist_tof_pid=new TH2F("TOF","TOF;|p|/q (GeV/c);1/#beta",1000,-5.0,5.0,1000,0.0,10.0);
	
	hist_phi_vs_eta =new TH2F("hist_phi_vs_eta","#phi vs. #eta ",500,-5.0,5.0,500,0.0,2.0*pi);
	hist_phi_vs_eta->GetXaxis()->SetTitle("#eta");
	hist_phi_vs_eta->GetYaxis()->SetTitle("#phi");
	
	hist_phi_vs_eta_TOF =new TH2F("hist_phi_vs_eta_TOF","#phi vs. #eta TOF hit",500,-5.0,5.0,500,0.0,2.0*pi);
	hist_phi_vs_eta_TOF->GetXaxis()->SetTitle("#eta");
	hist_phi_vs_eta_TOF->GetYaxis()->SetTitle("#phi");
	
	hist_eta_vs_pt =new TH2F("hist_eta_vs_pt","#eta vs. P_{T}/q",500,-5.0,5.0,500,-5.0,5.0);
	hist_eta_vs_pt->GetXaxis()->SetTitle("P_{T}/q (GeV)");
	hist_eta_vs_pt->GetYaxis()->SetTitle("#eta");
	
	for(int i =0;i<10;i++)
	{
		p_name =Form("mass2_pt_proton_corrected_eta_%d",i);
		p_title =Form("mass^{2} vs P_{T} for proton %f<#eta<%f;P_{T}/q (GeV) ;mass^{2} (GeV^{2})",(0.2*i)-1.0,(0.2*i)-0.8);
		mass2_pt_proton_corrected_eta[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,2.0);
		p_name =Form("mass2_pt_proton_uncorrected_eta_%d",i);
		p_title =Form("mass^{2} vs P_{T} for proton %f<#eta<%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(0.2*i)-1.0,(0.2*i)-0.8);
		mass2_pt_proton_uncorrected_eta[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,2.0);
		p_name =Form("mass2_pt_pion_uncorrected_eta_%d",i);
		p_title =Form("mass^{2} vs P_{T} for #pi %f<#eta<%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(0.2*i)-1.0,(0.2*i)-0.8);
		mass2_pt_pion_uncorrected_eta[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,0.5);
		p_name =Form("mass2_pt_kaon_uncorrected_eta_%d",i);
		p_title =Form("mass^{2} vs P_{T} for k %f<#eta<%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(0.2*i)-1.0,(0.2*i)-0.8);
		mass2_pt_kaon_uncorrected_eta[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,1.0);
		p_name =Form("mass2_inv_pt_proton_uncorrected_eta_%d",i);
		p_title =Form("mass^{2} vs P_{T}^{-1} for proton %f<#eta<%f;P_{T}/q (GeV) ;mass^{2} (GeV^{2})",(0.2*i)-1.0,(0.2*i)-0.8);
		mass2_inv_pt_proton_uncorrected_eta[i] =new TH2F(p_name.Data(),p_title.Data(),1000,-2.0,2.0,1000,0.0,2.0);
		p_name =Form("mass2_phi_proton_eta_%d",i);
		p_title =Form("mass^{2} vs #phi for proton %f<#eta<%f;#phi;mass^{2} (GeV^{2}/c^{4})",(0.2*i)-1.0,(0.2*i)-0.8);
		mass2_phi_proton_eta[i] = new TH2F(p_name.Data(),p_title.Data(),384,0.0,2.0*pi,2000,-2.0,2.0);
	}
	for(int i =0;i<14;i++)
	{
		p_name =Form("mass2_pt_proton_corrected_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for proton %f<V_{Z}<%f;P_{T}/q (GeV) ;mass^{2} (GeV^{2})",(10.0*i)-70.0,(10.0*i)-60.0);
		mass2_pt_proton_corrected_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,2.0);
		p_name =Form("mass2_pt_proton_uncorrected_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for proton %f<V_{Z}<%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(10.0*i)-70.0,(10.0*i)-60.0);
		mass2_pt_proton_uncorrected_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,2.0);
		p_name =Form("mass2_pt_pion_uncorrected_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for #pi %f<V_{Z}<%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(10.0*i)-70.0,(10.0*i)-60.0);
		mass2_pt_pion_uncorrected_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,0.5);
		p_name =Form("mass2_pt_kaon_uncorrected_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for k %f<V_{Z}<%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(10.0*i)-70.0,(10.0*i)-60.0);
		mass2_pt_kaon_uncorrected_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,1.0);
	}
	for(int i = 0;i<40;i++)
	{
		p_name =Form("mass2_pt_proton_corrected_mean_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for proton %f<<V_{Z,track}><%f;P_{T}/q (GeV) ;mass^{2} (GeV^{2})",(10.0*i)-200.0,(10.0*i)-190.0);
		mass2_pt_proton_corrected_mean_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,2.0);
		p_name =Form("mass2_pt_proton_uncorrected_mean_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for proton %f<<V_{Z,track}><%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(10.0*i)-200.0,(10.0*i)-190.0);
		mass2_pt_proton_uncorrected_mean_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,2.0);
		p_name =Form("mass2_pt_pion_uncorrected_mean_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for #pi %f<<V_{Z,track}><%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(10.0*i)-200.0,(10.0*i)-190.0);
		mass2_pt_pion_uncorrected_mean_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,0.5);
		p_name =Form("mass2_pt_kaon_uncorrected_mean_vz_%d",i);
		p_title =Form("mass^{2} vs P_{T} for k %f<<V_{Z,track}><%f uncorrected;P_{T}/q (GeV);mass^{2} (GeV^{2})",(10.0*i)-200.0,(10.0*i)-190.0);
		mass2_pt_kaon_uncorrected_mean_vz[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,1.0);
		p_name =Form("mass2_pt_proton_uncorrected_mean_vz_no_membrane_%d",i);
		p_title =Form("mass^{2} vs P_{T} for proton %f<<V_{Z,track}><%f uncorrected no membrane;P_{T}/q (GeV);mass^{2} (GeV^{2})",(10.0*i)-200.0,(10.0*i)-190.0);
		mass2_pt_proton_uncorrected_mean_vz_no_membrane[i] =new TH2F(p_name.Data(),p_title.Data(),500,-5.0,5.0,1000,0.0,2.0);
	}
	for(int i = 0;i<12;i++)
	{
		p_name =Form("mass2_phi_proton_%d",i);
		p_title =Form("mass^{2} vs #phi for proton bin %d;#phi;mass^{2} (GeV^{2}/c^{4})",i);
		mass2_phi_proton[i] =new TH2F(p_name.Data(),p_title.Data(),384,0.0,2.0*pi,2000,-2.0,2.0);
		p_name =Form("mass2_phi_pion_%d",i);
		p_title =Form("mass^{2} vs #phi for #pi bin %d;#phi;mass^{2} (GeV^{2}/c^{4})",i);
		mass2_phi_pion[i] =new TH2F(p_name.Data(),p_title.Data(),384,0.0,2.0*pi,2000,-0.5,0.5);
		p_name =Form("mass2_phi_kaon_%d",i);
		p_title =Form("mass^{2} vs #phi for kaon bin %d;#phi;mass^{2} (GeV^{2}/c^{4})",i);
		mass2_phi_kaon[i] =new TH2F(p_name.Data(),p_title.Data(),384,0.0,2.0*pi,2000,-1.0,1.0);
		p_name =Form("mass2_phi_proton_eta_pos_%d",i);
		p_title =Form("mass^{2} vs #phi for proton #eta>0.0 bin %d;#phi;mass^{2} (GeV^{2}/c^{4})",i);
		mass2_phi_proton_eta_pos[i] =new TH2F(p_name.Data(),p_title.Data(),384,0.0,2.0*pi,2000,-2.0,2.0);
		p_name =Form("mass2_phi_proton_eta_neg_%d",i);
		p_title =Form("mass^{2} vs #phi for proton #eta<0.0 bin %d;#phi;mass^{2} (GeV^{2}/c^{4})",i);
		mass2_phi_proton_eta_neg[i] =new TH2F(p_name.Data(),p_title.Data(),384,0.0,2.0*pi,2000,-2.0,2.0);
	}
	//pid histograms:
	
	for(int iTrackHisto=0; iTrackHisto<NTrackHistoFolders; iTrackHisto++)
	{
		p_name =Form("hdEdX_%d",iTrackHisto);
		p_title =Form("hdEdX_%d",iTrackHisto);
		// fTrackPdgToHistoIndex[ pdgTrackHisto[iTrackHisto] ] = iTrackHisto;
		fHistodEdXTracks[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 500, 0, 10, 1500, 0, 100);
		p_name =Form("hdEdXwithToF_%d",iTrackHisto);
		p_title =Form("hdEdXwithToF_%d",iTrackHisto);
		fHistodEdXwithToFTracks[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 500, 0, 10, 1500, 0, 100);
		p_name =Form("hTofPID_%d",iTrackHisto);
		p_title =Form("hTofPID_%d",iTrackHisto);
		fHistoTofPIDTracks[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 1000, -6, 6, 500, 0, 6);
		p_name =Form("hMomentum_%d",iTrackHisto);
		p_title =Form("hMomentum_%d",iTrackHisto);
		fHistoMomentumTracks[iTrackHisto] = new TH1D(p_name.Data(),p_title.Data(), 500, 0, 10);
		p_name =Form("hdEdXPull_%d",iTrackHisto);
		p_title =Form("hdEdXPull_%d",iTrackHisto);
		fHistodEdXPull[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 500, 0, 10, 600, -30, 30);
		p_name =Form("hdEdXZ_%d",iTrackHisto);
		p_title =Form("hdEdXZ_%d",iTrackHisto);
		fHistodEdXZ[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 200, -5, 5, 280, -1, 6);
		p_name =Form("dedx_of_tof_%d",iTrackHisto);
		p_title =Form("dedx_of_tof_%d",iTrackHisto);
		fHistodEdX_of_tofpid[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 500, 0, 10, 1500, 0, 30);
		p_name =Form("tof_of_dedx_%d",iTrackHisto);
		p_title =Form("tof_of_dedx_%d",iTrackHisto);
		fHistotof_of_dEdXpid[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 500, 0, 6, 500, 0, 6);
		p_name =Form("pid_accuracy_%d",iTrackHisto);
		p_title =Form("pid_accuracy_%d",iTrackHisto);
		fHisto_accuracy[iTrackHisto] =new TH2F(p_name.Data(),p_title.Data(),60,-3.0,3.0,1500,mass[iTrackHisto]*mass[iTrackHisto]*0.7,mass[iTrackHisto]*mass[iTrackHisto]*1.3);
		p_name =Form("mass_2_%d",iTrackHisto);
		p_title =Form("mass_2_%d",iTrackHisto);
		fHisto_mass2[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(), 500, 0, 6.0, 1500,mass[iTrackHisto]*mass[iTrackHisto]*0.7,mass[iTrackHisto]*mass[iTrackHisto]*1.3);
		p_name =Form("mass_2_1d_%d",iTrackHisto);
		p_title =Form("mass_2_1d_%d",iTrackHisto);
		fh1d_mass2[iTrackHisto] = new TH1D(p_name.Data(),p_title.Data(),5000,mass[iTrackHisto]*mass[iTrackHisto]*0.7,mass[iTrackHisto]*mass[iTrackHisto]*1.3);
		p_name =Form("acceptance_%d",iTrackHisto);
		p_title =Form("acceptance_%d",iTrackHisto);
		fHisto_particle_acceptance[iTrackHisto] = new TH2F(p_name.Data(),p_title.Data(),600,-3.0,3.0,500,0.0,6.0);
		
		// hist_p_v_gb[iTrackHisto]= new TH2F("hist_p_v_gb","momentum vs #gamma#beta",500,0,10,1500,0,150);
		// hist_p_v_gb[iTrackHisto]->GetXaxis()->SetTitle("|P| (GeV)");
		// hist_p_v_gb[iTrackHisto]->GetYaxis()->SetTitle("#gamma#beta");

		// p_v_gb_particle[iTrackHisto]= new TProfile("p_v_gb_particle","momentum_vs_#gamma#beta",50,0,10);
		// p_v_gb_particle[iTrackHisto]->GetXaxis()->SetTitle("|P| (GeV)");
		// p_v_gb_particle[iTrackHisto]->GetYaxis()->SetTitle("#gamma#beta");

		// hist_pt_m2[iTrackHisto] = new TH2F("hist_pt_m2","P_{T} vs m^{2}",100,0,10,1500,mass[iTrackHisto]*mass[iTrackHisto]*0.7,mass[iTrackHisto]*mass[iTrackHisto]*1.3);
		// hist_pt_m2[iTrackHisto]->GetXaxis()->SetTitle("P_{T} (GeV)");
		// hist_pt_m2[iTrackHisto]->GetYaxis()->SetTitle("(m)^{2} (GeV^{2})");

		// hist_y_m2[iTrackHisto] = new TH2F("hist_y_m2","y vs m^{2}",100,0,3,1500,mass[iTrackHisto]*mass[iTrackHisto]*0.7,mass[iTrackHisto]*mass[iTrackHisto]*1.3);
		// hist_y_m2[iTrackHisto]->GetXaxis()->SetTitle("y");
		// hist_y_m2[iTrackHisto]->GetYaxis()->SetTitle("(m)^{2} (GeV^{2})");

		// profile_pt_m2[iTrackHisto] = new TProfile("profile_pt_m2","P_{T} vs m^{2}",100,0,10);
		// profile_pt_m2[iTrackHisto]->GetXaxis()->SetTitle("P_{T} (GeV)");
		// profile_pt_m2[iTrackHisto]->GetYaxis()->SetTitle("(m)^{2} (GeV^{2})");	

		// profile_y_m2[iTrackHisto] = new TProfile("profile_y_m2","y vs m^{2}",100,0,3);
		// profile_y_m2[iTrackHisto]->GetXaxis()->SetTitle("y");
		// profile_y_m2[iTrackHisto]->GetYaxis()->SetTitle("(m)^{2} (GeV^{2})");
	}
  

	//tree
	PID_tree = new TTree("PID_tree","PID_tree");
	PID_tree->Branch("event_number",&event_number,"event_number/I");
	PID_tree->Branch("nHitsDedx",&nHitsDedx,"nHitsDedx/I");
	PID_tree->Branch("nHitsFit",&nHitsFit,"nHitsFit/I");
	PID_tree->Branch("vx_tpc",&vx_tpc,"vx_tpc/F");
	PID_tree->Branch("vy_tpc",&vy_tpc,"vy_tpc/F");
	PID_tree->Branch("vz_tpc",&vz_tpc,"vz_tpc/F");
	PID_tree->Branch("vz_tpc_error",&vz_tpc_error,"vz_tpc_error/F");
	PID_tree->Branch("vz_vpd",&vz_vpd,"vz_vpd/F");
	// PID_tree->Branch("pt",&pt,"pt/F");
	PID_tree->Branch("px",&px,"px/F");
	PID_tree->Branch("py",&py,"py/F");
	PID_tree->Branch("pz",&pz,"pz/F");
	PID_tree->Branch("dca_x",&dca_x,"dca_x/F");
	PID_tree->Branch("dca_y",&dca_y,"dca_y/F");
	PID_tree->Branch("dca_z",&dca_z,"dca_z/F");
	PID_tree->Branch("q",&q,"q/I");
	PID_tree->Branch("m2tof",&m2tof,"m2tof/F");
	PID_tree->Branch("tof",&tof,"tof/F");
	// PID_tree->Branch("gb2_prim",&gb2_prim,"gb2_prim/F");
	// PID_tree->Branch("gb2_dca",&gb2_dca,"gb2_dca/F");
	PID_tree->Branch("btofYLocal",&btofYLocal,"btofYLocal/F");
	PID_tree->Branch("btofZLocal",&btofZLocal,"btofZLocal/F");
	PID_tree->Branch("nsigma_He3",&nsigma_He3,"nsigma_He3/F");
	PID_tree->Branch("nsigma_He4",&nsigma_He4,"nsigma_He4/F");
	PID_tree->Branch("nsigma_tri",&nsigma_tri,"nsigma_tri/F");
	// PID_tree->Branch("nsigma_He3_old",&nsigma_He3_old,"nsigma_He3_old/F");
	// PID_tree->Branch("nsigma_He4_old",&nsigma_He4_old,"nsigma_He4_old/F");
	// PID_tree->Branch("nsigma_tri_old",&nsigma_tri_old,"nsigma_tri_old/F");
	PID_tree->Branch("magn",&magn,"magn/F");
	PID_tree->Branch("rescale",&rescale,"rescale/F");  
	PID_tree->Branch("dedx",&dedx,"dedx/F");
	PID_tree->Branch("dedx_error",&dedx_error,"dedx_error/F");
	
    cout << "End of Histograms" << endl;
    return kStOK;
}
//__________________________________________________________________________________
void Shift::Clear(Option_t *opt)
{
    StMaker::Clear();
}

//__________________________________________________________________________________
Int_t Shift::Finish() {
    cout << "Shift::Finish()\n";
    //===============================
    //  Write Histograms
    //===============================

    File->cd();
	// PID_tree->Write();
	
	hist_VyVx->Write();
	hist_Vz->Write();
	hist_vpd_vz->Write();
	// hist_vz_vpd_tpc->Write();
	
	hist_dEdx->Write();
	hist_dEdx_w_tof->Write();

	hist_tof_pid->Write();
	
	// hist_phi_vs_eta->Write();
	// hist_phi_vs_eta_TOF->Write();
	// hist_eta_vs_pt->Write();
	
	for(int i =0;i<10;i++)
	{
		// mass2_pt_proton_corrected_eta[i]->Write();
		mass2_pt_proton_uncorrected_eta[i]->Write();
		// mass2_pt_pion_uncorrected_eta[i]->Write();
		// mass2_pt_kaon_uncorrected_eta[i]->Write();
		mass2_inv_pt_proton_uncorrected_eta[i]->Write();
		mass2_phi_proton_eta[i]->Write();
	}
	for(int i =0;i<14;i++)
	{
		// mass2_pt_proton_corrected_vz[i]->Write();
		mass2_pt_proton_uncorrected_vz[i]->Write();
		// mass2_pt_pion_uncorrected_vz[i]->Write();
		// mass2_pt_kaon_uncorrected_vz[i]->Write();
	}
	
	for(int i =0;i<40;i++)
	{
		// mass2_pt_proton_corrected_mean_vz[i]->Write();
		mass2_pt_proton_uncorrected_mean_vz[i]->Write();
		mass2_pt_pion_uncorrected_mean_vz[i]->Write();
		mass2_pt_kaon_uncorrected_mean_vz[i]->Write();
		mass2_pt_proton_uncorrected_mean_vz_no_membrane[i]->Write();
	}
	for(int i =0;i<12;i++)
	{
		mass2_phi_proton[i]->Write();
		mass2_phi_pion[i]->Write();
		mass2_phi_kaon[i]->Write();
		mass2_phi_proton_eta_pos[i]->Write();
		mass2_phi_proton_eta_neg[i]->Write();
	}
	// for(int iTrackHisto=0; iTrackHisto<NTrackHistoFolders; iTrackHisto++)
	// {
		// fTrackPdgToHistoIndex[ pdgTrackHisto[iTrackHisto] ] = iTrackHisto;
		// fHistodEdXTracks[iTrackHisto]->Write();
		// fHistodEdXwithToFTracks[iTrackHisto]->Write();
		// fHistoTofPIDTracks[iTrackHisto]->Write();
		// fHistoMomentumTracks[iTrackHisto]->Write();
		// fHistodEdXPull[iTrackHisto]->Write();
		// fHistodEdXZ[iTrackHisto]->Write();
		// fHistodEdX_of_tofpid[iTrackHisto]->Write();
		// fHistotof_of_dEdXpid[iTrackHisto]->Write();
		// fHisto_accuracy[iTrackHisto]->Write();
		// fHisto_mass2[iTrackHisto]->Write();
		// fh1d_mass2[iTrackHisto]->Write();
		// fHisto_particle_acceptance[iTrackHisto]->Write();
	// }
	
    return kStOK;
}
//__________________________________________________________________________________
Int_t Shift::GetRunIndex( const Int_t run ) {return -1;}
//__________________________________________________________________________________
Int_t Shift::Make() {
    //Begining of Event loop  
    //------------------------------------------------------------------
    if(!mPicoDstMaker) {LOG_WARN << " No PicoDstMaker! Skip! " << endm;return kStWarn;}
    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {LOG_WARN << " No PicoDst! Skip! " << endm;return kStWarn;}
    picoEvent = (StPicoEvent*)mPicoDst->event();
    if( !picoEvent ){LOG_WARN << " No PicoEvent! Skip! " << endm; return kStWarn;}
    //------------------------------------------------------------------
    TVector3 pVertex = picoEvent->primaryVertex();
	TVector3 picoPVError = picoEvent->primaryVertexError();
	vx_tpc=pVertex.X();
	vy_tpc=pVertex.Y();
	vz_tpc =pVertex.Z();
	double vr= sqrt(pVertex.X()*pVertex.X()+pVertex.Y()*pVertex.Y());
	
	vz_tpc_error= picoPVError.z();
	StThreeVectorF pv(pVertex.x(),pVertex.y(),pVertex.x());
	Int_t vz_sign = 0; if(vz_tpc<0){vz_sign = 1;}
	magn = picoEvent->bField();
	vz_vpd = picoEvent->vzVpd();
	double zdcx = picoEvent->ZDCx();
	double bbcx = picoEvent->BBCx();
	hist_VyVx->Fill(vx_tpc,vy_tpc);
	//if(!isGoodTrigger(picoEvent)){return 0;}
    if(!isGoodEvent(picoEvent)) {return 0;}
    int run_number = picoEvent->runId();
    event_number = picoEvent->eventId();
    // int runindex = GetRunIndex(run_number);
	mRefMultCorrUtil->init(run_number);

    //------ remove bad run---------
	if(mRefMultCorrUtil->isBadRun( run_number )) {return kStOk;}
	
	//vz histograms:
	hist_Vz->Fill(vz_tpc);
	hist_vpd_vz->Fill(vz_vpd);
	hist_vz_vpd_tpc->Fill(vz_vpd,vz_tpc);
	
	const Int_t nTrack = mPicoDst->numberOfTracks();
	
	int vz_bin =abs((vz_tpc+70)/10);
	
	////////////////////////////////////////////////////PID loop
	for (Int_t itr=0;itr<nTrack;itr++)
	{
        const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);

        if(!ptrk)  continue;
        if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
        if(!isGoodTrack(ptrk))  continue;
		// Get PID parameters
		px = ptrk->pMom().X(); py = ptrk->pMom().Y(); pz = ptrk->pMom().Z();
		pt = ptrk->pPt();
		dca_x = ptrk->gDCA( pVertex ).X(); dca_y = ptrk->gDCA( pVertex ).Y(); dca_z = ptrk->gDCA( pVertex ).Z();
		StThreeVectorF track_dca(ptrk->gDCA( pVertex ).X(),ptrk->gDCA( pVertex ).Y(),ptrk->gDCA( pVertex ).Z());
		double trackP = ptrk->pPtot();
		StPicoPhysicalHelix helix       = ptrk->helix(magn);
		double eta = ptrk->pMom().Eta();
		double phi = ptrk->pMom().Phi(); if(phi < 0.0) phi += 2.0*pi; if(phi > 2.0*pi) phi -= 2.0*pi;
		q = ptrk->charge();
		dedx=ptrk->dEdx();
		dedx_error = ptrk->dEdxError();
		
		//calculate mean z:
		double exit_point = 0;
		bool membrane_passage = 0;
		if((-(vz_tpc+dca_z)+200)*pt/pz>0){exit_point = (-(vz_tpc+dca_z)+200)*pt/pz;}
		if((-(vz_tpc+dca_z)-200)*pt/pz>0){exit_point = (-(vz_tpc+dca_z)-200)*pt/pz;}
		if(exit_point > 200){exit_point = 200;}
		double max_z = vz_tpc +exit_point*(pz/pt);
		double mean_z = (vz_tpc +max_z)/2.0;
		if((vz_tpc<0 && 0<max_z) || (max_z<0&& 0<vz_tpc)){membrane_passage = 1;}
		int mean_z_bin =abs((mean_z+200)/10);
		hist_dEdx->Fill(trackP/q,dedx);
		hist_phi_vs_eta->Fill(eta,phi);
		hist_eta_vs_pt->Fill(pt/q,eta);
		
		rescale = 1.0;
		
		//////////////////get Tof info:
		m2tof = -1.e6;
		tof = -1.e6;
		gb2_prim = -1.e6;
		gb2_dca = -1.e6;
		btofYLocal= -1.e6;
		btofZLocal= -1.e6;
		double gb2 = -1.e6;
		bool isTofm2 = false;
		double C = 29.97924580; //speed of light in centimeters/nanosecond
		double betaTof =- 1.e6;
		if(ptrk->bTofPidTraitsIndex() > 0)
		{
			//cout<<"what about here?"<<endl;
			const StPicoBTofPidTraits* btofPid = mPicoDst->btofPidTraits(ptrk->bTofPidTraitsIndex());
			double betaTof2 = btofPid->btofBeta() * btofPid->btofBeta();
			betaTof = btofPid->btofBeta();
			TVector3 tbtofHisPos = btofPid->btofHitPos();  

			//StThreeVectorF btofHitPos = tofPid->btofHitPos();
			StThreeVectorF btofHitPos(tbtofHisPos.X(),tbtofHisPos.Y(),tbtofHisPos.Z());
			float L_prim = tofPathLength(&pv,&btofHitPos,helix.curvature());
			StThreeVectorF dca_vertex = pv+track_dca;
			float L_dca = tofPathLength(&dca_vertex,&btofHitPos,helix.curvature());
			tof = btofPid->btof();

			gb2_prim =L_prim*L_prim/(tof*tof*C*C-L_prim*L_prim);
			gb2_dca = L_dca*L_dca/(tof*tof*C*C-L_dca*L_dca);
			//cout<<betaTof2<<endl;
			gb2 =betaTof2/(1.0-betaTof2);
			if(fabs(betaTof2) > 1.e-6 && gb2>0.0)
			{
				m2tof = trackP*trackP/gb2;
				btofYLocal= btofPid->btofYLocal();
				btofZLocal= btofPid->btofZLocal();
				isTofm2 = true;
			}
			hist_dEdx_w_tof->Fill(trackP/q, ptrk->dEdx());
			hist_tof_pid->Fill(trackP/q, 1.0/(btofPid->btofBeta()));
		}
		
		double dEdXPull[NTrackHistoFolders] = {0};
		for(int iTrackHisto=0; iTrackHisto<NTrackHistoFolders; iTrackHisto++)
		{
			if(fabs(charge[iTrackHisto])==1.0){dEdXPull[iTrackHisto] = TMath::Log(dedx/TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(trackP/mass[iTrackHisto]),fdEdXMode)))/dedx_error;}
			if(fabs(charge[iTrackHisto])==2.0){dEdXPull[iTrackHisto] = TMath::Log(dedx/(4.0*Bichsel::Instance()->GetI70M(TMath::Log10(2.0*trackP/mass[iTrackHisto]),3.0)))/dedx_error;}
		}
		
		int pt_bin =0;
		if(pt<0.75){pt_bin = 0;}
		else if(pt<0.875){pt_bin = 1;}
		else if(pt<1.0){pt_bin = 2;}
		else if(pt<1.125){pt_bin = 3;}
		else if(pt<1.25){pt_bin = 4;}
		else if(pt<1.375){pt_bin = 5;}
		else if(pt<1.5){pt_bin = 6;}
		else if(pt<1.75){pt_bin = 7;}
		else if(pt<2.0){pt_bin = 8;}
		else if(pt<2.5){pt_bin = 9;}
		else if(pt<3.0){pt_bin = 10;}
		else {pt_bin = 11;}
		
		int eta_bin = abs((eta+1.0)*5.0);
		if(eta_bin<0){eta_bin = 0;} if(eta_bin>9){eta_bin = 9;}
		if(fabs(dEdXPull[4])<3.0&& gb2 !=-1.e6)//pion
		{
			mass2_pt_pion_uncorrected_eta[eta_bin]->Fill(pt/q,m2tof);
			mass2_pt_pion_uncorrected_vz[vz_bin]->Fill(pt/q,m2tof);
			mass2_pt_pion_uncorrected_mean_vz[mean_z_bin]->Fill(pt/q,m2tof);
			mass2_phi_pion[pt_bin]->Fill(phi,m2tof/q);
		}
		if(fabs(dEdXPull[6])<3.0&& gb2 !=-1.e6)//kaon
		{
			mass2_pt_kaon_uncorrected_eta[eta_bin]->Fill(pt/q,m2tof);
			mass2_pt_kaon_uncorrected_vz[vz_bin]->Fill(pt/q,m2tof);
			mass2_pt_kaon_uncorrected_mean_vz[mean_z_bin]->Fill(pt/q,m2tof);
			mass2_phi_kaon[pt_bin]->Fill(phi,m2tof/q);
		}
		if(fabs(dEdXPull[8])<3.0&& gb2 !=-1.e6)//proton
		{
			mass2_pt_proton_corrected_eta[eta_bin]->Fill(pt/q,m2tof);
			mass2_pt_proton_corrected_vz[vz_bin]->Fill(pt/q,m2tof);
			mass2_pt_proton_uncorrected_eta[eta_bin]->Fill(pt/q,m2tof);
			mass2_inv_pt_proton_uncorrected_eta[eta_bin]->Fill((1.0*q)/pt,m2tof);
			mass2_pt_proton_uncorrected_vz[vz_bin]->Fill(pt/q,m2tof);
			mass2_pt_proton_uncorrected_mean_vz[mean_z_bin]->Fill(pt/q,m2tof);
			if(!membrane_passage){mass2_pt_proton_uncorrected_mean_vz_no_membrane[mean_z_bin]->Fill(pt/q,m2tof);}
			mass2_phi_proton[pt_bin]->Fill(phi,m2tof/q);
			mass2_phi_proton_eta[eta_bin]->Fill(phi,m2tof/q);
			if(eta>0){mass2_phi_proton_eta_pos[pt_bin]->Fill(phi,m2tof/q);}
			if(eta<0){mass2_phi_proton_eta_neg[pt_bin]->Fill(phi,m2tof/q);}
			

		}
		
		if(fabs(dEdXPull[12])<4.0||fabs(dEdXPull[14])<4.0||fabs(dEdXPull[16])<4.0)
		{
			nsigma_tri=dEdXPull[12];
			nsigma_He3=dEdXPull[14];
			nsigma_He4=dEdXPull[16]; 

			PID_tree->Fill();
		}
		for(int iTrackHisto=0; iTrackHisto<NTrackHistoFolders; iTrackHisto++)
		{
			double mdq2 = (mass[iTrackHisto]/charge[iTrackHisto])*(mass[iTrackHisto]/charge[iTrackHisto]);
			if(fabs(dEdXPull[iTrackHisto])<3.0 &&(m2tof == -1.e6 || fabs(mdq2-m2tof)/mdq2<0.2))
			{
				fHistodEdXwithToFTracks[iTrackHisto]->Fill(trackP/q,dedx);
				fHistoMomentumTracks[iTrackHisto]->Fill(trackP);
				fHistodEdXPull[iTrackHisto]->Fill(trackP,dEdXPull[iTrackHisto]);
				if(m2tof!=-1.e6)
				{
					fHisto_accuracy[iTrackHisto]->Fill(dEdXPull[iTrackHisto],m2tof);
					fHisto_mass2[iTrackHisto]->Fill(trackP,m2tof);
				}
			}
			if(fabs(dEdXPull[iTrackHisto])<3.0)
			{
				fHistoTofPIDTracks[iTrackHisto]->Fill(trackP/q,dedx);
				if(m2tof!=-1.e6)
				{
					fHistotof_of_dEdXpid[iTrackHisto]->Fill(trackP,1.0/betaTof);
				}
			}
			if(fabs(mdq2-m2tof)/mdq2<0.2)
			{
				fHistodEdX_of_tofpid[iTrackHisto]->Fill(trackP,dedx);
			}
		}
	}
    return kStOK;
}
//end loop

//__________________________________________________________________________________
bool Shift::isGoodEvent(const StPicoEvent *event)
{
    Float_t vx=event->primaryVertex().X();
    Float_t vy=event->primaryVertex().Y();
    Float_t vz=event->primaryVertex().Z();
	if(vz < -70 || vz > 70) return false;
	if(( (vx)*(vx) + (vy)*(vy) ) > 4) return false;
    return true;
}

bool Shift::isGoodTrack(const StPicoTrack *ptrk) {
	const Float_t dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
	nHitsFit = ptrk->nHitsFit();
    const Int_t nHitsPoss = ptrk->nHitsMax();
    nHitsDedx = ptrk->nHitsDedx();
    const Float_t quality = (Float_t)nHitsFit/(Float_t)nHitsPoss;
	
	// if( fabs(dca)>3.0 ) return false;
	if( nHitsDedx < 5 ) return false;
	if( nHitsFit < 15 )  return false; //note, for excess_proton analysis, i use 10
	if( quality < 0.52 )  return false;
	return true;
}
Double_t Shift::Pileup_rejection(Double_t gRefmult, Int_t gNtofMatch, Double_t gvz,const StPicoEvent *event)
{
	int vz_bin =abs((gvz+145)/10); if(vz_bin<0||vz_bin>28){return -1.0;}
	double b0=10000000,b1=0,b2=0,b3=0,b4=0;
	double d0=0,d1=0,d2=0,d3=0,d4=0;
	double vz_corr_7p7[29]={1.00654,1.00053,0.998852,0.996999,0.996159,0.996542,0.996349,1.00047,1.00468,1.00304,0.99581,0.9983,0.999919,0.999749,
							1,1.00118,0.9992,0.999113,0.995557,1.00442,1.00548,0.999452,0.996127,0.996413,0.996018,0.996163,0.998376,0.998865,1.00318};
	double vz_corr_9p2_trig1[29]={
				0.981917, 0.981996, 0.980511, 0.97507, 0.974313, 0.973339, 0.964442, 0.975903, 0.980922, 0.982125, 0.97381, 0.979041, 0.987953, 0.994754, 1.0,
				1.00207,  1.01171,  1.01019,  1.01176, 1.02353,  1.02679,  1.02377,  1.02388,  1.02206,  1.03417,  1.03343, 1.04365,  1.05574,  1.07789};
	double vz_corr_9p2_trig2[29]={
				1.00894,  1.00431,  1.00199,  1.00086,  1.00139, 1.00079,  1.00023,  1.00587,  1.01,     1.00751,  0.999888, 1.00118,  1.00301, 1.00237, 1.0,
				0.997929, 0.995314, 0.990985, 0.98731, 0.993973, 0.994798, 0.989476, 0.985773, 0.986552, 0.986979, 0.989481, 0.993353, 1.00088, 1.01531};
	gRefmult+=1.0*r1->Uniform()-0.5;
	if(dataset == "7p7_Gev_2021")
	{
		gRefmult = vz_corr_7p7[vz_bin]*gRefmult;
		if(gvz<-87.0)// -145,-87 cm
		{
			b0 = 39.578630496797,b1 = 1.46561577132993,b2 = 0.006515367058115,b3 = -4.06391982010589e-05,b4 = 5.51203917383809e-08;
			d0 = -14.8817460248614,d1 = 0.764539480062978,d2 = 0.00368901349656326,d3 = -1.27602217700865e-05,d4 = 8.02618485000158e-10;
		}
		else if(gvz<-29.0)// -87,-29 cm
		{
			b0 = 26.1841414192908,b1 = 1.73354655107464,b2 = -0.00280668326418846,b3 = 1.22370803379957e-05,b4 = -3.15068617200212e-08;
			d0 = -13.1831127837376,d1 = 0.760227210117286,d2 = 0.00195873375843822,d3 = -2.69378951644624e-06,d4 = -1.05344843941749e-08;
		}
		else if(gvz<29.0)// -29,29 cm
		{
			b0 = 23.3635904884101,b1 = 1.58179764458174,b2 = -0.00100184372825271,b3 = 7.76378744751984e-07,b4 = -6.46469867000365e-09;
			d0 = -11.4340781454132,d1 = 0.72398407747444,d2 = 0.00121092416745035,d3 = 1.17875404059176e-07,d4 = -9.81658682040738e-09;
		}
		else if(gvz<87.0)// 29,87 cm
		{
			b0 = 29.4343991835005,b1 = 1.48353715105631,b2 = 0.00106271734149745,b3 = -9.07835076338586e-06,b4 = 6.7722581625238e-09;
			d0 = -9.97159163811459,d1 = 0.591000613390771,d2 = 0.00449768928484714,d3 = -1.71667412152202e-05,d4 = 1.6467383813372e-08;
		}
		else// 87,145 cm
		{
			b0 = 37.0772875081557,b1 = 1.53484162926915,b2 = 0.00471873506675937,b3 = -2.94958548877277e-05,b4 = 3.60887574265838e-08;
			d0 = -13.3927733032856,d1 = 0.704319390196747,d2 = 0.00485360248820988,d3 = -2.10416804123978e-05,d4 = 1.92342533435503e-08;
		}
	}
	else if(dataset =="9p2_GeV_2020")
	{
		if(event->isTrigger(780010)){gRefmult = vz_corr_9p2_trig1[vz_bin]*gRefmult;}
		else if(event->isTrigger(780020)){gRefmult = vz_corr_9p2_trig2[vz_bin]*gRefmult;}
		if(gvz<-87.0)// -145,-87 cm
		{
			b0=25.6055790979197, b1=2.02528136596901, b2=-0.0058370984051939, b3=2.59602314466234e-05, b4=-5.3014743584261e-08;
			d0=-17.7059596791057, d1=0.614538168662738, d2=0.00534180935164814, d3=-1.79582873880806e-05, d4=1.01623054170579e-08;
		}
		else if(gvz<-29.0)// -87,-29 cm
		{
			b0=23.0160060308621, b1=1.61885832757588, b2=-0.00275873189631398, b3=1.31262550392554e-05, b4=-2.94368020941846e-08;
			d0=-17.3591842617911, d1=0.796170989774258, d2=0.000670722514533827, d3=3.26258075150876e-06, d4=-1.60611460182112e-08;
		}
		else if(gvz<29.0)// -29,29 cm
		{
			b0=16.4277056306649, b1=1.71652229539398, b2=-0.00406847684302521, b3=1.65203560938885e-05, b4=-2.96250329214512e-08;
			d0=-15.7887025834219, d1=0.789786364309292, d2=-0.000637115144252616, d3=1.00019972792727e-05, d4=-2.45208851616324e-08;
		}
		else if(gvz<87.0)// 29,87 cm
		{
			b0=21.2024767158778, b1=1.70521848381614, b2=-0.00352260930859763, b3=1.60905730948817e-05, b4=-3.37443468806432e-08;
			d0=-17.1166088395929, d1=0.814739436616432, d2=0.000227197779215977, d3=6.55397838050604e-06, d4=-2.28812912596058e-08;
		}
		else// 87,145 cm
		{
			b0=26.0970905882739, b1=1.88889714311734, b2=-0.00195374948885512, b3=-6.14244087431038e-06, b4=1.99930095058841e-08;
			d0=-15.6624325989392, d1=0.52385751891358, d2=0.00794996911844969, d3=-4.09239155250494e-05, d4=6.40163739983216e-08;
		}
	}
	else if(dataset =="14p6_gev_2019") //https://drupal.star.bnl.gov/STAR/system/files/Centrality_Study_at_14p6_final_0.pdf
	{
		b0 = 36.4811873257854, b1 = 1.96363692967013, b2 = -0.00491528146300182, b3 = 1.45179464078414e-05, b4 = -1.82634741809226e-08;
		d0 = -16.176117733536, d1 = 0.780745107634961, d2 = -2.03347057620351e-05, d3 = 3.80646723724747e-06,d4 = -9.43403282145648e-09;
	}
	if(gRefmult<d0+
			  d1*gNtofMatch+
			  d2*gNtofMatch*gNtofMatch+
			  d3*gNtofMatch*gNtofMatch*gNtofMatch+
			  d4*gNtofMatch*gNtofMatch*gNtofMatch*gNtofMatch){return -1;}
	if(gRefmult>b0+
			  b1*gNtofMatch+
			  b2*gNtofMatch*gNtofMatch+
			  b3*gNtofMatch*gNtofMatch*gNtofMatch+
			  b4*gNtofMatch*gNtofMatch*gNtofMatch*gNtofMatch){return -1;}
	return gRefmult;
}
//bool Shift::isGoodTrackSooraj(const StPicoTrack *ptrk) {return true;}
bool Shift::isGoodTrigger(const StPicoEvent* event)
{
	if(dataset == "9p2_GeV_2020")
	{
		if(event->isTrigger(780010)|| event->isTrigger(780020)) return true; //9p2GeV minbias only
	}
	else if(dataset == "7p7_Gev_2021")
	{
		if(event->isTrigger(810010)|| event->isTrigger(810020)|| event->isTrigger(810030)|| event->isTrigger(810040)) return true; //7p7GeV minbias only
	}
	else if(dataset == "11p5_gev_2020")
	{
		if(event->isTrigger(710000)|| event->isTrigger(710010)|| event->isTrigger(710020)) return true; //11p5 GeV minbias only
	}
	else if(dataset == "14p6_gev_2019")
	{
		if(event->isTrigger(650000)) return true; //minbias only
	}
	else if(dataset == "17p3_gev_2021")
	{
		if(event->isTrigger(870010)) return true; //minbias only
	}
	else if(dataset == "19p6_Gev_2019")
	{
		if(event->isTrigger(640001) || event->isTrigger(640011) || event->isTrigger(640021) || event->isTrigger(640031) || event->isTrigger(640041) || event->isTrigger(640051)) return true;
	}
	else if(dataset == "27_Gev_2018")
	{
		if(event->isTrigger(1)||event->isTrigger(610001)||event->isTrigger(610011)||event->isTrigger(610021)||event->isTrigger(610031) ||event->isTrigger(610041)||event->isTrigger(610051)) return true;//27 gev
	}
	else if(dataset == "4p5_gev_FXT_2020")
	{
		if(event->isTrigger(740000)|| event->isTrigger(740010)) return true; //epde-or-bbce-or-vpde-tof1 trigger, since there is no minbias trigger
		return true;
	}
	return false;
}
bool Shift::isProton(const StPicoTrack *ptrk)
{

	if(fabs(ptrk->nSigmaProton())<3.0) return true;
	return false;
}
bool Shift::isPion(const StPicoTrack *ptrk)
{
	double trackP =ptrk->pPtot();
	double Beta = -1.0;
	// TMath::Abs(ptrk->nSigmaPion())<3.0
	if(ptrk->isTofTrack())
	{	
		StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits(ptrk->bTofPidTraitsIndex());
		if(trait){Beta = trait->btofBeta();}
		return TMath::Abs(TMath::Sqrt((0.13957/trackP)*(0.13957/trackP)+1.0)-(1.0/Beta))<0.03 && trackP<1.0;
	}
	return false;
}
int Shift::EP_group(int run_number, int cent)
{
	int ep_group_num = 0;
	if(dataset == "27_Gev_2018")
	{
		int event_cuts[14]= {19131037, 19135016, 19137041, 19139063, 19140030, 19141030, 19144012, 19144033,
						   19145034, 19147021, 19147048, 19155057, 19158020, 19268002};
		for(int i = 0;i<14;i++)
		{
			if(run_number<=event_cuts[13-i]){ ep_group_num=(13-i)*9+cent;}
		}
	}
	else
	{
		ep_group_num = GetRunIndex(run_number)*9+cent;
	}
	return ep_group_num;
}
		
 

//__________________________________________________________________________________
