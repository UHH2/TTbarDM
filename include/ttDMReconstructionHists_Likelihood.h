#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "TF1.h"

class ttDMReconstructionHists_Likelihood: public uhh2::Hists {
 public:
   ttDMReconstructionHists_Likelihood(uhh2::Context & ctx, const std::string & dirname);
   
   virtual void fill(const uhh2::Event & ev) override;
  
private:
   TH1F *hist_chi2;
   TH1F *hist_pxrec_pxgen;
   TH1F *hist_pyrec_pygen;
   TH1F *hist_pzrec_pzgen;
   TH1F *hist_DM_MET,  *hist_DM_MET_chi2, *hist_DM_MET_semilept, *hist_DM_MET_chi2_semilept;
   TH1F *hist_MET_Neutrino_px; 
   TH1F *hist_MET_Neutrino_py; 
   TH1F *hist_MET_Neutrino_DM_px;
   TH1F *hist_MET_Neutrino_DM_py;
   TH1F *hist_MET_RecNeutrino_px; 
   TH1F *hist_MET_RecNeutrino_py; 
   TH1F *hist_MET_RecNeutrino_px_abs; 
   TH1F *hist_MET_RecNeutrino_py_abs; 
   TH1F *hist_MET_RecNeutrino_DM_px;
   TH1F *hist_MET_RecNeutrino_DM_py;
   TH1F *hist_DM_MET_gen;
   
   uhh2::Event::Handle<double> h_likelihood;
   uhh2::Event::Handle<LorentzVector> h_recneutrino;
   uhh2::Event::Handle<TTbarGen> h_ttbargen;
   uhh2::Event::Handle<std::vector<TopJet>> h_heptopjets_WP3;
   uhh2::Event::Handle<Jet> h_bjets;

   TH1F *hist_DeltaR_genb_nextjet;
   TH1F *hist_DeltaR_genb_nextjet_lowchi2;
   TH1F *hist_DeltaR_genb_nextjet_highchi2;

   TH1F *hist_DMMET_GenDMpT;
   TH1F *hist_DMMETgen_GenDMpT;
   TH1F *hist_DMMET_GenDMpT_semilept;
   TH1F *hist_DMMETgen_GenDMpT_only;
   TH1F *hist_DMMET_GenDMpT_genmet;
   TH1F *hist_n_nu;

   TH1F *hist_pxrec_pxgen_100;
   TH1F *hist_pyrec_pygen_100;
   TH1F *hist_DMMET_GenDMpT_vec;
   TH1F *hist_DMMET_GenDMpT_px;
   TH1F *hist_DMMET_GenDMpT_py;

   TH1F *hist_neutrino_pT;
   TH1F *hist_neutrino_phi;
   TH1F *hist_neutrino_pT_gen;
   TH1F *hist_neutrino_phi_gen;
   TH1F *hist_neutrino_gen_pT;
   TH1F *hist_neutrino_gen_phi;
   TH1F *hist_neutrino_pT_120;
   TH1F *hist_neutrino_phi_120;
   TH1F *hist_neutrino_pT_120_220;
   TH1F *hist_neutrino_phi_120_220;
   TH1F *hist_neutrino_pT_220_320;
   TH1F *hist_neutrino_phi_220_320;
   TH1F *hist_neutrino_pT_320;
   TH1F *hist_neutrino_phi_320;

   TH2F *hist_DeltaX_DeltaY;
   TH2F *hist_pX_pY_gen;
   TH2F *hist_pX_pY_rec;

   TH1F *hist_neutrino_pT_met160_320;
   TH1F *hist_neutrino_phi_met160_320;
   TH1F *hist_neutrino_pT_met320_500;
   TH1F *hist_neutrino_phi_met320_500;
   TH1F *hist_neutrino_pT_500;
   TH1F *hist_neutrino_phi_500;
   
   TH1F *hist_neutrino_gen_pT_reweighted;
   TH1F *hist_neutrino_pT_120_reweighted;
   TH1F *hist_neutrino_pT_120_220_reweighted;
   TH1F *hist_neutrino_pT_220_320_reweighted;
   TH1F *hist_neutrino_pT_320_reweighted;

   TH2F *hist_2D_ptmiss_pT;
   TH2F *hist_2D_METx_px;
   TH2F *hist_2D_METy_py;

   TH1F *hist_chi2_METpT_components;
   TH1F *hist_pTmisspT;
   TH1F *hist_pTmisspT_fabs;
   TH1F *hist_pTmisspT_quotient;
   TH1F *hist_Metxpx;
   TH1F *hist_Metypy;
   TH1F *hist_Metxpx_quad;
   TH1F *hist_Metypy_quad;
   TH1F *hist_Metxpx_Metypy_quad;
      
   TH1F *hist_unc_px;
   TH1F *hist_unc_py;
   TH1F *hist_relunc_px;
   TH1F *hist_relunc_py;
   TH1F *hist_pTttbar;
   TH2F* hist_2d_DMMET_MET;
   TH2F* hist_2d_DMMET_MT2W;
   TH2F* hist_2d_DMMET_deltaphijetmet12;
   TH2F* hist_2d_DMMET_deltaphijet123;
   TH2F* hist_2d_DMMET_deltphilep;
   TH2F* hist_2d_DMMET_mtlep;

   TH1F *hist_deltaphi_tj_met;
   TH1F *hist_deltaphi_tj_lep;
   TH1F *hist_deltaphi_neutrino_lep;
   TH1F *hist_deltaphi_tj_neutrino;
   TH1F *hist_deltaphi_neutrino_met;
   TH1F *hist_deltaphi_thad_tlep;
   TH1F *hist_deltaphi_thad_MET;
   
   

};  

