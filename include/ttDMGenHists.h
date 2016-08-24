#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/common/include/TTbarGen.h"

#include <string>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>

class ttDMGenHists: public uhh2::Hists {
 public:
  ttDMGenHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  
 private:
   double weight;

   TH1F* top_pt;
   TH1F* top_eta;
   TH1F* top_phi;
   TH1F* antitop_pt;
   TH1F* antitop_eta;
   TH1F* antitop_phi;
   TH1F* bothtops_pt;
   TH1F* bothtops_eta;
   TH1F* bothtops_phi;
   TH1F* neutrino_top_pt;
   TH1F* neutrino_top_eta;
   TH1F* neutrino_top_phi;
   TH1F* neutrino_antitop_pt;
   TH1F* neutrino_antitop_eta;
   TH1F* neutrino_antitop_phi;
   TH1F* neutrino_pt;
   TH1F* neutrino_eta;
   TH1F* neutrino_phi;
   TH1F* neutrino_px;
   TH1F* neutrino_py;
   TH1F* neutrino_pz;
   TH1F* DM_pt;
   TH1F* DM_eta;
   TH1F* DM_phi;
   TH1F* DeltaR_top_DM;
   TH1F* DeltaR_antitop_DM;
   TH1F* DeltaR_neutrino_DM;
   TH1F* DeltaR_top_neutrino;
   TH1F* DeltaR_bothtops_neutrino;
   TH1F* DeltaR_antitop_neutrino;
   TH1F* InvMass_b_lep;
   TH1F* InvMass_neutrino_lep;
   TH1F* InvMass_W_closestJet;
   TH1F* lepton_pt;
   TH1F* lepton_eta;
   TH1F* lepton_phi;
   TH1F* DeltaR_lepton_neutrino;
   TH1F* Deltaeta_lepton_neutrino;
   TH1F* Deltaphi_lepton_neutrino;
   TH1F* Deltaeta_lepton_neutrino_2;
   TH1F* Deltaphi_lepton_neutrino_2;
   TH1F* DeltaPx_lepton_neutrino;
   TH1F* DeltaPy_lepton_neutrino;
   TH1F* DeltaPz_lepton_neutrino;
   TH1F* DeltaPt_lepton_neutrino;
   TH1F* DeltaP_lepton_neutrino;
   TH1F* DeltaE_lepton_neutrino;
   TH1F* Deltax_lepton_neutrino;
   TH1F* Deltay_lepton_neutrino;
   TH1F* Deltaz_lepton_neutrino;
   TH1F* DeltaPx_lepton_neutrino_2;
   TH1F* DeltaPy_lepton_neutrino_2;
   TH1F* DeltaPz_lepton_neutrino_2;
   TH1F* DeltaPt_lepton_neutrino_2;
   TH1F* Deltax_lepton_neutrino_2;
   TH1F* Deltay_lepton_neutrino_2;
   TH1F* Deltaz_lepton_neutrino_2;
   TH1F* Resolution_MET;
   TH1F* Gen_MET;

   // TH3F* h_likelihood;
   // TH3F* h_likelihood_bin1;
   // TH3F* h_likelihood_bin2;
   // TH3F* h_likelihood_bin2_noweights;
 
   TH2F* h_dphi_deta;
   TH2F* h_dW_deta;
   TH2F* h_dW_dphi;
   TH1F* DeltaR_lepton_genb;
   uhh2::Event::Handle<TTbarGen> h_ttbargen;
   
   TH1F* DM_mass;
   TH1F* DMparticles_mass;
   TH1F* mediator_mass;

   TH1F* MET_Neutrino_px;
   TH1F* MET_Neutrino_py;
};
