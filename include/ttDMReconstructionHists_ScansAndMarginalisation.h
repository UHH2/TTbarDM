#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "TF1.h"

class ttDMReconstructionHists_ScansAndMarginalisation: public uhh2::Hists {
 public:
   ttDMReconstructionHists_ScansAndMarginalisation(uhh2::Context & ctx, const std::string & dirname);
   
   virtual void fill(const uhh2::Event & ev) override;
   Double_t MWParametrization(Double_t mw);
   Double_t DPhiParametrization(Double_t dphi);
   Double_t DEtaParametrization(Double_t dphi);
   
 private:
   double weight;
   uhh2::Event::Handle<TTbarGen> h_ttbargen;
   int k =0;

   TH2F* h_px_py_likelihood_1;
   TH2F* h_px_pz_likelihood_1;
   TH2F* h_py_pz_likelihood_1;
   TH2F* h_px_py_likelihood_2;
   TH2F* h_px_pz_likelihood_2;
   TH2F* h_py_pz_likelihood_2;
   TH2F* h_px_py_likelihood_3;
   TH2F* h_px_pz_likelihood_3;
   TH2F* h_py_pz_likelihood_3;
   TH2F* h_px_py_likelihood_4;
   TH2F* h_px_pz_likelihood_4;
   TH2F* h_py_pz_likelihood_4;
   TH2F* h_px_py_likelihood_5;
   TH2F* h_px_pz_likelihood_5;
   TH2F* h_py_pz_likelihood_5;
   TH2F* h_px_py_likelihood_6;
   TH2F* h_px_pz_likelihood_6;
   TH2F* h_py_pz_likelihood_6;
   TH2F* h_px_py_likelihood_7;
   TH2F* h_px_pz_likelihood_7;
   TH2F* h_py_pz_likelihood_7;
   TH2F* h_px_py_likelihood_8;
   TH2F* h_px_pz_likelihood_8;
   TH2F* h_py_pz_likelihood_8;
   TH2F* h_px_py_likelihood_9;
   TH2F* h_px_pz_likelihood_9;
   TH2F* h_py_pz_likelihood_9;
   TH2F* h_px_py_likelihood_10;
   TH2F* h_px_pz_likelihood_10;
   TH2F* h_py_pz_likelihood_10;
   TH1F* neutrino_px_1; 
   TH1F* neutrino_py_1; 
   TH1F* neutrino_pz_1; 
   TH1F* neutrino_px_2; 
   TH1F* neutrino_py_2; 
   TH1F* neutrino_pz_2; 
   TH1F* neutrino_px_3; 
   TH1F* neutrino_py_3; 
   TH1F* neutrino_pz_3; 
   TH1F* neutrino_px_4; 
   TH1F* neutrino_py_4; 
   TH1F* neutrino_pz_4; 
   TH1F* neutrino_px_5; 
   TH1F* neutrino_py_5; 
   TH1F* neutrino_pz_5; 
   TH1F* neutrino_px_6; 
   TH1F* neutrino_py_6; 
   TH1F* neutrino_pz_6; 
   TH1F* neutrino_px_7; 
   TH1F* neutrino_py_7; 
   TH1F* neutrino_pz_7; 
   TH1F* neutrino_px_8; 
   TH1F* neutrino_py_8; 
   TH1F* neutrino_pz_8; 
   TH1F* neutrino_px_9; 
   TH1F* neutrino_py_9; 
   TH1F* neutrino_pz_9; 
   TH1F* neutrino_px_10;
   TH1F* neutrino_py_10;
   TH1F* neutrino_pz_10;

   TH1F* met_px_1; 
   TH1F* met_py_1; 
   TH1F* met_px_2; 
   TH1F* met_py_2; 
   TH1F* met_px_3; 
   TH1F* met_py_3; 
   TH1F* met_px_4; 
   TH1F* met_py_4; 
   TH1F* met_px_5; 
   TH1F* met_py_5; 
   TH1F* met_px_6; 
   TH1F* met_py_6; 
   TH1F* met_px_7; 
   TH1F* met_py_7; 
   TH1F* met_px_8; 
   TH1F* met_py_8; 
   TH1F* met_px_9; 
   TH1F* met_py_9; 
   TH1F* met_px_10;
   TH1F* met_py_10;

};
