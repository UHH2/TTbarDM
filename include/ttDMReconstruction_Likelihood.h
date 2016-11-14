#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"
#include "TFile.h"
#include "TH3F.h"

class ttDMReconstruction_Likelihood: public uhh2::AnalysisModule {
public:
   
   explicit ttDMReconstruction_Likelihood(uhh2::Context & ctx);
   
   virtual bool process(uhh2::Event & event);

private:
  
   int NFitPar;
   Particle lepton;
   uhh2::Event::Handle<double> h_likelihood;
   uhh2::Event::Handle<LorentzVector> h_recneutrino;
   uhh2::Event::Handle<TTbarGen>  h_ttbargen;
   uhh2::Event::Handle<Jet> h_bjet;
   uhh2::Event::Handle<std::vector<Muon>> h_muons;
   
   
   TFile *f;
   TH3F *h;
   };

