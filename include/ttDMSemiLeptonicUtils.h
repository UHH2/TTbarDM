#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

class TopJetLeptonDeltaRCleaner : public uhh2::AnalysisModule {
 public:
  explicit TopJetLeptonDeltaRCleaner(float mindr=0.8): minDR_(mindr) {}
  virtual bool process(uhh2::Event&) override;

 private:
  float minDR_;
};
class HepTTLeptonDeltaRCleaner : public uhh2::AnalysisModule {
 public:
   explicit HepTTLeptonDeltaRCleaner(uhh2::Context & ctx, float mindr=1.5): minDR_(mindr){
      h_heptopjets_WP3=ctx.get_handle<std::vector<TopJet>>("h_heptopjets_WP3_wobtag");
   }
   virtual bool process(uhh2::Event&) override;

 private:
   float minDR_;
   TString collection_; 
   uhh2::Event::Handle<std::vector<TopJet>> h_heptopjets_WP3;
};

double CalculateMT2W(const uhh2::Event &);
