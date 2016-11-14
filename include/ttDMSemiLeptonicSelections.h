#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"

#include <string>
#include <vector>

namespace uhh2 {

  class HTlepCut : public Selection {
   public:
    explicit HTlepCut(float, float max_htlep=infinity);
    virtual bool passes(const Event&) override;

   private:
    float min_htlep_, max_htlep_;
  };
  /////

  class MTlepCut : public Selection {
  public:
    explicit MTlepCut(float min_mtlep, float max_mtlep=infinity);
    virtual bool passes(const Event &) override;

  private:
    float min_mtlep_, max_mtlep_;
  };
  /////

  class METCut : public Selection {
   public:
    explicit METCut(float, float max_met=infinity);
    virtual bool passes(const Event&) override;

   private:
    float min_met_, max_met_;
  };
  /////

  class NJetCut : public Selection {
   public:
    explicit NJetCut(int, int nmax=999, float ptmin=0., float etamax=infinity);
    virtual bool passes(const Event&) override;

   private:
    int nmin, nmax;
    float ptmin, etamax;
  };
  /////

  class METJetDPhiCut : public Selection {
   public:
    explicit METJetDPhiCut(float, int);
    virtual bool passes(const Event&) override;

   private:
    float min_dphi_, maxjetindex_;
  };
  /////

   class METLeptonDPhiCut : public Selection {
   public:
      explicit METLeptonDPhiCut(float min_dphi);
      virtual bool passes(const Event&) override;
      
   private:
      float min_dphi_;
   };
   /////

  class MT2WCut : public Selection {
   public:
    explicit MT2WCut(float);
    virtual bool passes(const Event&) override;

   private:
    float min_mt2w_;
  };
  /////

  class TwoDCut : public Selection {
   public:
     explicit TwoDCut(Context &ctx, float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {
       h_muons=ctx.get_handle<std::vector<Muon>>("h_muons_medium");
       h_electrons=ctx.get_handle<std::vector<Electron>>("h_electrons_tight");      
    }
    virtual bool passes(const Event&) override;

   private:
    float min_deltaR_, min_pTrel_;
     uhh2::Event::Handle<std::vector<Muon>> h_muons;
     uhh2::Event::Handle<std::vector<Electron>> h_electrons;
  };
  /////

  class TriangularCuts : public Selection {
   public:
    explicit TriangularCuts(float, float);
    virtual bool passes(const Event&) override;

   private:
    float a_, b_;
  };
  /////

  class TopTagEventSelection: public Selection {
   public:
    explicit TopTagEventSelection(const TopJetId& tjet_id=CMSTopTag(), float minDR_jet_ttag=1.2);
    virtual bool passes(const Event&) override;

   private:
    TopJetId topjetID_;
    std::unique_ptr<Selection> topjet1_sel_;
    float minDR_jet_toptag_;
  };
  /////

   class Type2TopTagEventSelection: public Selection {
   public:
      explicit Type2TopTagEventSelection(const TopJetId& tjet_id=CMSTopTag(), float minDR_jet_ttag=0.8, float maxDR_jet_ttag=2.0, float mininvmass_jet_ttag=140., float maxinvmass_jet_ttag=250.);
    virtual bool passes(const Event&) override;

   private:
      TopJetId topjetID_;
      std::unique_ptr<Selection> topjet1_sel_;
      float minDR_jet_toptag_;
      float maxDR_jet_toptag_;  
      float mininvmass_jet_ttag_;
      float maxinvmass_jet_ttag_;
   };
  /////


  class LeptonicTopPtCut: public Selection {
   public:
    explicit LeptonicTopPtCut(Context&, float, float, const std::string& hyps="TTbarReconstruction", const std::string& disc="Chi2");
    virtual bool passes(const Event&) override;

   private:
    float tlep_pt_min_, tlep_pt_max_;
    Event::Handle<std::vector<ReconstructionHypothesis>> h_hyps_;
    std::string disc_name_;
  };
  /////

  class HypothesisDiscriminatorCut: public Selection {
   public:
    explicit HypothesisDiscriminatorCut(Context&, float, float, const std::string& disc="Chi2", const std::string& hyps="TTbarReconstruction");
    virtual bool passes(const Event&) override;

   private:
    float disc_min_, disc_max_;
    Event::Handle<std::vector<ReconstructionHypothesis>> h_hyps_;
    std::string disc_name_;
  };
  /////

  class LikelihoodSelection: public Selection {
   public:
     explicit LikelihoodSelection(uhh2::Context& ctx, float lmax );
     virtual bool passes(const Event&) override;

  private:
     float lmax_;
     Event::Handle<double> h_likelihood_;
  };   


   class DeltaPhiMetNeutrino: public Selection {
   public:
     explicit DeltaPhiMetNeutrino(uhh2::Context& ctx, float deltaphimin);
     virtual bool passes(const Event&) override;

  private:
     float deltaphimin_;
     Event::Handle<LorentzVector> h_neutrino_;
  }; 

  class DeltaPhiTaggedJetNeutrino: public Selection {
   public:
     explicit DeltaPhiTaggedJetNeutrino(uhh2::Context& ctx, float deltaphimax);
     virtual bool passes(const Event&) override;

  private:
     float deltaphimax_;
     Event::Handle<LorentzVector> h_neutrino_;
     Event::Handle<std::vector<TopJet>> h_taggedjet_;
  }; 

   class NeutrinopTSelection: public Selection {
   public:
      explicit NeutrinopTSelection(uhh2::Context& ctx, float pTmin);
      virtual bool passes(const Event&) override;
      
   private:
      float pTmin_;
      Event::Handle<LorentzVector> h_neutrino_;
   };

class DeltaPhiTaggedJetTopLep: public Selection {
   public:
     explicit DeltaPhiTaggedJetTopLep(uhh2::Context& ctx, float deltaphimax);
     virtual bool passes(const Event&) override;

  private:
   float deltaphimax_;
   Event::Handle<LorentzVector> h_neutrino_;
   Event::Handle<std::vector<TopJet>> h_taggedjet_;
   uhh2::Event::Handle<Jet> h_b_jets_;
};  

class ttbarpTSel: public Selection {
   public:
     explicit ttbarpTSel(uhh2::Context& ctx, float pTmin);
     virtual bool passes(const Event&) override;

  private:
   float pTmin_;
   Event::Handle<LorentzVector> h_neutrino_;
   Event::Handle<std::vector<TopJet>> h_taggedjet_;
   uhh2::Event::Handle<Jet> h_b_jets_;
};  

   class DMMETSelection: public Selection {
   public:
      explicit DMMETSelection(uhh2::Context& ctx, float DMMETmin);
      virtual bool passes(const Event&) override;
      
   private:
      float DMMETmin_;
      Event::Handle<LorentzVector> h_neutrino_;
   };

   /////
}
