#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/common/include/HypothesisHists.h"

#include "UHH2/TTbarDM/include/ttDMSemiLeptonicSelections.h"
#include "UHH2/TTbarDM/include/ttDMPostSelectionHists.h"

using namespace uhh2;

class ttDMPostSelectionModule: public AnalysisModule {
 public:
  explicit ttDMPostSelectionModule(Context&);
  virtual bool process(Event&) override;

 private:
  std::unique_ptr<AnalysisModule> ttgenprod;
  Event::Handle<TTbarGen> h_ttbargen;

  //Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;

  JetId btagAK4_wp;
  Event::Handle<int> h_flag_toptagevent;

  // selections
  std::unique_ptr<Selection> btagAK4_sel;
  std::unique_ptr<Selection> leptoppt_sel;
  //std::unique_ptr<Selection> chi2_sel;
  std::unique_ptr<Selection> met_sel;
  std::unique_ptr<Selection> dphijet1_sel;
  std::unique_ptr<Selection> dphijet2_sel;

  // hists
  std::unique_ptr<Hists> hi_input;
  //std::unique_ptr<Hists> hi_input__hyp;
  std::unique_ptr<Hists> hi_leptoppt;
  //std::unique_ptr<Hists> hi_leptoppt__hyp;
  //std::unique_ptr<Hists> hi_chi2;
  //std::unique_ptr<Hists> hi_chi2__hyp;
  std::unique_ptr<Hists> met_h;
  std::unique_ptr<Hists> dphijet1_h;
  std::unique_ptr<Hists> dphijet2_h;
  std::unique_ptr<Hists> hi_t0b0;
  //std::unique_ptr<Hists> hi_t0b0__hyp;
  std::unique_ptr<Hists> hi_t0b1;
  //std::unique_ptr<Hists> hi_t0b1__hyp;
  std::unique_ptr<Hists> hi_t1;
  //std::unique_ptr<Hists> hi_t1__hyp;
};

ttDMPostSelectionModule::ttDMPostSelectionModule(Context& ctx){

  //bool muon(false), elec(false);
  const std::string channel(ctx.get("channel", ""));
  // if(channel == "muon") muon = true;
  // else if(channel == "electron") elec = true;
  // else throw std::runtime_error("undefined argument for 'channel' key in xml file (must be 'muon' or 'electron'): "+channel);

  // ttbar GEN
  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // ttbar RECO hypotheses
  //h_ttbar_hyps = ctx.declare_event_input<std::vector<ReconstructionHypothesis>>("TTbarReconstruction");

  // b-tagging
  btagAK4_wp = CSVBTag(CSVBTag::WP_MEDIUM);
  btagAK4_sel.reset(new NJetSelection(1, -1, btagAK4_wp));

  // top-tagging flag (from ttDMSelection ntuple)
  h_flag_toptagevent = ctx.declare_event_input<int>("flag_toptagevent");

  // SELECTION
  // if(elec) leptoppt_sel.reset(new LeptonicTopPtCut(ctx, 140., infinity, "TTbarReconstruction", "Chi2"));
  // else if(muon) leptoppt_sel.reset(new AndSelection(ctx));
  leptoppt_sel.reset(new AndSelection(ctx));

  //chi2_sel.reset(new HypothesisDiscriminatorCut(ctx, 0., 50., "TTbarReconstruction", "Chi2"));

  //ttDM SELECTION
  met_sel.reset(new METCut(320., std::numeric_limits<double>::infinity()));
  dphijet1_sel.reset(new METJetDPhiCut(1.2, 0));
  dphijet2_sel.reset(new METJetDPhiCut(1.2, 1));

  // HISTS
  hi_input.reset(new ttDMPostSelectionHists(ctx, "input"));
  //hi_input__hyp.reset(new HypothesisHists(ctx, "input__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  hi_leptoppt.reset(new ttDMPostSelectionHists(ctx, "leptoppt"));
  //hi_leptoppt__hyp.reset(new HypothesisHists(ctx, "leptoppt__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  //hi_chi2.reset(new ttDMPostSelectionHists(ctx, "chi2"));
  //hi_chi2__hyp.reset(new HypothesisHists(ctx, "chi2__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  met_h.reset(new ttDMPostSelectionHists(ctx, "met"));
  dphijet1_h.reset(new ttDMPostSelectionHists(ctx, "dphijet1"));
  dphijet2_h.reset(new ttDMPostSelectionHists(ctx, "dphijet2"));

  hi_t0b0.reset(new ttDMPostSelectionHists(ctx, "t0b0"));
  //hi_t0b0__hyp.reset(new HypothesisHists(ctx, "t0b0__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  hi_t0b1.reset(new ttDMPostSelectionHists(ctx, "t0b1"));
  //hi_t0b1__hyp.reset(new HypothesisHists(ctx, "t0b1__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  hi_t1.reset(new ttDMPostSelectionHists(ctx, "t1"));
  //hi_t1__hyp.reset(new HypothesisHists(ctx, "t1__hyp_chi2min", "TTbarReconstruction", "Chi2"));
}

bool ttDMPostSelectionModule::process(Event& event) {

  hi_input->fill(event);
  //hi_input__hyp->fill(event);

  //// LEPTONIC-TOP pt selection
  bool pass_leptoppt = leptoppt_sel->passes(event);
  if(!pass_leptoppt) return false;
  hi_leptoppt->fill(event);
  //hi_leptoppt__hyp->fill(event);
  ////

  //// CHI2 selection
  // bool pass_chi2 = chi2_sel->passes(event);
  // if(!pass_chi2) return false;
  // hi_chi2->fill(event);
  // hi_chi2__hyp->fill(event);
  ////

  //// MET Post-Selection
  bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;
  met_h->fill(event);

  //// DeltaPhi Lead Jet Selection
  bool pass_dphijet1 = dphijet1_sel->passes(event);
  if(!pass_dphijet1) return false;
  dphijet1_h->fill(event);

  //// DeltaPhi Sub-Lead Jet Selection
  bool pass_dphijet2 = dphijet2_sel->passes(event);
  if(!pass_dphijet2) return false;
  dphijet2_h->fill(event);

  bool btag(btagAK4_sel->passes(event));
  bool toptag(event.get(h_flag_toptagevent));

  if(!toptag){

    if(!btag){

      hi_t0b0->fill(event);
      //hi_t0b0__hyp->fill(event);
    }
    else {

      hi_t0b1->fill(event);
      //hi_t0b1__hyp->fill(event);
    }
  }
  else {

    hi_t1->fill(event);
    //hi_t1__hyp->fill(event);
  }

  return false;
}

UHH2_REGISTER_ANALYSIS_MODULE(ttDMPostSelectionModule)
