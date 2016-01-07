#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/LumiSelection.h"

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
  bool is_mc;
  bool lumisel;
  bool mcpileupreweight;
  bool pvfilter;

  JetId btagAK4_wp;
  Event::Handle<int> h_flag_toptagevent;
  Event::Handle<int> h_flag_heptoptagevent;

  // selections
  std::unique_ptr<Selection> btagAK4_sel;
  //std::unique_ptr<Selection> leptoppt_sel;
  //std::unique_ptr<Selection> chi2_sel;
  std::unique_ptr<Selection> met_sel;
  std::unique_ptr<Selection> mtlep_sel;
  std::unique_ptr<Selection> jetmetdphi_sel;
  std::unique_ptr<Selection> mt2w_sel;
  std::unique_ptr<AnalysisModule> pu_reweight;
  std::unique_ptr<AnalysisModule> primaryvertex_filter;
  std::unique_ptr<Selection> lumi_selection;

  // hists
  std::unique_ptr<Hists> hi_input;
  //std::unique_ptr<Hists> hi_input__hyp;
  std::unique_ptr<Hists> hi_filter;
  //std::unique_ptr<Hists> hi_leptoppt;
  //std::unique_ptr<Hists> hi_leptoppt__hyp;
  //std::unique_ptr<Hists> hi_chi2;
  //std::unique_ptr<Hists> hi_chi2__hyp;
  std::unique_ptr<Hists> met_h;
  std::unique_ptr<Hists> mtlep_h;
  std::unique_ptr<Hists> jetmetdphi_h;
  std::unique_ptr<Hists> mt2w_h;
  std::unique_ptr<Hists> hi_t0b0;
  //std::unique_ptr<Hists> hi_t0b0__hyp;
  std::unique_ptr<Hists> hi_t0b1;
  //std::unique_ptr<Hists> hi_t0b1__hyp;
  std::unique_ptr<Hists> hi_t1;
  //std::unique_ptr<Hists> hi_t1__hyp;
  std::unique_ptr<Hists> hi_cmstag;
  std::unique_ptr<Hists> hi_heptag;
  std::unique_ptr<Hists> hi_bothtag;
};

ttDMPostSelectionModule::ttDMPostSelectionModule(Context& ctx){

  is_mc = ctx.get("dataset_type") == "MC";
  mcpileupreweight = ctx.get("pileup_directory_25ns") != "";
  pvfilter = string2bool(ctx.get("dopvfilter","true"));
  lumisel = ctx.get("lumi_file") != "";

  PrimaryVertexId pvid = StandardPrimaryVertexId();
  if(is_mc) {
    if(mcpileupreweight) pu_reweight.reset(new MCPileupReweight(ctx));
  } else {
    if(lumisel) lumi_selection.reset(new LumiSelection(ctx));
    if(pvfilter) primaryvertex_filter.reset(new PrimaryVertexCleaner(pvid));
  }
  
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
  h_flag_heptoptagevent = ctx.declare_event_input<int>("flag_heptoptagevent");

  // SELECTION
  // if(elec) leptoppt_sel.reset(new LeptonicTopPtCut(ctx, 140., infinity, "TTbarReconstruction", "Chi2"));
  // else if(muon) leptoppt_sel.reset(new AndSelection(ctx));
  //leptoppt_sel.reset(new AndSelection(ctx));

  //chi2_sel.reset(new HypothesisDiscriminatorCut(ctx, 0., 50., "TTbarReconstruction", "Chi2"));

  //ttDM SELECTION
  met_sel.reset(new METCut(320., std::numeric_limits<double>::infinity()));
  mtlep_sel.reset(new MTlepCut(160., std::numeric_limits<double>::infinity()));
  jetmetdphi_sel.reset(new METJetDPhiCut(1.2, 1));
  mt2w_sel.reset(new MT2WCut(200));

  // HISTS
  hi_input.reset(new ttDMPostSelectionHists(ctx, "input"));
  //hi_input__hyp.reset(new HypothesisHists(ctx, "input__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  hi_filter.reset(new ttDMPostSelectionHists(ctx, "filter"));
  //hi_leptoppt.reset(new ttDMPostSelectionHists(ctx, "leptoppt"));
  //hi_leptoppt__hyp.reset(new HypothesisHists(ctx, "leptoppt__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  //hi_chi2.reset(new ttDMPostSelectionHists(ctx, "chi2"));
  //hi_chi2__hyp.reset(new HypothesisHists(ctx, "chi2__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  met_h.reset(new ttDMPostSelectionHists(ctx, "met"));
  mtlep_h.reset(new ttDMPostSelectionHists(ctx, "mtlep"));
  jetmetdphi_h.reset(new ttDMPostSelectionHists(ctx, "jetmetdphi"));
  mt2w_h.reset(new ttDMPostSelectionHists(ctx, "mt2w"));

  hi_t0b0.reset(new ttDMPostSelectionHists(ctx, "t0b0"));
  //hi_t0b0__hyp.reset(new HypothesisHists(ctx, "t0b0__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  hi_t0b1.reset(new ttDMPostSelectionHists(ctx, "t0b1"));
  //hi_t0b1__hyp.reset(new HypothesisHists(ctx, "t0b1__hyp_chi2min", "TTbarReconstruction", "Chi2"));

  hi_t1.reset(new ttDMPostSelectionHists(ctx, "t1"));
  //hi_t1__hyp.reset(new HypothesisHists(ctx, "t1__hyp_chi2min", "TTbarReconstruction", "Chi2"));
  hi_cmstag.reset(new ttDMPostSelectionHists(ctx, "cmstag"));
  hi_heptag.reset(new ttDMPostSelectionHists(ctx, "heptag"));
  hi_bothtag.reset(new ttDMPostSelectionHists(ctx, "bothtag"));
}

bool ttDMPostSelectionModule::process(Event& event) {

  hi_input->fill(event);
  if(lumisel && !is_mc) if(!lumi_selection->passes(event)) return false;
  if(mcpileupreweight && is_mc) pu_reweight->process(event);
  if(pvfilter && !is_mc) primaryvertex_filter->process(event);

  hi_filter->fill(event);
  //hi_input__hyp->fill(event);

  //// LEPTONIC-TOP pt selection
  // bool pass_leptoppt = leptoppt_sel->passes(event);
  // if(!pass_leptoppt) return false;
  // hi_leptoppt->fill(event);
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

  /* MT_lep selection */
  bool pass_mtlep = mtlep_sel->passes(event);
  if(!pass_mtlep) return false;
  mtlep_h->fill(event);
  ////

  //// DeltaPhi Jet MET Selection
  bool pass_jetmetdphi = jetmetdphi_sel->passes(event);
  if(!pass_jetmetdphi) return false;
  jetmetdphi_h->fill(event);

  //// MT2W Selection
  bool pass_mt2w = mt2w_sel->passes(event);
  if(!pass_mt2w) return false;
  mt2w_h->fill(event);

  bool btag(btagAK4_sel->passes(event));
  bool cmstoptag(event.get(h_flag_toptagevent));
  bool heptoptag(event.get(h_flag_heptoptagevent));
  bool toptag = cmstoptag || heptoptag;

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
    if (cmstoptag) hi_cmstag->fill(event);
    if (heptoptag) hi_heptag->fill(event);
    if (cmstoptag && heptoptag) hi_bothtag->fill(event);
    //hi_t1__hyp->fill(event);
  }

  return false;
}

UHH2_REGISTER_ANALYSIS_MODULE(ttDMPostSelectionModule)
