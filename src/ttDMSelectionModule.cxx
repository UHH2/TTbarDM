#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/common/include/HypothesisHists.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"

#include "UHH2/TTbarDM/include/ttDMSemiLeptonicSelections.h"
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicUtils.h"
#include "UHH2/TTbarDM/include/ttDMSelectionHists.h"

#include "UHH2/TTbarDM/include/Mt2Com_bisect.h"

/** \brief module to produce "Selection" ntuples for the ttDM Analysis
 *
 *  -- GOALS:
 *   * complete object reconstruction (pt/eta cuts, IDs, jet-lepton cleaning, JER smearing)
 *   * apply (most of) the kinematic cuts for the lepton+jets SR. current cutflow:
 *     * HLT
 *     * ==1 lepton (w/ pt+eta+ID cuts)
 *     * >=2 AK4 jets w/ pt> 50 |eta|<2.4
 *     * >=1 AK4 jets w/ pt>200 |eta|<2.4
 *     * MET > 160 GeV (Final cut will be 320)
 *   * perform ttbar kinematical reconstruction (hyps stored in output ntuple)
 *
 * -- ITEMS TO BE IMPLEMENTED:
 *   * primary vertex selection
 *   * JER smearing for TopJet collection
 *   * update 2D cut values (after validation)
 *
 */
using namespace uhh2;

class ttDMSelectionModule: public AnalysisModule {
 public:
  explicit ttDMSelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

 private:
  Event::Handle<int> h_flag_toptagevent;
  Event::Handle<int> h_flag_heptoptagevent;

  bool is_mc;
  bool lumisel;
  bool mcpileupreweight;
  bool metfilters;
  bool pvfilter;

  std::unique_ptr<AnalysisModule> ht_calculator;
  std::unique_ptr<Selection> lumi_selection;
  std::unique_ptr<AndSelection> metfilters_selection;
  std::unique_ptr<AnalysisModule> primaryvertex_filter;
  std::unique_ptr<AnalysisModule> pu_reweight;
  
  // cleaners
  std::unique_ptr<MuonCleaner> muo_cleaner;
  std::unique_ptr<ElectronCleaner> ele_cleaner;
  std::unique_ptr<JetCorrector> jet_corrector;
  std::unique_ptr<JetResolutionSmearer> jetER_smearer;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner;
  std::unique_ptr<JetCleaner> jet_cleaner1;
  std::unique_ptr<JetCleaner> jet_cleaner2;
  std::unique_ptr<TopJetCorrector> topjet_corrector;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner;
  std::unique_ptr<TopJetCleaner> topjet_cleaner;

  // selections
  std::unique_ptr<Selection> trigger_sel;
  std::unique_ptr<AndSelection> lep1_sel;
  std::unique_ptr<Selection> jet2_sel;
  std::unique_ptr<Selection> jet1_sel;
  std::unique_ptr<Selection> met_sel;
  std::unique_ptr<Selection> toptagevent_sel;
  std::unique_ptr<Selection> heptoptagevent_sel;

  // ttbar reconstruction
  std::unique_ptr<AnalysisModule> ttgenprod;
  std::unique_ptr<AnalysisModule> reco_primlep;
  std::unique_ptr<AnalysisModule> ttbar_reco_toptag0, ttbar_reco_toptag1;

  // hists
  std::unique_ptr<Hists> input_h;
  std::unique_ptr<Hists> filter_h;
  std::unique_ptr<Hists> trigger_h;
  std::unique_ptr<Hists> lep1_h;
  std::unique_ptr<Hists> jet2_h;
  std::unique_ptr<Hists> jet1_h;
  std::unique_ptr<Hists> met_h;
  std::unique_ptr<Hists> toptagevent_h;
  std::unique_ptr<Hists> heptoptagevent_h;
};

ttDMSelectionModule::ttDMSelectionModule(Context & ctx){

  is_mc = ctx.get("dataset_type") == "MC";
  lumisel = ctx.get("lumi_file") != "";
  mcpileupreweight = ctx.get("pileup_directory_25ns") != "";
  metfilters = string2bool(ctx.get("dometfilters","true"));
  pvfilter = string2bool(ctx.get("dopvfilter","true"));

  PrimaryVertexId pvid = StandardPrimaryVertexId();
  if(is_mc){
    if(mcpileupreweight) pu_reweight.reset(new MCPileupReweight(ctx));
  } else{
    if(lumisel) lumi_selection.reset(new LumiSelection(ctx));
    if(metfilters){
      metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
      metfilters_selection->add<TriggerSelection>("CSCTightHalo", "Flag_CSCTightHaloFilter");
      metfilters_selection->add<TriggerSelection>("eeBadSc", "Flag_eeBadScFilter");
      if(pvfilter) metfilters_selection->add<NPVSelection>("1 good PV",1,-1,pvid);
    }
    if(pvfilter) primaryvertex_filter.reset(new PrimaryVertexCleaner(pvid));
  }
  ht_calculator.reset(new HTCalculator(ctx));

  //// OBJ CLEANING
  muo_cleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDMedium(), PtEtaCut(45., 2.1))));
  ele_cleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_MVAnotrig_Spring15_25ns_loose, PtEtaCut(50., 2.5))));
  if (is_mc) {
    jet_corrector.reset(new JetCorrector(JERFiles::Summer15_25ns_L123_AK4PFchs_MC));
    jetlepton_cleaner.reset(new JetLeptonCleaner(JERFiles::Summer15_25ns_L123_AK4PFchs_MC));
    topjet_corrector.reset(new TopJetCorrector(JERFiles::Summer15_25ns_L123_AK8PFchs_MC));
  } else {
    jet_corrector.reset(new JetCorrector(JERFiles::Summer15_25ns_L123_AK4PFchs_DATA));
    jetlepton_cleaner.reset(new JetLeptonCleaner(JERFiles::Summer15_25ns_L123_AK4PFchs_DATA));
    topjet_corrector.reset(new TopJetCorrector(JERFiles::Summer15_25ns_L123_AK8PFchs_DATA));
  }
  jetER_smearer.reset(new JetResolutionSmearer(ctx));
  jetlepton_cleaner->set_drmax(.4);
  jet_cleaner1.reset(new JetCleaner(25., std::numeric_limits<double>::infinity()));
  jet_cleaner2.reset(new JetCleaner(30., 2.4));
  topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));
  topjet_cleaner.reset(new TopJetCleaner(TopJetId(PtEtaCut(400., 2.4))));

  //// EVENT SELECTION
  bool muon(false), elec(false);
  const std::string channel(ctx.get("channel", ""));
  const std::string triggername(ctx.get("triggername", "NotSet"));
  if(channel == "muon") muon = true;
  else if(channel == "electron") elec = true;
  else throw std::runtime_error("undefined argument for 'channel' key in xml file (must be 'muon' or 'electron'): "+channel);

  lep1_sel.reset(new AndSelection(ctx));
  if(muon){
    lep1_sel->add<NMuonSelection>("muoN == 1", 1, 1);
    lep1_sel->add<NElectronSelection>("eleN == 0", 0, 0);

    if(triggername != "NotSet") trigger_sel = make_unique<TriggerSelection>(triggername);
    else trigger_sel = make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v*");
  }
  else if(elec){
    lep1_sel->add<NMuonSelection>("muoN == 0", 0, 0);
    lep1_sel->add<NElectronSelection>("eleN == 1", 1, 1);

    if(triggername != "NotSet") trigger_sel = make_unique<TriggerSelection>(triggername);
    else {
      trigger_sel = make_unique<TriggerSelection>("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1");
    }
  }

  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut( 50., 2.4))));
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(200., 2.4))));
  met_sel.reset(new METCut(160., std::numeric_limits<double>::infinity()));

  const TopJetId topjetID = AndId<TopJet>(CMSTopTag(), Tau32());
  const float minDR_topjet_jet(1.2);
  const TopJetId heptopjetID = AndId<TopJet>(HEPTopTag(), Tau21());

  toptagevent_sel.reset(new TopTagEventSelection(topjetID, minDR_topjet_jet));
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");

  heptoptagevent_sel.reset(new TopTagEventSelection(heptopjetID, minDR_topjet_jet));
  h_flag_heptoptagevent = ctx.declare_event_output<int>("flag_heptoptagevent");

  //// TTBAR KINEMATICAL RECO
  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  reco_primlep.reset(new PrimaryLepton(ctx));

  //// HISTS
  input_h.reset(new ttDMSelectionHists(ctx, "input"));
  filter_h.reset(new ttDMSelectionHists(ctx, "filter"));
  trigger_h.reset(new ttDMSelectionHists(ctx, "trigger"));
  lep1_h.reset(new ttDMSelectionHists(ctx, "lep1"));
  jet2_h.reset(new ttDMSelectionHists(ctx, "jet2"));
  jet1_h.reset(new ttDMSelectionHists(ctx, "jet1"));
  met_h.reset(new ttDMSelectionHists(ctx, "met"));
  toptagevent_h.reset(new ttDMSelectionHists(ctx, "toptagevent"));
  heptoptagevent_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent"));

}

bool ttDMSelectionModule::process(Event & event){

  input_h->fill(event);

  if(lumisel && !is_mc) if(!lumi_selection->passes(event)) return false;
  if(metfilters && !is_mc) if(!metfilters_selection->passes(event)) return false;
  if(mcpileupreweight && is_mc) pu_reweight->process(event);
  if(pvfilter && !is_mc) primaryvertex_filter->process(event);
  ht_calculator->process(event);

  filter_h->fill(event);
  
  //// HLT selection
  bool pass_trigger = trigger_sel->passes(event);
  if(!pass_trigger) return false;
  trigger_h->fill(event);

  //// LEPTON selection
  muo_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);

  ele_cleaner->process(event);
  sort_by_pt<Electron>(*event.electrons);

  bool pass_lep1 = lep1_sel->passes(event);
  if(!pass_lep1) return false;
  lep1_h->fill(event);

  //// JET selection
  jet_corrector->process(event);
  if (is_mc) jetER_smearer->process(event);
  jetlepton_cleaner->process(event);

  /* lepton-2Dcut boolean */
  jet_cleaner1->process(event); // jets w/ pt>25 GeV for lepton-2Dcut
  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);

  topjet_corrector->process(event);
  topjetlepton_cleaner->process(event);
  topjet_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  /* 2nd AK4 jet selection */
  bool pass_jet2 = jet2_sel->passes(event);
  if(!pass_jet2) return false;
  jet2_h->fill(event);

  /* 1st AK4 jet selection */
  bool pass_jet1 = jet1_sel->passes(event);
  if(!pass_jet1) return false;
  jet1_h->fill(event);

  //// MET selection
  bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;
  met_h->fill(event);

  //// TOPTAG-EVENT selection
  bool pass_toptagevent = toptagevent_sel->passes(event);
  if(pass_toptagevent) toptagevent_h->fill(event);

  //// HEPTOPTAG-EVENT selection
  bool pass_heptoptagevent = heptoptagevent_sel->passes(event);
  if(pass_heptoptagevent) heptoptagevent_h->fill(event);

  /* add flag_toptagevent to output ntuple */
  event.set(h_flag_toptagevent, int(pass_toptagevent));
  event.set(h_flag_heptoptagevent, int(pass_heptoptagevent));

  //// TTBAR KIN RECO
  if (is_mc) ttgenprod->process(event);

  reco_primlep->process(event);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ttDMSelectionModule)
