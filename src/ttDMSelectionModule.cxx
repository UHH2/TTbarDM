#include <iostream>
#include <memory>
#include "TError.h"
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
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/CollectionProducer.h"
#include "UHH2/common/include/HandleSelection.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/ElectronHists.h"

#include "UHH2/TTbarDM/include/ttDMSemiLeptonicSelections.h"
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicUtils.h"
#include "UHH2/TTbarDM/include/ttDMSelectionHists.h"
#include "UHH2/TTbarDM/include/ttDMGenHists.h"
#include "UHH2/TTbarDM/include/ttDMReconstruction_Likelihood.h"
#include "UHH2/TTbarDM/include/ttDMReconstructionHists_Likelihood.h"
#include "UHH2/TTbarDM/include/ttDMReconstructionHists_ScansAndMarginalisation.h"

#include "UHH2/TTbarDM/include/Mt2Com_bisect.h"

/** \brief module to produce "Selection" ntuples for the ttDM Analysis
 *
 *  -- GOALS:
 *   * complete object reconstruction (pt/eta cuts, IDs, jet-lepton cleaning, JER smearing)
 *   * apply (most of) the kinematic cuts for the lepton+jets SR. current cutflow:
 *     * HLT
 *     * ==1 lepton (w/ pt+eta+ID cuts)
 *     * >=3 AK4 jets w/ pt> 50 |eta|<2.4
 *     * >=1 AK4 jets w/ pt>200 |eta|<2.4
 *     * MET > 160 GeV (Final cut will be 320)
 *     * LEP Transverse Mass > 160
 *     * Need to add MT2W Razor cut
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
   
   Event::Handle<bool> h_flag_twodcut;

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
   std::unique_ptr<AnalysisModule> lumi_weight; 
   std::unique_ptr<AnalysisModule> ttgenprod;
   std::unique_ptr<AnalysisModule> muo_med_noniso_SF;
   //std::unique_ptr<AnalysisModule> muo_triggerSF;
    std::unique_ptr<AnalysisModule> reco_primlep;
   std::unique_ptr<AnalysisModule> collectionprod_muonmed;
   std::unique_ptr<AnalysisModule> collectionprod_eletight;
   // cleaners
   std::unique_ptr<MuonCleaner> muo_cleaner;
   std::unique_ptr<ElectronCleaner> ele_cleaner;
   std::unique_ptr<JetCorrector> jet_corrector;
   std::unique_ptr<JetResolutionSmearer> jetER_smearer;
   std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner;
   std::unique_ptr<JetCleaner> jet_cleaner0;
   std::unique_ptr<JetCleaner> jet_cleaner1;
   std::unique_ptr<JetCleaner> jet_cleaner2;
   std::unique_ptr<TopJetCorrector> topjet_corrector;
   //std::unique_ptr<TopJetResolutionSmearer> topjetER_smearer;
   std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner;
   std::unique_ptr<TopJetCleaner> topjet_cleaner;
         
   // selections
   std::unique_ptr<Selection> trigger_sel; 
   std::unique_ptr<AndSelection> lep1_sel;
   std::unique_ptr<Selection> jet2_sel;
   std::unique_ptr<Selection> jet1_sel;
   std::unique_ptr<Selection> met_sel;
   std::unique_ptr<Selection> pv_sel;
   std::unique_ptr<Selection> twodcut_sel;  
   
   //neutrino reconstruction
   std::unique_ptr<AnalysisModule> ttbar_DM_reco_likelihood;
   std::unique_ptr<Hists> likelihood_hists;

   // hists
   std::unique_ptr<Hists> lumihists;
   std::unique_ptr<Hists> filter_h;
   std::unique_ptr<Hists> input_h,electrons_input_h,jets_input_h,muons_input_h,events_input_h,topjets_input_h;
   std::unique_ptr<Hists> trigger_h, electrons_trigger_h,jets_trigger_h,muons_trigger_h,events_trigger_h,topjets_trigger_h;
   std::unique_ptr<Hists> lep1_h,electrons_lep1_h,jets_lep1_h,muons_lep1_h,events_lep1_h,topjets_lep1_h;
   std::unique_ptr<Hists> jet1_h, electrons_jet1_h,jets_jet1_h,muons_jet1_h,events_jet1_h,topjets_jet1_h;
   std::unique_ptr<Hists> jet2_h, electrons_jet2_h,jets_jet2_h,muons_jet2_h,events_jet2_h,topjets_jet2_h;
   std::unique_ptr<Hists> met_h, electrons_met_h,jets_met_h,muons_met_h,events_met_h,topjets_met_h, genhists_met_h;//scans_h;
};

ttDMSelectionModule::ttDMSelectionModule(Context & ctx){
   
   h_flag_twodcut = ctx.declare_event_output<bool>("flag_twodcut");
   
   is_mc = ctx.get("dataset_type") == "MC";
   lumisel = ctx.get("lumi_file") != ""; 
   mcpileupreweight = ctx.get("pileup_directory") != "";
   metfilters = string2bool(ctx.get("dometfilters","true")); 
   pvfilter = string2bool(ctx.get("dopvfilter","true"));
   
   PrimaryVertexId pvid = StandardPrimaryVertexId();
   if(pvfilter) 
      {
         primaryvertex_filter.reset(new PrimaryVertexCleaner(pvid)); 
         pv_sel.reset(new NPVSelection(1,-1,pvid));
      }
   if(is_mc){
      if(mcpileupreweight) pu_reweight.reset(new MCPileupReweight(ctx));
      lumi_weight.reset(new MCLumiWeight(ctx));
   } else{
      if(lumisel) lumi_selection.reset(new LumiSelection(ctx));
   }
   if(metfilters){
      metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
      metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
      metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
      metfilters_selection->add<TriggerSelection>("CSCTightHalo2015Filter", "Flag_CSCTightHalo2015Filter");
      metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
      metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
      metfilters_selection->add<TriggerSelection>("chargedHadronTrackResolutionFilter", "Flag_chargedHadronTrackResolutionFilter"); 
      metfilters_selection->add<TriggerSelection>("muonBadTrackFilter", "Flag_muonBadTrackFilter");
   }
   
   ht_calculator.reset(new HTCalculator(ctx));
   muo_med_noniso_SF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_7_6_3/src/UHH2/common/data/MuonID_Z_RunCD_Reco76X_Feb15.root","MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1", 1., "id"));
   //muo_triggerSF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_7_4_15_patch1/src/UHH2/common/data/SingleMuonTrigger_Z_RunD_Reco74X_Nov20.root", "Mu45_eta2p1_PtEtaBins", 1., "trg")); 
   //do we have to veto tight or medium muons? NO SF for Iso applied is approx. 1

   ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
   reco_primlep.reset(new PrimaryLepton(ctx));
   
   //neutrino reconstruction
   ttbar_DM_reco_likelihood.reset(new ttDMReconstruction_Likelihood(ctx));
   likelihood_hists.reset(new ttDMReconstructionHists_Likelihood(ctx, "hists_likelihood"));

   //// OBJ CLEANING
   //muo_cleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDMedium(),PtEtaCut(47., 2.1)))); 
   //ele_cleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_Spring15_25ns_tight_noIso, PtEtaCut(50., 2.5))));
   muo_cleaner.reset(new MuonCleaner    (AndId<Muon>    (PtEtaCut  (10., 2.1), MuonIDLoose()))); 
   ele_cleaner.reset(new ElectronCleaner(AndId<Electron>(PtEtaSCCut(10., 2.5), ElectronID_MVAnotrig_Spring15_25ns_loose)));
   
   if (is_mc) {
      jet_corrector.reset(new JetCorrector(ctx, JERFiles::Fall15_25ns_L123_AK4PFchs_MC));
      jetlepton_cleaner.reset(new JetLeptonCleaner(ctx, JERFiles::Fall15_25ns_L123_AK4PFchs_MC));
      topjet_corrector.reset(new TopJetCorrector(ctx, JERFiles::Fall15_25ns_L123_AK8PFchs_MC));
      jetER_smearer.reset(new JetResolutionSmearer(ctx));
   } else {
      jet_corrector.reset(new JetCorrector(ctx, JERFiles::Fall15_25ns_L123_AK4PFchs_DATA));
      jetlepton_cleaner.reset(new JetLeptonCleaner(ctx, JERFiles::Fall15_25ns_L123_AK4PFchs_DATA));
      topjet_corrector.reset(new TopJetCorrector(ctx, JERFiles::Fall15_25ns_L123_AK8PFchs_DATA));
   }
   jetlepton_cleaner->set_drmax(.4);
   //jet_cleaner1.reset(new JetCleaner(ctx, AndId<Jet>(JetPFID(JetPFID::WP_LOOSE),PtEtaCut(25., std::numeric_limits<double>::infinity())))); 
   //jet_cleaner2.reset(new JetCleaner(ctx,AndId<Jet>(JetPFID(JetPFID::WP_LOOSE),PtEtaCut(30., 2.4))));
   jet_cleaner1.reset(new JetCleaner(ctx, AndId<Jet>(JetPFID(JetPFID::WP_LOOSE),PtEtaCut(15., 3.0)))); 
   jet_cleaner2.reset(new JetCleaner(ctx,AndId<Jet>(JetPFID(JetPFID::WP_LOOSE),PtEtaCut(15., 2.4))));

   topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));
   topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(400., 2.4))));
   jet_cleaner0.reset(new JetCleaner(ctx, AndId<Jet>(JetPFID(JetPFID::WP_LOOSE),PtEtaCut(0., 5.0))));
   ////
        
   //// EVENT SELECTION
   bool muon(false), elec(false);
   const std::string channel(ctx.get("channel", ""));
   const std::string triggername(ctx.get("triggername", "NotSet"));
   if(channel == "muon") muon = true;
   else if(channel == "electron") elec = true;
   else throw std::runtime_error("undefined argument for 'channel' key in xml file (must be 'muon' or 'electron'): "+channel);
   
   const MuonId muonid = AndId<Muon>(PtEtaCut  (47., 2.1), MuonIDMedium());
   const ElectronId electronid = AndId<Electron>(ElectronID_Spring15_25ns_tight_noIso,PtEtaSCCut(50., 2.5));
   
   lep1_sel.reset(new AndSelection(ctx));
   if(muon){
      lep1_sel->add<NMuonSelection>("muo==1", 1, 1, muonid);
      lep1_sel->add<NElectronSelection>("ele=0",0, 0, electronid);
      
      if(triggername != "NotSet") trigger_sel.reset(new TriggerSelection(triggername));
      else trigger_sel.reset(new TriggerSelection("HLT_Mu45_eta2p1_v*"));
   }
   else if(elec){
      lep1_sel->add<NMuonSelection>("muo==0", 0, 0, muonid);
      lep1_sel->add<NElectronSelection>("elle==1",1, 1, electronid);

      if(triggername != "NotSet") trigger_sel.reset(new TriggerSelection(triggername));
      else {
         if (!is_mc) trigger_sel.reset(new TriggerSelection("HLT_Ele15_IsoVVVL_PFHT350_PFMET70_v*"));
         else trigger_sel.reset (new TriggerSelection("HLT_Ele15_IsoVVVL_PFHT400_PFMET70_v*"));
      }
   }
   collectionprod_muonmed.reset(new CollectionProducer<Muon>(ctx, "muons", "h_muons_medium", muonid)); 
   collectionprod_eletight.reset(new CollectionProducer<Electron>(ctx, "electrons", "h_electrons_tight", electronid)); 
   // jet2_sel.reset(new NJetSelection(3, -1, JetId(PtEtaCut( 50., 2.4))));    
   // jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(200., 2.4))));
   jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut( 50., 2.4))));    
   jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(50., 2.4))));
   
   met_sel.reset(new METCut(50., std::numeric_limits<double>::infinity()));  
   twodcut_sel.reset(new TwoDCut(ctx,.4, 20.));
      
   //// HISTS
   lumihists.reset(new LuminosityHists(ctx, "lumi"));
   input_h.reset(new ttDMSelectionHists(ctx, "input"));
   filter_h.reset(new ttDMSelectionHists(ctx, "filter"));
   trigger_h.reset(new ttDMSelectionHists(ctx, "trigger"));
   lep1_h.reset(new ttDMSelectionHists(ctx, "lep1"));
   jet2_h.reset(new ttDMSelectionHists(ctx, "jet2"));
   jet1_h.reset(new ttDMSelectionHists(ctx, "jet1"));
   met_h.reset(new ttDMSelectionHists(ctx, "met"));
     
   electrons_input_h.reset(new ElectronHists(ctx,"electrons_input"));
   jets_input_h.reset(new JetHists(ctx,"jets_input"));
   muons_input_h.reset(new MuonHists(ctx,"muons_input"));
   events_input_h.reset(new EventHists(ctx,"events_input"));
   topjets_input_h.reset(new TopJetHists(ctx,"topjets_input",4,"patJetsHepTopTagCHSPacked_daughters"));
   
   electrons_trigger_h.reset(new ElectronHists(ctx,"electrons_trigger"));
   jets_trigger_h.reset(new JetHists(ctx,"jets_trigger"));
   muons_trigger_h.reset(new MuonHists(ctx,"muons_trigger"));
   events_trigger_h.reset(new EventHists(ctx,"events_trigger"));
   topjets_trigger_h.reset(new TopJetHists(ctx,"topjets_trigger",4,"patJetsHepTopTagCHSPacked_daughters"));
   
   electrons_lep1_h.reset(new ElectronHists(ctx,"electrons_lep1"));
   jets_lep1_h.reset(new JetHists(ctx,"jets_lep1"));
   muons_lep1_h.reset(new MuonHists(ctx,"muons_lep1"));
   events_lep1_h.reset(new EventHists(ctx,"events_lep1"));
   topjets_lep1_h.reset(new TopJetHists(ctx,"topjets_lep1",4,"patJetsHepTopTagCHSPacked_daughters"));
   
   electrons_jet1_h.reset(new ElectronHists(ctx,"electrons_jet1"));
   jets_jet1_h.reset(new JetHists(ctx,"jets_jet1"));
   muons_jet1_h.reset(new MuonHists(ctx,"muons_jet1"));
   events_jet1_h.reset(new EventHists(ctx,"events_jet1"));
   topjets_jet1_h.reset(new TopJetHists(ctx,"topjets_jet1",4,"patJetsHepTopTagCHSPacked_daughters"));
   
   electrons_jet2_h.reset(new ElectronHists(ctx,"electrons_jet2"));
   jets_jet2_h.reset(new JetHists(ctx,"jets_jet2"));
   muons_jet2_h.reset(new MuonHists(ctx,"muons_jet2"));
   events_jet2_h.reset(new EventHists(ctx,"events_jet2"));
   topjets_jet2_h.reset(new TopJetHists(ctx,"topjets_jet2",4,"patJetsHepTopTagCHSPacked_daughters"));
   
   electrons_met_h.reset(new ElectronHists(ctx,"electrons_met"));
   jets_met_h.reset(new JetHists(ctx,"jets_met"));
   muons_met_h.reset(new MuonHists(ctx,"muons_met"));
   events_met_h.reset(new EventHists(ctx,"events_met"));
   topjets_met_h.reset(new TopJetHists(ctx,"topjets_met",4,"patJetsHepTopTagCHSPacked_daughters"));
   genhists_met_h.reset(new ttDMGenHists(ctx,"genhists_met"));
   //scans_h.reset(new ttDMReconstructionHists_ScansAndMarginalisation(ctx,"scans"));
}

bool ttDMSelectionModule::process(Event & event){
   if (is_mc) ttgenprod->process(event);
   input_h->fill(event);
   if(lumisel && !is_mc) if(!lumi_selection->passes(event)) return false;
   if(metfilters) if(!metfilters_selection->passes(event)) return false;
   if(mcpileupreweight && is_mc) pu_reweight->process(event);
   if(is_mc) lumi_weight->process(event); 
   if(pvfilter) 
      {
         primaryvertex_filter->process(event);
         bool pass_goodpv = pv_sel->passes(event);
         if (!pass_goodpv) return false; 
      }
   ht_calculator->process(event);
   filter_h->fill(event); 
   //// HLT selection
   bool pass_trigger = trigger_sel->passes(event);
   if(!pass_trigger) return false;
   
   trigger_h->fill(event);
   electrons_trigger_h->fill(event);
   jets_trigger_h->fill(event);
   muons_trigger_h->fill(event);
   events_trigger_h->fill(event);
   topjets_trigger_h->fill(event);

   //// LEPTON selection
   muo_cleaner->process(event);
   sort_by_pt<Muon>(*event.muons);
   ele_cleaner->process(event);
   sort_by_pt<Electron>(*event.electrons);
   reco_primlep->process(event);
   
   jet_corrector->process(event);
   jet_cleaner0->process(event);
   if (is_mc) jetER_smearer->process(event);
   jetlepton_cleaner->process(event);
  
   bool pass_lep1 = lep1_sel->passes(event);
   if(!pass_lep1) return false;
   collectionprod_muonmed->process(event);
   collectionprod_eletight->process(event);

   jet_cleaner1->process(event); // jets w/ pt>15 GeV for lepton-2Dcut
   bool pass_twodcut = twodcut_sel->passes(event); //ACHTUNG DARF NUR EIN MUON!!
   event.set(h_flag_twodcut, bool(pass_twodcut));
   
   jet_cleaner2->process(event);
   sort_by_pt<Jet>(*event.jets);
   topjet_corrector->process(event);
   //topjetER_smearer->process(event);
   topjetlepton_cleaner->process(event);
   topjet_cleaner->process(event);
   sort_by_pt<TopJet>(*event.topjets);
   
   muo_med_noniso_SF->process(event);
   //muo_triggerSF->process(event);
   
   lep1_h->fill(event);
   electrons_lep1_h->fill(event);
   jets_lep1_h->fill(event);
   muons_lep1_h->fill(event);
   events_lep1_h->fill(event);
   topjets_lep1_h->fill(event);
   
   if (event.isRealData) lumihists->fill(event);
   
   //// JET selection
   /* 2nd AK4 jet selection */
   bool pass_jet2 = jet2_sel->passes(event);
   if(!pass_jet2) return false;
   
   jet2_h->fill(event);
   electrons_jet2_h->fill(event);
   jets_jet2_h->fill(event);
   muons_jet2_h->fill(event);
   events_jet2_h->fill(event);
   topjets_jet2_h->fill(event);
   
   /* 1st AK4 jet selection */
   bool pass_jet1 = jet1_sel->passes(event);
   if(!pass_jet1) return false;
   
   jet1_h->fill(event);
   electrons_jet1_h->fill(event);
   jets_jet1_h->fill(event);
   muons_jet1_h->fill(event);
   events_jet1_h->fill(event);
   topjets_jet1_h->fill(event);
   
   //// MET selection
   bool pass_met = met_sel->passes(event);
   if(!pass_met) return false;
   
   met_h->fill(event);
   electrons_met_h->fill(event);
   jets_met_h->fill(event);
   muons_met_h->fill(event);
   events_met_h->fill(event);
   topjets_met_h->fill(event);
   if(!event.isRealData) genhists_met_h->fill(event);
   // neutrino reconstruction
   ttbar_DM_reco_likelihood->process(event);
   likelihood_hists->fill(event);
   //if(!event.isRealData) scans_h->fill(event);
   return true;
}



UHH2_REGISTER_ANALYSIS_MODULE(ttDMSelectionModule)

