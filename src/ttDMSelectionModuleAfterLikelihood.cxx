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
#include "UHH2/common/include/TTbarGenHists.h"
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
#include "UHH2/common/include/PartonHT.h"
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicSelections.h"
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicUtils.h"
#include "UHH2/TTbarDM/include/ttDMSelectionHists.h"
#include "UHH2/TTbarDM/include/ttDMGenHists.h"
#include "UHH2/TTbarDM/include/ttDMReconstruction_Likelihood.h"
#include "UHH2/TTbarDM/include/ttDMReconstructionHists_Likelihood.h"

#include "UHH2/TTbarDM/include/Mt2Com_bisect.h"

//to check: changed jet eta to 2.4 from 2.5

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

class ttDMSelectionModuleAfterLikelihood: public AnalysisModule {
 public:
  explicit ttDMSelectionModuleAfterLikelihood(Context & ctx);
  virtual bool process(Event & event) override;

 private:
   Event::Handle<bool> h_flag_twodcut;
   Event::Handle<double> h_likelihood;
   Event::Handle<LorentzVector> h_recneutrino;
   Event::Handle<Jet> h_bjets;
   Event::Handle<TTbarGen> h_ttbargen;
   Event::Handle<float> h_isoweight;
   Event::Handle<float> h_isoweight_down;
   Event::Handle<float> h_isoweight_up;

   Event::Handle<float> h_tagweight;
   Event::Handle<float> h_tagweight_down;
   Event::Handle<float> h_tagweight_up;
   
   bool is_mc;
   bool lumisel;
   bool mcpileupreweight;
   bool metfilters;
   bool pvfilter;
     
   bool isSemiLept;
   bool isDiLept;
   bool isOther;
   
   bool isResolved=false;
   bool isBoosted=false;
   bool IsWJetsInclusive=false;
   bool IsWJetsHT =false;

   //analysis modules
   std::unique_ptr<AnalysisModule> ht_calculator, partonht_calculator;
   
   std::unique_ptr<AnalysisModule> primaryvertex_filter;
   std::unique_ptr<AnalysisModule> pu_reweight;
   std::unique_ptr<AnalysisModule> lumi_weight; 
   std::unique_ptr<AnalysisModule> muo_med_noniso_SF;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_pteta;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_pteta;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_WP3_wobtag;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_WP3_wobtag;
   //std::unique_ptr<AnalysisModule> muo_triggerSF;
   std::unique_ptr<AnalysisModule> ttgenprod;
   std::unique_ptr<AnalysisModule> hepttleptoncleaner;
   std::unique_ptr<AnalysisModule> jetmuon_cleaner;
   std::unique_ptr<AnalysisModule> jetmuon_cleaner_resolved, muon_cleaner_resolved;
   
   // cleaners
   std::unique_ptr<MuonCleaner> muon_cleaner_iso;
   std::unique_ptr<GenericTopJetCorrector> heptopjet_corrector;
   std::unique_ptr<GenericTopJetCorrector> heptopjet_correctorL2L3;
   std::unique_ptr<JetCleaner> jet_cleaner_resolved;
   std::unique_ptr<ElectronCleaner> ele_cleaner_resolved;
   // selections
   std::unique_ptr<Selection> heptoptagevent_WP3_wobtag_sel;
   std::unique_ptr<Selection> pv_sel;
   std::unique_ptr<Selection> btag_sel;
   std::unique_ptr<AndSelection> lep1_sel, lepVeto_sel;
   std::unique_ptr<Selection> mtlep100_sel;
   std::unique_ptr<Selection> mtlep160_sel;
   std::unique_ptr<Selection> lumi_selection;
   std::unique_ptr<Selection> selection_Met160;
   std::unique_ptr<Selection> selection_Met250;
   std::unique_ptr<AndSelection> lep1_resolved_sel;
   std::unique_ptr<AndSelection> metfilters_selection;
   std::unique_ptr<Selection> btagtight_sel;
   std::unique_ptr<Selection> jet_resolved_sel;
   std::unique_ptr<Selection> jetmetdphijet2_sel;
   std::unique_ptr<Selection> MT2W_sel;
   std::unique_ptr<Selection> triggerMu50_sel, triggerTrkMu50_sel, partonht_sel;
   // hists
   std::unique_ptr<Hists> lumihists;
   std::unique_ptr<Hists> filter_h;
   std::unique_ptr<Hists> input_h;
   std::unique_ptr<Hists> met_h, electrons_met_h,jets_met_h,muons_met_h,events_met_h,topjets_met_h, genhists_met_h,likelihood_met_h,ttbargenhists_met_h;
   std::unique_ptr<Hists> nonisolep_h, electrons_nonisolep_h,jets_nonisolep_h,muons_nonisolep_h,events_nonisolep_h,topjets_nonisolep_h,genhists_nonisolep_h;
   std::unique_ptr<Hists> isolep_h, electrons_isolep_h,jets_isolep_h,muons_isolep_h,events_isolep_h,topjets_isolep_h;
   std::unique_ptr<Hists> semileptcontrolregion_h,electrons_semileptcontrolregion_h,jets_semileptcontrolregion_h,muons_semileptcontrolregion_h,events_semileptcontrolregion_h,topjets_semileptcontrolregion_h, likelihood_semileptcontrolregion_h;
   std::unique_ptr<TopJetHists> topjets_isolep_tagged_WP3_wobtag_h,topjets_semileptcontrolregion_tagged_h;
   std::unique_ptr<Hists> preselection_h,electrons_preselection_h, jets_preselection_h,topjets_preselection_h,events_preselection_h,muons_preselection_h, genhists_preselection_h,likelihood_preselection_h,ttbargenhists_preselection_h;
   std::unique_ptr<Hists> preselection_Met160_h,electrons_preselection_Met160_h, jets_preselection_Met160_h,topjets_preselection_Met160_h,events_preselection_Met160_h,muons_preselection_Met160_h, genhists_preselection_Met160_h,likelihood_preselection_Met160_h,ttbargenhists_preselection_Met160_h;
   std::unique_ptr<Hists> leptonveto_h,electrons_leptonveto_h, jets_leptonveto_h,topjets_leptonveto_h,events_leptonveto_h,muons_leptonveto_h, genhists_leptonveto_h,likelihood_leptonveto_h,ttbargenhists_leptonveto_h;
   std::unique_ptr<Hists> dileptcontrolregion_h,electrons_dileptcontrolregion_h, jets_dileptcontrolregion_h,topjets_dileptcontrolregion_h,events_dileptcontrolregion_h,muons_dileptcontrolregion_h, genhists_dileptcontrolregion_h,ttbargenhists_dileptcontrolregion_h,likelihood_dileptcontrolregion_h;
   std::unique_ptr<Hists> electrons_signalregion_noniso_h,jets_signalregion_noniso_h,topjets_signalregion_noniso_h,muons_signalregion_noniso_h,events_signalregion_noniso_h,likelihood_signalregion_noniso_h, signalregion_noniso_h;
   std::unique_ptr<Hists> signalregion_iso_h,electrons_signalregion_iso_h,jets_signalregion_iso_h,muons_signalregion_iso_h,events_signalregion_iso_h,topjets_signalregion_iso_h,likelihood_signalregion_iso_h;
   std::unique_ptr<Hists> signalregion_iso_Met160_h,electrons_signalregion_iso_Met160_h,jets_signalregion_iso_Met160_h,muons_signalregion_iso_Met160_h,events_signalregion_iso_Met160_h,topjets_signalregion_iso_Met160_h,likelihood_signalregion_iso_Met160_h;
   std::unique_ptr<Hists> btag_h,electrons_btag_h,jets_btag_h,muons_btag_h,events_btag_h,topjets_btag_h,likelihood_btag_h;
   std::unique_ptr<Hists> mtlepsel_h,electrons_mtlepsel_h,jets_mtlepsel_h,muons_mtlepsel_h,events_mtlepsel_h,topjets_mtlepsel_h,likelihood_mtlepsel_h;
   std::unique_ptr<Hists> Wjetscontrolregion_h,electrons_Wjetscontrolregion_h,jets_Wjetscontrolregion_h,muons_Wjetscontrolregion_h,events_Wjetscontrolregion_h,topjets_Wjetscontrolregion_h,likelihood_Wjetscontrolregion_h;
   std::unique_ptr<Hists> resolved_h, combination_h,resolved_preselection_h;
   std::unique_ptr<BTagMCEfficiencyHists> h_btag_SF;
   std::unique_ptr<Selection> DeltaRMuonJet_sel, jet2_sel;
   std::unique_ptr<AndSelection> lepVeto_resolved_sel;
   std::unique_ptr<AnalysisModule> muonTRK_SF, muonID_SF, muonHLT_SF, muonHLT2_SF, muonIso_SF, toptag_SF, apply_btag_SF, elecGsf_SF, elecID_SF;

   //For usage of run dependent muonHLT Effs
   double lumi_tot;
   double lumi1;
   double lumi2;
};

ttDMSelectionModuleAfterLikelihood::ttDMSelectionModuleAfterLikelihood(Context & ctx){
   
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
      metfilters_selection->add<TriggerSelection>("globalTightHalo2016Filter", "Flag_globalTightHalo2016Filter");
      metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
      if (!is_mc) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
      metfilters_selection->add<TriggerSelection>("chargedHadronTrackResolutionFilter", "Flag_chargedHadronTrackResolutionFilter"); 
      metfilters_selection->add<TriggerSelection>("muonBadTrackFilter", "Flag_muonBadTrackFilter");
   }
   
   ht_calculator.reset(new HTCalculator(ctx));
   partonht_calculator.reset(new PartonHT(ctx.get_handle<double>("parton_ht")));
      
   if (!is_mc){
   triggerMu50_sel.reset(new TriggerSelection("HLT_Mu50_v*")); 
   triggerTrkMu50_sel.reset(new TriggerSelection("HLT_TkMu50_v*")); 
   }
   
   //SF
   muonTRK_SF.reset(new MCMuonTrkScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_8_0_20/src/UHH2/common/data/general_eff_aeta_dr030e030_corr_ratio.txt", 0.0, "TRK"));
   muonIso_SF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_8_0_20/src/UHH2/common/data/MuonIso_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1", 1.0, "ISO"));
   muonID_SF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_8_0_20/src/UHH2/common/data/MuonID_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1", 1.0, "ID"));
   muonHLT2_SF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_8_0_20/src/UHH2/common/data/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root", "Mu50_PtEtaBins_Run273158_to_274093/efficienciesDATA/", 0.0, "HLT",true));
   muonHLT_SF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_8_0_20/src/UHH2/common/data/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root", "Mu50_OR_TkMu50_PtEtaBins_Run274954_to_276097/efficienciesDATA/", 0.0, "HLT",true));
   lumi_tot = string2double(ctx.get("target_lumi"));
   lumi1 = 622.;//0.622/fb in 2016 data
   lumi2 = lumi_tot - lumi1;
   toptag_SF.reset(new MCTopTagScaleFactor(ctx, "nominal", "h_heptopjets_WP3_wobtag", 1.5));
   elecID_SF.reset(new MCElecScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_8_0_20/src/UHH2/common/data/egammaEffi.txt_SF2D_looseId.root", 0.0, "ID"));
   elecGsf_SF.reset(new MCElecScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_8_0_20/src/UHH2/common/data/egammaEffi.txt_SF2D_GSFtrack.root", 0.0, "Gsf"));

   //// OBJ CLEANING
   muon_cleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(53., 2.4),MuonIso(0.15))));
   muon_cleaner_resolved.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), PtEtaCut(30., 2.1),MuonIso(0.15))));
   jetmuon_cleaner.reset(new JetMuonOverlapRemoval(0.4));
   jetmuon_cleaner_resolved.reset(new JetMuonOverlapRemoval(0.5));
      ////
   ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
   h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
   h_isoweight=ctx.get_handle<float>("weight_sfmu_ISO");
   h_isoweight_down=ctx.get_handle<float>("weight_sfmu_ISO_down");
   h_isoweight_up=ctx.get_handle<float>("weight_sfmu_ISO_up");
                                      
   h_tagweight=ctx.get_handle<float>("weight_toptag");
   h_tagweight_down=ctx.get_handle<float>("weight_toptag_down");
   h_tagweight_up=ctx.get_handle<float>("weight_toptag_up");
   
   isSemiLept = ctx.get("dataset_version") == "TTbar_SemiLept";
   isDiLept = ctx.get("dataset_version") == "TTbar_Dilept";
   isOther = ctx.get("dataset_version") == "TTbar_other";

   IsWJetsInclusive = ctx.get("dataset_version") =="WJetsToLNu";
   IsWJetsHT = ctx.get("dataset_version") =="WJets_LNu_HT-100To200" || ctx.get("dataset_version") =="WJets_LNu_HT-200To400" || ctx.get("dataset_version") =="WJets_LNu_HT-400To600" || ctx.get("dataset_version") =="WJets_LNu_HT-600To800" || ctx.get("dataset_version") =="WJets_LNu_HT-800To1200";
   //// EVENT SELECTION
   bool muon(false), elec(false);
   const std::string channel(ctx.get("channel", ""));
   if(channel == "muon") muon = true;
   else if(channel == "electron") elec = true;
   else throw std::runtime_error("undefined argument for 'channel' key in xml file (must be 'muon' or 'electron'): "+channel);
   
   const MuonId muonid = AndId<Muon>(PtEtaCut  (53., 2.4), MuonIDMedium_ICHEP());
  
   lep1_sel.reset(new AndSelection(ctx));
   if(muon){
      lep1_sel->add<NMuonSelection>("muoN == 1", 1, 1, muonid);
      //lep1_sel->add<NElectronSelection>("eleN == 0", 0, 0);
   }
   else if(elec){
      lep1_sel->add<NMuonSelection>("muoN == 0", 0, 0);
      lep1_sel->add<NElectronSelection>("eleN == 1", 1, 1);
   }
   
   lepVeto_sel.reset(new AndSelection(ctx));
   if(muon){
      lepVeto_sel->add<NMuonSelection>("muoN == 1", 1, 1);  
      lepVeto_sel->add<NMuonSelection>("muoN == 0", 1, 1, muonid);
      lepVeto_sel->add<NElectronSelection>("eleN == 0", 0, 0);
   }
   else if(elec){
      lepVeto_sel->add<NMuonSelection>("muoN == 0", 0, 0);
      lepVeto_sel->add<NElectronSelection>("eleN == 1", 1, 1);
   }
   JetId Btag = CSVBTag(0.8);  
   btag_sel.reset(new NJetSelection(1, -1, Btag));
   
   mtlep100_sel.reset(new MTlepCut(100., std::numeric_limits<double>::infinity()));
   selection_Met160.reset(new METCut(160., std::numeric_limits<double>::infinity()));
   partonht_sel.reset(new PartonHTCut(ctx,100.));

   h_flag_twodcut = ctx.get_handle<bool>("flag_twodcut");
   h_likelihood = ctx.get_handle<double>("likelihood");
   h_recneutrino = ctx.get_handle<LorentzVector>("rec_neutrino");
   h_bjets = ctx.get_handle<Jet>("bjet");
   
   const TopJetId HTTTopJetId_pteta = PtEtaCut(150., 2.4); 
   const TopJetId HTTTopJetId_WP3_wobtag = AndId<TopJet>(HEPTopTagV2(85,280,0.47), Tau32groomed(0.97));
   heptoptagevent_WP3_wobtag_sel.reset(new HandleSelection<int>(ctx, "n_heptopjets_WP3_wobtag", 1, 999));
   //heptt collection korrigieren
   if (is_mc){
   heptopjet_corrector.reset(new GenericTopJetCorrector(ctx, JERFiles::Spring16_25ns_L123_AK8PFchs_MC,"patJetsHepTopTagCHSPacked_daughters"));
   heptopjet_correctorL2L3.reset(new GenericTopJetCorrector(ctx, JERFiles::Spring16_25ns_L23_AK8PFchs_MC,"h_heptopjets_pteta"));
   }
   else{
      heptopjet_corrector.reset(new GenericTopJetCorrector(ctx, JERFiles::Spring16_25ns_L123_AK8PFchs_DATA,"patJetsHepTopTagCHSPacked_daughters"));
      heptopjet_correctorL2L3.reset(new GenericTopJetCorrector(ctx, JERFiles::Spring16_25ns_L23_AK8PFchs_DATA,"h_heptopjets_pteta"));
   }
   //mit collections aufpassen!! HIER ALLE RICHTIg? WERDEN AUCH IN HISTS >> VERWENDET
   collectionprod_heptt_pteta.reset(new CollectionProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "h_heptopjets_pteta", TopJetId(HTTTopJetId_pteta))); 
   collectionsizeprod_heptt_pteta.reset(new CollectionSizeProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "n_heptopjets_pteta", TopJetId(HTTTopJetId_pteta)));  
   collectionprod_heptt_WP3_wobtag.reset(new CollectionProducer<TopJet>(ctx, "h_heptopjets_pteta", "h_heptopjets_WP3_wobtag", TopJetId(HTTTopJetId_WP3_wobtag))); 
   hepttleptoncleaner.reset(new HepTTLeptonDeltaRCleaner(ctx, 1.5));
   collectionsizeprod_heptt_WP3_wobtag.reset(new CollectionSizeProducer<TopJet>(ctx, "h_heptopjets_WP3_wobtag", "n_heptopjets_WP3_wobtag", TopJetId(HTTTopJetId_WP3_wobtag)));  
   jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut( 50., 2.4))));

   //reolved analysis: tightmuon, 
   const MuonId muonid_tight = AndId<Muon>(PtEtaCut  (30., 2.1), MuonIDTight(), MuonIso(0.15));
   const MuonId muonid_tight_iso = AndId<Muon>(PtEtaCut  (10., 2.1), MuonIDLoose(), MuonIso(0.25));
   lep1_resolved_sel.reset(new AndSelection(ctx));
   lep1_resolved_sel->add<NMuonSelection>("muoN == 1", 1, 1, muonid_tight);
   lep1_resolved_sel->add<NElectronSelection>("eleN == 0", 0, 0);
   mtlep160_sel.reset(new MTlepCut(160., std::numeric_limits<double>::infinity()));
   jet_cleaner_resolved.reset(new JetCleaner(ctx, AndId<Jet>(JetPFID(JetPFID::WP_LOOSE),PtEtaCut(30., 2.4)))); 
   jetmuon_cleaner_resolved.reset(new JetMuonOverlapRemoval(0.5));
   JetId Btag_tight = CSVBTag(0.8);
   btagtight_sel.reset(new NJetSelection(1, -1, Btag_tight));
   jet_resolved_sel.reset(new NJetSelection(3, -1, JetId(PtEtaCut( 30., 2.4))));
   MT2W_sel.reset(new MT2WCut(200.));
   jetmetdphijet2_sel.reset(new METJetDPhiCut(1.2, 2));
   lepVeto_resolved_sel.reset(new AndSelection(ctx));
   lepVeto_resolved_sel->add<NMuonSelection>("muoN == 1", 1, 1, muonid_tight_iso);  
   lepVeto_resolved_sel->add<NMuonSelection>("muoN == 0", 1, 1, muonid_tight);
   lepVeto_resolved_sel->add<NElectronSelection>("eleN == 0", 0, 0);
   ele_cleaner_resolved.reset(new ElectronCleaner(AndId<Electron>(PtEtaSCCut(10., 2.5), ElectronID_Spring16_loose)));
  
   //// HISTS
   lumihists.reset(new LuminosityHists(ctx, "lumi"));
   input_h.reset(new ttDMSelectionHists(ctx, "input"));
   filter_h.reset(new ttDMSelectionHists(ctx, "filter"));
  
   likelihood_met_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met"));
   met_h.reset(new ttDMSelectionHists(ctx, "met"));
   electrons_met_h.reset(new ElectronHists(ctx,"electrons_met"));
   jets_met_h.reset(new JetHists(ctx,"jets_met"));
   muons_met_h.reset(new MuonHists(ctx,"muons_met"));
   events_met_h.reset(new EventHists(ctx,"events_met"));
   topjets_met_h.reset(new TopJetHists(ctx,"topjets_met",4,"patJetsHepTopTagCHSPacked_daughters"));
   genhists_met_h.reset(new ttDMGenHists(ctx,"genhists_met"));
   ttbargenhists_met_h.reset(new TTbarGenHists(ctx,"ttbargenhists_met"));

   nonisolep_h.reset(new ttDMSelectionHists(ctx, "nonisolep"));
   electrons_nonisolep_h.reset(new ElectronHists(ctx,"electrons_nonisolep"));
   jets_nonisolep_h.reset(new JetHists(ctx,"jets_nonisolep"));
   muons_nonisolep_h.reset(new MuonHists(ctx,"muons_nonisolep"));
   events_nonisolep_h.reset(new EventHists(ctx,"events_nonisolep"));
   topjets_nonisolep_h.reset(new TopJetHists(ctx,"topjets_nonisolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   genhists_nonisolep_h.reset(new  ttDMGenHists(ctx, "genhists_nonisolep_h"));

   isolep_h.reset(new ttDMSelectionHists(ctx, "isolep"));
   electrons_isolep_h.reset(new ElectronHists(ctx,"electrons_isolep"));
   jets_isolep_h.reset(new JetHists(ctx,"jets_isolep"));
   muons_isolep_h.reset(new MuonHists(ctx,"muons_isolep"));
   events_isolep_h.reset(new EventHists(ctx,"events_isolep"));
   topjets_isolep_h.reset(new TopJetHists(ctx,"topjets_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_tagged_WP3_wobtag_h.reset(new TopJetHists(ctx,"topjets_tagged_WP3_wobtag_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_tagged_WP3_wobtag_h->set_TopJetId(HTTTopJetId_WP3_wobtag);
   
   semileptcontrolregion_h.reset(new ttDMSelectionHists(ctx,"semileptcontrolregion"));
   electrons_semileptcontrolregion_h.reset(new ElectronHists(ctx,"electrons_semileptcontrolregion"));
   jets_semileptcontrolregion_h.reset(new JetHists(ctx,"jets_semileptcontrolregion"));
   muons_semileptcontrolregion_h.reset(new MuonHists(ctx,"muons_semileptcontrolregion"));
   events_semileptcontrolregion_h.reset(new EventHists(ctx,"events_semileptcontrolregion"));
   topjets_semileptcontrolregion_h.reset(new TopJetHists(ctx,"topjets_semileptcontrolregion", 4,"h_heptopjets_WP3_wobtag"));
   likelihood_semileptcontrolregion_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_semileptcontrolregion"));

   Wjetscontrolregion_h.reset(new ttDMSelectionHists(ctx,"Wjetscontrolregion"));
   electrons_Wjetscontrolregion_h.reset(new ElectronHists(ctx,"electrons_Wjetscontrolregion"));
   jets_Wjetscontrolregion_h.reset(new JetHists(ctx,"jets_Wjetscontrolregion"));
   muons_Wjetscontrolregion_h.reset(new MuonHists(ctx,"muons_Wjetscontrolregion"));
   events_Wjetscontrolregion_h.reset(new EventHists(ctx,"events_Wjetscontrolregion"));
   topjets_Wjetscontrolregion_h.reset(new TopJetHists(ctx,"topjets_Wjetscontrolregion", 4,"patJetsHepTopTagCHSPacked_daughters"));
   likelihood_Wjetscontrolregion_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_Wjetscontrolregion"));

   signalregion_iso_h.reset(new ttDMSelectionHists(ctx,"signalregion_iso"));
   electrons_signalregion_iso_h.reset(new ElectronHists(ctx,"electrons_signalregion_iso"));
   jets_signalregion_iso_h.reset(new JetHists(ctx,"jets_signalregion_iso"));
   muons_signalregion_iso_h.reset(new MuonHists(ctx,"muons_signalregion_iso"));
   events_signalregion_iso_h.reset(new EventHists(ctx,"events_signalregion_iso"));
   topjets_signalregion_iso_h.reset(new TopJetHists(ctx,"topjets_signalregion_iso", 4,"patJetsHepTopTagCHSPacked_daughters"));
   likelihood_signalregion_iso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_signalregion_iso"));

   signalregion_iso_Met160_h.reset(new ttDMSelectionHists(ctx,"signalregion_iso_Met160"));
   electrons_signalregion_iso_Met160_h.reset(new ElectronHists(ctx,"electrons_signalregion_iso_Met160"));
   jets_signalregion_iso_Met160_h.reset(new JetHists(ctx,"jets_signalregion_iso_Met160"));
   muons_signalregion_iso_Met160_h.reset(new MuonHists(ctx,"muons_signalregion_iso_Met160"));
   events_signalregion_iso_Met160_h.reset(new EventHists(ctx,"events_signalregion_iso_Met160"));
   topjets_signalregion_iso_Met160_h.reset(new TopJetHists(ctx,"topjets_signalregion_iso_Met160", 4,"patJetsHepTopTagCHSPacked_daughters"));
   likelihood_signalregion_iso_Met160_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_signalregion_iso_Met160"));

   btag_h.reset(new ttDMSelectionHists(ctx,"btag"));
   electrons_btag_h.reset(new ElectronHists(ctx,"electrons_btag"));
   jets_btag_h.reset(new JetHists(ctx,"jets_btag"));
   muons_btag_h.reset(new MuonHists(ctx,"muons_btag"));
   events_btag_h.reset(new EventHists(ctx,"events_btag"));
   topjets_btag_h.reset(new TopJetHists(ctx,"topjets_btag", 4,"patJetsHepTopTagCHSPacked_daughters"));
   likelihood_btag_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btag"));

   mtlepsel_h.reset(new ttDMSelectionHists(ctx,"mtlepsel"));
   electrons_mtlepsel_h.reset(new ElectronHists(ctx,"electrons_mtlepsel"));
   jets_mtlepsel_h.reset(new JetHists(ctx,"jets_mtlepsel"));
   muons_mtlepsel_h.reset(new MuonHists(ctx,"muons_mtlepsel"));
   events_mtlepsel_h.reset(new EventHists(ctx,"events_mtlepsel"));
   topjets_mtlepsel_h.reset(new TopJetHists(ctx,"topjets_mtlepsel", 4,"patJetsHepTopTagCHSPacked_daughters"));
   likelihood_mtlepsel_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_mtlepsel"));

   likelihood_preselection_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_preselection"));
   preselection_h.reset(new ttDMSelectionHists(ctx,"preselection"));
   electrons_preselection_h.reset(new ElectronHists(ctx,"electrons_preselection"));
   jets_preselection_h.reset(new JetHists(ctx,"jets_preselection"));
   topjets_preselection_h.reset(new TopJetHists(ctx,"topjets_preselection", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_preselection_h.reset(new MuonHists(ctx,"muons_preselection"));
   events_preselection_h.reset(new EventHists(ctx,"events_preselection"));
   genhists_preselection_h.reset(new ttDMGenHists(ctx,"genhists_preselection"));
   ttbargenhists_preselection_h.reset(new TTbarGenHists(ctx, "ttbargenhists_preselection"));

   likelihood_preselection_Met160_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_preselection_Met160"));
   preselection_Met160_h.reset(new ttDMSelectionHists(ctx,"preselection_Met160"));
   electrons_preselection_Met160_h.reset(new ElectronHists(ctx,"electrons_preselection_Met160"));
   jets_preselection_Met160_h.reset(new JetHists(ctx,"jets_preselection_Met160"));
   topjets_preselection_Met160_h.reset(new TopJetHists(ctx,"topjets_preselection_Met160", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_preselection_Met160_h.reset(new MuonHists(ctx,"muons_preselection_Met160"));
   events_preselection_Met160_h.reset(new EventHists(ctx,"events_preselection_Met160"));
   genhists_preselection_Met160_h.reset(new ttDMGenHists(ctx,"genhists_preselection_Met160"));
   ttbargenhists_preselection_Met160_h.reset(new TTbarGenHists(ctx, "ttbargenhists_preselection_Met160"));

   leptonveto_h.reset(new ttDMSelectionHists(ctx,"leptonveto"));
   electrons_leptonveto_h.reset(new ElectronHists(ctx,"electrons_leptonveto"));
   jets_leptonveto_h.reset(new JetHists(ctx,"jets_leptonveto"));
   topjets_leptonveto_h.reset(new TopJetHists(ctx,"topjets_leptonveto", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_leptonveto_h.reset(new MuonHists(ctx,"muons_leptonveto"));
   events_leptonveto_h.reset(new EventHists(ctx,"events_leptonveto"));
   genhists_leptonveto_h.reset(new ttDMGenHists(ctx,"genhists_leptonveto"));
   likelihood_leptonveto_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_leptonveto"));
   ttbargenhists_leptonveto_h.reset(new TTbarGenHists(ctx, "ttbargenhists_leptonveto"));

   dileptcontrolregion_h.reset(new ttDMSelectionHists(ctx,"dileptcontrolregion"));
   electrons_dileptcontrolregion_h.reset(new ElectronHists(ctx,"electrons_dileptcontrolregion"));
   jets_dileptcontrolregion_h.reset(new JetHists(ctx,"jets_dileptcontrolregion"));
   topjets_dileptcontrolregion_h.reset(new TopJetHists(ctx,"topjets_dileptcontrolregion", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_dileptcontrolregion_h.reset(new MuonHists(ctx,"muons_dileptcontrolregion"));
   events_dileptcontrolregion_h.reset(new EventHists(ctx,"events_dileptcontrolregion"));
   genhists_dileptcontrolregion_h.reset(new ttDMGenHists(ctx,"genhists_dileptcontrolregion"));
   likelihood_dileptcontrolregion_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_dileptcontrolregion"));
   ttbargenhists_dileptcontrolregion_h.reset(new TTbarGenHists(ctx, "ttbargenhists_dileptcontrolregion"));

   electrons_signalregion_noniso_h.reset(new ElectronHists(ctx,"electrons_signalregion_noniso"));
   jets_signalregion_noniso_h.reset(new JetHists(ctx,"jets_signalregion_noniso"));
   topjets_signalregion_noniso_h.reset(new TopJetHists(ctx,"topjets_signalregion_noniso", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_signalregion_noniso_h.reset(new MuonHists(ctx,"muons_signalregion_noniso"));
   events_signalregion_noniso_h.reset(new EventHists(ctx,"events_signalregion_noniso"));
   likelihood_signalregion_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_signalregion_noniso"));
   signalregion_noniso_h.reset(new ttDMSelectionHists(ctx,"signalregion_noniso_h"));
  
   combination_h.reset(new ttDMSelectionHists(ctx,"combination"));
   resolved_h.reset(new ttDMSelectionHists(ctx,"resolved"));
   resolved_preselection_h.reset(new ttDMSelectionHists(ctx,"resolved_preselection"));
   DeltaRMuonJet_sel.reset(new DeltaRMuonJet(0.2));

   h_btag_SF.reset(new BTagMCEfficiencyHists(ctx,"btag_SF",CSVBTag::WP_MEDIUM));
   apply_btag_SF.reset(new MCBTagScaleFactor(ctx,CSVBTag::WP_MEDIUM,"jets","central","mujets","incl","MCBtagEfficiencies"));
}

bool ttDMSelectionModuleAfterLikelihood::process(Event & event){
   event.set(h_isoweight, 1.0);
   event.set(h_isoweight_down, 1.0);
   event.set(h_isoweight_up, 1.0);

   event.set(h_tagweight, 1.0);
   event.set(h_tagweight_down, 1.0);
   event.set(h_tagweight_up, 1.0);

   if (is_mc) {
      double w = event.weight;
      double phi_rew = (1.01+0.28*TMath::Cos(1.06*event.met->phi()-0.28));
      event.weight=w*phi_rew;
   }
   isResolved=false;
   isBoosted=false;
   if (is_mc) partonht_calculator->process(event);
   if (IsWJetsInclusive) {
      bool pass_partonht= partonht_sel->passes(event);
      if (!pass_partonht) return false;
   }      
   if (IsWJetsHT) {
      bool pass_partonht= partonht_sel->passes(event);
      if (pass_partonht) return false;
   }
   input_h->fill(event);
   
 
   if (is_mc) {
      ttgenprod->process(event);
      TTbarGen ttbargen = event.get(h_ttbargen);
      if (isSemiLept && !ttbargen.IsSemiLeptonicDecay()) return false;
      if (isDiLept && (ttbargen.DecayChannel() < 4 || ttbargen.DecayChannel() >9)) return false;
      if (isOther && !(ttbargen.DecayChannel()==0 || ttbargen.DecayChannel() == 10)) return false;
   }
  
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
   
   if (!is_mc){
      bool pass_Mu50 = triggerMu50_sel->passes(event);
      if (event.run<274954) {
         if (!pass_Mu50) return false;
      }
      else {
         bool pass_TrkMu50 = triggerTrkMu50_sel->passes(event);
         if (!(pass_Mu50||pass_TrkMu50)) return false;
      }
   }

   jet_cleaner_resolved->process(event);
   
   //DATEN MUESSEN NOCH KORRIGIERT WERDEN!
    //correct topjets with L1L2L3
    heptopjet_corrector->process(event);
    //apply eta and pt requirement on topjets
    collectionprod_heptt_pteta->process(event);
    collectionsizeprod_heptt_pteta->process(event);
  
    //correct topjets with L2L3 only (HTTV2 requires cut on L2L3-corrected mass) 
    heptopjet_correctorL2L3->process(event);
    collectionprod_heptt_WP3_wobtag->process(event); 
    //hepttleptoncleaner->process(event); //be careful with lepton veto
    collectionsizeprod_heptt_WP3_wobtag->process(event);
    bool pass_HEPTT_WP3_wobtag = heptoptagevent_WP3_wobtag_sel->passes(event);
    //histograms filled with uncorrected mass
    
   //re-correct topjets with L1L2L3
    heptopjet_corrector->process(event);
    
    std::unique_ptr< std::vector<Jet> >  uncleaned_jets   (new std::vector<Jet>   (*event.jets));
    std::unique_ptr< std::vector<Muon> > uncleaned_muons(new std::vector<Muon>(*event.muons));
    std::unique_ptr< std::vector<Electron> > uncleaned_electrons(new std::vector<Electron>(*event.electrons));
    
    ele_cleaner_resolved->process(event); 
    bool pass_lepVeto_resolved = lepVeto_resolved_sel->passes(event);
    muon_cleaner_resolved->process(event); 
    bool lepresolved = lep1_resolved_sel->passes(event);
    bool pass_Met160 = selection_Met160->passes(event);
    if (pass_lepVeto_resolved && lepresolved)
       {
          jetmuon_cleaner_resolved->process(event);
          bool jet_resolved = jet_resolved_sel->passes(event);
          if (jet_resolved){
             bool btagtight = btagtight_sel->passes(event);
             if (btagtight && pass_Met160){
                resolved_preselection_h->fill(event);
                bool jetmetdphijet2= jetmetdphijet2_sel->passes(event);
                bool mtlep160 = mtlep160_sel->passes(event);
                bool MT2W = MT2W_sel->passes(event);
                if (mtlep160 && jetmetdphijet2 && MT2W){
                   isResolved=true;
                   resolved_h->fill(event);
                }
             }
          }
       }
    if (isResolved && !pass_HEPTT_WP3_wobtag) combination_h->fill(event);
    
    event.jets->clear();
    event.jets->reserve(uncleaned_jets->size());
    for(const auto & j : *uncleaned_jets) event.jets->push_back(j); 
    sort_by_pt<Jet>(*event.jets);
    
    event.muons->clear();
    event.muons->reserve(uncleaned_muons->size());
    for(const auto & j : *uncleaned_muons) event.muons->push_back(j); 
    sort_by_pt<Muon>(*event.muons);
    
    event.electrons->clear();
    event.electrons->reserve(uncleaned_electrons->size());
    for(const auto & j : *uncleaned_electrons) event.electrons->push_back(j); 
    sort_by_pt<Electron>(*event.electrons);
    
    bool jet2 = jet2_sel->passes(event);
    if (!jet2) return false;
    //two jets
    bool pass_lep1 = lep1_sel->passes(event);
    if (!pass_lep1) return false;   
    
    //due to unknown efficiency in 2016 data, skip events with: muon pt>500 GeV in the endcap
    bool good_muon(0);
    Muon lep1 = event.muons->at(0);   //so richtig???
    good_muon = (!(lep1.pt()> 500 && fabs(lep1.eta())> 1.2));
    if(!good_muon) return false;
    if (event.isRealData) lumihists->fill(event);
    
    elecGsf_SF->process(event);
    elecID_SF->process(event);

    muonID_SF->process(event);
    muonIso_SF->process(event);
    muonTRK_SF->process(event);
    double w_wo_HLT = event.weight;
    // muon-HLT eff
    muonHLT2_SF->process(event);
    double w1 = event.weight;
    event.weight = w_wo_HLT;
    muonHLT_SF->process(event);
    double w2 = event.weight;
    double w = (lumi1*w1+lumi2*w2)/(lumi_tot);
    event.weight = w;
    
    met_h->fill(event);  //MET >50GeV
    electrons_met_h->fill(event);
    jets_met_h->fill(event);
    muons_met_h->fill(event);
    events_met_h->fill(event);
    topjets_met_h->fill(event);
    if(!event.isRealData) genhists_met_h->fill(event);
    likelihood_met_h->fill(event); 
    ttbargenhists_met_h->fill(event);
    
    bool pass_twodcut = event.get(h_flag_twodcut);
    if(pass_twodcut)                                               //pre-selection control plots 
       {
          likelihood_preselection_h->fill(event);        
          preselection_h->fill(event);
          electrons_preselection_h->fill(event);
          muons_preselection_h->fill(event);
          events_preselection_h->fill(event);
          jets_preselection_h->fill(event);
          topjets_preselection_h->fill(event);
          genhists_preselection_h->fill(event);
          ttbargenhists_preselection_h->fill(event);
      }
    if(pass_twodcut &&  pass_Met160)                                               //pre-selection control plots 
       {
          likelihood_preselection_Met160_h->fill(event);        
          preselection_Met160_h->fill(event);
          electrons_preselection_Met160_h->fill(event);
          muons_preselection_Met160_h->fill(event);
          events_preselection_Met160_h->fill(event);
          jets_preselection_Met160_h->fill(event);
          topjets_preselection_Met160_h->fill(event);
          genhists_preselection_Met160_h->fill(event);
          ttbargenhists_preselection_Met160_h->fill(event);
       }
    
    //veto additional leptons
    bool pass_lepVeto = lepVeto_sel->passes(event);
    if (pass_lepVeto){
       likelihood_leptonveto_h->fill(event);        
       leptonveto_h->fill(event);
       electrons_leptonveto_h->fill(event);
       muons_leptonveto_h->fill(event);
       events_leptonveto_h->fill(event);
       jets_leptonveto_h->fill(event);
       topjets_leptonveto_h->fill(event);
       genhists_leptonveto_h->fill(event);
       ttbargenhists_leptonveto_h->fill(event);
    }
    
    //---------------------------------------------------- control regions --------------------------------------  

    bool pass_btagsel = btag_sel->passes(event);
    bool pass_mtlep100 = mtlep100_sel->passes(event);
    
    // lepton + jets control region
    if (pass_twodcut && pass_btagsel && !pass_mtlep100 && pass_lepVeto){              //muon Iso fehlt
       apply_btag_SF->process(event); 
       semileptcontrolregion_h->fill(event);
       electrons_semileptcontrolregion_h->fill(event);
       jets_semileptcontrolregion_h->fill(event);
       muons_semileptcontrolregion_h->fill(event);
       events_semileptcontrolregion_h->fill(event);
       topjets_semileptcontrolregion_h->fill(event); 
       likelihood_semileptcontrolregion_h->fill(event);
    }
    //dilepton control region
    if (pass_twodcut && !pass_lepVeto && pass_btagsel && pass_mtlep100) {           //muon Iso fehlt, ele SF
       apply_btag_SF->process(event);
       likelihood_dileptcontrolregion_h->fill(event);        
       dileptcontrolregion_h->fill(event);
       electrons_dileptcontrolregion_h->fill(event);
       muons_dileptcontrolregion_h->fill(event);
       events_dileptcontrolregion_h->fill(event);
       jets_dileptcontrolregion_h->fill(event);
       topjets_dileptcontrolregion_h->fill(event);
       genhists_dileptcontrolregion_h->fill(event);
       ttbargenhists_dileptcontrolregion_h->fill(event);
    }  
    //W+jets control region
    if (pass_twodcut && !pass_btagsel && pass_mtlep100 && pass_lepVeto){             //muon Iso fehlt
       apply_btag_SF->process(event);
       Wjetscontrolregion_h->fill(event);
       electrons_Wjetscontrolregion_h->fill(event);
       jets_Wjetscontrolregion_h->fill(event);
       muons_Wjetscontrolregion_h->fill(event);
       events_Wjetscontrolregion_h->fill(event);
       topjets_Wjetscontrolregion_h->fill(event); 
       likelihood_Wjetscontrolregion_h->fill(event);
    }
    
    //---------------------------------------------------- signal selections --------------------------------------  

     if (!pass_lepVeto) return false;
     if(!pass_btagsel) return false;
     apply_btag_SF->process(event);
     btag_h->fill(event);
     electrons_btag_h->fill(event);
     jets_btag_h->fill(event);
     muons_btag_h->fill(event);
     events_btag_h->fill(event);
     topjets_btag_h->fill(event); 
     likelihood_btag_h->fill(event);
     
     if (!pass_mtlep100) return false;
     mtlepsel_h->fill(event);
     electrons_mtlepsel_h->fill(event);
     jets_mtlepsel_h->fill(event);
     muons_mtlepsel_h->fill(event);
     events_mtlepsel_h->fill(event);
     topjets_mtlepsel_h->fill(event); 
     likelihood_mtlepsel_h->fill(event);
     
   //---------------------------------------------------- boosted analysis -------------------------------------- 

    std::unique_ptr< std::vector<Muon> > cleaned_muons_noIso(new std::vector<Muon>(*event.muons)); //only muons, to do: add electrons
    muon_cleaner_iso->process(event);
   
    sort_by_pt<Muon>(*event.muons);
    pass_lep1 =  lep1_sel->passes(event);
    if (pass_lep1)    
      {
         bool pass_DeltaRMuonJet_sel = DeltaRMuonJet_sel-> passes(event);
         if (!pass_DeltaRMuonJet_sel) 
            {
               isolep_h->fill(event);                       
               electrons_isolep_h->fill(event);
               jets_isolep_h->fill(event);
               muons_isolep_h->fill(event);
               events_isolep_h->fill(event);
               topjets_isolep_h->fill(event);
            }
         //cleaning 
         jetmuon_cleaner->process(event);
         bool jet2 = jet2_sel->passes(event);
         if (!jet2) return false;
         
         pass_btagsel = btag_sel->passes(event);
         if(!pass_btagsel) return false;
         if (!pass_HEPTT_WP3_wobtag) return false;
         //std::cout<<"found a top tag"<<std::endl;
         toptag_SF->process(event);
         
         //h_btag_SF->fill(event);
         
         signalregion_iso_h->fill(event);
         electrons_signalregion_iso_h->fill(event);
         jets_signalregion_iso_h->fill(event);
         muons_signalregion_iso_h->fill(event);
         events_signalregion_iso_h->fill(event);
         topjets_signalregion_iso_h->fill(event); 
         likelihood_signalregion_iso_h->fill(event);
         isBoosted =true;
                  
         if (pass_Met160){
            signalregion_iso_Met160_h->fill(event);
            electrons_signalregion_iso_Met160_h->fill(event);
            jets_signalregion_iso_Met160_h->fill(event);
            muons_signalregion_iso_Met160_h->fill(event);
            events_signalregion_iso_Met160_h->fill(event);
            topjets_signalregion_iso_Met160_h->fill(event); 
            likelihood_signalregion_iso_Met160_h->fill(event);
         }
         
         return true;
      }
  
    event.muons->clear();
    event.muons->reserve(cleaned_muons_noIso->size());
    for(auto & muon : *cleaned_muons_noIso) event.muons->push_back(muon); 
    sort_by_pt<Muon>(*event.muons);
   
    nonisolep_h->fill(event);
    genhists_nonisolep_h->fill(event);
    electrons_nonisolep_h->fill(event);
    jets_nonisolep_h->fill(event);
    muons_nonisolep_h->fill(event);
    events_nonisolep_h->fill(event);
    topjets_nonisolep_h->fill(event);
    
    if(!pass_twodcut) return false;
    //h_btag_SF->fill(event);
    signalregion_noniso_h->fill(event);
    electrons_signalregion_noniso_h->fill(event);
    muons_signalregion_noniso_h->fill(event);
    events_signalregion_noniso_h->fill(event);
    jets_signalregion_noniso_h->fill(event);
    topjets_signalregion_noniso_h->fill(event);
    likelihood_signalregion_noniso_h->fill(event);        
    
   return true;
}



UHH2_REGISTER_ANALYSIS_MODULE(ttDMSelectionModuleAfterLikelihood)

