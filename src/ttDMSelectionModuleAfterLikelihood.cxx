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

#include "UHH2/TTbarDM/include/ttDMSemiLeptonicSelections.h"
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicUtils.h"
#include "UHH2/TTbarDM/include/ttDMSelectionHists.h"
#include "UHH2/TTbarDM/include/ttDMGenHists.h"
#include "UHH2/TTbarDM/include/ttDMReconstruction_Likelihood.h"
#include "UHH2/TTbarDM/include/ttDMReconstructionHists_Likelihood.h"

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

class ttDMSelectionModuleAfterLikelihood: public AnalysisModule {
 public:
  explicit ttDMSelectionModuleAfterLikelihood(Context & ctx);
  virtual bool process(Event & event) override;

 private:
   Event::Handle<int> h_flag_heptoptagevent;
   Event::Handle<bool> h_flag_twodcut;
   Event::Handle<double> h_likelihood;
   Event::Handle<LorentzVector> h_recneutrino;
   Event::Handle<Jet> h_bjets;
  
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

   // cleaners
   std::unique_ptr<MuonCleaner> muon_cleaner_iso;
   std::unique_ptr<GenericTopJetCorrector> heptopjet_corrector;
   std::unique_ptr<GenericTopJetCorrector> heptopjet_correctorL2L3;
      
   // selections
   std::unique_ptr<Selection> heptoptagevent_WP1_sel;
   std::unique_ptr<Selection> heptoptagevent_WP3_sel;
   std::unique_ptr<Selection> heptoptagevent_WP1_wobtag_sel;
   std::unique_ptr<Selection> heptoptagevent_WP3_wobtag_sel;
   std::unique_ptr<Selection> pv_sel;
   std::unique_ptr<Selection> MET_sel;
   std::unique_ptr<Selection> MET100_sel;
   std::unique_ptr<Selection> MET160_sel,btag_sel;
   std::unique_ptr<AndSelection> lep1_sel, lepVeto_sel;
   std::unique_ptr<Selection> jetmetdphijet2_sel;
   std::unique_ptr<Selection> DeltaPhiMetNeutrino_sel2;
   std::unique_ptr<Selection> DeltaPhiMetNeutrino_sel15;
   std::unique_ptr<Hists> mtlep_noniso_h;
   std::unique_ptr<Hists> likelihood_noiso_h;
   std::unique_ptr<AnalysisModule> muo_med_noniso_SF;

   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_pteta;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_pteta;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_WP3;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_WP3;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_WP1;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_WP1;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_WP3_wobtag;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_WP3_wobtag;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_WP1_wobtag;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_WP1_wobtag;
   //std::unique_ptr<AnalysisModule> muo_triggerSF;
   std::unique_ptr<AnalysisModule> ttgenprod;
   std::unique_ptr<AnalysisModule> hepttleptoncleaner;
   
   // hists
   std::unique_ptr<Hists> lumihists;
   std::unique_ptr<Hists> filter_h;
   std::unique_ptr<Hists> input_h;
   std::unique_ptr<Hists> met_h, electrons_met_h,jets_met_h,muons_met_h,events_met_h,topjets_met_h, genhists_met_h;
   std::unique_ptr<Hists> nonisolep_h, electrons_nonisolep_h,jets_nonisolep_h,muons_nonisolep_h,events_nonisolep_h,topjets_nonisolep_h;
   std::unique_ptr<Hists> isolep_h, electrons_isolep_h,jets_isolep_h,muons_isolep_h,events_isolep_h,topjets_isolep_h;
   std::unique_ptr<TopJetHists> topjets_isolep_tagged_WP1_h, topjets_isolep_tagged_WP3_h, topjets_isolep_notag_WP1_h, topjets_isolep_notag_WP3_h;
   std::unique_ptr<TopJetHists> topjets_isolep_tagged_WP1_wobtag_h, topjets_isolep_tagged_WP3_wobtag_h;
   std::unique_ptr<Hists> heptoptagevent_WP3_h,electrons_heptoptagevent_WP3_h,jets_heptoptagevent_WP3_h,muons_heptoptagevent_WP3_h,events_heptoptagevent_WP3_h,topjets_heptoptagevent_WP3_h;
   std::unique_ptr<Hists> heptoptagevent_WP3_wobtag_h,electrons_heptoptagevent_WP3_wobtag_h,jets_heptoptagevent_WP3_wobtag_h,muons_heptoptagevent_WP3_wobtag_h,events_heptoptagevent_WP3_wobtag_h,topjets_heptoptagevent_WP3_wobtag_h;
   std::unique_ptr<Hists> mediumMET_h,electrons_mediumMET_h,jets_mediumMET_h,topjets_mediumMET_h,muons_mediumMET_h,events_mediumMET_h, likelihood_mediumMET_h;
   std::unique_ptr<TopJetHists>  topjets_heptoptagevent_WP3_tagged_h;
   std::unique_ptr<TopJetHists>  topjets_heptoptagevent_WP3_wobtag_tagged_h;
   std::unique_ptr<Hists> twodcut_h,electrons_twodcut_h, jets_twodcut_h,topjets_twodcut_h,events_twodcut_h,muons_twodcut_h;
   std::unique_ptr<Hists> preselection_h,electrons_preselection_h, jets_preselection_h,topjets_preselection_h,events_preselection_h,muons_preselection_h, genhists_preselection_h;
   std::unique_ptr<Hists> leptonveto_h,electrons_leptonveto_h, jets_leptonveto_h,topjets_leptonveto_h,events_leptonveto_h,muons_leptonveto_h, genhists_leptonveto_h;
   std::unique_ptr<Hists> dileptcontrolregion_h,electrons_dileptcontrolregion_h, jets_dileptcontrolregion_h,topjets_dileptcontrolregion_h,events_dileptcontrolregion_h,muons_dileptcontrolregion_h, genhists_dileptcontrolregion_h;
   std::unique_ptr<Hists> heptoptagevent_WP1_h;
   std::unique_ptr<Hists> heptoptagevent_WP1_wobtag_h;
   std::unique_ptr<Hists> jetmetdphi_noniso_h;
   std::unique_ptr<Hists> likelihood_preselection_h;
   std::unique_ptr<Hists> likelihood_leptonveto_h;
   std::unique_ptr<Hists> likelihood_dileptcontrolregion_h;
   std::unique_ptr<Hists> likelihood_met_h;
   std::unique_ptr<Hists> likelihood_twodcut_h;
   std::unique_ptr<Selection> jetmetdphi_sel;
   std::unique_ptr<Selection> mtlep_sel;
   std::unique_ptr<Selection> thirdjet_sel, mtlep100_sel, mtlep150_sel;
   std::unique_ptr<Selection> jet123metdphi_sel;
   std::unique_ptr<Selection> DeltaPhiMetNeutrino_sel;
   std::unique_ptr<Selection> DMMETSel;
   std::unique_ptr<Hists> genhists_nonisolep_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_h,likelihood_heptoptagevent_WP1_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_wobtag_h,likelihood_heptoptagevent_WP1_wobtag_h;
 
   std::unique_ptr<Hists> thirdjet_h,electrons_thirdjet_h,jets_thirdjet_h,topjets_thirdjet_h,muons_thirdjet_h,events_thirdjet_h,likelihood_thirdjet_h;
   std::unique_ptr<Hists> electrons_mtlep_h,jets_mtlep_h,topjets_mtlep_h,muons_mtlep_h,events_mtlep_h,likelihood_mtlep_h, mtlep_h;
   std::unique_ptr<Hists> jetmetdphijet2_h,electrons_jetmetdphijet2_h,jets_jetmetdphijet2_h,topjets_jetmetdphijet2_h,muons_jetmetdphijet2_h,events_jetmetdphijet2_h,likelihood_jetmetdphijet2_h;
   std::unique_ptr<Hists> electrons_btagsel_h,jets_btagsel_h,topjets_btagsel_h,muons_btagsel_h,events_btagsel_h,ttbargenhists_btagsel_h;

   std::unique_ptr<Hists> likelihood_thirdjet_mtlep_h;
   std::unique_ptr<Hists> thirdjet_mtlep_h;
 
   std::unique_ptr<Hists> jet123metdphi_h;
   std::unique_ptr<Hists> hists_likelihood_jet123metdphi_h;
   
   std::unique_ptr<Hists> ttbargenhists_preselection_h;
   std::unique_ptr<Hists> ttbargenhists_leptonveto_h;
   std::unique_ptr<Hists> ttbargenhists_dileptcontrolregion_h;
   std::unique_ptr<Hists> ttbargenhists_heptoptagevent_WP3_h;
   std::unique_ptr<Hists> ttbargenhists_met_h;
   uhh2::Event::Handle<TTbarGen> h_ttbargen;
   
   std::unique_ptr<Hists> btagsel_h;
   std::unique_ptr<Hists> likelihood_btagsel;
   std::unique_ptr<Hists> met160_h;
   std::unique_ptr<Hists> likelihood_met160;
   std::unique_ptr<Hists> likelihood_met160_btag;
   std::unique_ptr<Hists> jets_met160_btag_h,topjets_met160_btag_h,muons_met160_btag_h,events_met160_btag_h,met160_btag_h;
  
   std::unique_ptr<Hists> jets_met50_btag_h;
   std::unique_ptr<Hists> topjets_met50_btag_h;
   std::unique_ptr<Hists> muons_met50_btag_h;
   std::unique_ptr<Hists> events_met50_btag_h;
   std::unique_ptr<Hists> met50_btag_h;
   std::unique_ptr<Hists> likelihood_met50_btag;

   std::unique_ptr<Hists> jets_met100_btag_h;
   std::unique_ptr<Hists> topjets_met100_btag_h;
   std::unique_ptr<Hists> muons_met100_btag_h;
   std::unique_ptr<Hists> events_met100_btag_h;
   std::unique_ptr<Hists> met100_btag_h;
   std::unique_ptr<Hists> likelihood_met100_btag;
   
   std::unique_ptr<Hists> btag_noniso_h;
   std::unique_ptr<Hists> electrons_btag_noniso_h;
   std::unique_ptr<Hists> jets_btag_noniso_h;
   std::unique_ptr<Hists> topjets_btag_noniso_h;
   std::unique_ptr<Hists> muons_btag_noniso_h;
   std::unique_ptr<Hists> events_btag_noniso_h;
   
   std::unique_ptr<Hists> met100_btag_mtlep100_h;
   std::unique_ptr<Hists> likelihood_met100_btag_mtlep100;
   std::unique_ptr<Hists> HEPTT_WP3_noniso_h;
   std::unique_ptr<Hists> likelihood_WP3_noniso_h;
   std::unique_ptr<Hists> likelihood_btag_noniso_h;
   std::unique_ptr<Hists> nomet_noDeltaPhi_h;
   std::unique_ptr<Hists> likelihood_nomet_noDeltaPhi_h;
   std::unique_ptr<Hists> nobtag_h;
   std::unique_ptr<Hists> likelihood_nobtag_h;
   std::unique_ptr<Hists> nothirdjet_h;
   std::unique_ptr<Hists> likelihood_nothirdjet_h;
   std::unique_ptr<Hists> noDeltaPhiJets_h;
   std::unique_ptr<Hists> likelihood_noDeltaPhiJets_h;
   std::unique_ptr<Hists> nomet_h;
   std::unique_ptr<Hists> likelihood_nomet_h;
   std::unique_ptr<Hists> nomtlep;
   std::unique_ptr<Hists> likelihood_nomtlep;
   std::unique_ptr<Hists> noDeltaPhiMetNeutrino;
   std::unique_ptr<Hists> likelihood_noDeltaPhiMetNeutrino;

   std::unique_ptr<Hists> noDeltaPhiMetNeutrino_met160;
   std::unique_ptr<Hists> likelihood_noDeltaPhiMetNeutrino_met160;
   std::unique_ptr<Hists> DeltaPhiMetNeutrinotight;
   std::unique_ptr<Hists> likelihood_DeltaPhiMetNeutrinotight;
   std::unique_ptr<Hists> noDeltaPhiMetNeutrino_DMMET;
   std::unique_ptr<Hists> likelihood_noDeltaPhiMetNeutrino_DMMET;
   std::unique_ptr<Hists> noDeltaPhiMetNeutrino_MT2W;
   std::unique_ptr<Hists> likelihood_noDeltaPhiMetNeutrino_MT2W;
   std::unique_ptr<Hists> noDeltaPhiMetNeutrino_DeltaPhi;
   std::unique_ptr<Hists> likelihood_noDeltaPhiMetNeutrino_DeltaPhi;
   
   std::unique_ptr<Hists> heptoptagevent_WP3_wobtag_noDeltaPhiMET_h;
   std::unique_ptr<Hists> electrons_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h;
   std::unique_ptr<Hists> jets_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h;
   std::unique_ptr<Hists> muons_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h;
   std::unique_ptr<Hists> events_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h;
   std::unique_ptr<Hists> topjets_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h;

  std::unique_ptr<Hists>  heptoptagevent_WP3_wobtag_nobtag_h;
   std::unique_ptr<Hists> electrons_heptoptagevent_WP3_wobtag_nobtag_h;
   std::unique_ptr<Hists> jets_heptoptagevent_WP3_wobtag_nobtag_h;
   std::unique_ptr<Hists> muons_heptoptagevent_WP3_wobtag_nobtag_h;
   std::unique_ptr<Hists> events_heptoptagevent_WP3_wobtag_nobtag_h;
   std::unique_ptr<Hists> topjets_heptoptagevent_WP3_wobtag_nobtag_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_wobtag_nobtag_h;
   
   std::unique_ptr<Hists> heptoptagevent_WP3_wobtag_nomtlep_h;
   std::unique_ptr<Hists> electrons_heptoptagevent_WP3_wobtag_nomtlep_h;
   std::unique_ptr<Hists> jets_heptoptagevent_WP3_wobtag_nomtlep_h;
   std::unique_ptr<Hists> muons_heptoptagevent_WP3_wobtag_nomtlep_h;
   std::unique_ptr<Hists> events_heptoptagevent_WP3_wobtag_nomtlep_h;
   std::unique_ptr<Hists> topjets_heptoptagevent_WP3_wobtag_nomtlep_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_wobtag_nomtlep_h;
   std::unique_ptr<Selection>  DeltaPhiTaggedJetTopLep_sel, MT2W_sel,tightbtag_sel;
   std::unique_ptr<Hists> HEPTT_WP3_noniso_mtlep_h,likelihood_WP3_noniso_mtlep_h;
   std::unique_ptr<Hists> mtlep200_h;
   std::unique_ptr<Hists> likelihood_mtlep200;
   std::unique_ptr<Hists> btagtight_h;
   std::unique_ptr<Hists> likelihood_btagtight_h;
   bool isSemiLept;
   bool isDiLept;
   bool isOther;

   std::unique_ptr<Hists> WJetsControlRegion_noniso_h,electrons_WJetsControlRegion_noniso_h,jets_WJetsControlRegion_noniso_h, topjets_WJetsControlRegion_noniso_h, muons_WJetsControlRegion_noniso_h, events_WJetsControlRegion_noniso_h, likelihood_WJetsControlRegion_noniso_h;

   std::unique_ptr<Hists> DileptControlRegion_noniso_h, electrons_DileptControlRegion_noniso_h, jets_DileptControlRegion_noniso_h,topjets_DileptControlRegion_noniso_h, muons_DileptControlRegion_noniso_h, events_DileptControlRegion_noniso_h, likelihood_DileptControlRegion_noniso_h;

   std::unique_ptr<Hists>  SemileptControlRegion_noniso_h,electrons_SemileptControlRegion_noniso_h,jets_SemileptControlRegion_noniso_h,topjets_SemileptControlRegion_noniso_h,muons_SemileptControlRegion_noniso_h,events_SemileptControlRegion_noniso_h,likelihood_SemileptControlRegion_noniso_h;
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
      metfilters_selection->add<TriggerSelection>("CSCTightHalo2015Filter", "Flag_CSCTightHalo2015Filter");
      metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
      metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
      metfilters_selection->add<TriggerSelection>("chargedHadronTrackResolutionFilter", "Flag_chargedHadronTrackResolutionFilter"); 
      metfilters_selection->add<TriggerSelection>("muonBadTrackFilter", "Flag_muonBadTrackFilter");
   }
   
   ht_calculator.reset(new HTCalculator(ctx));
   muo_med_noniso_SF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_7_6_3/src/UHH2/common/data/MuonID_Z_RunCD_Reco76X_Feb15.root","MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1", 1., "id"));
   // muo_triggerSF.reset(new MCMuonScaleFactor(ctx, "/afs/desy.de/user/m/mameyer/xxl-af-cms/CMSSW_7_4_15_patch1/src/UHH2/common/data/SingleMuonTrigger_Z_RunD_Reco74X_Nov20.root", "Mu45_eta2p1_PtEtaBins", 1., "trg")); 
   //do we have to veto tight or medium muons? NO SF for Iso applied is approx. 1

   //Eff(b)=3%
   
   JetId subBtag = CSVBTag(0.46f);
   const TopJetId HTTTopJetId_pteta = PtEtaCut(150., 2.4); 
   const TopJetId HTTTopJetId_WP3 = AndId<TopJet>(HEPTopTagV2(85,280,0.47, subBtag), Tau32(0.97)); 
   const TopJetId HTTTopJetId_WP3_wobtag = AndId<TopJet>(HEPTopTagV2(85,280,0.47), Tau32(0.97));
   const TopJetId HTTTopJetId_WP1= AndId<TopJet>(HEPTopTagV2(110,185,0.2, subBtag), Tau32(0.93));
   const TopJetId HTTTopJetId_WP1_wobtag= AndId<TopJet>(HEPTopTagV2(110,185,0.2), Tau32(0.93));
   //// OBJ CLEANING
   muon_cleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuonIDMedium(), PtEtaCut(47., 2.5),MuonIso(0.15))));
      ////
   ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
   h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

   isSemiLept = ctx.get("dataset_version") == "TTbar_SemiLept";
   isDiLept = ctx.get("dataset_version") == "TTbar_Dilept";
   isOther = ctx.get("dataset_version") == "TTbar_other";
   
   

   //// EVENT SELECTION
   bool muon(false), elec(false);
   const std::string channel(ctx.get("channel", ""));
   if(channel == "muon") muon = true;
   else if(channel == "electron") elec = true;
   else throw std::runtime_error("undefined argument for 'channel' key in xml file (must be 'muon' or 'electron'): "+channel);
   
   const MuonId muonid = AndId<Muon>(PtEtaCut  (47., 2.1), MuonIDMedium());
 
   lep1_sel.reset(new AndSelection(ctx));
   if(muon){
      lep1_sel->add<NMuonSelection>("muoN == 1", 1, 1, muonid);
      lep1_sel->add<NElectronSelection>("eleN == 0", 0, 0);
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
   JetId BtagTight = CSVBTag(0.935);
   btag_sel.reset(new NJetSelection(1, -1, Btag));
   tightbtag_sel.reset(new NJetSelection(1, -1, BtagTight));

   DMMETSel.reset(new DMMETSelection(ctx,250));
   MET_sel.reset(new METCut(320., std::numeric_limits<double>::infinity()));
   MET100_sel.reset(new METCut(100., std::numeric_limits<double>::infinity()));
   MET160_sel.reset(new METCut(160., std::numeric_limits<double>::infinity()));
   mtlep_sel.reset(new MTlepCut(200., std::numeric_limits<double>::infinity()));
   mtlep100_sel.reset(new MTlepCut(100., std::numeric_limits<double>::infinity()));
   mtlep150_sel.reset(new MTlepCut(150., std::numeric_limits<double>::infinity()));
   jetmetdphi_sel.reset(new METJetDPhiCut(1.2, 2));
   thirdjet_sel.reset(new NJetSelection(3, -1, JetId(PtEtaCut( 50., 2.4))));
   jet123metdphi_sel.reset(new METJetDPhiCut(1.4,3));
   jetmetdphijet2_sel.reset(new METJetDPhiCut(1.2, 2));
   DeltaPhiMetNeutrino_sel.reset(new DeltaPhiMetNeutrino(ctx,1));
   DeltaPhiMetNeutrino_sel15.reset(new DeltaPhiMetNeutrino(ctx,1.5));
   DeltaPhiMetNeutrino_sel2.reset(new DeltaPhiMetNeutrino(ctx,2));
   //h_flag_heptoptagevent = ctx.declare_event_output<int>("flag_heptoptagevent");
   h_flag_twodcut = ctx.get_handle<bool>("flag_twodcut");
   h_likelihood = ctx.get_handle<double>("likelihood");
   h_recneutrino = ctx.get_handle<LorentzVector>("rec_neutrino");
   h_bjets = ctx.get_handle<Jet>("bjet");
   MT2W_sel.reset(new MT2WCut(200.));
   DeltaPhiTaggedJetTopLep_sel.reset(new DeltaPhiTaggedJetTopLep(ctx,1.0));
   
   heptoptagevent_WP3_sel.reset(new HandleSelection<int>(ctx, "n_heptopjets_WP3", 1, 999));
   heptoptagevent_WP1_sel.reset(new HandleSelection<int>(ctx, "n_heptopjets_WP1", 1, 999));
   heptoptagevent_WP3_wobtag_sel.reset(new HandleSelection<int>(ctx, "n_heptopjets_WP3_wobtag", 1, 999));
   heptoptagevent_WP1_wobtag_sel.reset(new HandleSelection<int>(ctx, "n_heptopjets_WP1_wobtag", 1, 999));
   
   //heptt collection korrigieren
   if (is_mc){
   heptopjet_corrector.reset(new GenericTopJetCorrector(ctx, JERFiles::Fall15_25ns_L123_AK8PFchs_MC,"patJetsHepTopTagCHSPacked_daughters"));
   heptopjet_correctorL2L3.reset(new GenericTopJetCorrector(ctx, JERFiles::Fall15_25ns_L23_AK8PFchs_MC,"h_heptopjets_pteta"));
   }
   else{
      heptopjet_corrector.reset(new GenericTopJetCorrector(ctx, JERFiles::Fall15_25ns_L123_AK8PFchs_DATA,"patJetsHepTopTagCHSPacked_daughters"));
      heptopjet_correctorL2L3.reset(new GenericTopJetCorrector(ctx, JERFiles::Fall15_25ns_L23_AK8PFchs_DATA,"h_heptopjets_pteta"));
   }
   //mit collections aufpassen!! HIER ALLE RICHTIg? WERDEN AUCH IN HISTS >> VERWENDET

   collectionprod_heptt_pteta.reset(new CollectionProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "h_heptopjets_pteta", TopJetId(HTTTopJetId_pteta))); 
   collectionsizeprod_heptt_pteta.reset(new CollectionSizeProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "n_heptopjets_pteta", TopJetId(HTTTopJetId_pteta)));  

   collectionprod_heptt_WP3_wobtag.reset(new CollectionProducer<TopJet>(ctx, "h_heptopjets_pteta", "h_heptopjets_WP3_wobtag", TopJetId(HTTTopJetId_WP3_wobtag))); 
   hepttleptoncleaner.reset(new HepTTLeptonDeltaRCleaner(ctx, 1.5));
   collectionsizeprod_heptt_WP3_wobtag.reset(new CollectionSizeProducer<TopJet>(ctx, "h_heptopjets_WP3_wobtag", "n_heptopjets_WP3_wobtag", TopJetId(HTTTopJetId_WP3_wobtag)));  
   
   collectionprod_heptt_WP1_wobtag.reset(new CollectionProducer<TopJet>(ctx, "h_heptopjets_pteta", "h_heptopjets_WP1_wobtag", TopJetId(HTTTopJetId_WP1_wobtag))); 
   collectionsizeprod_heptt_WP1_wobtag.reset(new CollectionSizeProducer<TopJet>(ctx, "h_heptopjets_pteta", "n_heptopjets_WP1_wobtag", TopJetId(HTTTopJetId_WP1_wobtag))); 
   
   collectionprod_heptt_WP3.reset(new CollectionProducer<TopJet>(ctx, "h_heptopjets_pteta", "h_heptopjets_WP3", TopJetId(HTTTopJetId_WP3))); 
   collectionsizeprod_heptt_WP3.reset(new CollectionSizeProducer<TopJet>(ctx, "h_heptopjets_pteta", "n_heptopjets_WP3", TopJetId(HTTTopJetId_WP3)));  
   
   collectionprod_heptt_WP1.reset(new CollectionProducer<TopJet>(ctx, "h_heptopjets_pteta", "h_heptopjets_WP1", TopJetId(HTTTopJetId_WP1))); 
   collectionsizeprod_heptt_WP1.reset(new CollectionSizeProducer<TopJet>(ctx, "h_heptopjets_pteta", "n_heptopjets_WP1", TopJetId(HTTTopJetId_WP1))); 
   
   //// HISTS
   lumihists.reset(new LuminosityHists(ctx, "lumi"));
   input_h.reset(new ttDMSelectionHists(ctx, "input"));
   filter_h.reset(new ttDMSelectionHists(ctx, "filter"));
   met_h.reset(new ttDMSelectionHists(ctx, "met"));
   
   electrons_met_h.reset(new ElectronHists(ctx,"electrons_met"));
   jets_met_h.reset(new JetHists(ctx,"jets_met"));
   muons_met_h.reset(new MuonHists(ctx,"muons_met"));
   events_met_h.reset(new EventHists(ctx,"events_met"));
   topjets_met_h.reset(new TopJetHists(ctx,"topjets_met",4,"patJetsHepTopTagCHSPacked_daughters"));
   genhists_met_h.reset(new ttDMGenHists(ctx,"genhists_met"));

   nonisolep_h.reset(new ttDMSelectionHists(ctx, "nonisolep"));
   electrons_nonisolep_h.reset(new ElectronHists(ctx,"electrons_nonisolep"));
   jets_nonisolep_h.reset(new JetHists(ctx,"jets_nonisolep"));
   muons_nonisolep_h.reset(new MuonHists(ctx,"muons_nonisolep"));
   events_nonisolep_h.reset(new EventHists(ctx,"events_nonisolep"));
   topjets_nonisolep_h.reset(new TopJetHists(ctx,"topjets_nonisolep",4,"patJetsHepTopTagCHSPacked_daughters"));

   isolep_h.reset(new ttDMSelectionHists(ctx, "isolep"));
   electrons_isolep_h.reset(new ElectronHists(ctx,"electrons_isolep"));
   jets_isolep_h.reset(new JetHists(ctx,"jets_isolep"));
   muons_isolep_h.reset(new MuonHists(ctx,"muons_isolep"));
   events_isolep_h.reset(new EventHists(ctx,"events_isolep"));
   topjets_isolep_h.reset(new TopJetHists(ctx,"topjets_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_tagged_WP1_h.reset(new TopJetHists(ctx,"topjets_tagged_WP1_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_tagged_WP1_h->set_TopJetId(HTTTopJetId_WP1);
   topjets_isolep_tagged_WP3_h.reset(new TopJetHists(ctx,"topjets_tagged_WP3_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_tagged_WP3_h->set_TopJetId(HTTTopJetId_WP3);
   topjets_isolep_notag_WP1_h.reset(new TopJetHists(ctx,"topjets_notag_WP1_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_notag_WP1_h->set_TopJetId(HTTTopJetId_WP1);
   topjets_isolep_notag_WP3_h.reset(new TopJetHists(ctx,"topjets_notag_WP3_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_notag_WP3_h->set_TopJetId(HTTTopJetId_WP3);
   topjets_isolep_tagged_WP1_wobtag_h.reset(new TopJetHists(ctx,"topjets_tagged_WP1_wobtag_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_tagged_WP1_wobtag_h->set_TopJetId(HTTTopJetId_WP1_wobtag);
   topjets_isolep_tagged_WP3_wobtag_h.reset(new TopJetHists(ctx,"topjets_tagged_WP3_wobtag_isolep",4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_isolep_tagged_WP3_wobtag_h->set_TopJetId(HTTTopJetId_WP3_wobtag);
   
   heptoptagevent_WP3_h.reset(new ttDMSelectionHists(ctx,"heptoptagevent_WP3"));
   electrons_heptoptagevent_WP3_h.reset(new ElectronHists(ctx,"electrons_heptoptagevent_WP3"));
   jets_heptoptagevent_WP3_h.reset(new JetHists(ctx,"jets_heptoptagevent_WP3"));
   muons_heptoptagevent_WP3_h.reset(new MuonHists(ctx,"muons_heptoptagevent_WP3"));
   events_heptoptagevent_WP3_h.reset(new EventHists(ctx,"events_heptoptagevent_WP3"));
   topjets_heptoptagevent_WP3_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3", 4,"patJetsHepTopTagCHSPacked_daughters"));
  
   heptoptagevent_WP3_wobtag_h.reset(new ttDMSelectionHists(ctx,"heptoptagevent_WP3_wobtag"));
   electrons_heptoptagevent_WP3_wobtag_h.reset(new ElectronHists(ctx,"electrons_heptoptagevent_WP3_wobtag"));
   jets_heptoptagevent_WP3_wobtag_h.reset(new JetHists(ctx,"jets_heptoptagevent_WP3_wobtag"));
   muons_heptoptagevent_WP3_wobtag_h.reset(new MuonHists(ctx,"muons_heptoptagevent_WP3_wobtag"));
   events_heptoptagevent_WP3_wobtag_h.reset(new EventHists(ctx,"events_heptoptagevent_WP3_wobtag"));
   topjets_heptoptagevent_WP3_wobtag_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_wobtag", 4,"h_heptopjets_WP3_wobtag"));
  
   heptoptagevent_WP3_wobtag_nomtlep_h.reset(new ttDMSelectionHists(ctx,"heptoptagevent_WP3_wobtag_nomtlep"));
   electrons_heptoptagevent_WP3_wobtag_nomtlep_h.reset(new ElectronHists(ctx,"electrons_heptoptagevent_WP3_wobtag_nomtlep"));
   jets_heptoptagevent_WP3_wobtag_nomtlep_h.reset(new JetHists(ctx,"jets_heptoptagevent_WP3_wobtag_nomtlep"));
   muons_heptoptagevent_WP3_wobtag_nomtlep_h.reset(new MuonHists(ctx,"muons_heptoptagevent_WP3_wobtag_nomtlep"));
   events_heptoptagevent_WP3_wobtag_nomtlep_h.reset(new EventHists(ctx,"events_heptoptagevent_WP3_wobtag_nomtlep"));
   topjets_heptoptagevent_WP3_wobtag_nomtlep_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_wobtag_nomtlep", 4,"patJetsHepTopTagCHSPacked_daughters"));
   
   heptoptagevent_WP3_wobtag_nobtag_h.reset(new ttDMSelectionHists(ctx,"heptoptagevent_WP3_wobtag_nobtag"));
   electrons_heptoptagevent_WP3_wobtag_nobtag_h.reset(new ElectronHists(ctx,"electrons_heptoptagevent_WP3_wobtag_nobtag"));
   jets_heptoptagevent_WP3_wobtag_nobtag_h.reset(new JetHists(ctx,"jets_heptoptagevent_WP3_wobtag_nobtag"));
   muons_heptoptagevent_WP3_wobtag_nobtag_h.reset(new MuonHists(ctx,"muons_heptoptagevent_WP3_wobtag_nobtag"));
   events_heptoptagevent_WP3_wobtag_nobtag_h.reset(new EventHists(ctx,"events_heptoptagevent_WP3_wobtag_nobtag"));
   topjets_heptoptagevent_WP3_wobtag_nobtag_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_wobtag_nobtag", 4,"patJetsHepTopTagCHSPacked_daughters"));
   
   heptoptagevent_WP3_wobtag_noDeltaPhiMET_h.reset(new ttDMSelectionHists(ctx,"heptoptagevent_WP3_wobtag_noDeltaPhiMET"));
   electrons_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h.reset(new ElectronHists(ctx,"electrons_heptoptagevent_WP3_wobtag_noDeltaPhiMET"));
   jets_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h.reset(new JetHists(ctx,"jets_heptoptagevent_WP3_wobtag_noDeltaPhiMET"));
   muons_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h.reset(new MuonHists(ctx,"muons_heptoptagevent_WP3_wobtag_noDeltaPhiMET"));
   events_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h.reset(new EventHists(ctx,"events_heptoptagevent_WP3_wobtag_noDeltaPhiMET"));
   topjets_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_wobtag_noDeltaPhiMET", 4,"patJetsHepTopTagCHSPacked_daughters"));
  
   
   twodcut_h.reset(new ttDMSelectionHists(ctx,"twodcut"));
   electrons_twodcut_h.reset(new ElectronHists(ctx,"electrons_twodcut"));
   jets_twodcut_h.reset(new JetHists(ctx,"jets_twodcut"));
   topjets_twodcut_h.reset(new TopJetHists(ctx,"topjets_twodcut", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_twodcut_h.reset(new MuonHists(ctx,"muons_twodcut"));
   events_twodcut_h.reset(new EventHists(ctx,"events_twodcut"));

   preselection_h.reset(new ttDMSelectionHists(ctx,"preselection"));
   electrons_preselection_h.reset(new ElectronHists(ctx,"electrons_preselection"));
   jets_preselection_h.reset(new JetHists(ctx,"jets_preselection"));
   topjets_preselection_h.reset(new TopJetHists(ctx,"topjets_preselection", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_preselection_h.reset(new MuonHists(ctx,"muons_preselection"));
   events_preselection_h.reset(new EventHists(ctx,"events_preselection"));
   genhists_preselection_h.reset(new ttDMGenHists(ctx,"genhists_preselection"));

   leptonveto_h.reset(new ttDMSelectionHists(ctx,"leptonveto"));
   electrons_leptonveto_h.reset(new ElectronHists(ctx,"electrons_leptonveto"));
   jets_leptonveto_h.reset(new JetHists(ctx,"jets_leptonveto"));
   topjets_leptonveto_h.reset(new TopJetHists(ctx,"topjets_leptonveto", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_leptonveto_h.reset(new MuonHists(ctx,"muons_leptonveto"));
   events_leptonveto_h.reset(new EventHists(ctx,"events_leptonveto"));
   genhists_leptonveto_h.reset(new ttDMGenHists(ctx,"genhists_leptonveto"));

   dileptcontrolregion_h.reset(new ttDMSelectionHists(ctx,"dileptcontrolregion"));
   electrons_dileptcontrolregion_h.reset(new ElectronHists(ctx,"electrons_dileptcontrolregion"));
   jets_dileptcontrolregion_h.reset(new JetHists(ctx,"jets_dileptcontrolregion"));
   topjets_dileptcontrolregion_h.reset(new TopJetHists(ctx,"topjets_dileptcontrolregion", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_dileptcontrolregion_h.reset(new MuonHists(ctx,"muons_dileptcontrolregion"));
   events_dileptcontrolregion_h.reset(new EventHists(ctx,"events_dileptcontrolregion"));
   genhists_dileptcontrolregion_h.reset(new ttDMGenHists(ctx,"genhists_dileptcontrolregion"));


   thirdjet_h.reset(new ttDMSelectionHists(ctx,"thirdjet"));
   electrons_thirdjet_h.reset(new ElectronHists(ctx,"electrons_thirdjet"));
   jets_thirdjet_h.reset(new JetHists(ctx,"jets_thirdjet"));
   topjets_thirdjet_h.reset(new TopJetHists(ctx,"topjets_thirdjet", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_thirdjet_h.reset(new MuonHists(ctx,"muons_thirdjet"));
   events_thirdjet_h.reset(new EventHists(ctx,"events_thirdjet"));

   electrons_btagsel_h.reset(new ElectronHists(ctx,"electrons_btagsel"));
   jets_btagsel_h.reset(new JetHists(ctx,"jets_btagsel"));
   topjets_btagsel_h.reset(new TopJetHists(ctx,"topjets_btagsel", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_btagsel_h.reset(new MuonHists(ctx,"muons_btagsel"));
   events_btagsel_h.reset(new EventHists(ctx,"events_btagsel"));

   jets_met160_btag_h.reset(new JetHists(ctx,"jets_met160_btag"));
   topjets_met160_btag_h.reset(new TopJetHists(ctx,"topjets_met160_btag", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_met160_btag_h.reset(new MuonHists(ctx,"muons_met160_btag"));
   events_met160_btag_h.reset(new EventHists(ctx,"events_met160_btag"));

   jets_met50_btag_h.reset(new JetHists(ctx,"jets_met50_btag"));
   topjets_met50_btag_h.reset(new TopJetHists(ctx,"topjets_met50_btag", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_met50_btag_h.reset(new MuonHists(ctx,"muons_met50_btag"));
   events_met50_btag_h.reset(new EventHists(ctx,"events_met50_btag"));
   met50_btag_h.reset(new ttDMSelectionHists(ctx, "met50_btag"));
   likelihood_met50_btag.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met50_btag"));

   jets_met100_btag_h.reset(new JetHists(ctx,"jets_met100_btag"));
   topjets_met100_btag_h.reset(new TopJetHists(ctx,"topjets_met100_btag", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_met100_btag_h.reset(new MuonHists(ctx,"muons_met100_btag"));
   events_met100_btag_h.reset(new EventHists(ctx,"events_met100_btag"));
   met100_btag_h.reset(new ttDMSelectionHists(ctx, "met100_btag"));
   likelihood_met100_btag.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met100_btag"));

   electrons_mtlep_h.reset(new ElectronHists(ctx,"electrons_mtlep"));
   jets_mtlep_h.reset(new JetHists(ctx,"jets_mtlep"));
   topjets_mtlep_h.reset(new TopJetHists(ctx,"topjets_mtlep", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_mtlep_h.reset(new MuonHists(ctx,"muons_mtlep"));
   events_mtlep_h.reset(new EventHists(ctx,"events_mtlep"));

   btag_noniso_h.reset(new ttDMSelectionHists(ctx,"btag_noniso"));
   electrons_btag_noniso_h.reset(new ElectronHists(ctx,"electrons_btag_noniso"));
   jets_btag_noniso_h.reset(new JetHists(ctx,"jets_btag_noniso"));
   topjets_btag_noniso_h.reset(new TopJetHists(ctx,"topjets_btag_noniso", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_btag_noniso_h.reset(new MuonHists(ctx,"muons_btag_noniso"));
   events_btag_noniso_h.reset(new EventHists(ctx,"events_btag_noniso"));
             
   likelihood_preselection_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_preselection"));
   likelihood_leptonveto_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_leptonveto"));
   likelihood_dileptcontrolregion_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_dileptcontrolregion"));
   likelihood_met_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met"));
   likelihood_twodcut_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_twodcut"));
   likelihood_noiso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noiso"));
   likelihood_heptoptagevent_WP3_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3"));
   likelihood_heptoptagevent_WP1_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP1"));
   likelihood_heptoptagevent_WP3_wobtag_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3_wobtag"));
   likelihood_heptoptagevent_WP3_wobtag_nomtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3_wobtag_nomtlep"));
   likelihood_heptoptagevent_WP3_wobtag_nobtag_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3_wobtag_nobtag"));
   likelihood_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3_wobtag_noDeltaPhiMET"));
   likelihood_heptoptagevent_WP1_wobtag_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP1_wobtag"));
   heptoptagevent_WP1_wobtag_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP1_wobtag"));
   heptoptagevent_WP1_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP1"));

   likelihood_thirdjet_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_thirdjet"));
   likelihood_mtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_mtlep"));
   likelihood_mediumMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_mediumMET"));
   mtlep_noniso_h.reset(new ttDMSelectionHists(ctx, "mtlep_noniso"));
   jetmetdphi_noniso_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi_noniso"));
   genhists_nonisolep_h.reset(new  ttDMGenHists(ctx, "genhists_nonisolep_h"));

   likelihood_thirdjet_mtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_thirdjet_mtlep"));
   thirdjet_mtlep_h.reset(new ttDMSelectionHists(ctx, "thirdjet_mtlep"));

   jet123metdphi_h.reset(new ttDMSelectionHists(ctx, "jet123metdphi_h"));
   hists_likelihood_jet123metdphi_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_jet123metdphi_h"));
  
   ttbargenhists_preselection_h.reset(new TTbarGenHists(ctx, "ttbargenhists_preselection"));
   ttbargenhists_dileptcontrolregion_h.reset(new TTbarGenHists(ctx, "ttbargenhists_dileptcontrolregion"));
   ttbargenhists_leptonveto_h.reset(new TTbarGenHists(ctx, "ttbargenhists_leptonveto"));
   ttbargenhists_heptoptagevent_WP3_h.reset(new TTbarGenHists(ctx, "ttbargenhists_heptoptagevent_WP3"));
   ttbargenhists_met_h.reset(new TTbarGenHists(ctx,"ttbargenhists_met"));
   ttbargenhists_btagsel_h.reset(new TTbarGenHists(ctx,"ttbargenhists_btagsel"));

   btagsel_h.reset(new ttDMSelectionHists(ctx, "btagsel"));
   likelihood_btagsel.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btagsel"));

   met160_h.reset(new ttDMSelectionHists(ctx, "met160"));
   likelihood_met160.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met160"));
   met160_btag_h.reset(new ttDMSelectionHists(ctx, "met160_btag"));
   likelihood_met160_btag.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met160_btag"));

   met100_btag_mtlep100_h.reset(new ttDMSelectionHists(ctx, "met100_btag_mtlep100"));
   likelihood_met100_btag_mtlep100.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met100_btag_mtlep100"));
   
   HEPTT_WP3_noniso_h.reset(new ttDMSelectionHists(ctx, "HEPTT_WP3_noniso"));
   likelihood_WP3_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_WP3_noniso"));
   HEPTT_WP3_noniso_mtlep_h.reset(new ttDMSelectionHists(ctx, "HEPTT_WP3_noniso_mtlep"));
   likelihood_WP3_noniso_mtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_WP3_noniso_mtlep"));
   likelihood_btag_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btag_noniso"));

    mediumMET_h.reset(new ttDMSelectionHists(ctx,"mediumMET"));
    electrons_mediumMET_h.reset(new ElectronHists(ctx,"electrons_mediumMET"));
    jets_mediumMET_h.reset(new JetHists(ctx,"jets_mediumMET"));
    topjets_mediumMET_h.reset(new TopJetHists(ctx,"topjets_mediumMET", 4, "patJetsHepTopTagCHSPacked_daughters"));
    muons_mediumMET_h.reset(new MuonHists(ctx,"muons_mediumMET"));
    events_mediumMET_h.reset(new EventHists(ctx,"events_mediumMET"));  

    nobtag_h.reset(new ttDMSelectionHists(ctx,"nobtag"));
    likelihood_nobtag_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_nobtag"));
    nothirdjet_h.reset(new ttDMSelectionHists(ctx,"nothirdjet"));
    likelihood_nothirdjet_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_nothirdjet"));
    noDeltaPhiJets_h.reset(new ttDMSelectionHists(ctx,"noDeltaPhiJets"));
    likelihood_noDeltaPhiJets_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noDeltaPhiJets"));
    nomet_h.reset(new ttDMSelectionHists(ctx,"nomet"));
    likelihood_nomet_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_nomet"));
    nomtlep.reset(new ttDMSelectionHists(ctx,"nomtlep"));
    likelihood_nomtlep.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_nomtlep"));
    noDeltaPhiMetNeutrino.reset(new ttDMSelectionHists(ctx,"noDeltaPhiMetNeutrino"));
    likelihood_noDeltaPhiMetNeutrino.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noDeltaPhiMetNeutrino"));

    noDeltaPhiMetNeutrino_met160.reset(new ttDMSelectionHists(ctx,"noDeltaPhiMetNeutrino_met160"));
    likelihood_noDeltaPhiMetNeutrino_met160.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noDeltaPhiMetNeutrino_met160"));
    DeltaPhiMetNeutrinotight.reset(new ttDMSelectionHists(ctx,"DeltaPhiMetNeutrinotight"));
    likelihood_DeltaPhiMetNeutrinotight.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_DeltaPhiMetNeutrinotight"));
    noDeltaPhiMetNeutrino_DMMET.reset(new ttDMSelectionHists(ctx,"noDeltaPhiMetNeutrino_DMMET"));
    likelihood_noDeltaPhiMetNeutrino_DMMET.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noDeltaPhiMetNeutrino_DMMET"));
    noDeltaPhiMetNeutrino_MT2W.reset(new ttDMSelectionHists(ctx,"noDeltaPhiMetNeutrino_MT2W"));
    likelihood_noDeltaPhiMetNeutrino_MT2W.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noDeltaPhiMetNeutrino_MT2W"));
    noDeltaPhiMetNeutrino_DeltaPhi.reset(new ttDMSelectionHists(ctx,"noDeltaPhiMetNeutrino_DeltaPhi"));
    likelihood_noDeltaPhiMetNeutrino_DeltaPhi.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noDeltaPhiMetNeutrino_DeltaPhi"));
    nomet_noDeltaPhi_h.reset(new ttDMSelectionHists(ctx,"nomet_noDeltaPhi"));
    likelihood_nomet_noDeltaPhi_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_nomet_noDeltaPhi"));
    mtlep_h.reset(new ttDMSelectionHists(ctx,"mtlep_h"));
    mtlep200_h.reset(new ttDMSelectionHists(ctx,"mtlep200"));
    likelihood_mtlep200.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_mtlep200"));
    btagtight_h.reset(new ttDMSelectionHists(ctx,"btagtight"));
    likelihood_btagtight_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btagtight"));

    WJetsControlRegion_noniso_h.reset(new ttDMSelectionHists(ctx,"WJetsControlRegion_noniso"));
    electrons_WJetsControlRegion_noniso_h.reset(new ElectronHists(ctx,"electrons_WJetsControlRegion_noniso"));
    jets_WJetsControlRegion_noniso_h.reset(new JetHists(ctx,"jets_WJetsControlRegion_noniso"));
    topjets_WJetsControlRegion_noniso_h.reset(new TopJetHists(ctx,"topjets_WJetsControlRegion_noniso", 4, "patJetsHepTopTagCHSPacked_daughters"));
    muons_WJetsControlRegion_noniso_h.reset(new MuonHists(ctx,"muons_WJetsControlRegion_noniso"));
    events_WJetsControlRegion_noniso_h.reset(new EventHists(ctx,"events_WJetsControlRegion_noniso"));
    likelihood_WJetsControlRegion_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_WJetsControlRegion_noniso"));

    DileptControlRegion_noniso_h.reset(new ttDMSelectionHists(ctx,"DileptControlRegion_noniso"));
    electrons_DileptControlRegion_noniso_h.reset(new ElectronHists(ctx,"electrons_DileptControlRegion_noniso"));
    jets_DileptControlRegion_noniso_h.reset(new JetHists(ctx,"jets_DileptControlRegion_noniso"));
    topjets_DileptControlRegion_noniso_h.reset(new TopJetHists(ctx,"topjets_DileptControlRegion_noniso", 4, "patJetsHepTopTagCHSPacked_daughters"));
    muons_DileptControlRegion_noniso_h.reset(new MuonHists(ctx,"muons_DileptControlRegion_noniso"));
    events_DileptControlRegion_noniso_h.reset(new EventHists(ctx,"events_DileptControlRegion_noniso"));
    likelihood_DileptControlRegion_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_DileptControlRegion_noniso"));

    SemileptControlRegion_noniso_h.reset(new ttDMSelectionHists(ctx,"SemileptControlRegion_noniso"));
    electrons_SemileptControlRegion_noniso_h.reset(new ElectronHists(ctx,"electrons_SemileptControlRegion_noniso"));
    jets_SemileptControlRegion_noniso_h.reset(new JetHists(ctx,"jets_SemileptControlRegion_noniso"));
    topjets_SemileptControlRegion_noniso_h.reset(new TopJetHists(ctx,"topjets_SemileptControlRegion_noniso", 4, "patJetsHepTopTagCHSPacked_daughters"));
    muons_SemileptControlRegion_noniso_h.reset(new MuonHists(ctx,"muons_SemileptControlRegion_noniso"));
    events_SemileptControlRegion_noniso_h.reset(new EventHists(ctx,"events_SemileptControlRegion_noniso"));
    likelihood_SemileptControlRegion_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_SemileptControlRegion_noniso"));
}

bool ttDMSelectionModuleAfterLikelihood::process(Event & event){
 
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
   muo_med_noniso_SF->process(event); //spaeter
   filter_h->fill(event); 
   bool pass_lep1 = lep1_sel->passes(event);
   if (!pass_lep1) return false;   
   if (event.isRealData) lumihists->fill(event);
 
   met_h->fill(event);  //MET >50GeV
   electrons_met_h->fill(event);
   jets_met_h->fill(event);
   muons_met_h->fill(event);
   events_met_h->fill(event);
   topjets_met_h->fill(event);
   if(!event.isRealData) genhists_met_h->fill(event);
   likelihood_met_h->fill(event); 
   ttbargenhists_met_h->fill(event);
    bool pass_met160 = MET160_sel->passes(event);
   bool pass_twodcut = event.get(h_flag_twodcut);
   if(pass_twodcut && pass_met160)                                               //pre-selection control plots 
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
   //----------------------------------------------------boosted analysis-------------------------------------- 
    //DATEN MUESSEN NOCH KORRIGIERT WERDEN!
    //correct topjets with L1L2L3
    heptopjet_corrector->process(event);
    //apply eta and pt requirement on topjets
    collectionprod_heptt_pteta->process(event);
    collectionsizeprod_heptt_pteta->process(event);
  
    //correct topjets with L2L3 only (HTTV2 requires cut on L2L3-corrected mass) 
    heptopjet_correctorL2L3->process(event);
  
    collectionprod_heptt_WP3->process(event); 
    collectionsizeprod_heptt_WP3->process(event);
    collectionprod_heptt_WP1->process(event); 
    collectionsizeprod_heptt_WP1->process(event);
    collectionprod_heptt_WP3_wobtag->process(event); 
    //hepttleptoncleaner->process(event); //be careful with lepton veto
    collectionsizeprod_heptt_WP3_wobtag->process(event);
    collectionprod_heptt_WP1_wobtag->process(event); 
    collectionsizeprod_heptt_WP1_wobtag->process(event);
 
    //event.set(h_flag_heptoptagevent, int(pass_HEPTT)); 
    bool pass_HEPTT_WP3 = heptoptagevent_WP3_sel->passes(event);
    bool pass_HEPTT_WP1 = heptoptagevent_WP1_sel->passes(event);
    bool pass_HEPTT_WP3_wobtag = heptoptagevent_WP3_wobtag_sel->passes(event);
    bool pass_HEPTT_WP1_wobtag = heptoptagevent_WP1_wobtag_sel->passes(event);
    
   //re-correct topjets with L1L2L3
    heptopjet_corrector->process(event);
    std::unique_ptr< std::vector<Muon> > cleaned_muons_noIso(new std::vector<Muon>(*event.muons)); //only muons, to do: add electrons
    muon_cleaner_iso->process(event);
    sort_by_pt<Muon>(*event.muons);
    pass_lep1 =  lep1_sel->passes(event);
    if (pass_lep1)    
      {
         isolep_h->fill(event);                       
         electrons_isolep_h->fill(event);
         jets_isolep_h->fill(event);
         muons_isolep_h->fill(event);
         events_isolep_h->fill(event);
         topjets_isolep_h->fill(event);
         topjets_isolep_tagged_WP1_h->fill(event);  //histograms filled with uncorrected mass
         topjets_isolep_tagged_WP3_h->fill(event);
        
         bool pass_btagsel = btag_sel->passes(event);
         bool pass_thirdjet_sel = thirdjet_sel->passes(event);
         bool pass_jet123metdphi_sel=false;
         if (pass_thirdjet_sel) pass_jet123metdphi_sel=jet123metdphi_sel->passes(event);
         bool pass_DeltaPhiMetNeutrino_sel=DeltaPhiMetNeutrino_sel->passes(event);
         bool pass_DeltaPhiMetNeutrino_sel2=DeltaPhiMetNeutrino_sel2->passes(event);
        
         bool pass_met100 = MET100_sel->passes(event);
         bool pass_mtlep100 = mtlep100_sel->passes(event);
         bool pass_met160= MET160_sel->passes(event);

         if (pass_HEPTT_WP1_wobtag && pass_lepVeto)  {
            heptoptagevent_WP1_wobtag_h->fill(event);
            likelihood_heptoptagevent_WP1_wobtag_h->fill(event); 
         }
         if (pass_HEPTT_WP1 && pass_btagsel && pass_thirdjet_sel && pass_jet123metdphi_sel && pass_DeltaPhiMetNeutrino_sel && pass_met100 && pass_mtlep100 && pass_lepVeto)  {
            heptoptagevent_WP1_h->fill(event);
            likelihood_heptoptagevent_WP1_h->fill(event); 
         }
         
         else topjets_isolep_notag_WP1_h->fill(event);  //histograms filled with uncorrected mass
        
         if (pass_HEPTT_WP3_wobtag && pass_btagsel && !pass_mtlep100 && pass_lepVeto)
            {
               heptoptagevent_WP3_wobtag_h->fill(event);
               electrons_heptoptagevent_WP3_wobtag_h->fill(event);
               jets_heptoptagevent_WP3_wobtag_h->fill(event);
               muons_heptoptagevent_WP3_wobtag_h->fill(event);
               events_heptoptagevent_WP3_wobtag_h->fill(event);
               topjets_heptoptagevent_WP3_wobtag_h->fill(event); 
               likelihood_heptoptagevent_WP3_wobtag_h->fill(event);
            }
        
          if (!pass_lepVeto && pass_HEPTT_WP3_wobtag && pass_btagsel && pass_mtlep100) {
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

         //NO DeltaPhiMetNeutrino
         if (pass_HEPTT_WP3_wobtag && pass_btagsel && pass_mtlep100 && pass_lepVeto)
            {
               heptoptagevent_WP3_wobtag_noDeltaPhiMET_h->fill(event);
               electrons_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h->fill(event);
               jets_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h->fill(event);
               muons_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h->fill(event);
               events_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h->fill(event);
               topjets_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h->fill(event); 
               likelihood_heptoptagevent_WP3_wobtag_noDeltaPhiMET_h->fill(event);
               
               bool pass_DMMETSelection = DMMETSel->passes(event);
               bool pass_DeltaPhi = DeltaPhiTaggedJetTopLep_sel->passes(event);
               bool pass_Mt2W = MT2W_sel->passes(event);
               if (pass_DMMETSelection){
                  noDeltaPhiMetNeutrino_DMMET->fill(event);
                  likelihood_noDeltaPhiMetNeutrino_DMMET->fill(event);   
               }
               if (pass_Mt2W){
                  noDeltaPhiMetNeutrino_MT2W->fill(event);
                  likelihood_noDeltaPhiMetNeutrino_MT2W->fill(event);   
               }
               if (pass_DeltaPhi){
                  noDeltaPhiMetNeutrino_DeltaPhi->fill(event);
                  likelihood_noDeltaPhiMetNeutrino_DeltaPhi->fill(event);  
               }
               bool pass_met160= MET160_sel->passes(event);
               if (pass_met160)
                  {
                     met160_h->fill(event);  
                     likelihood_met160->fill(event);
                  } 
               bool pass_mtlep200=mtlep_sel->passes(event);
               if(pass_mtlep200){
                  mtlep200_h->fill(event);  
                  likelihood_mtlep200->fill(event);
               }
               bool pass_btagtight= tightbtag_sel->passes(event);
               if (pass_btagtight){
                  btagtight_h->fill(event);  
                  likelihood_btagtight_h->fill(event);  
               }

            }
         //no btag 
         if (pass_HEPTT_WP3_wobtag &&  !pass_btagsel && pass_mtlep100 && pass_lepVeto)
            {
               heptoptagevent_WP3_wobtag_nobtag_h->fill(event);
               electrons_heptoptagevent_WP3_wobtag_nobtag_h->fill(event);
               jets_heptoptagevent_WP3_wobtag_nobtag_h->fill(event);
               muons_heptoptagevent_WP3_wobtag_nobtag_h->fill(event);
               events_heptoptagevent_WP3_wobtag_nobtag_h->fill(event);
               topjets_heptoptagevent_WP3_wobtag_nobtag_h->fill(event); 
               likelihood_heptoptagevent_WP3_wobtag_nobtag_h->fill(event);
            }
         //no mtlep
         if (pass_HEPTT_WP3_wobtag && pass_btagsel && pass_DeltaPhiMetNeutrino_sel && pass_lepVeto)
            {
               heptoptagevent_WP3_wobtag_nomtlep_h->fill(event);
               electrons_heptoptagevent_WP3_wobtag_nomtlep_h->fill(event);
               jets_heptoptagevent_WP3_wobtag_nomtlep_h->fill(event);
               muons_heptoptagevent_WP3_wobtag_nomtlep_h->fill(event);
               events_heptoptagevent_WP3_wobtag_nomtlep_h->fill(event);
               topjets_heptoptagevent_WP3_wobtag_nomtlep_h->fill(event); 
               likelihood_heptoptagevent_WP3_wobtag_nomtlep_h->fill(event);
            }


        
         if (!pass_HEPTT_WP3) {
            topjets_isolep_notag_WP3_h->fill(event);
            return false;
         }
         if (!pass_lepVeto) return false;
         heptoptagevent_WP3_h->fill(event);
         electrons_heptoptagevent_WP3_h->fill(event);
         jets_heptoptagevent_WP3_h->fill(event);
         muons_heptoptagevent_WP3_h->fill(event);
         events_heptoptagevent_WP3_h->fill(event);
         topjets_heptoptagevent_WP3_h->fill(event); 
         // topjets_heptoptagevent_WP3_tagged_h->fill(event); 
         likelihood_heptoptagevent_WP3_h->fill(event); 
         ttbargenhists_heptoptagevent_WP3_h->fill(event);                  
         
         //     remove cuts to check selection
         // test with met100, otheriwse very low stat
         //no btag
          if (pass_thirdjet_sel && pass_jet123metdphi_sel && pass_DeltaPhiMetNeutrino_sel && pass_met100 && pass_mtlep100){
             nobtag_h->fill(event);
             likelihood_nobtag_h->fill(event);
          }
         // no third jet, no DeltaPhiJets
          if (pass_btagsel && pass_DeltaPhiMetNeutrino_sel && pass_mtlep100){
             nothirdjet_h->fill(event);
             likelihood_nothirdjet_h->fill(event);
          }
          //no DeltaPhiJets
          if (pass_thirdjet_sel && pass_btagsel && pass_DeltaPhiMetNeutrino_sel && pass_met100 && pass_mtlep100){
             noDeltaPhiJets_h->fill(event);
             likelihood_noDeltaPhiJets_h->fill(event);
          }
        
          //no met cut
          if (pass_thirdjet_sel && pass_jet123metdphi_sel && pass_DeltaPhiMetNeutrino_sel &&  pass_mtlep100 && pass_btagsel){
             nomet_h->fill(event);
             likelihood_nomet_h->fill(event);
          }
          if (pass_thirdjet_sel &&  pass_DeltaPhiMetNeutrino_sel &&  pass_mtlep100 && pass_btagsel){
             nomet_noDeltaPhi_h->fill(event);
             likelihood_nomet_noDeltaPhi_h->fill(event);
          }
          //no mtcut
          if (pass_thirdjet_sel && pass_jet123metdphi_sel && pass_DeltaPhiMetNeutrino_sel && pass_met100 && pass_btagsel){
             nomtlep->fill(event);
             likelihood_nomtlep->fill(event);
          }
          // no DeltaPhiMetNeutrino
          if (pass_thirdjet_sel && pass_jet123metdphi_sel && pass_btagsel && pass_met100 && pass_mtlep100){
             noDeltaPhiMetNeutrino->fill(event);
             likelihood_noDeltaPhiMetNeutrino->fill(event);
             if (pass_met160){
                noDeltaPhiMetNeutrino_met160->fill(event);
                likelihood_noDeltaPhiMetNeutrino_met160->fill(event);
             }
             if (pass_DeltaPhiMetNeutrino_sel2){
                DeltaPhiMetNeutrinotight->fill(event);
                likelihood_DeltaPhiMetNeutrinotight->fill(event);
             }
            
          }
         
          //bool pass_btagsel = btag_sel->passes(event);
         if (!pass_btagsel) return false;
         btagsel_h->fill(event);
         likelihood_btagsel->fill(event);
         electrons_btagsel_h->fill(event);
         jets_btagsel_h->fill(event);
         muons_btagsel_h->fill(event);
         events_btagsel_h->fill(event);
         topjets_btagsel_h->fill(event); 
         ttbargenhists_btagsel_h->fill(event); 
         
         pass_thirdjet_sel = thirdjet_sel->passes(event);
         if (!pass_thirdjet_sel)return false;
         
         pass_jet123metdphi_sel=jet123metdphi_sel->passes(event);
         if (!pass_jet123metdphi_sel) return false;
         jet123metdphi_h->fill(event);
         hists_likelihood_jet123metdphi_h->fill(event);
         
         //bool pass_DeltaPhiMetNeutrino_sel=DeltaPhiMetNeutrino_sel->passes(event);
         // if (!pass_DeltaPhiMetNeutrino_sel) return false;
         // met50_btag_h->fill(event);
         // likelihood_met50_btag->fill(event);
         // muons_met50_btag_h->fill(event);
         // jets_met50_btag_h->fill(event);
         // events_met50_btag_h->fill(event);
         // topjets_met50_btag_h->fill(event);
         
         // bool pass_met100 = MET100_sel->passes(event);
         // if (!pass_met100) return false; 
         // met100_btag_h->fill(event);
         // likelihood_met100_btag->fill(event);
         // muons_met100_btag_h->fill(event);
         // jets_met100_btag_h->fill(event);
         // events_met100_btag_h->fill(event);
         // topjets_met100_btag_h->fill(event);
        
         //bool pass_mtlep100 = mtlep100_sel->passes(event);
         if (!pass_mtlep100) return false;
         met100_btag_mtlep100_h->fill(event);
         likelihood_met100_btag_mtlep100->fill(event);
      
         //bool pass_met160= MET160_sel->passes(event);
         // if (pass_met160)
         //    {
         //       met160_h->fill(event);  //FSP bestes fuer DMMET
         //       likelihood_met160->fill(event);
         //    } 
         
         return false;
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
   twodcut_h->fill(event);
   electrons_twodcut_h->fill(event);
   muons_twodcut_h->fill(event);
   events_twodcut_h->fill(event);
   jets_twodcut_h->fill(event);
   topjets_twodcut_h->fill(event);
   likelihood_twodcut_h->fill(event);
 
   //bool pass_jetmetdphijet2_sel = jetmetdphijet2_sel->passes(event);
  bool pass_met100 = MET100_sel->passes(event);
   if (pass_met100){
      mediumMET_h->fill(event);
      electrons_mediumMET_h->fill(event);
      muons_mediumMET_h->fill(event);
      events_mediumMET_h->fill(event);
      jets_mediumMET_h->fill(event);
      topjets_mediumMET_h->fill(event);
      likelihood_mediumMET_h->fill(event);  
   }
   
   bool pass_btag_sel = btag_sel->passes(event);
   bool pass_thirdjet_sel = thirdjet_sel->passes(event);
   bool pass_DeltaPhiMetNeutrino_sel=DeltaPhiMetNeutrino_sel2->passes(event);
   bool pass_mtlep100 = mtlep100_sel->passes(event);
   
   //control regions
   //btag veto
   if (!pass_btag_sel && pass_lepVeto &&pass_mtlep100){
      WJetsControlRegion_noniso_h->fill(event);
      electrons_WJetsControlRegion_noniso_h->fill(event);
      muons_WJetsControlRegion_noniso_h->fill(event);
      events_WJetsControlRegion_noniso_h->fill(event);
      jets_WJetsControlRegion_noniso_h->fill(event);
      topjets_WJetsControlRegion_noniso_h->fill(event);
      likelihood_WJetsControlRegion_noniso_h->fill(event);
   }
   //lepton veto 
   if (pass_btag_sel && !pass_lepVeto && pass_mtlep100){
      DileptControlRegion_noniso_h->fill(event);
      electrons_DileptControlRegion_noniso_h->fill(event);
      muons_DileptControlRegion_noniso_h->fill(event);
      events_DileptControlRegion_noniso_h->fill(event);
      jets_DileptControlRegion_noniso_h->fill(event);
      topjets_DileptControlRegion_noniso_h->fill(event);
      likelihood_DileptControlRegion_noniso_h->fill(event);
   }

   //mTsel
   if (pass_btag_sel && pass_lepVeto && !pass_mtlep100){
      SemileptControlRegion_noniso_h->fill(event);
      electrons_SemileptControlRegion_noniso_h->fill(event);
      muons_SemileptControlRegion_noniso_h->fill(event);
      events_SemileptControlRegion_noniso_h->fill(event);
      jets_SemileptControlRegion_noniso_h->fill(event);
      topjets_SemileptControlRegion_noniso_h->fill(event);
      likelihood_SemileptControlRegion_noniso_h->fill(event);
   }

   //signal selections
   if (!pass_lepVeto) return false;
   if (!pass_btag_sel) return false;
   btag_noniso_h->fill(event);
   electrons_btag_noniso_h->fill(event);
   muons_btag_noniso_h->fill(event);
   events_btag_noniso_h->fill(event);
   jets_btag_noniso_h->fill(event);
   topjets_btag_noniso_h->fill(event);
   likelihood_btag_noniso_h->fill(event);
 
 
   
   if (pass_DeltaPhiMetNeutrino_sel)
      {
         met50_btag_h->fill(event);
         likelihood_met50_btag->fill(event);
         muons_met50_btag_h->fill(event);
         jets_met50_btag_h->fill(event);
         events_met50_btag_h->fill(event);
         topjets_met50_btag_h->fill(event);
      }
   
   if (pass_met100){
      met100_btag_h->fill(event);
      likelihood_met100_btag->fill(event);
      muons_met100_btag_h->fill(event);
      jets_met100_btag_h->fill(event);
      events_met100_btag_h->fill(event);
      topjets_met100_btag_h->fill(event);
   }
   if (pass_mtlep100){
      mtlep_h->fill(event);
      electrons_mtlep_h->fill(event);
      muons_mtlep_h->fill(event);
      events_mtlep_h->fill(event);
      jets_mtlep_h->fill(event);
      topjets_mtlep_h->fill(event);
      likelihood_mtlep_h->fill(event);        
   }

   if (pass_HEPTT_WP3_wobtag){
      HEPTT_WP3_noniso_h->fill(event);
      likelihood_WP3_noniso_h->fill(event);
   }

   if (pass_HEPTT_WP3_wobtag &&  pass_mtlep100){
      HEPTT_WP3_noniso_mtlep_h->fill(event);
      likelihood_WP3_noniso_mtlep_h->fill(event);

   }
 
   bool pass_mtlep_sel = mtlep_sel->passes(event); 
        
   if (pass_thirdjet_sel){
      bool pass_jet123metdphi_sel=jet123metdphi_sel->passes(event);
      if (pass_jet123metdphi_sel){
         thirdjet_h->fill(event);
         electrons_thirdjet_h->fill(event);
         muons_thirdjet_h->fill(event);
         events_thirdjet_h->fill(event);
         jets_thirdjet_h->fill(event);
         topjets_thirdjet_h->fill(event);
         likelihood_thirdjet_h->fill(event);    
         if (pass_mtlep_sel){
            likelihood_thirdjet_mtlep_h->fill(event);
            thirdjet_mtlep_h->fill(event);
         }
      }
   }
 
   
   return true;
}



UHH2_REGISTER_ANALYSIS_MODULE(ttDMSelectionModuleAfterLikelihood)

