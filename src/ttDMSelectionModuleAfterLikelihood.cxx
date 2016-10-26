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
   std::unique_ptr<Selection> highMET_sel;
   std::unique_ptr<Selection> MET_sel;
   std::unique_ptr<Selection> MET100_sel;
   std::unique_ptr<Selection> MET80_sel;
   std::unique_ptr<Selection> MET160_sel,btag_sel;
   std::unique_ptr<AndSelection> jetsel1208050;
   std::unique_ptr<AndSelection> lep1_sel, lepVeto_sel;

   std::unique_ptr<Hists> highMET_h;
   std::unique_ptr<Hists> mtlep_noniso_h;
   std::unique_ptr<Hists> highMET_noniso_h;
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
   
   // hists
   std::unique_ptr<Hists> lumihists;
   std::unique_ptr<Hists> filter_h;
   std::unique_ptr<Hists> input_h;
   std::unique_ptr<Hists> met_h, electrons_met_h,jets_met_h,muons_met_h,events_met_h,topjets_met_h, genhists_met_h;
   std::unique_ptr<Hists> met100_h, electrons_met100_h,jets_met100_h,muons_met100_h,events_met100_h,topjets_met100_h, genhists_met100_h;
   std::unique_ptr<Hists> nonisolep_h, electrons_nonisolep_h,jets_nonisolep_h,muons_nonisolep_h,events_nonisolep_h,topjets_nonisolep_h;
   std::unique_ptr<Hists> isolep_h, electrons_isolep_h,jets_isolep_h,muons_isolep_h,events_isolep_h,topjets_isolep_h;
   std::unique_ptr<TopJetHists> topjets_isolep_tagged_WP1_h, topjets_isolep_tagged_WP3_h, topjets_isolep_notag_WP1_h, topjets_isolep_notag_WP3_h;
   std::unique_ptr<TopJetHists> topjets_isolep_tagged_WP1_wobtag_h, topjets_isolep_tagged_WP3_wobtag_h;
   std::unique_ptr<Hists> heptoptagevent_WP3_h,electrons_heptoptagevent_WP3_h,jets_heptoptagevent_WP3_h,muons_heptoptagevent_WP3_h,events_heptoptagevent_WP3_h,topjets_heptoptagevent_WP3_h;
   std::unique_ptr<Hists> heptoptagevent_WP3_wobtag_h,electrons_heptoptagevent_WP3_wobtag_h,jets_heptoptagevent_WP3_wobtag_h,muons_heptoptagevent_WP3_wobtag_h,events_heptoptagevent_WP3_wobtag_h,topjets_heptoptagevent_WP3_wobtag_h;
   std::unique_ptr<TopJetHists>  topjets_heptoptagevent_WP3_tagged_h;
   std::unique_ptr<TopJetHists>  topjets_heptoptagevent_WP3_wobtag_tagged_h;
   std::unique_ptr<Hists> twodcut_h,electrons_twodcut_h, jets_twodcut_h,topjets_twodcut_h,events_twodcut_h,muons_twodcut_h;
   std::unique_ptr<Hists> metlepdphi_noiso_h,electrons_metlepdphi_noiso_h, jets_metlepdphi_noiso_h,topjets_metlepdphi_noiso_h,events_metlepdphi_noiso_h,muons_metlepdphi_noiso_h;
   std::unique_ptr<Hists> preselection_h,electrons_preselection_h, jets_preselection_h,topjets_preselection_h,events_preselection_h,muons_preselection_h, genhists_preselection_h;
   std::unique_ptr<Hists> leptonveto_h,electrons_leptonveto_h, jets_leptonveto_h,topjets_leptonveto_h,events_leptonveto_h,muons_leptonveto_h, genhists_leptonveto_h;
   std::unique_ptr<Hists> heptoptagevent_WP1_h;
   std::unique_ptr<Hists> heptoptagevent_WP1_wobtag_h;
   std::unique_ptr<Hists> jetmetdphi_h;
   std::unique_ptr<Hists> jetmetdphi_noniso_h;
   std::unique_ptr<Hists> likelihood_preselection_h;
   std::unique_ptr<Hists> likelihood_leptonveto_h;
   std::unique_ptr<Hists> likelihood_met_h;
   std::unique_ptr<Hists> likelihood_twodcut_h;
   std::unique_ptr<Hists> likelihood_metlepdphi_h;
   std::unique_ptr<Hists> likelihood_jetmetdphi_h;
   std::unique_ptr<Hists> likelihood_cut1_noiso_h;
   std::unique_ptr<Hists> likelihood_cut2_noiso_h;
   std::unique_ptr<Hists> metlepdphi_h,jetmetdphi_mediumMET_h;
   std::unique_ptr<Hists> likelihood_metlepdphi_noiso_h;
   std::unique_ptr<Selection> jetmetdphi_sel, jetmetdphijet2_sel;
   std::unique_ptr<Selection> mtlep_sel,mediumMET_sel;
   std::unique_ptr<Selection> cut1_likelihood;
   std::unique_ptr<Selection> cut2_likelihood;
   std::unique_ptr<Selection> highpTJetSel;
   std::unique_ptr<Selection> leadingjet_sel;
   std::unique_ptr<Selection> HTlep_sel;
   std::unique_ptr<Selection> thirdjet_sel, mtlep100_sel, mtlep150_sel;
   std::unique_ptr<Selection> metleptondphi_sel;
   std::unique_ptr<Selection> jet123metdphi_sel;
   std::unique_ptr<Selection> MT2WCut_sel;
   std::unique_ptr<Selection> DeltaPhiMetNeutrino_sel;
   std::unique_ptr<Selection> DeltaPhiTaggedjet_Neutrino_sel;
   std::unique_ptr<Selection> DeltaPhiTaggedJetTopLep_sel;
   std::unique_ptr<Selection> neutrinopT_sel;
   std::unique_ptr<Selection> ttbarpT_sel;
   
   std::unique_ptr<Hists> genhists_nonisolep_h;
   std::unique_ptr<Hists> likelihood_metlepdphi_mtlep_h, likelihood_heptoptagevent_WP3_h,likelihood_heptoptagevent_WP1_h, likelihood_met100_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_wobtag_h,likelihood_heptoptagevent_WP1_wobtag_h;
   std::unique_ptr<Hists> hists_metlepdphi_mtlep;
   std::unique_ptr<Hists> cut2_noiso_mtlep;
   std::unique_ptr<Hists> likelihood_cut2_noiso_mtlep_h;
  
   std::unique_ptr<Hists> met50_btag_ttbarpTSel;
   std::unique_ptr<Hists> likelihood_met50_btag_ttbarpTSel;

   std::unique_ptr<Hists> mediumMET_h,electrons_mediumMET_h,jets_mediumMET_h,topjets_mediumMET_h,muons_mediumMET_h,events_mediumMET_h, likelihood_mediumMET_h;
   std::unique_ptr<Hists> thirdjet_h,electrons_thirdjet_h,jets_thirdjet_h,topjets_thirdjet_h,muons_thirdjet_h,events_thirdjet_h,likelihood_thirdjet_h;
   std::unique_ptr<Hists> mtlep_h,electrons_mtlep_h,jets_mtlep_h,topjets_mtlep_h,muons_mtlep_h,events_mtlep_h,likelihood_mtlep_h;
   std::unique_ptr<Hists> jetmetdphijet2_h,electrons_jetmetdphijet2_h,jets_jetmetdphijet2_h,topjets_jetmetdphijet2_h,muons_jetmetdphijet2_h,events_jetmetdphijet2_h,likelihood_jetmetdphijet2_h;
   std::unique_ptr<Hists> electrons_btagsel_h,jets_btagsel_h,topjets_btagsel_h,muons_btagsel_h,events_btagsel_h,ttbargenhists_btagsel_h;

   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_wobtag_highMET_h;
   std::unique_ptr<Hists> heptoptagevent_WP3_wobtag_highMET_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP1_wobtag_highMET_h;
   std::unique_ptr<Hists> heptoptagevent_WP1_wobtag_highMET_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP1_highMET_h;
   std::unique_ptr<Hists> heptoptagevent_WP1_highMET_h;

   std::unique_ptr<Hists> likelihood_thirdjet_mtlep_h;
   std::unique_ptr<Hists> thirdjet_mtlep_h;
   std::unique_ptr<Hists> likelihood_thirdjet_jetmetdphi_h;
   std::unique_ptr<Hists> thirdjet_jetmetdphi_h;

   std::unique_ptr<Hists> electrons_thirdjet_jetmetdphi_h;
   std::unique_ptr<Hists> jets_thirdjet_jetmetdphi_h;
   std::unique_ptr<Hists> topjets_thirdjet_jetmetdphi_h;
   std::unique_ptr<Hists> muons_thirdjet_jetmetdphi_h;
   std::unique_ptr<Hists> events_thirdjet_jetmetdphi_h;
   std::unique_ptr<Hists> likelihood_jetmetdphi_mediumMET_h;
   std::unique_ptr<Hists> jet123metdphi_h;
   std::unique_ptr<Hists> hists_likelihood_jet123metdphi_h;
   
   std::unique_ptr<Hists> jet123metdphi_MT2WCut_h;
   std::unique_ptr<Hists> hists_likelihood_jet123metdphi_MT2WCut_h;
   std::unique_ptr<Hists> metlepdphi_MT2WCut_h;
   std::unique_ptr<Hists> likelihood_metlepdphi_MT2WCut_h;
   std::unique_ptr<Hists> ttbargenhists_preselection_h;
   std::unique_ptr<Hists> ttbargenhists_leptonveto_h;
   std::unique_ptr<Hists> ttbargenhists_heptoptagevent_WP3_h;
   std::unique_ptr<Hists> ttbargenhists_jetmetdphi_mediumMET_h;
   std::unique_ptr<Hists> ttbargenhists_metlepdphi_mtlep;
   std::unique_ptr<Hists> ttbargenhists_metlepdphi_MT2WCut_h;
   std::unique_ptr<Hists> jet123metdphi_MT2WCut_MET240_h;
   std::unique_ptr<Hists> ttbargenhists_met_h;
   std::unique_ptr<Hists> likelihood_jet123metdphi_MT2WCut_MET240_h;
   std::unique_ptr<Hists> thirdjet_MT2WCut_MetSel_h;
   std::unique_ptr<Hists> likelihood_thirdjet_MT2WCut_MetSel_h;
   std::unique_ptr<Hists> thirdjet_MT2WCut_h;
   std::unique_ptr<Hists> likelihood_thirdjet_MT2WCut_h;
   uhh2::Event::Handle<TTbarGen> h_ttbargen;
   
   std::unique_ptr<Hists> electrons_jet123metdphi_MT2WCut_h;
   std::unique_ptr<Hists> jets_jet123metdphi_MT2WCut_h;
   std::unique_ptr<Hists> topjets_jet123metdphi_MT2WCut_h;
   std::unique_ptr<Hists> muons_jet123metdphi_MT2WCut_h;
   std::unique_ptr<Hists> events_jet123metdphi_MT2WCut_h;
   std::unique_ptr<Hists> jetsel1208050_h;
   std::unique_ptr<Hists> likelihood_jetsel1208050_h;
   std::unique_ptr<Hists> met80_h;
   std::unique_ptr<Hists> likelihood_met80;
   std::unique_ptr<Hists> btagsel_h;
   std::unique_ptr<Hists> likelihood_btagsel;
   std::unique_ptr<Hists> met160_h;
   std::unique_ptr<Hists> likelihood_met160;
   std::unique_ptr<Hists> met160_btag_h;
   std::unique_ptr<Hists> likelihood_met160_btag;
   std::unique_ptr<Hists> hists_likelihood_DeltaPhiMetNeutrino;
   std::unique_ptr<Hists> hists_DeltaPhiMetNeutrino;
   std::unique_ptr<Hists> heptoptagevent_WP3_highMET_h;
   std::unique_ptr<Hists> likelihood_heptoptagevent_WP3_highMET_h;
   std::unique_ptr<Hists> DeltaPhiTaggedjet_Neutrino_mtlep150;
   std::unique_ptr<Hists> likelihood_DeltaPhiTaggedjet_Neutrino_mtlep150;
   std::unique_ptr<Hists> electrons_btagsel_MT2WCut_h,jets_btagsel_MT2WCut_h,topjets_btagsel_MT2WCut_h,muons_btagsel_MT2WCut_h,events_btagsel_MT2WCut_h,ttbargenhists_btagsel_MT2WCut_h,btagsel_MT2WCut_h,likelihood_btagsel_MT2WCut;
   std::unique_ptr<Hists> electrons_btagsel_mtlep_h,jets_btagsel_mtlep_h,topjets_btagsel_mtlep_h,muons_btagsel_mtlep_h,events_btagsel_mtlep_h,ttbargenhists_btagsel_mtlep_h,btagsel_mtlep_h,likelihood_btagsel_mtlep;
   std::unique_ptr<Hists> electrons_btagsel_MT2WCut_metsel_h,jets_btagsel_MT2WCut_metsel_h,topjets_btagsel_MT2WCut_metsel_h,muons_btagsel_MT2WCut_metsel_h,events_btagsel_MT2WCut_metsel_h,ttbargenhists_btagsel_MT2WCut_metsel_h,btagsel_MT2WCut_metsel_h,likelihood_btagsel_MT2WCut_metsel;
   std::unique_ptr<Hists> jets_met160_btag_h,topjets_met160_btag_h,muons_met160_btag_h,events_met160_btag_h;
  
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
   
   std::unique_ptr<Hists> met160_btag_metlep100_h;
   std::unique_ptr<Hists> likelihood_met160_btag_metlep100;
   std::unique_ptr<Hists> met160_btag_metlepdphi_h;
   std::unique_ptr<Hists> likelihood_met160_btag_metlepdphi;
   std::unique_ptr<Hists> DeltaPhiTaggedjet_Neutrino;
   std::unique_ptr<Hists> likelihood_DeltaPhiTaggedjet_Neutrino;
   std::unique_ptr<Hists> met100_btag_mtlep100_h;
   std::unique_ptr<Hists> likelihood_met100_btag_mtlep100;
   std::unique_ptr<Hists> met100_btag_neutrinopT_h;
   std::unique_ptr<Hists> likelihood_met100_btag_neutrinopT;
   std::unique_ptr<Hists> met100_btag_DeltaPhiTaggedJetTopLep_h;
   std::unique_ptr<Hists> likelihood_met100_btag_DeltaPhiTaggedJetTopLep;
   std::unique_ptr<Hists> HEPTT_WP3_noniso_h;
   std::unique_ptr<Hists> likelihood_WP3_noniso_h;
   std::unique_ptr<Hists> likelihood_btag_noniso_h;
   std::unique_ptr<Hists> met100_btag_neutrinopT_ttbarpTSel_h;
   std::unique_ptr<Hists> likelihood_met100_btag_neutrinopT_ttbarpTSel;

   bool isSemiLept;
   bool isDiLept;
   bool isOther;
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
   
   
   highpTJetSel.reset(new NJetSelection(2, -1, JetId(PtEtaCut( 100., 2.4))));
   HTlep_sel.reset(new HTlepCut(300, std::numeric_limits<double>::infinity()));
   
   jetsel1208050.reset(new AndSelection(ctx));
   jetsel1208050->add<NJetSelection>("",1, -1, JetId(PtEtaCut( 120., 2.4)));
   jetsel1208050->add<NJetSelection>("",2, -1, JetId(PtEtaCut( 80., 2.4)));
   jetsel1208050->add<NJetSelection>("",3, -1, JetId(PtEtaCut( 50., 2.4)));

   JetId Btag = CSVBTag(0.8);
   btag_sel.reset(new NJetSelection(1, -1, Btag));

   leadingjet_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut( 150., 2.4))));
   highMET_sel.reset(new METCut(150., std::numeric_limits<double>::infinity()));
   mediumMET_sel.reset(new METCut(260., std::numeric_limits<double>::infinity()));
   MET_sel.reset(new METCut(320., std::numeric_limits<double>::infinity()));
   MET100_sel.reset(new METCut(100., std::numeric_limits<double>::infinity()));
   MET160_sel.reset(new METCut(160., std::numeric_limits<double>::infinity()));
   MET80_sel.reset(new METCut(80., std::numeric_limits<double>::infinity()));
   mtlep_sel.reset(new MTlepCut(200., std::numeric_limits<double>::infinity()));
   mtlep100_sel.reset(new MTlepCut(100., std::numeric_limits<double>::infinity()));
   mtlep150_sel.reset(new MTlepCut(150., std::numeric_limits<double>::infinity()));
   jetmetdphi_sel.reset(new METJetDPhiCut(1.2, 2));
   jetmetdphijet2_sel.reset(new METJetDPhiCut(1.2, 2));
   metleptondphi_sel.reset(new METLeptonDPhiCut(1.1));
   thirdjet_sel.reset(new NJetSelection(3, -1, JetId(PtEtaCut( 50., 2.4))));
   ttbarpT_sel.reset(new ttbarpTSel(ctx,350));
   jet123metdphi_sel.reset(new METJetDPhiCut(1.4,3));
   
   cut1_likelihood.reset(new LikelihoodSelection(ctx,22));
   cut2_likelihood.reset(new LikelihoodSelection(ctx,14));
    
   DeltaPhiMetNeutrino_sel.reset(new DeltaPhiMetNeutrino(ctx,1));
   DeltaPhiTaggedjet_Neutrino_sel.reset(new  DeltaPhiTaggedJetNeutrino(ctx,1.5));
   neutrinopT_sel.reset(new NeutrinopTSelection(ctx, 150));
   DeltaPhiTaggedJetTopLep_sel.reset(new DeltaPhiTaggedJetTopLep(ctx,1.5));
   //h_flag_heptoptagevent = ctx.declare_event_output<int>("flag_heptoptagevent");
   h_flag_twodcut = ctx.get_handle<bool>("flag_twodcut");
   h_likelihood = ctx.get_handle<double>("likelihood");
   h_recneutrino = ctx.get_handle<LorentzVector>("rec_neutrino");
   h_bjets = ctx.get_handle<Jet>("bjet");

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
   collectionprod_heptt_pteta.reset(new CollectionProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "h_heptopjets_pteta", TopJetId(HTTTopJetId_pteta))); 
   collectionsizeprod_heptt_pteta.reset(new CollectionSizeProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "n_heptopjets_pteta", TopJetId(HTTTopJetId_pteta)));  

   collectionprod_heptt_WP3_wobtag.reset(new CollectionProducer<TopJet>(ctx, "h_heptopjets_pteta", "h_heptopjets_WP3_wobtag", TopJetId(HTTTopJetId_WP3_wobtag))); 
   collectionsizeprod_heptt_WP3_wobtag.reset(new CollectionSizeProducer<TopJet>(ctx, "h_heptopjets_pteta", "n_heptopjets_WP3_wobtag", TopJetId(HTTTopJetId_WP3_wobtag)));  
   
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
   met100_h.reset(new ttDMSelectionHists(ctx, "met100"));

   electrons_met_h.reset(new ElectronHists(ctx,"electrons_met"));
   jets_met_h.reset(new JetHists(ctx,"jets_met"));
   muons_met_h.reset(new MuonHists(ctx,"muons_met"));
   events_met_h.reset(new EventHists(ctx,"events_met"));
   topjets_met_h.reset(new TopJetHists(ctx,"topjets_met",4,"patJetsHepTopTagCHSPacked_daughters"));
   genhists_met_h.reset(new ttDMGenHists(ctx,"genhists_met"));

   electrons_met100_h.reset(new ElectronHists(ctx,"electrons_met100"));
   jets_met100_h.reset(new JetHists(ctx,"jets_met100"));
   muons_met100_h.reset(new MuonHists(ctx,"muons_met100"));
   events_met100_h.reset(new EventHists(ctx,"events_met100"));
   topjets_met100_h.reset(new TopJetHists(ctx,"topjets_met100",4,"patJetsHepTopTagCHSPacked_daughters"));
   genhists_met100_h.reset(new ttDMGenHists(ctx,"genhists_met100"));

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
   topjets_heptoptagevent_WP3_tagged_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_tagged", 4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_heptoptagevent_WP3_tagged_h->set_TopJetId(HTTTopJetId_WP3);

   heptoptagevent_WP3_wobtag_h.reset(new ttDMSelectionHists(ctx,"heptoptagevent_WP3_wobtag"));
   electrons_heptoptagevent_WP3_wobtag_h.reset(new ElectronHists(ctx,"electrons_heptoptagevent_WP3_wobtag"));
   jets_heptoptagevent_WP3_wobtag_h.reset(new JetHists(ctx,"jets_heptoptagevent_WP3_wobtag"));
   muons_heptoptagevent_WP3_wobtag_h.reset(new MuonHists(ctx,"muons_heptoptagevent_WP3_wobtag"));
   events_heptoptagevent_WP3_wobtag_h.reset(new EventHists(ctx,"events_heptoptagevent_WP3_wobtag"));
   topjets_heptoptagevent_WP3_wobtag_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_wobtag", 4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_heptoptagevent_WP3_wobtag_tagged_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_wobtag_tagged", 4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_heptoptagevent_WP3_wobtag_tagged_h->set_TopJetId(HTTTopJetId_WP3_wobtag);
   
   twodcut_h.reset(new ttDMSelectionHists(ctx,"twodcut"));
   electrons_twodcut_h.reset(new ElectronHists(ctx,"electrons_twodcut"));
   jets_twodcut_h.reset(new JetHists(ctx,"jets_twodcut"));
   topjets_twodcut_h.reset(new TopJetHists(ctx,"topjets_twodcut", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_twodcut_h.reset(new MuonHists(ctx,"muons_twodcut"));
   events_twodcut_h.reset(new EventHists(ctx,"events_twodcut"));

   metlepdphi_noiso_h.reset(new ttDMSelectionHists(ctx,"metlepdphi_noiso"));
   electrons_metlepdphi_noiso_h.reset(new ElectronHists(ctx,"electrons_metlepdphi_noiso"));
   jets_metlepdphi_noiso_h.reset(new JetHists(ctx,"jets_metlepdphi_noiso"));
   topjets_metlepdphi_noiso_h.reset(new TopJetHists(ctx,"topjets_metlepdphi_noiso", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_metlepdphi_noiso_h.reset(new MuonHists(ctx,"muons_metlepdphi_noiso"));
   events_metlepdphi_noiso_h.reset(new EventHists(ctx,"events_metlepdphi_noiso"));

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

   mediumMET_h.reset(new ttDMSelectionHists(ctx,"mediumMET"));
   electrons_mediumMET_h.reset(new ElectronHists(ctx,"electrons_mediumMET"));
   jets_mediumMET_h.reset(new JetHists(ctx,"jets_mediumMET"));
   topjets_mediumMET_h.reset(new TopJetHists(ctx,"topjets_mediumMET", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_mediumMET_h.reset(new MuonHists(ctx,"muons_mediumMET"));
   events_mediumMET_h.reset(new EventHists(ctx,"events_mediumMET"));

   thirdjet_h.reset(new ttDMSelectionHists(ctx,"thirdjet"));
   electrons_thirdjet_h.reset(new ElectronHists(ctx,"electrons_thirdjet"));
   jets_thirdjet_h.reset(new JetHists(ctx,"jets_thirdjet"));
   topjets_thirdjet_h.reset(new TopJetHists(ctx,"topjets_thirdjet", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_thirdjet_h.reset(new MuonHists(ctx,"muons_thirdjet"));
   events_thirdjet_h.reset(new EventHists(ctx,"events_thirdjet"));

   electrons_thirdjet_jetmetdphi_h.reset(new ElectronHists(ctx,"electrons_thirdjet_jetmetdphi"));
   jets_thirdjet_jetmetdphi_h.reset(new JetHists(ctx,"jets_thirdjet_jetmetdphi"));
   topjets_thirdjet_jetmetdphi_h.reset(new TopJetHists(ctx,"topjets_thirdjet_jetmetdphi", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_thirdjet_jetmetdphi_h.reset(new MuonHists(ctx,"muons_thirdjet_jetmetdphi"));
   events_thirdjet_jetmetdphi_h.reset(new EventHists(ctx,"events_thirdjet_jetmetdphi"));

   electrons_btagsel_h.reset(new ElectronHists(ctx,"electrons_btagsel"));
   jets_btagsel_h.reset(new JetHists(ctx,"jets_btagsel"));
   topjets_btagsel_h.reset(new TopJetHists(ctx,"topjets_btagsel", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_btagsel_h.reset(new MuonHists(ctx,"muons_btagsel"));
   events_btagsel_h.reset(new EventHists(ctx,"events_btagsel"));

   electrons_btagsel_MT2WCut_h.reset(new ElectronHists(ctx,"electrons_btagsel_MT2WCut"));
   jets_btagsel_MT2WCut_h.reset(new JetHists(ctx,"jets_btagsel_MT2WCut"));
   topjets_btagsel_MT2WCut_h.reset(new TopJetHists(ctx,"topjets_btagsel_MT2WCut", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_btagsel_MT2WCut_h.reset(new MuonHists(ctx,"muons_btagsel_MT2WCut"));
   events_btagsel_MT2WCut_h.reset(new EventHists(ctx,"events_btagsel_MT2WCut"));
   ttbargenhists_btagsel_MT2WCut_h.reset(new TTbarGenHists(ctx,"ttbargenhists_btagsel_MT2WCut"));
   btagsel_MT2WCut_h.reset(new ttDMSelectionHists(ctx, "btagsel_MT2WCut"));
   likelihood_btagsel_MT2WCut.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btagsel_MT2WCut"));

   electrons_btagsel_mtlep_h.reset(new ElectronHists(ctx,"electrons_btagsel_mtlep"));
   jets_btagsel_mtlep_h.reset(new JetHists(ctx,"jets_btagsel_mtlep"));
   topjets_btagsel_mtlep_h.reset(new TopJetHists(ctx,"topjets_btagsel_mtlep", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_btagsel_mtlep_h.reset(new MuonHists(ctx,"muons_btagsel_mtlep"));
   events_btagsel_mtlep_h.reset(new EventHists(ctx,"events_btagsel_mtlep"));
   ttbargenhists_btagsel_mtlep_h.reset(new TTbarGenHists(ctx,"ttbargenhists_btagsel_mtlep"));
   btagsel_mtlep_h.reset(new ttDMSelectionHists(ctx, "btagsel_mtlep"));
   likelihood_btagsel_mtlep.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btagsel_mtlep"));

   electrons_btagsel_MT2WCut_metsel_h.reset(new ElectronHists(ctx,"electrons_btagsel_MT2WCut_metsel"));
   jets_btagsel_MT2WCut_metsel_h.reset(new JetHists(ctx,"jets_btagsel_MT2WCut_metsel"));
   topjets_btagsel_MT2WCut_metsel_h.reset(new TopJetHists(ctx,"topjets_btagsel_MT2WCut_metsel", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_btagsel_MT2WCut_metsel_h.reset(new MuonHists(ctx,"muons_btagsel_MT2WCut_metsel"));
   events_btagsel_MT2WCut_metsel_h.reset(new EventHists(ctx,"events_btagsel_MT2WCut_metsel"));
   ttbargenhists_btagsel_MT2WCut_metsel_h.reset(new TTbarGenHists(ctx,"ttbargenhists_btagsel_MT2WCut_metsel"));
   btagsel_MT2WCut_metsel_h.reset(new ttDMSelectionHists(ctx, "btagsel_MT2WCut_metsel"));
   likelihood_btagsel_MT2WCut_metsel.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btagsel_MT2WCut_metsel"));

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

   mtlep_h.reset(new ttDMSelectionHists(ctx,"mtlep"));
   electrons_mtlep_h.reset(new ElectronHists(ctx,"electrons_mtlep"));
   jets_mtlep_h.reset(new JetHists(ctx,"jets_mtlep"));
   topjets_mtlep_h.reset(new TopJetHists(ctx,"topjets_mtlep", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_mtlep_h.reset(new MuonHists(ctx,"muons_mtlep"));
   events_mtlep_h.reset(new EventHists(ctx,"events_mtlep"));

   jetmetdphijet2_h.reset(new ttDMSelectionHists(ctx,"jetmetdphijet2"));
   electrons_jetmetdphijet2_h.reset(new ElectronHists(ctx,"electrons_jetmetdphijet2"));
   jets_jetmetdphijet2_h.reset(new JetHists(ctx,"jets_jetmetdphijet2"));
   topjets_jetmetdphijet2_h.reset(new TopJetHists(ctx,"topjets_jetmetdphijet2", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_jetmetdphijet2_h.reset(new MuonHists(ctx,"muons_jetmetdphijet2"));
   events_jetmetdphijet2_h.reset(new EventHists(ctx,"events_jetmetdphijet2"));

   btag_noniso_h.reset(new ttDMSelectionHists(ctx,"btag_noniso"));
   electrons_btag_noniso_h.reset(new ElectronHists(ctx,"electrons_btag_noniso"));
   jets_btag_noniso_h.reset(new JetHists(ctx,"jets_btag_noniso"));
   topjets_btag_noniso_h.reset(new TopJetHists(ctx,"topjets_btag_noniso", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_btag_noniso_h.reset(new MuonHists(ctx,"muons_btag_noniso"));
   events_btag_noniso_h.reset(new EventHists(ctx,"events_btag_noniso"));
  
   electrons_jet123metdphi_MT2WCut_h.reset(new ElectronHists(ctx,"electrons_jet123metdphi_MT2WCut"));
   jets_jet123metdphi_MT2WCut_h.reset(new JetHists(ctx,"jets_jet123metdphi_MT2WCut"));
   topjets_jet123metdphi_MT2WCut_h.reset(new TopJetHists(ctx,"topjets_jet123metdphi_MT2WCut", 4, "patJetsHepTopTagCHSPacked_daughters"));
   muons_jet123metdphi_MT2WCut_h.reset(new MuonHists(ctx,"muons_jet123metdphi_MT2WCut"));
   events_jet123metdphi_MT2WCut_h.reset(new EventHists(ctx,"events_jet123metdphi_MT2WCut"));

   likelihood_metlepdphi_mtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_metlepdphi_mtlep"));
   likelihood_cut2_noiso_mtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_cut2_noiso_mtlep"));
             
   likelihood_preselection_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_preselection"));
   likelihood_leptonveto_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_leptonveto"));
   likelihood_met_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met"));
   likelihood_met100_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met100"));
   likelihood_twodcut_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_twodcut"));
   likelihood_metlepdphi_noiso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_metlepdphi_noiso"));
   likelihood_jetmetdphi_h.reset(new ttDMReconstructionHists_Likelihood(ctx,"likelihood_jetmetdphi"));
   likelihood_noiso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_noiso"));
   likelihood_cut1_noiso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_cut1_noiso"));
   likelihood_cut2_noiso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_cut2_noiso"));
   likelihood_metlepdphi_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_metlepdphi_h"));
   likelihood_heptoptagevent_WP3_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3"));
   likelihood_heptoptagevent_WP1_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP1"));
   likelihood_heptoptagevent_WP3_wobtag_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3_wobtag"));
   likelihood_heptoptagevent_WP1_wobtag_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP1_wobtag"));
   hists_metlepdphi_mtlep.reset(new ttDMSelectionHists(ctx, "hists_metlepdphi_mtlep"));
   heptoptagevent_WP1_wobtag_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP1_wobtag"));
   heptoptagevent_WP1_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP1"));
    
   likelihood_thirdjet_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_thirdjet"));
   likelihood_mediumMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_mediumMET"));
   likelihood_jetmetdphijet2_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_jetmetdphijet2"));
   likelihood_mtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_mtlep"));

   highMET_h.reset(new ttDMSelectionHists(ctx, "highMET"));
   jetmetdphi_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi"));
   mtlep_noniso_h.reset(new ttDMSelectionHists(ctx, "mtlep_noniso"));
   highMET_noniso_h.reset(new ttDMSelectionHists(ctx, "highMET_noniso"));
   jetmetdphi_noniso_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi_noniso"));
   metlepdphi_h.reset(new ttDMSelectionHists(ctx, "metlepdphi_h"));
   jetmetdphi_mediumMET_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi_mediumMET_h"));
   cut2_noiso_mtlep.reset(new ttDMSelectionHists(ctx, "cut2_noiso_mtlep"));
   
   genhists_nonisolep_h.reset(new  ttDMGenHists(ctx, "genhists_nonisolep_h"));
 
   likelihood_heptoptagevent_WP3_wobtag_highMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3_wobtag_highMET"));
   heptoptagevent_WP3_wobtag_highMET_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP3_wobtag_highMET"));
   likelihood_heptoptagevent_WP1_wobtag_highMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP1_wobtag_highMET"));
   heptoptagevent_WP1_wobtag_highMET_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP1_wobtag_highMET"));
   likelihood_heptoptagevent_WP1_highMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP1_highMET"));
   heptoptagevent_WP1_highMET_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP1_highMET"));
   likelihood_jetmetdphi_mediumMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_jetmetdphi_mediumMET"));
   likelihood_thirdjet_mtlep_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_thirdjet_mtlep"));
   thirdjet_mtlep_h.reset(new ttDMSelectionHists(ctx, "thirdjet_mtlep"));
   likelihood_thirdjet_jetmetdphi_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_thirdjet_jetmetdphi"));
   thirdjet_jetmetdphi_h.reset(new ttDMSelectionHists(ctx, "thirdjet_jetmetdphi"));
   jet123metdphi_h.reset(new ttDMSelectionHists(ctx, "jet123metdphi_h"));
   hists_likelihood_jet123metdphi_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_jet123metdphi_h"));
   jet123metdphi_MT2WCut_h.reset(new ttDMSelectionHists(ctx, "jet123metdphi_MT2WCut"));
   hists_likelihood_jet123metdphi_MT2WCut_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_jet123metdphi_MT2WCut"));
   metlepdphi_MT2WCut_h.reset(new ttDMSelectionHists(ctx, "metlepdphi_MT2WCut"));
   likelihood_metlepdphi_MT2WCut_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_metlepdphi_MT2WCut"));
   jet123metdphi_MT2WCut_MET240_h.reset(new ttDMSelectionHists(ctx, "jet123metdphi_MT2WCut_MET240"));
   likelihood_jet123metdphi_MT2WCut_MET240_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_jet123metdphi_MT2WCut_MET240"));
   thirdjet_MT2WCut_MetSel_h.reset(new ttDMSelectionHists(ctx, "thirdjet_MT2WCut_MetSel"));
   likelihood_thirdjet_MT2WCut_MetSel_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_thirdjet_MT2WCut_MetSel_h"));
   thirdjet_MT2WCut_h.reset(new ttDMSelectionHists(ctx, "thirdjet_MT2WCut_h"));
   likelihood_thirdjet_MT2WCut_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_thirdjet_MT2WCut_h"));
   MT2WCut_sel.reset(new MT2WCut(200));

   ttbargenhists_preselection_h.reset(new TTbarGenHists(ctx, "ttbargenhists_preselection"));
   ttbargenhists_leptonveto_h.reset(new TTbarGenHists(ctx, "ttbargenhists_leptonveto"));
   ttbargenhists_heptoptagevent_WP3_h.reset(new TTbarGenHists(ctx, "ttbargenhists_heptoptagevent_WP3"));
   ttbargenhists_jetmetdphi_mediumMET_h.reset(new TTbarGenHists(ctx, "ttbargenhists_jetmetdphi_mediumMET"));
   ttbargenhists_metlepdphi_mtlep.reset(new TTbarGenHists(ctx, "ttbargenhists_metlepdphi_mtlep"));
   ttbargenhists_metlepdphi_MT2WCut_h.reset(new TTbarGenHists(ctx, "ttbargenhists_metlepdphi_MT2WCut"));
   ttbargenhists_met_h.reset(new TTbarGenHists(ctx,"ttbargenhists_met"));
   ttbargenhists_btagsel_h.reset(new TTbarGenHists(ctx,"ttbargenhists_btagsel"));

   jetsel1208050_h.reset(new ttDMSelectionHists(ctx, "jetsel1208050"));
   likelihood_jetsel1208050_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_jetsel1208050"));

   met80_h.reset(new ttDMSelectionHists(ctx, "met80"));
   likelihood_met80.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met80"));

   btagsel_h.reset(new ttDMSelectionHists(ctx, "btagsel"));
   likelihood_btagsel.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btagsel"));

   heptoptagevent_WP3_highMET_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP3_highMET"));
   likelihood_heptoptagevent_WP3_highMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_heptoptagevent_WP3_highMET"));

   met160_h.reset(new ttDMSelectionHists(ctx, "met160"));
   likelihood_met160.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met160"));
   met160_btag_h.reset(new ttDMSelectionHists(ctx, "met160_btag"));
   likelihood_met160_btag.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met160_btag"));

   hists_likelihood_DeltaPhiMetNeutrino.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_DeltaPhiMetNeutrino"));
   hists_DeltaPhiMetNeutrino.reset(new ttDMSelectionHists(ctx, "DeltaPhiMetNeutrino"));

   met160_btag_metlep100_h.reset(new ttDMSelectionHists(ctx, "met160_btag_metlep100"));
   likelihood_met160_btag_metlep100.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met160_btag_metlep100"));
   met160_btag_metlepdphi_h.reset(new ttDMSelectionHists(ctx, "met160_btag_metlepdphi"));
   likelihood_met160_btag_metlepdphi.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met160_btag_metlepdphi"));
 
   DeltaPhiTaggedjet_Neutrino.reset(new ttDMSelectionHists(ctx, "DeltaPhiTaggedjet_Neutrino"));
   likelihood_DeltaPhiTaggedjet_Neutrino.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_DeltaPhiTaggedjet_Neutrino"));

    met100_btag_mtlep100_h.reset(new ttDMSelectionHists(ctx, "met100_btag_mtlep100"));
    likelihood_met100_btag_mtlep100.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met100_btag_mtlep100"));

    met100_btag_neutrinopT_h.reset(new ttDMSelectionHists(ctx, "met100_btag_neutrinopT"));
    likelihood_met100_btag_neutrinopT.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met100_btag_neutrinopT"));
    
    met100_btag_DeltaPhiTaggedJetTopLep_h.reset(new ttDMSelectionHists(ctx, "met100_btag_DeltaPhiTaggedJetTopLep"));
    likelihood_met100_btag_DeltaPhiTaggedJetTopLep.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met100_btag_DeltaPhiTaggedJetTopLep"));

    DeltaPhiTaggedjet_Neutrino_mtlep150.reset(new ttDMSelectionHists(ctx, "DeltaPhiTaggedjet_Neutrino_mtlep150"));
    likelihood_DeltaPhiTaggedjet_Neutrino_mtlep150.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_DeltaPhiTaggedjet_Neutrino_mtlep150"));

    HEPTT_WP3_noniso_h.reset(new ttDMSelectionHists(ctx, "HEPTT_WP3_noniso"));
    likelihood_WP3_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_WP3_noniso"));
    likelihood_btag_noniso_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_btag_noniso"));

    met50_btag_ttbarpTSel.reset(new ttDMSelectionHists(ctx, "met50_btag_ttbarpTSel"));
    likelihood_met50_btag_ttbarpTSel.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met50_btag_ttbarpTSel"));

    met100_btag_neutrinopT_ttbarpTSel_h.reset(new ttDMSelectionHists(ctx, "met100_btag_neutrinopT_ttbarpTSel"));
    likelihood_met100_btag_neutrinopT_ttbarpTSel.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_met100_btag_neutrinopT_ttbarpTSel"));
}

bool ttDMSelectionModuleAfterLikelihood::process(Event & event){
   
   input_h->fill(event);
   
   if (is_mc) {
      ttgenprod->process(event);
      TTbarGen ttbargen = event.get(h_ttbargen);
      if (isSemiLept && !ttbargen.IsSemiLeptonicDecay()) return false;
      if (isDiLept && (ttbargen.DecayChannel() < 4 || ttbargen.DecayChannel() >9)) return false;
      if (isOther && ttbargen.DecayChannel()!=0 && ttbargen.DecayChannel() != 10) return false;
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
 
   //MET >100GeV
   bool pass_met100 = MET100_sel->passes(event);
   if (pass_met100){
      met100_h->fill(event);  //MET >100GeV
      electrons_met100_h->fill(event);
      jets_met100_h->fill(event);
      muons_met100_h->fill(event);
      events_met100_h->fill(event);
      topjets_met100_h->fill(event);
      if(!event.isRealData) genhists_met100_h->fill(event);
      likelihood_met100_h->fill(event);
   }
   

   bool pass_met160 = MET160_sel->passes(event);
   //MET >160GeV
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
   if (!pass_lepVeto) return false;
   
   likelihood_leptonveto_h->fill(event);        
   leptonveto_h->fill(event);
   electrons_leptonveto_h->fill(event);
   muons_leptonveto_h->fill(event);
   events_leptonveto_h->fill(event);
   jets_leptonveto_h->fill(event);
   topjets_leptonveto_h->fill(event);
   genhists_leptonveto_h->fill(event);
   ttbargenhists_leptonveto_h->fill(event);
    
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
         
         bool pass_highMET_sel = MET_sel->passes(event);
         if (pass_HEPTT_WP3_wobtag)
            {
               heptoptagevent_WP3_wobtag_h->fill(event);
               electrons_heptoptagevent_WP3_wobtag_h->fill(event);
               jets_heptoptagevent_WP3_wobtag_h->fill(event);
               muons_heptoptagevent_WP3_wobtag_h->fill(event);
               events_heptoptagevent_WP3_wobtag_h->fill(event);
               topjets_heptoptagevent_WP3_wobtag_h->fill(event); 
               topjets_heptoptagevent_WP3_wobtag_tagged_h->fill(event); 
               likelihood_heptoptagevent_WP3_wobtag_h->fill(event);
               if (pass_highMET_sel){
                  likelihood_heptoptagevent_WP3_wobtag_highMET_h->fill(event);
                  heptoptagevent_WP3_wobtag_highMET_h->fill(event);
               }
            }
         
         if (pass_HEPTT_WP3)
            {
               heptoptagevent_WP3_h->fill(event);
               electrons_heptoptagevent_WP3_h->fill(event);
               jets_heptoptagevent_WP3_h->fill(event);
               muons_heptoptagevent_WP3_h->fill(event);
               events_heptoptagevent_WP3_h->fill(event);
               topjets_heptoptagevent_WP3_h->fill(event); 
               topjets_heptoptagevent_WP3_tagged_h->fill(event); 
               likelihood_heptoptagevent_WP3_h->fill(event); 
               ttbargenhists_heptoptagevent_WP3_h->fill(event);                           
               
               if (pass_highMET_sel){
                  heptoptagevent_WP3_highMET_h->fill(event);
                  likelihood_heptoptagevent_WP3_highMET_h->fill(event);
               }


               bool pass_metlepdphi= metleptondphi_sel->passes(event);
               bool pass_mtlep_sel = mtlep_sel->passes(event);
               bool pass_jetmetdphi_sel = jetmetdphi_sel->passes(event);
               //      bool pass_jetmetdphi_sel2 = jetmetdphijet2_sel->passes(event);
               bool pass_thirdjet_sel = thirdjet_sel->passes(event);
               bool pass_MT2WCut = MT2WCut_sel->passes(event);
               bool pass_btagsel = btag_sel->passes(event);
               if (pass_btagsel){
                  btagsel_h->fill(event);
                  likelihood_btagsel->fill(event);
                  electrons_btagsel_h->fill(event);
                  jets_btagsel_h->fill(event);
                  muons_btagsel_h->fill(event);
                  events_btagsel_h->fill(event);
                  topjets_btagsel_h->fill(event); 
                  ttbargenhists_btagsel_h->fill(event); 
                  
                  if (pass_MT2WCut)
                     {
                        btagsel_MT2WCut_h->fill(event);
                        likelihood_btagsel_MT2WCut->fill(event);
                        electrons_btagsel_MT2WCut_h->fill(event);
                        jets_btagsel_MT2WCut_h->fill(event);
                        muons_btagsel_MT2WCut_h->fill(event);
                        events_btagsel_MT2WCut_h->fill(event);
                        topjets_btagsel_MT2WCut_h->fill(event); 
                        ttbargenhists_btagsel_MT2WCut_h->fill(event);  
                        if (pass_highMET_sel){
                           btagsel_MT2WCut_metsel_h->fill(event);
                           likelihood_btagsel_MT2WCut_metsel->fill(event);
                           electrons_btagsel_MT2WCut_metsel_h->fill(event);
                           jets_btagsel_MT2WCut_metsel_h->fill(event);
                           muons_btagsel_MT2WCut_metsel_h->fill(event);
                           events_btagsel_MT2WCut_metsel_h->fill(event);
                           topjets_btagsel_MT2WCut_metsel_h->fill(event); 
                           ttbargenhists_btagsel_MT2WCut_metsel_h->fill(event); 
                        }

                     }
                  if (pass_mtlep_sel)
                     {
                        btagsel_mtlep_h->fill(event);          
                        likelihood_btagsel_mtlep->fill(event);
                        electrons_btagsel_mtlep_h->fill(event);
                        jets_btagsel_mtlep_h->fill(event);
                        muons_btagsel_mtlep_h->fill(event);
                        events_btagsel_mtlep_h->fill(event);
                        topjets_btagsel_mtlep_h->fill(event); 
                        ttbargenhists_btagsel_mtlep_h->fill(event);  
                        bool pass_met100 = MET100_sel->passes(event);
                        if (pass_met100){
                           //jetmetdphi_mediumMET_h->fill(event);
                           //likelihood_jetmetdphi_mediumMET_h->fill(event);
                           if (pass_jetmetdphi_sel){
                              jetmetdphi_mediumMET_h->fill(event);
                              likelihood_jetmetdphi_mediumMET_h->fill(event);
                           }
                        }
                        
                     }

               }               

               if (pass_mtlep_sel) mtlep_h->fill(event);
               if (pass_jetmetdphi_sel) 
                  {
                     // jetmetdphi_h->fill(event);
                     // if(mediumMET_sel->passes(event)) {
                     //    jetmetdphi_mediumMET_h->fill(event);
                     //    likelihood_jetmetdphi_mediumMET_h->fill(event);
                     //    ttbargenhists_jetmetdphi_mediumMET_h->fill(event);
                     // }
                     likelihood_jetmetdphi_h->fill(event);        
                  }
               if (pass_metlepdphi) 
                  {
                     metlepdphi_h->fill(event);
                     likelihood_metlepdphi_h->fill(event);
                     ttbargenhists_metlepdphi_mtlep->fill(event);
                     // likelihood_metlepdphi_mtlep_h->fill(event);
                     // hists_metlepdphi_mtlep->fill(event);
                     
                     // metlepdphi_MT2WCut_h->fill(event);
                     // likelihood_metlepdphi_MT2WCut_h->fill(event);
                     // ttbargenhists_metlepdphi_MT2WCut_h->fill(event);
                     
                  }
               if (pass_thirdjet_sel){
                  bool pass_jet123metdphi_sel=jet123metdphi_sel->passes(event);
                  bool pass_DeltaPhiMetNeutrino_sel=DeltaPhiMetNeutrino_sel->passes(event);
                  bool pass_met100 = MET100_sel->passes(event);
                  bool pass_mtlep100 = mtlep100_sel->passes(event);
                  if (pass_jet123metdphi_sel){
                     jet123metdphi_h->fill(event);
                     hists_likelihood_jet123metdphi_h->fill(event);
                     if (pass_DeltaPhiMetNeutrino_sel && pass_btagsel){
                        met50_btag_h->fill(event);
                        likelihood_met50_btag->fill(event);
                        muons_met50_btag_h->fill(event);
                        jets_met50_btag_h->fill(event);
                        events_met50_btag_h->fill(event);
                        topjets_met50_btag_h->fill(event);
                        bool pass_ttbarpTSel=ttbarpT_sel->passes(event);
                        if (pass_ttbarpTSel){
                           met50_btag_ttbarpTSel->fill(event);
                           likelihood_met50_btag_ttbarpTSel->fill(event);
                        }
                        
                        bool pass_DeltaPhiTaggedjet_Neutrino =DeltaPhiTaggedjet_Neutrino_sel->passes(event);
                        if (pass_DeltaPhiTaggedjet_Neutrino && pass_mtlep100){
                           DeltaPhiTaggedjet_Neutrino->fill(event);  
                           likelihood_DeltaPhiTaggedjet_Neutrino->fill(event);
                           bool pass_mtlep150 = mtlep150_sel->passes(event);
                           if (pass_mtlep150)
                              {
                                 DeltaPhiTaggedjet_Neutrino_mtlep150->fill(event);  
                                 likelihood_DeltaPhiTaggedjet_Neutrino_mtlep150->fill(event);
                              }
                        }
                     }
                     if (pass_met100 && pass_DeltaPhiMetNeutrino_sel && pass_btagsel){
                        met100_btag_h->fill(event);
                        likelihood_met100_btag->fill(event);
                        muons_met100_btag_h->fill(event);
                        jets_met100_btag_h->fill(event);
                        events_met100_btag_h->fill(event);
                        topjets_met100_btag_h->fill(event);
                        if (pass_mtlep100){
                           met100_btag_mtlep100_h->fill(event);
                           likelihood_met100_btag_mtlep100->fill(event);
                        }
                        bool pass_neutrinopT_sel=neutrinopT_sel->passes(event);
                        if (pass_neutrinopT_sel){
                           met100_btag_neutrinopT_h->fill(event);
                           likelihood_met100_btag_neutrinopT->fill(event);
                           bool pass_ttbarpTSel=ttbarpT_sel->passes(event);
                           if (pass_ttbarpTSel){
                              met100_btag_neutrinopT_ttbarpTSel_h->fill(event);
                              likelihood_met100_btag_neutrinopT_ttbarpTSel->fill(event);
                           }
                        }
                        bool pass_DeltaPhiTaggedJetTopLep_sel=DeltaPhiTaggedJetTopLep_sel->passes(event);
                        if (pass_DeltaPhiTaggedJetTopLep_sel){
                           met100_btag_DeltaPhiTaggedJetTopLep_h->fill(event);
                           likelihood_met100_btag_DeltaPhiTaggedJetTopLep->fill(event);
                           //cut hat nicht ricgtig funktioniert
                        }
                     }
                     bool pass_met160= MET160_sel->passes(event);
                     if (pass_met160)
                        {
                           met160_h->fill(event);  //FSP bestes fuer DMMET
                           likelihood_met160->fill(event);
                            bool pass_DeltaPhiMetNeutrino_sel=DeltaPhiMetNeutrino_sel->passes(event);
                            if (pass_DeltaPhiMetNeutrino_sel){
                               hists_likelihood_DeltaPhiMetNeutrino->fill(event);
                               hists_DeltaPhiMetNeutrino->fill(event);                     
                               if (pass_btagsel)
                                  {
                                     met160_btag_h->fill(event);
                                     likelihood_met160_btag->fill(event);
                                     muons_met160_btag_h->fill(event);
                                     jets_met160_btag_h->fill(event);
                                     events_met160_btag_h->fill(event);
                                     topjets_met160_btag_h->fill(event);
                                     if (pass_metlepdphi)
                                        {
                                           met160_btag_metlepdphi_h->fill(event);
                                           likelihood_met160_btag_metlepdphi->fill(event);
                                        }
                                     
                                     if (pass_mtlep100){
                                        met160_btag_metlep100_h->fill(event);
                                        likelihood_met160_btag_metlep100->fill(event);
                                     }

                                  }
                            }
                            
                           // if (pass_btagsel)
                           //    {
                           //       met160_btag_h->fill(event);
                           //       likelihood_met160_btag->fill(event);
                           //    }
                        }

                     
                     bool pass_MetSel = MET_sel->passes(event);  
                     if (pass_MetSel)
                        {
                           jet123metdphi_MT2WCut_h->fill(event);
                           hists_likelihood_jet123metdphi_MT2WCut_h->fill(event);
                           electrons_jet123metdphi_MT2WCut_h->fill(event);
                           jets_jet123metdphi_MT2WCut_h->fill(event);
                           muons_jet123metdphi_MT2WCut_h->fill(event);
                           events_jet123metdphi_MT2WCut_h->fill(event);
                           topjets_jet123metdphi_MT2WCut_h->fill(event);
                          
                           if (pass_MetSel){
                              jet123metdphi_MT2WCut_MET240_h->fill(event);
                              likelihood_jet123metdphi_MT2WCut_MET240_h->fill(event);
                           }
                        }
                  }
               }
               bool pass_jetsel1208050 = jetsel1208050->passes(event);
               if (pass_jetsel1208050)
                  {
                     jetsel1208050_h->fill(event);
                     likelihood_jetsel1208050_h->fill(event);
                     bool pass_met80= MET80_sel->passes(event);
                     if (pass_met80)
                        {
                           met80_h->fill(event);
                           likelihood_met80->fill(event);
                        }
                  }

            }
         else topjets_isolep_notag_WP3_h->fill(event);
         
         if (pass_HEPTT_WP1_wobtag)  {
            heptoptagevent_WP1_wobtag_h->fill(event);
            likelihood_heptoptagevent_WP1_wobtag_h->fill(event); 
            if (pass_highMET_sel){
               likelihood_heptoptagevent_WP1_wobtag_highMET_h->fill(event);
               heptoptagevent_WP1_wobtag_highMET_h->fill(event);
            }
         }
         if (pass_HEPTT_WP1)  {
            heptoptagevent_WP1_h->fill(event);
            likelihood_heptoptagevent_WP1_h->fill(event); 
            if (pass_highMET_sel){
               likelihood_heptoptagevent_WP1_highMET_h->fill(event);
               heptoptagevent_WP1_highMET_h->fill(event);
            }
         }
         else topjets_isolep_notag_WP1_h->fill(event);  //histograms filled with uncorrected mass
         
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

   bool pass_btag_sel = btag_sel->passes(event);
   if (!pass_btag_sel) return false;
   btag_noniso_h->fill(event);
   electrons_btag_noniso_h->fill(event);
   muons_btag_noniso_h->fill(event);
   events_btag_noniso_h->fill(event);
   jets_btag_noniso_h->fill(event);
   topjets_btag_noniso_h->fill(event);
   likelihood_btag_noniso_h->fill(event);

   if (pass_HEPTT_WP3){
      HEPTT_WP3_noniso_h->fill(event);
      likelihood_WP3_noniso_h->fill(event);

   }

   //   bool pass_medMET_sel = mediumMET_sel->passes(event);
   bool pass_thirdjet_sel = thirdjet_sel->passes(event);
   bool pass_mtlep_sel = mtlep_sel->passes(event); 
   bool pass_jetmetdphijet2_sel = jetmetdphijet2_sel->passes(event);
   bool pass_jetmetdphi_sel = jetmetdphi_sel->passes(event);
   bool pass_metlepdphi= metleptondphi_sel->passes(event);
   bool pass_MetSel = MET_sel->passes(event);
   bool pass_MT2WCut = MT2WCut_sel->passes(event);
 
   if (pass_jetmetdphijet2_sel){
      mediumMET_h->fill(event);
      electrons_mediumMET_h->fill(event);
      muons_mediumMET_h->fill(event);
      events_mediumMET_h->fill(event);
      jets_mediumMET_h->fill(event);
      topjets_mediumMET_h->fill(event);
      likelihood_mediumMET_h->fill(event);  
   }
   
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
         if (pass_MT2WCut){
            thirdjet_MT2WCut_h->fill(event);
            likelihood_thirdjet_MT2WCut_h->fill(event); 
            if (pass_MetSel){
               thirdjet_MT2WCut_MetSel_h->fill(event);
               likelihood_thirdjet_MT2WCut_MetSel_h->fill(event);
            }
            
         }

         if (pass_mtlep_sel){
            likelihood_thirdjet_mtlep_h->fill(event);
            thirdjet_mtlep_h->fill(event);
         }
         if (pass_jetmetdphi_sel&&pass_metlepdphi){
            likelihood_thirdjet_jetmetdphi_h->fill(event); 
            thirdjet_jetmetdphi_h->fill(event);
            electrons_thirdjet_jetmetdphi_h->fill(event);
            muons_thirdjet_jetmetdphi_h->fill(event);
            events_thirdjet_jetmetdphi_h->fill(event);
            jets_thirdjet_jetmetdphi_h->fill(event);
            topjets_thirdjet_jetmetdphi_h->fill(event);
         }
      }
   }
   if(pass_metlepdphi){
      metlepdphi_noiso_h->fill(event);
      electrons_metlepdphi_noiso_h->fill(event);
      muons_metlepdphi_noiso_h->fill(event);
      events_metlepdphi_noiso_h->fill(event);
      jets_metlepdphi_noiso_h->fill(event);
      topjets_metlepdphi_noiso_h->fill(event);
      likelihood_metlepdphi_noiso_h->fill(event);
   }

  if (pass_mtlep_sel){
      mtlep_h->fill(event);
      electrons_mtlep_h->fill(event);
      muons_mtlep_h->fill(event);
      events_mtlep_h->fill(event);
      jets_mtlep_h->fill(event);
      topjets_mtlep_h->fill(event);
      likelihood_mtlep_h->fill(event);        
  }
  
  if (pass_jetmetdphijet2_sel){
     jetmetdphijet2_h->fill(event);
     electrons_jetmetdphijet2_h->fill(event);
     muons_jetmetdphijet2_h->fill(event);
     events_jetmetdphijet2_h->fill(event);
     jets_jetmetdphijet2_h->fill(event);
     topjets_jetmetdphijet2_h->fill(event);
     likelihood_jetmetdphijet2_h->fill(event);
         bool likelihood_cut1_nosiso = cut1_likelihood->passes(event);
         if (likelihood_cut1_nosiso) likelihood_cut1_noiso_h->fill(event);
         bool likelihood_cut2_nosiso = cut2_likelihood->passes(event);
         if (likelihood_cut2_nosiso) 
            {
               likelihood_cut2_noiso_h->fill(event);
               bool pass_mtlep_sel = mtlep_sel->passes(event);
               if (pass_mtlep_sel){
                  likelihood_cut2_noiso_mtlep_h->fill(event);
                  cut2_noiso_mtlep->fill(event);
               }
            }
      }
   return true;
}



UHH2_REGISTER_ANALYSIS_MODULE(ttDMSelectionModuleAfterLikelihood)

