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
      
   // selections
   std::unique_ptr<Selection> heptoptagevent_WP1_sel;
   std::unique_ptr<Selection> heptoptagevent_WP3_sel;   
   std::unique_ptr<Selection> pv_sel;
   std::unique_ptr<Selection> highMET_sel;
   std::unique_ptr<AndSelection> lep1_sel;

   std::unique_ptr<Hists> highMET_h;
   std::unique_ptr<Hists> mtlep_noniso_h;
   std::unique_ptr<Hists> highMET_noniso_h;
   std::unique_ptr<Hists> likelihood_hists_noiso;
   std::unique_ptr<AnalysisModule> muo_med_noniso_SF;

   std::unique_ptr<AnalysisModule> collectionprod_heptt_WP3;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_WP3;
   std::unique_ptr<AnalysisModule> collectionprod_heptt_WP1;
   std::unique_ptr<AnalysisModule> collectionsizeprod_heptt_WP1;
   
   //std::unique_ptr<AnalysisModule> muo_triggerSF;
   std::unique_ptr<AnalysisModule> ttgenprod;
   
   // hists
   std::unique_ptr<Hists> lumihists;
   std::unique_ptr<Hists> filter_h;
   std::unique_ptr<Hists> input_h;
   std::unique_ptr<Hists> met_h, electrons_met_h,jets_met_h,muons_met_h,events_met_h,topjets_met_h, genhists_met_h;
   std::unique_ptr<Hists> nonisolep_h, electrons_nonisolep_h,jets_nonisolep_h,muons_nonisolep_h,events_nonisolep_h,topjets_nonisolep_h;
   std::unique_ptr<Hists> isolep_h, electrons_isolep_h,jets_isolep_h,muons_isolep_h,events_isolep_h,topjets_isolep_h;
   std::unique_ptr<TopJetHists> topjets_isolep_tagged_WP1_h, topjets_isolep_tagged_WP3_h;
   std::unique_ptr<Hists> heptoptagevent_WP3_h,electrons_heptoptagevent_WP3_h,jets_heptoptagevent_WP3_h,muons_heptoptagevent_WP3_h,events_heptoptagevent_WP3_h,topjets_heptoptagevent_WP3_h;
   std::unique_ptr<TopJetHists>  topjets_heptoptagevent_WP3_tagged_h;
   std::unique_ptr<Hists> twodcut_h,electrons_twodcut_h, jets_twodcut_h,topjets_twodcut_h,events_twodcut_h,muons_twodcut_h;
   std::unique_ptr<Hists> preselection_h,electrons_preselection_h, jets_preselection_h,topjets_preselection_h,events_preselection_h,muons_preselection_h;
   std::unique_ptr<Hists> heptoptagevent_WP1_h;
   std::unique_ptr<Hists> jetmetdphi_h;
   std::unique_ptr<Hists> jetmetdphi_noniso_h;
   std::unique_ptr<Hists> likelihood_hists_preselection;
   std::unique_ptr<Hists> likelihood_hists_met;
   std::unique_ptr<Hists> likelihood_hists_twodcut;
   std::unique_ptr<Hists> likelihood_hists_jetmetdphi_highMET_h;
   std::unique_ptr<Hists> likelihood_hists_highMET;
   std::unique_ptr<Hists> likelihood_hists_cut1_noiso;
   std::unique_ptr<Hists> likelihood_hists_cut2_noiso;
   std::unique_ptr<Hists> jetmetdphi_highMET_h,jetmetdphi_mediumMET_h;
  
   std::unique_ptr<Selection> jetmetdphi_sel, jetmetdphijet2_sel;
   std::unique_ptr<Selection> mtlep_sel,mediumMET_sel;
   std::unique_ptr<Selection> cut1_likelihood;
   std::unique_ptr<Selection> cut2_likelihood;
   std::unique_ptr<Selection> highpTJetSel;
   std::unique_ptr<Selection> leadingjet_sel;
   std::unique_ptr<Selection> HTlep_sel;
   std::unique_ptr<Selection> thirdjet_sel;

   std::unique_ptr<Hists> genhists_nonisolep_h;
   std::unique_ptr<Hists> likelihood_hists_jetmetdphi_highMET_mtlep, likelihood_hists_heptoptagevent_WP3,likelihood_hists_heptoptagevent_WP1;
   std::unique_ptr<Hists> hists_jetmetdphi_highMET_mtlep;
   std::unique_ptr<Hists> cut2_noiso_mtlep;
   std::unique_ptr<Hists> likelihood_hists_cut2_noiso_mtlep;

   std::unique_ptr<Hists> mediumMET_h,electrons_mediumMET_h,jets_mediumMET_h,topjets_mediumMET_h,muons_mediumMET_h,events_mediumMET_h, likelihood_hists_mediumMET;
   std::unique_ptr<Hists> thirdjet_h,electrons_thirdjet_h,jets_thirdjet_h,topjets_thirdjet_h,muons_thirdjet_h,events_thirdjet_h,likelihood_hists_thirdjet;
   std::unique_ptr<Hists> mtlep_h,electrons_mtlep_h,jets_mtlep_h,topjets_mtlep_h,muons_mtlep_h,events_mtlep_h,likelihood_hists_mtlep;
   std::unique_ptr<Hists> jetmetdphijet2_h,electrons_jetmetdphijet2_h,jets_jetmetdphijet2_h,topjets_jetmetdphijet2_h,muons_jetmetdphijet2_h,events_jetmetdphijet2_h,likelihood_hists_jetmetdphijet2;

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
   //const TopJetId HTTTopJetId_WP3 = AndId<TopJet>(PtEtaCut(150., 2.4),HEPTopTagV2(85,280,0.47, subBtag), Tau32(0.97)); 
   const TopJetId HTTTopJetId_WP3 = AndId<TopJet>(PtEtaCut(150., 2.4),HEPTopTagV2(120,330,0.25,-0.76,0.24), Tau32(0.71));   //old HepTT, as long as subjet b-tag does not work
   //const TopJetId HTTTopJetId_WP1= AndId<TopJet>(PtEtaCut(150., 2.4),HEPTopTagV2(110,185,0.2, subBtag), Tau32(0.93));
   const TopJetId HTTTopJetId_WP1= AndId<TopJet>(PtEtaCut(150., 2.4),HEPTopTagV2(80,300,0.5,-0.67,0.71), Tau32(0.91));//old HepTT, as long as subjet b-tag does not work
   //// OBJ CLEANING
   muon_cleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuonIDMedium(), PtEtaCut(47., 2.5),MuonIso(0.15))));
   ////
   ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
   
   //// EVENT SELECTION
   bool muon(false), elec(false);
   const std::string channel(ctx.get("channel", ""));
   if(channel == "muon") muon = true;
   else if(channel == "electron") elec = true;
   else throw std::runtime_error("undefined argument for 'channel' key in xml file (must be 'muon' or 'electron'): "+channel);
   lep1_sel.reset(new AndSelection(ctx));
   if(muon){
      lep1_sel->add<NMuonSelection>("muoN == 1", 1, 1);
      lep1_sel->add<NElectronSelection>("eleN == 0", 0, 0);
   }
   else if(elec){
      lep1_sel->add<NMuonSelection>("muoN == 0", 0, 0);
      lep1_sel->add<NElectronSelection>("eleN == 1", 1, 1);
   }
   
   highpTJetSel.reset(new NJetSelection(2, -1, JetId(PtEtaCut( 100., 2.4))));
   HTlep_sel.reset(new HTlepCut(300, std::numeric_limits<double>::infinity()));
   
   leadingjet_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut( 150., 2.4))));
   highMET_sel.reset(new METCut(320., std::numeric_limits<double>::infinity()));
   mediumMET_sel.reset(new METCut(260., std::numeric_limits<double>::infinity()));
   mtlep_sel.reset(new MTlepCut(140., std::numeric_limits<double>::infinity()));
   jetmetdphi_sel.reset(new METJetDPhiCut(2.0, 1));
   jetmetdphijet2_sel.reset(new METJetDPhiCut(1.3, 2));
   thirdjet_sel.reset(new NJetSelection(3, -1, JetId(PtEtaCut( 50., 2.4))));


   cut1_likelihood.reset(new LikelihoodSelection(ctx,22));
   cut2_likelihood.reset(new LikelihoodSelection(ctx,14));
    
   //h_flag_heptoptagevent = ctx.declare_event_output<int>("flag_heptoptagevent");
   h_flag_twodcut = ctx.get_handle<bool>("flag_twodcut");
   h_likelihood = ctx.get_handle<double>("likelihood");
   h_recneutrino = ctx.get_handle<LorentzVector>("rec_neutrino");

   heptoptagevent_WP3_sel.reset(new HandleSelection<int>(ctx, "n_heptopjets_WP3", 1, 999));
   heptoptagevent_WP1_sel.reset(new HandleSelection<int>(ctx, "n_heptopjets_WP1", 1, 999));
   
   //heptt collection korrigieren
   heptopjet_corrector.reset(new GenericTopJetCorrector(ctx, JERFiles::Fall15_25ns_L123_AK8PFchs_MC,"patJetsHepTopTagCHSPacked_daughters"));
   
   collectionprod_heptt_WP3.reset(new CollectionProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "h_heptopjets_WP3", TopJetId(HTTTopJetId_WP3))); 
   collectionsizeprod_heptt_WP3.reset(new CollectionSizeProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "n_heptopjets_WP3", TopJetId(HTTTopJetId_WP3)));  
   
   collectionprod_heptt_WP1.reset(new CollectionProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "h_heptopjets_WP1", TopJetId(HTTTopJetId_WP1))); 
   collectionsizeprod_heptt_WP1.reset(new CollectionSizeProducer<TopJet>(ctx, "patJetsHepTopTagCHSPacked_daughters", "n_heptopjets_WP1", TopJetId(HTTTopJetId_WP1))); 
   
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
   
   heptoptagevent_WP3_h.reset(new ttDMSelectionHists(ctx,"heptoptagevent_WP3"));
   electrons_heptoptagevent_WP3_h.reset(new ElectronHists(ctx,"electrons_heptoptagevent_WP3"));
   jets_heptoptagevent_WP3_h.reset(new JetHists(ctx,"jets_heptoptagevent_WP3"));
   muons_heptoptagevent_WP3_h.reset(new MuonHists(ctx,"muons_heptoptagevent_WP3"));
   events_heptoptagevent_WP3_h.reset(new EventHists(ctx,"events_heptoptagevent_WP3"));
   topjets_heptoptagevent_WP3_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3", 4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_heptoptagevent_WP3_tagged_h.reset(new TopJetHists(ctx,"topjets_heptoptagevent_WP3_tagged", 4,"patJetsHepTopTagCHSPacked_daughters"));
   topjets_heptoptagevent_WP3_tagged_h->set_TopJetId(HTTTopJetId_WP3);
   
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

   likelihood_hists_jetmetdphi_highMET_mtlep.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_jetmetdphi_highMET_mtlep"));
   likelihood_hists_cut2_noiso_mtlep.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_cut2_noiso_mtlep"));
                  
   likelihood_hists_preselection.reset(new ttDMReconstructionHists_Likelihood(ctx, "hists_likelihood_preselection"));
   likelihood_hists_met.reset(new ttDMReconstructionHists_Likelihood(ctx, "hists_likelihood_met"));
   likelihood_hists_twodcut.reset(new ttDMReconstructionHists_Likelihood(ctx, "hists_likelihood_twodcut"));
   likelihood_hists_highMET.reset(new ttDMReconstructionHists_Likelihood(ctx,"likelihood_hists_highMET"));
   likelihood_hists_noiso.reset(new ttDMReconstructionHists_Likelihood(ctx, "hists_likelihood_noiso"));
   likelihood_hists_cut1_noiso.reset(new ttDMReconstructionHists_Likelihood(ctx, "hists_likelihood_cut1_noiso"));
   likelihood_hists_cut2_noiso.reset(new ttDMReconstructionHists_Likelihood(ctx, "hists_likelihood_cut2_noiso"));
   likelihood_hists_jetmetdphi_highMET_h.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_jetmetdphi_highMET_h"));
   likelihood_hists_heptoptagevent_WP3.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_heptoptagevent_WP3"));
   likelihood_hists_heptoptagevent_WP1.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_heptoptagevent_WP1"));
   hists_jetmetdphi_highMET_mtlep.reset(new ttDMSelectionHists(ctx, "hists_jetmetdphi_highMET_mtlep"));
   heptoptagevent_WP1_h.reset(new ttDMSelectionHists(ctx, "heptoptagevent_WP1"));
   
   likelihood_hists_thirdjet.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_thirdjet"));
   likelihood_hists_mediumMET.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_mediumMET"));
   likelihood_hists_jetmetdphijet2.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_jetmetdphijet2"));
   likelihood_hists_mtlep.reset(new ttDMReconstructionHists_Likelihood(ctx, "likelihood_hists_mtlep"));

   highMET_h.reset(new ttDMSelectionHists(ctx, "highMET"));
   jetmetdphi_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi"));
   mtlep_noniso_h.reset(new ttDMSelectionHists(ctx, "mtlep_noniso"));
   highMET_noniso_h.reset(new ttDMSelectionHists(ctx, "highMET_noniso"));
   jetmetdphi_noniso_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi_noniso"));
   jetmetdphi_highMET_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi_highMET_h"));
   jetmetdphi_mediumMET_h.reset(new ttDMSelectionHists(ctx, "jetmetdphi_mediumMET_h"));
   cut2_noiso_mtlep.reset(new ttDMSelectionHists(ctx, "cut2_noiso_mtlep"));
   
   genhists_nonisolep_h.reset(new  ttDMGenHists(ctx, "genhists_nonisolep_h"));
 
   
 }

bool ttDMSelectionModuleAfterLikelihood::process(Event & event){
   
   input_h->fill(event);
   
   if (is_mc) ttgenprod->process(event);

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
   muo_med_noniso_SF->process(event);
   filter_h->fill(event); 
   
   bool pass_lep1 = lep1_sel->passes(event);
   if (event.isRealData) lumihists->fill(event);
      
   met_h->fill(event);
   electrons_met_h->fill(event);
   jets_met_h->fill(event);
   muons_met_h->fill(event);
   events_met_h->fill(event);
   topjets_met_h->fill(event);
   if(!event.isRealData) genhists_met_h->fill(event);
   likelihood_hists_met->fill(event); 
 
   bool pass_twodcut = event.get(h_flag_twodcut);
   
   if(pass_twodcut)                                               //pre-selection control plots 
       {
          likelihood_hists_preselection->fill(event);        
          preselection_h->fill(event);
          electrons_preselection_h->fill(event);
          muons_preselection_h->fill(event);
          events_preselection_h->fill(event);
          jets_preselection_h->fill(event);
          topjets_preselection_h->fill(event);
       }

   //----------------------------------------------------boosted analysis-------------------------------------- 
   
   heptopjet_corrector->process(event);
   collectionprod_heptt_WP3->process(event); 
   collectionsizeprod_heptt_WP3->process(event);
   collectionprod_heptt_WP1->process(event); 
   collectionsizeprod_heptt_WP1->process(event);
 
   //event.set(h_flag_heptoptagevent, int(pass_HEPTT)); 
   bool pass_HEPTT_WP3 = heptoptagevent_WP3_sel->passes(event);
   bool pass_HEPTT_WP1 = heptoptagevent_WP1_sel->passes(event);
   
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
         topjets_isolep_tagged_WP1_h->fill(event);
         topjets_isolep_tagged_WP3_h->fill(event);
         if (pass_HEPTT_WP3)
            {
               heptoptagevent_WP3_h->fill(event);
               electrons_heptoptagevent_WP3_h->fill(event);
               jets_heptoptagevent_WP3_h->fill(event);
               muons_heptoptagevent_WP3_h->fill(event);
               events_heptoptagevent_WP3_h->fill(event);
               topjets_heptoptagevent_WP3_h->fill(event); 
               topjets_heptoptagevent_WP3_tagged_h->fill(event); 
               likelihood_hists_heptoptagevent_WP3->fill(event); 

               bool pass_highMET_sel = highMET_sel->passes(event);
               bool pass_mtlep_sel = mtlep_sel->passes(event);
               bool pass_jetmetdphi_sel = jetmetdphi_sel->passes(event);
               bool pass_jetmetdphi_sel2 = jetmetdphijet2_sel->passes(event);
               if (pass_mtlep_sel) mtlep_h->fill(event);
               if (pass_jetmetdphi_sel && pass_jetmetdphi_sel2) 
                  {
                     jetmetdphi_h->fill(event);
                     if(mediumMET_sel->passes(event)) jetmetdphi_mediumMET_h->fill(event);;
                  }
               if (pass_highMET_sel) {
                  highMET_h->fill(event);
                  likelihood_hists_highMET->fill(event);        
                  if (pass_jetmetdphi_sel && pass_jetmetdphi_sel2) 
                     {
                        jetmetdphi_highMET_h->fill(event);
                        likelihood_hists_jetmetdphi_highMET_h->fill(event);
                        bool likelihood_cut2 = cut2_likelihood->passes(event);
                        if (likelihood_cut2) {
                           likelihood_hists_jetmetdphi_highMET_mtlep->fill(event);
                           hists_jetmetdphi_highMET_mtlep->fill(event);
                        }
                        
                     }
               }
            }
         if (pass_HEPTT_WP1)  {
            heptoptagevent_WP1_h->fill(event);
            likelihood_hists_heptoptagevent_WP1->fill(event); 
         }
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
   likelihood_hists_twodcut->fill(event);
   

   bool pass_medMET_sel = mediumMET_sel->passes(event);
   bool pass_thirdjet_sel = thirdjet_sel->passes(event);
   bool pass_mtlep_sel = mtlep_sel->passes(event); 
   bool pass_jetmetdphijet2_sel = jetmetdphijet2_sel->passes(event);

   if (pass_medMET_sel){
      mediumMET_h->fill(event);
      electrons_mediumMET_h->fill(event);
      muons_mediumMET_h->fill(event);
      events_mediumMET_h->fill(event);
      jets_mediumMET_h->fill(event);
      topjets_mediumMET_h->fill(event);
      likelihood_hists_mediumMET->fill(event);  
   }
   
   if (pass_thirdjet_sel){
      thirdjet_h->fill(event);
      electrons_thirdjet_h->fill(event);
      muons_thirdjet_h->fill(event);
      events_thirdjet_h->fill(event);
      jets_thirdjet_h->fill(event);
      topjets_thirdjet_h->fill(event);
      likelihood_hists_thirdjet->fill(event);        
   }

  if (pass_mtlep_sel){
      mtlep_h->fill(event);
      electrons_mtlep_h->fill(event);
      muons_mtlep_h->fill(event);
      events_mtlep_h->fill(event);
      jets_mtlep_h->fill(event);
      topjets_mtlep_h->fill(event);
      likelihood_hists_mtlep->fill(event);        
  }
  
  if (pass_jetmetdphijet2_sel){
     jetmetdphijet2_h->fill(event);
     electrons_jetmetdphijet2_h->fill(event);
     muons_jetmetdphijet2_h->fill(event);
     events_jetmetdphijet2_h->fill(event);
     jets_jetmetdphijet2_h->fill(event);
     topjets_jetmetdphijet2_h->fill(event);
     likelihood_hists_jetmetdphijet2->fill(event);       
   }


   bool pass_highMET_sel_noniso = highMET_sel->passes(event);
   bool pass_mtlep_sel_noniso = mtlep_sel->passes(event);
   bool pass_jetmetdphi_sel_noniso = jetmetdphi_sel->passes(event);
   bool pass_jetmetdphi_sel2_noniso = jetmetdphijet2_sel->passes(event);
   if (pass_mtlep_sel_noniso) mtlep_noniso_h->fill(event);
   if (pass_jetmetdphi_sel_noniso && pass_jetmetdphi_sel2_noniso) jetmetdphi_noniso_h->fill(event);
   if (pass_highMET_sel_noniso && pass_jetmetdphi_sel_noniso && pass_jetmetdphi_sel2_noniso) 
      {
         highMET_noniso_h->fill(event);
         likelihood_hists_noiso->fill(event);
         bool likelihood_cut1_nosiso = cut1_likelihood->passes(event);
         if (likelihood_cut1_nosiso) likelihood_hists_cut1_noiso->fill(event);
         bool likelihood_cut2_nosiso = cut2_likelihood->passes(event);
         if (likelihood_cut2_nosiso) 
            {
               likelihood_hists_cut2_noiso->fill(event);
               bool pass_mtlep_sel = mtlep_sel->passes(event);
               if (pass_mtlep_sel){
                  likelihood_hists_cut2_noiso_mtlep->fill(event);
                  cut2_noiso_mtlep->fill(event);
               }
            }
      }
   return true;
}



UHH2_REGISTER_ANALYSIS_MODULE(ttDMSelectionModuleAfterLikelihood)

