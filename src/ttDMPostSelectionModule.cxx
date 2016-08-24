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
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/Utils.h"

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
   std::unique_ptr<AnalysisModule> PUReweight; 
   std::unique_ptr<AnalysisModule> lumi_weight;

   //Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;

   JetId btagAK4_wp;
   Event::Handle<int> h_flag_toptagevent;
   Event::Handle<int> h_flag_type2toptagevent;
   Event::Handle<int> h_flag_heptoptagevent;
   
   bool muon =false;
   bool elec= false;

   //cleaner
   std::unique_ptr<MuonCleaner> muon_cleaner_iso;
   std::unique_ptr<ElectronCleaner> ele_cleaner_iso;
   
   // selections
   std::unique_ptr<Selection> btagAK4_sel;
   std::unique_ptr<Selection> leptoppt_sel;
   //std::unique_ptr<Selection> chi2_sel;
   std::unique_ptr<Selection> met_sel;
   std::unique_ptr<Selection> mtlep_sel;
   std::unique_ptr<Selection> jetmetdphi_sel;
   std::unique_ptr<Selection> mt2w_sel;
   std::unique_ptr<AndSelection> lepiso_sel;
   std::unique_ptr<Selection> jet_sel;
   
   // hists
   std::unique_ptr<Hists> hi_input, hi_electrons_input,  hi_muons_input, hi_jets_input,hi_events_input,hi_topjets_input;
   //std::unique_ptr<Hists> hi_input__hyp;
   std::unique_ptr<Hists> hi_leptoppt, hi_electrons_leptoppt,  hi_muons_leptoppt, hi_jets_leptoppt,hi_events_leptoppt,hi_topjets_leptoppt;
   //std::unique_ptr<Hists> hi_leptoppt__hyp;
   //std::unique_ptr<Hists> hi_chi2;
   //std::unique_ptr<Hists> hi_chi2__hyp;
   std::unique_ptr<Hists> met_h, hi_electrons_met,  hi_muons_met, hi_jets_met,hi_events_met,hi_topjets_met;
   std::unique_ptr<Hists> mtlep_h, hi_electrons_mtlep,  hi_muons_mtlep, hi_jets_mtlep,hi_events_mtlep,hi_topjets_mtlep;
   std::unique_ptr<Hists> jetmetdphi_h, hi_electrons_jetmetdphi,  hi_muons_jetmetdphi, hi_jets_jetmetdphi,hi_events_jetmetdphi,hi_topjets_jetmetdphi;
   std::unique_ptr<Hists> mt2w_h, hi_electrons_mt2w,  hi_muons_mt2w, hi_jets_mt2w,hi_events_mt2w,hi_topjets_mt2w;
   std::unique_ptr<Hists> hi_t0b0, hi_electrons_t0b0,  hi_muons_t0b0, hi_jets_t0b0,hi_events_t0b0,hi_topjets_t0b0;
  //std::unique_ptr<Hists> hi_t0b0__hyp;
   std::unique_ptr<Hists> hi_t0b1, hi_electrons_t0b1,  hi_muons_t0b1, hi_jets_t0b1,hi_events_t0b1,hi_topjets_t0b1;
  //std::unique_ptr<Hists> hi_t0b1__hyp;
   std::unique_ptr<Hists> hi_t1, hi_electrons_t1,  hi_muons_t1, hi_jets_t1,hi_events_t1,hi_topjets_t1;
  //std::unique_ptr<Hists> hi_t1__hyp;
   std::unique_ptr<Hists> lumihists;
   std::unique_ptr<Hists> boosted_h, resolved_h,boosted_heptoptag_h, resolved_heptoptag_h,nonIsolated_h,toptag_h,type2toptag_h;
};

ttDMPostSelectionModule::ttDMPostSelectionModule(Context& ctx){
   
   const std::string channel(ctx.get("channel", ""));
   if(channel == "muon") muon = true;
   else if(channel == "electron") elec = true;
   else throw std::runtime_error("undefined argument for 'channel' key in xml file (must be 'muon' or 'electron'): "+channel);

   //cleaner
   muon_cleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIso(0.12),PtEtaCut(50., 2.1))));
   ele_cleaner_iso.reset(new ElectronCleaner(AndId<Electron>(ElectronID_Spring15_50ns_tight, PtEtaCut(50., 2.5))));
      
   if (ctx.get("dataset_version").find("Data") ==std::string::npos)  
      {
         lumi_weight.reset(new MCLumiWeight(ctx));
         PUReweight.reset(new MCPileupReweight(ctx));
      }
   
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
   h_flag_type2toptagevent = ctx.declare_event_input<int>("flag_type2toptagevent");
   h_flag_heptoptagevent = ctx.declare_event_input<int>("flag_heptoptagevent");
   
   // SELECTION
   // if(elec) leptoppt_sel.reset(new LeptonicTopPtCut(ctx, 140., infinity, "TTbarReconstruction", "Chi2"));
   // else if(muon) leptoppt_sel.reset(new AndSelection(ctx));
   leptoppt_sel.reset(new AndSelection(ctx));
   
   //chi2_sel.reset(new HypothesisDiscriminatorCut(ctx, 0., 50., "TTbarReconstruction", "Chi2"));
   
   //ttDM SELECTION
   met_sel.reset(new METCut(320., std::numeric_limits<double>::infinity()));
   mtlep_sel.reset(new MTlepCut(160., std::numeric_limits<double>::infinity()));
   jetmetdphi_sel.reset(new METJetDPhiCut(1.2, 1));
   mt2w_sel.reset(new MT2WCut(200));

   lepiso_sel.reset(new AndSelection(ctx));
   if (elec) lepiso_sel->add<NElectronSelection>("eleN == 1", 1, 1);
   else if (muon) lepiso_sel->add<NMuonSelection>("muonN == 1", 1, 1);
  
   jet_sel.reset(new NJetSelection(3, -1, JetId(PtEtaCut( 50., 2.4))));     

   // HISTS
   hi_input.reset(new ttDMPostSelectionHists(ctx, "input"));
   //hi_input__hyp.reset(new HypothesisHists(ctx, "input__hyp_chi2min", "TTbarReconstruction", "Chi2"));

   hi_leptoppt.reset(new ttDMPostSelectionHists(ctx, "leptoppt"));
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
   
   boosted_h.reset(new ttDMPostSelectionHists(ctx,"boosted"));
   resolved_h.reset(new ttDMPostSelectionHists(ctx,"resolved"));
   boosted_heptoptag_h.reset(new ttDMPostSelectionHists(ctx,"boosted_heptoptag"));
   resolved_heptoptag_h.reset(new ttDMPostSelectionHists(ctx,"resolved_heptoptag"));
   
   lumihists.reset(new LuminosityHists(ctx, "lumi"));
   
   hi_electrons_input.reset(new ElectronHists(ctx, "electrons_input"));
   hi_muons_input.reset(new MuonHists(ctx, "muons_input"));
   hi_events_input.reset(new EventHists(ctx, "events_input"));
   hi_jets_input.reset(new JetHists(ctx, "jets_input"));
   hi_topjets_input.reset(new TopJetHists(ctx, "topjets_input"));
   
   hi_electrons_leptoppt.reset(new ElectronHists(ctx, "electrons_leptoppt"));
   hi_muons_leptoppt.reset(new MuonHists(ctx, "muons_leptoppt"));
   hi_events_leptoppt.reset(new EventHists(ctx, "events_leptoppt"));
   hi_jets_leptoppt.reset(new JetHists(ctx, "jets_leptoppt"));
   hi_topjets_leptoppt.reset(new TopJetHists(ctx, "topjets_leptoppt"));
   
   hi_electrons_met.reset(new ElectronHists(ctx, "electrons_met"));
   hi_muons_met.reset(new MuonHists(ctx, "muons_met"));
   hi_events_met.reset(new EventHists(ctx, "events_met"));
   hi_jets_met.reset(new JetHists(ctx, "jets_met"));
   hi_topjets_met.reset(new TopJetHists(ctx, "topjets_met"));
   
   hi_electrons_mtlep.reset(new ElectronHists(ctx, "electrons_mtlep"));
   hi_muons_mtlep.reset(new MuonHists(ctx, "muons_mtlep"));
   hi_events_mtlep.reset(new EventHists(ctx, "events_mtlep"));
   hi_jets_mtlep.reset(new JetHists(ctx, "jets_mtlep"));
   hi_topjets_mtlep.reset(new TopJetHists(ctx, "topjets_mtlep"));
   
   hi_electrons_jetmetdphi.reset(new ElectronHists(ctx, "electrons_jetmetdphi"));
   hi_muons_jetmetdphi.reset(new MuonHists(ctx, "muons_jetmetdphi"));
   hi_events_jetmetdphi.reset(new EventHists(ctx, "events_jetmetdphi"));
   hi_jets_jetmetdphi.reset(new JetHists(ctx, "jets_jetmetdphi"));
   hi_topjets_jetmetdphi.reset(new TopJetHists(ctx, "topjets_jetmetdphi"));
   
   hi_electrons_mt2w.reset(new ElectronHists(ctx, "electrons_mt2w"));
   hi_muons_mt2w.reset(new MuonHists(ctx, "muons_mt2w"));
   hi_events_mt2w.reset(new EventHists(ctx, "events_mt2w"));
   hi_jets_mt2w.reset(new JetHists(ctx, "jets_mt2w"));
   hi_topjets_mt2w.reset(new TopJetHists(ctx, "topjets_mt2w"));
   
   hi_electrons_t0b0.reset(new ElectronHists(ctx, "electrons_t0b0"));
   hi_muons_t0b0.reset(new MuonHists(ctx, "muons_t0b0"));
   hi_events_t0b0.reset(new EventHists(ctx, "events_t0b0"));
   hi_jets_t0b0.reset(new JetHists(ctx, "jets_t0b0"));
   hi_topjets_t0b0.reset(new TopJetHists(ctx, "topjets_t0b0")); 
   
   hi_electrons_t0b1.reset(new ElectronHists(ctx, "electrons_t0b1"));
   hi_muons_t0b1.reset(new MuonHists(ctx, "muons_t0b1"));
   hi_events_t0b1.reset(new EventHists(ctx, "events_t0b1"));
   hi_jets_t0b1.reset(new JetHists(ctx, "jets_t0b1"));
   hi_topjets_t0b1.reset(new TopJetHists(ctx, "topjets_t0b1")); 
   
   hi_electrons_t1.reset(new ElectronHists(ctx, "electrons_t1"));
   hi_muons_t1.reset(new MuonHists(ctx, "muons_t1"));
   hi_events_t1.reset(new EventHists(ctx, "events_t1"));
   hi_jets_t1.reset(new JetHists(ctx, "jets_t1"));
   hi_topjets_t1.reset(new TopJetHists(ctx, "topjets_t1")); 
   
   nonIsolated_h.reset(new ttDMPostSelectionHists(ctx, "nonIsolated"));
   toptag_h.reset(new ttDMPostSelectionHists(ctx, "toptag"));
   type2toptag_h.reset(new ttDMPostSelectionHists(ctx, "type2toptag"));
}

bool ttDMPostSelectionModule::process(Event& event) {

     if (event.isRealData) lumihists->fill(event);
     else
        {
           lumi_weight->process(event); 
           PUReweight ->process(event);
        }
   
     hi_input->fill(event);
     hi_electrons_input->fill(event);
     hi_events_input->fill(event);
     hi_muons_input->fill(event);
     hi_topjets_input->fill(event);
     hi_jets_input->fill(event);
   
     //
     //sensitivity study: non isolated lepton or require one hep top tag or a cms toptag
     //if not toptagged then type two top tagged
  
     bool toptag(event.get(h_flag_toptagevent));
     bool type2toptag(event.get(h_flag_type2toptagevent));
     //    bool heptoptag(event.get(h_flag_heptoptagevent));

     
     std::vector<Muon>* all_muons;
     std::vector<Electron>* all_electrons;
     if (muon) 
        {
         all_muons = (new std::vector<Muon>(*event.muons));
         muon_cleaner_iso->process(event);
        }
     else if (elec) 
        {
           all_electrons = (new std::vector<Electron>(*event.electrons));
         ele_cleaner_iso->process(event);
        }
     bool pass_lepiso =  lepiso_sel->passes(event); //require exactly one lepton
     bool pass_jetsel = jet_sel->passes(event);
     // if (pass_lepiso && !heptoptag)
     //    {
     //       resolved_heptoptag_h->fill(event);
     //    }
     // else boosted_heptoptag_h ->fill(event);
     
     if (pass_lepiso && !toptag && !type2toptag && pass_jetsel) 
        {
           resolved_h ->fill(event);
        }
     if (pass_lepiso && !toptag && !type2toptag) 
        {
           // TO DO: && !not type-two //hep top tagger //cms top tagger ohne n-subjettiness 
           
           if (muon) {
              event.muons->clear();
              event.muons->reserve(all_muons->size());
              for(auto & muon : *all_muons) event.muons->push_back(muon); 
           }
           else if (elec) {
              event.electrons->clear();
              event.electrons->reserve(all_electrons->size());
              for(auto & ele : *all_electrons) event.electrons->push_back(ele); 
           }
           if(!pass_lepiso) nonIsolated_h->fill(event);
           if(toptag) toptag_h->fill(event);
           if(type2toptag) type2toptag_h->fill(event);
           
           boosted_h->fill(event);
        }
     if (muon) {
        event.muons->clear();
        event.muons->reserve(all_muons->size());
        for(auto & muon : *all_muons) event.muons->push_back(muon); 
     }
     else if (elec) {
        event.electrons->clear();
        event.electrons->reserve(all_electrons->size());
        for(auto & ele : *all_electrons) event.electrons->push_back(ele); 
     }
     
     //hi_input__hyp->fill(event);
     
     // //// LEPTONIC-TOP pt selection
     // bool pass_leptoppt = leptoppt_sel->passes(event);
     // if(!pass_leptoppt) return false;
     // hi_leptoppt->fill(event);
     // //hi_leptoppt__hyp->fill(event);
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
     hi_electrons_met->fill(event);
     hi_events_met->fill(event);
     hi_muons_met->fill(event);
     hi_topjets_met->fill(event);
     hi_jets_met->fill(event);
     
     
     /* MT_lep selection */
     bool pass_mtlep = mtlep_sel->passes(event);
     if(!pass_mtlep) return false;
     mtlep_h->fill(event);
     hi_electrons_mtlep->fill(event);
     hi_events_mtlep->fill(event);
     hi_muons_mtlep->fill(event);
     hi_topjets_mtlep->fill(event);
     hi_jets_mtlep->fill(event);
     ////
     
     //// DeltaPhi Jet MET Selection
     bool pass_jetmetdphi = jetmetdphi_sel->passes(event);
     if(!pass_jetmetdphi) return false;
     jetmetdphi_h->fill(event);
     hi_electrons_jetmetdphi->fill(event);
     hi_events_jetmetdphi->fill(event);
     hi_muons_jetmetdphi->fill(event);
     hi_topjets_jetmetdphi->fill(event);
     hi_jets_jetmetdphi->fill(event);
     
     //// MT2W Selection
     bool pass_mt2w = mt2w_sel->passes(event);
     if(!pass_mt2w) return false;
     mt2w_h->fill(event);
     hi_electrons_mt2w->fill(event);
     hi_events_mt2w->fill(event);
     hi_muons_mt2w->fill(event);
     hi_topjets_mt2w->fill(event);
     hi_jets_mt2w->fill(event);
     
     bool btag(btagAK4_sel->passes(event));
     
     if(!toptag){
        
        if(!btag){
           
           hi_t0b0->fill(event);
           hi_electrons_t0b0->fill(event);
           hi_events_t0b0->fill(event);
           hi_muons_t0b0->fill(event);
           hi_topjets_t0b0->fill(event);
           hi_jets_t0b0->fill(event);
           //hi_t0b0__hyp->fill(event);
        }
        else {
           
           hi_t0b1->fill(event);
           hi_electrons_t0b1->fill(event);
           hi_events_t0b1->fill(event);
           hi_muons_t0b1->fill(event);
           hi_topjets_t0b1->fill(event);
           hi_jets_t0b1->fill(event);
           //hi_t0b1__hyp->fill(event);
        }
     }
     else {
        
        hi_t1->fill(event);
        hi_electrons_t1->fill(event);
        hi_events_t1->fill(event);
        hi_muons_t1->fill(event);
        hi_topjets_t1->fill(event);
        hi_jets_t1->fill(event);
        //   //hi_t1__hyp->fill(event);
     }
     
     return false;
}

UHH2_REGISTER_ANALYSIS_MODULE(ttDMPostSelectionModule)
