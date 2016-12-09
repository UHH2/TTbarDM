
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicSelections.h"
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicUtils.h"

#include <iostream>
#include <memory>

#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/common/include/Utils.h"

uhh2::HTlepCut::HTlepCut(float min_htlep, float max_htlep):
  min_htlep_(min_htlep), max_htlep_(max_htlep) {}

bool uhh2::HTlepCut::passes(const uhh2::Event & event){

  assert(event.muons || event.electrons);
  assert(event.met);

  float plep_pt(0.);

  if(event.electrons){
    for(const auto & ele : *event.electrons){
      if(ele.pt() > plep_pt) plep_pt = ele.pt();
    }
  }

  if(event.muons) {
    for(const auto & mu : *event.muons){
      if(mu.pt() > plep_pt) plep_pt = mu.pt();
    }
  }

  float htlep = plep_pt + event.met->pt();
  return (htlep > min_htlep_) && (htlep < max_htlep_);
}
////////////////////////////////////////////////////////

uhh2::MTlepCut::MTlepCut(float min_mtlep, float max_mtlep):
  min_mtlep_(min_mtlep), max_mtlep_(max_mtlep) {}

bool uhh2::MTlepCut::passes(const uhh2::Event & event){

  assert(event.muons || event.electrons);
  assert(event.met);

  float plep_pt(0.);
  LorentzVector lep;

  if(event.electrons){
    for(const auto & ele : *event.electrons){
      if(ele.pt() > plep_pt) lep = ele.v4();
    }
  }

  if(event.muons) {
    for(const auto & mu : *event.muons){
      if(mu.pt() > plep_pt) lep = mu.v4();
    }
  }

  double deltaphi = uhh2::deltaPhi(*event.met, lep);
  double mtlep = sqrt(2*event.met->pt()*lep.pt()*(1-cos(deltaphi)));

  return (mtlep > min_mtlep_) && (mtlep < max_mtlep_);
}
////////////////////////////////////////////////////////

uhh2::METCut::METCut(float min_met, float max_met):
  min_met_(min_met), max_met_(max_met) {}

bool uhh2::METCut::passes(const uhh2::Event & event){

  assert(event.met);

  float MET = event.met->pt();
  return (MET > min_met_) && (MET < max_met_);
}
////////////////////////////////////////////////////////

uhh2::NJetCut::NJetCut(int nmin_, int nmax_, float ptmin_, float etamax_):
  nmin(nmin_), nmax(nmax_), ptmin(ptmin_), etamax(etamax_) {}

bool uhh2::NJetCut::passes(const uhh2::Event & event){

  int njet(0);
  for(auto & jet : *event.jets){
    if(jet.pt() > ptmin && fabs(jet.eta()) < etamax) ++njet;
  }

  return (njet >= nmin) && (njet <= nmax);
}
////////////////////////////////////////////////////////

uhh2::METJetDPhiCut::METJetDPhiCut(float min_dphi, int maxjetindex):
  min_dphi_(min_dphi), maxjetindex_(maxjetindex) {}

bool uhh2::METJetDPhiCut::passes(const uhh2::Event & event){
  assert(event.met);
  assert(event.jets->size() >= maxjetindex_);
  //assert(event.jets);
  //if (event.jets->size() < maxjetindex_) maxjetindex_ =event.jets->size();
  double deltaphi=0;
  double mindeltaphi = 9999;
  for(size_t i=0; i<maxjetindex_; i++){
     deltaphi = uhh2::deltaPhi(*event.met, event.jets->at(i));
     if (deltaphi < mindeltaphi) mindeltaphi=deltaphi;
  }
  return (mindeltaphi > min_dphi_);
}
////////////////////////////////////////////////////////

uhh2::METLeptonDPhiCut::METLeptonDPhiCut(float min_dphi):
   min_dphi_(min_dphi) {}

bool uhh2::METLeptonDPhiCut::passes(const uhh2::Event & event){
  
   assert(event.met);
   assert(event.electrons);
   assert(event.muons);
   
   double deltaphi;
   if (event.muons->size()>0) deltaphi = uhh2::deltaPhi(*event.met, event.muons->at(0));
   else deltaphi = uhh2::deltaPhi(*event.met, event.electrons->at(0));
   
   return (deltaphi > min_dphi_);
}
////////////////////////////////////////////////////////

uhh2::MT2WCut::MT2WCut(float min_mt2w):
  min_mt2w_(min_mt2w) {}

bool uhh2::MT2WCut::passes(const uhh2::Event & event){

  assert(event.met);
  assert(event.jets->size() > 1);

  double mt2w = CalculateMT2W(event);

  return (mt2w > min_mt2w_);
}
////////////////////////////////////////////////////////

bool uhh2::TwoDCut::passes(const uhh2::Event & event){

   //assert(event.muons && event.electrons && event.jets);
  assert(event.jets);
  std::vector<Muon> muons_med=event.get(h_muons);
  std::vector<Electron> electrons=event.get(h_electrons);
  // std::cout<<"electron size: "<<electrons.size()<<std::endl;
  // if (electrons.size()>0) std::cout<<"pt: "<<electrons.at(0).pt()<<std::endl;
  // if (electrons.size()>0)std::cout<<"eta: "<<electrons.at(0).eta()<<std::endl;
  // if (electrons.size()>0)std::cout<<ElectronID_Spring15_25ns_tight_noIso(electrons.at(0),event)<<std::endl;
  // std::cout<<"Muon size: im 2D cut "<<muons_med.size()<<std::endl;
  // for (unsigned int i=0; i<muons_med.size(); i++) 
  //    {
  //       if (muons_med.size()>0) std::cout<<"pt: "<<muons_med.at(i).pt()<<std::endl;
  //       if (muons_med.size()>0)std::cout<<"eta: "<<muons_med.at(i).eta()<<std::endl;
  //       std::cout<<"medium? "<< MuonIDMedium()(muons_med.at(i),event)<<std::endl;
  //       if (!MuonIDMedium()(muons_med.at(i),event)) std::cout<<"NOT MEDIUM MUON!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  //    }
  if((muons_med.size()+electrons.size()) != 1){
     std::cout << "\n @@@ WARNING -- TwoDCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
     
     return false;
  }

  float drmin=0;
  float ptrel=0;
  if(muons_med.size()) std::tie(drmin, ptrel) = drmin_pTrel(muons_med.at(0), *event.jets);
  else std::tie(drmin, ptrel) = drmin_pTrel(electrons.at(0), *event.jets);

  return (drmin > min_deltaR_) || (ptrel > min_pTrel_);
}
////////////////////////////////////////////////////////

uhh2::TriangularCuts::TriangularCuts(float a, float b): a_(a), b_(b) {

  if(!b_) std::runtime_error("TriangularCuts -- incorrect initialization (parameter 'b' is null)");
}

bool uhh2::TriangularCuts::passes(const uhh2::Event & event){

  assert(event.muons || event.electrons);
  assert(event.jets && event.met);

  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- TriangularCuts::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  if(!event.jets->size()){
    std::cout << "\n @@@ WARNING -- TriangularCuts::passes -- unexpected number of jets in the event (==0). returning 'false'\n";
    return false;
  }

  // charged lepton
  const Particle* lep1(0);
  if(event.muons->size()) lep1 = &event.muons->at(0);
  else lep1 = &event.electrons->at(0);

  // 1st entry in jet collection (should be the pt-leading jet)
  const Particle* jet1 = &event.jets->at(0);

  // MET-lepton triangular cut
  bool pass_tc_lep = fabs(fabs(deltaPhi(*event.met, *lep1)) - a_) < a_/b_ * event.met->pt();

  // MET-jet triangular cut
  bool pass_tc_jet = fabs(fabs(deltaPhi(*event.met, *jet1)) - a_) < a_/b_ * event.met->pt();

  return pass_tc_lep && pass_tc_jet;
}
////////////////////////////////////////////////////////

uhh2::TopTagEventSelection::TopTagEventSelection(const TopJetId& tjetID, float minDR_jet_ttag):
  topjetID_(tjetID), minDR_jet_toptag_(minDR_jet_ttag) {

  topjet1_sel_.reset(new NTopJetSelection(1, -1, topjetID_));
}

bool uhh2::TopTagEventSelection::passes(const uhh2::Event & event){

  if(!topjet1_sel_->passes(event)) return false;

  for(auto & topjet : * event.topjets){
    if(!topjetID_(topjet, event)) continue;

    for(auto & jet : * event.jets)
      if(deltaR(jet, topjet) > minDR_jet_toptag_) return true;
  }

  return false;
}
////////////////////////////////////////////////////////
uhh2::Type2TopTagEventSelection::Type2TopTagEventSelection(const TopJetId& tjetID, float minDR_jet_ttag, float maxDR_jet_ttag, float mininvmass_jet_ttag, float maxinvmass_jet_ttag):
   topjetID_(tjetID), minDR_jet_toptag_(minDR_jet_ttag), maxDR_jet_toptag_(maxDR_jet_ttag), mininvmass_jet_ttag_(mininvmass_jet_ttag), maxinvmass_jet_ttag_(maxinvmass_jet_ttag){
   topjet1_sel_.reset(new NTopJetSelection(1, -1, topjetID_));
}

bool uhh2::Type2TopTagEventSelection::passes(const uhh2::Event & event){

  if(!topjet1_sel_->passes(event)) return false;

  for(auto & topjet : * event.topjets){
    if(!topjetID_(topjet, event)) continue;

    for(auto & jet : * event.jets)
       if(deltaR(jet, topjet) > minDR_jet_toptag_ && deltaR(jet, topjet) < maxDR_jet_toptag_ && (topjet.v4()+jet.v4()).M() > mininvmass_jet_ttag_ && (topjet.v4()+jet.v4()).M() < maxinvmass_jet_ttag_) 
          {
             return true;
          }
  }
  return false;
}
////////////////////////////////////////////////////////

uhh2::LeptonicTopPtCut::LeptonicTopPtCut(uhh2::Context& ctx, float pt_min, float pt_max, const std::string& hyps_name, const std::string& disc_name):
  tlep_pt_min_(pt_min), tlep_pt_max_(pt_max), h_hyps_(ctx.get_handle<std::vector<ReconstructionHypothesis>>(hyps_name)), disc_name_(disc_name) {}

bool uhh2::LeptonicTopPtCut::passes(const uhh2::Event& event){

  std::vector<ReconstructionHypothesis> hyps = event.get(h_hyps_);
  const ReconstructionHypothesis* hyp = get_best_hypothesis(hyps, disc_name_);
  if(!hyp) std::runtime_error("LeptonicTopPtCut -- best hypothesis not found (discriminator="+disc_name_+")");

  float tlep_pt = hyp->toplep_v4().Pt();

  return (tlep_pt > tlep_pt_min_) && (tlep_pt < tlep_pt_max_);
}
////////////////////////////////////////////////////////

uhh2::HypothesisDiscriminatorCut::HypothesisDiscriminatorCut(uhh2::Context& ctx, float disc_min, float disc_max, const std::string& hyps_name, const std::string& disc_name):
  disc_min_(disc_min), disc_max_(disc_max), h_hyps_(ctx.get_handle<std::vector<ReconstructionHypothesis>>(hyps_name)), disc_name_(disc_name) {}

bool uhh2::HypothesisDiscriminatorCut::passes(const uhh2::Event & event){

  std::vector<ReconstructionHypothesis> hyps = event.get(h_hyps_);
  const ReconstructionHypothesis* hyp = get_best_hypothesis(hyps, disc_name_);
  if(!hyp) std::runtime_error("HypothesisDiscriminatorCut -- best hypothesis not found (discriminator="+disc_name_+")");

  float disc_val = hyp->discriminator(disc_name_);

  return (disc_val > disc_min_) && (disc_val < disc_max_);
}
////////////////////////////////////////////////////////

uhh2::LikelihoodSelection::LikelihoodSelection(uhh2::Context& ctx, float lmax ):
   lmax_(lmax), h_likelihood_(ctx.get_handle<double>("likelihood")) {}

bool uhh2::LikelihoodSelection::passes(const uhh2::Event & event){
   double likelihood = event.get(h_likelihood_);
   return ((likelihood < lmax_) && (likelihood > 0)) ;
}


uhh2::DeltaPhiMetNeutrino::DeltaPhiMetNeutrino(uhh2::Context& ctx, float deltaphimin ):
   deltaphimin_(deltaphimin), h_neutrino_(ctx.get_handle<LorentzVector>("rec_neutrino")) {}

bool uhh2::DeltaPhiMetNeutrino::passes(const uhh2::Event & event){
   LorentzVector neutrino = event.get(h_neutrino_);
   double deltaphi_neutrino_met = uhh2::deltaPhi(neutrino,*event.met);

   return (deltaphi_neutrino_met >deltaphimin_); 
}

uhh2::DeltaPhiTaggedJetNeutrino::DeltaPhiTaggedJetNeutrino(uhh2::Context& ctx, float deltaphimax ):
   deltaphimax_(deltaphimax), h_neutrino_(ctx.get_handle<LorentzVector>("rec_neutrino")), h_taggedjet_(ctx.get_handle<std::vector<TopJet>>("h_heptopjets_WP3_wobtag")) {}

bool uhh2::DeltaPhiTaggedJetNeutrino::passes(const uhh2::Event & event){
   LorentzVector neutrino = event.get(h_neutrino_);
   TopJet tj = event.get(h_taggedjet_).at(0);
   double deltaphi_neutrino_tj = uhh2::deltaPhi(neutrino,tj);

   return (deltaphi_neutrino_tj <deltaphimax_); 
}

uhh2::DeltaPhiTaggedJetTopLep::DeltaPhiTaggedJetTopLep(uhh2::Context& ctx, float  deltaphimax):
   deltaphimax_(deltaphimax), h_neutrino_(ctx.get_handle<LorentzVector>("rec_neutrino")), h_taggedjet_(ctx.get_handle<std::vector<TopJet>>("h_heptopjets_WP3_wobtag")), h_b_jets_(ctx.get_handle<Jet>("bjet")) {}

bool uhh2::DeltaPhiTaggedJetTopLep::passes(const uhh2::Event & event){
   LorentzVector neutrino = event.get(h_neutrino_);
   TopJet tj = event.get(h_taggedjet_).at(0);
   Jet b_jet = event.get(h_b_jets_);
   Muon muon = event.muons->at(0);

   LorentzVector toplep = muon.v4()+neutrino+b_jet.v4();
   double deltaphi_toplep_tj = uhh2::deltaPhi(toplep,tj);

   return (deltaphi_toplep_tj <deltaphimax_);
}

uhh2::NeutrinopTSelection::NeutrinopTSelection(uhh2::Context& ctx, float pTmin ):
   pTmin_(pTmin), h_neutrino_(ctx.get_handle<LorentzVector>("rec_neutrino")) {}

bool uhh2::NeutrinopTSelection::passes(const uhh2::Event & event){
   LorentzVector neutrino = event.get(h_neutrino_);
   return (pTmin_ < neutrino.pt());
}

uhh2::ttbarpTSel::ttbarpTSel(uhh2::Context& ctx, float pTmin ):
   pTmin_(pTmin), h_neutrino_(ctx.get_handle<LorentzVector>("rec_neutrino")), h_taggedjet_(ctx.get_handle<std::vector<TopJet>>("h_heptopjets_WP3_wobtag")),h_b_jets_(ctx.get_handle<Jet>("bjet"))  {}

bool uhh2::ttbarpTSel::passes(const uhh2::Event & event){
   LorentzVector neutrino = event.get(h_neutrino_);
   TopJet tj = event.get(h_taggedjet_).at(0);
   Muon lepton =event.muons->at(0);
   Jet b_jet = event.get(h_b_jets_);
   
   LorentzVector ttbar = tj.v4()+neutrino+lepton.v4()+b_jet.v4();
   return (ttbar.pt() >pTmin_); 
}

uhh2::DMMETSelection::DMMETSelection(uhh2::Context& ctx, float DMMETmin ):
   DMMETmin_(DMMETmin), h_neutrino_(ctx.get_handle<LorentzVector>("rec_neutrino")){}

bool uhh2::DMMETSelection::passes(const uhh2::Event & event){
   LorentzVector neutrino = event.get(h_neutrino_);
   double DM_MET = std::sqrt((event.met->v4().Px()-neutrino.Px())*(event.met->v4().Px()-neutrino.Px())+(event.met->v4().Py()-neutrino.Py())*(event.met->v4().Py()-neutrino.Py()));

   return (DM_MET > DMMETmin_);
}

uhh2::DeltaRMuonJet::DeltaRMuonJet(double deltaRmax):
   deltaRmax_(deltaRmax){}

bool uhh2::DeltaRMuonJet::passes(const uhh2::Event & event){
   
   Muon lepton =event.muons->at(0);
   
   float minDR_jet(-1.), pTrel_jet(-1.);
   std::tie(minDR_jet, pTrel_jet) = drmin_pTrel(lepton, *event.jets);
   return (minDR_jet< deltaRmax_);
}
