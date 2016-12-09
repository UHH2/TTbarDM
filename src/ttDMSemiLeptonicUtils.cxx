#include "UHH2/TTbarDM/include/ttDMSemiLeptonicUtils.h"
#include "UHH2/TTbarDM/include/Mt2Com_bisect.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/core/include/LorentzVector.h"

bool TopJetLeptonDeltaRCleaner::process(uhh2::Event & event){

  assert(event.topjets);
  std::vector<TopJet> cleaned_topjets;

  for(const auto & tjet : *event.topjets){
    bool skip_tjet(false);

    if(event.muons){
      for(const auto & muo : *event.muons)
        if(uhh2::deltaR(tjet, muo) < minDR_) skip_tjet = true;
    }

    if(skip_tjet) continue;

    if(event.electrons){
      for(const auto & ele : *event.electrons)
        if(uhh2::deltaR(tjet, ele) < minDR_) skip_tjet = true;
    }

    if(!skip_tjet) cleaned_topjets.push_back(tjet);
  }

  event.topjets->clear();
  event.topjets->reserve(cleaned_topjets.size());
  for(auto & j : cleaned_topjets) event.topjets->push_back(j);

  return true;
}

bool HepTTLeptonDeltaRCleaner::process(uhh2::Event & event){

   std::vector<TopJet> heptopjets_WP3;
   if (event.is_valid(h_heptopjets_WP3)) heptopjets_WP3=event.get(h_heptopjets_WP3);
   else return false;
   
   std::vector<TopJet> cleaned_topjets;
   
   for(const auto & tjet : heptopjets_WP3){
    bool skip_tjet(false);

    if(event.muons){
      for(const auto & muo : *event.muons)
        if(uhh2::deltaR(tjet, muo) < minDR_) skip_tjet = true;
    }

    if(skip_tjet) continue;

    if(event.electrons){
      for(const auto & ele : *event.electrons)
        if(uhh2::deltaR(tjet, ele) < minDR_) skip_tjet = true;
    }

    if(!skip_tjet) cleaned_topjets.push_back(tjet);
  }
  heptopjets_WP3.clear();
  heptopjets_WP3.reserve(cleaned_topjets.size());
  for(auto & j : cleaned_topjets) heptopjets_WP3.push_back(j);
  event.set(h_heptopjets_WP3, heptopjets_WP3);
  return true;
}

double CalculateMT2W(const uhh2::Event & event){

  assert(event.met);
  assert(event.jets->size()>=2);

  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> bjets;
  float plep_pt(0.);
  TLorentzVector lep;
  TVector2 met;

  for(const auto & jet : *event.jets){
     if (jet.pt()<30) continue; //changed from 50 to 30
    TLorentzVector tjet;
    tjet.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
    jets.push_back(tjet);
    if (CSVBTag(0.8)(jet, event)) bjets.push_back(tjet); //changed
  }

  if(event.electrons){
    for(const auto & ele : *event.electrons){
      if(ele.pt() > plep_pt) lep.SetPtEtaPhiE(ele.pt(), ele.eta(), ele.phi(), ele.energy());
    }
  }

  if(event.muons) {
    for(const auto & mu : *event.muons){
      if(mu.pt() > plep_pt) lep.SetPtEtaPhiE(mu.pt(), mu.eta(), mu.phi(), mu.energy());
    }
  }

  met.SetMagPhi(event.met->pt(), event.met->phi());
  double mt2w = Mt2Com_bisect().calculateMT2w(jets, bjets, lep, met, std::string("MT2w"));
  return mt2w;

}
