#include "UHH2/TTbarDM/include/ttDMGenHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"
#include "Math/LorentzVector.h"
#include "TFile.h"
#include "TVector2.h"
using namespace uhh2;


ttDMGenHists::ttDMGenHists(Context & ctx, const std::string & dirname): Hists(ctx, dirname){

   h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

   // tops
   top_pt = book<TH1F>("top_pt", "p_{T} top [GeV], semi-lep. top decays", 50, 0, 1000);
   top_eta = book<TH1F>("top_eta", "#eta top, semi-lep. top decays", 50, -5, 5);
   top_phi = book<TH1F>("top_phi", "#phi top, semi-lep. top decays", 50, -4, 4);
   antitop_pt =  book<TH1F>("antitop_pt", "p_{T} antitop [GeV], semi-lep. top decays", 50, 0, 1000);
   antitop_eta = book<TH1F>("antitop_eta", "#eta antitop, semi-lep. top decays", 50, -5, 5);
   antitop_phi = book<TH1F>("antitop_phi", "#phi antitop, semi-lep. top decays", 50, -4, 4);
   bothtops_pt = book<TH1F>("bothtops_pt ", "p_{T} tops [GeV], semi-lep. top decays", 50, 0, 1000);
   bothtops_eta = book<TH1F>("bothtops_eta", "#eta tops, semi-lep. top decays", 50, -5, 5);
   bothtops_phi = book<TH1F>("bothtops_phi", "#phi tops, semi-lep. top decays", 50, -4, 4);
   
   // neutrinos
   neutrino_pt = book<TH1F>("neutrino_pt", "p_{T,#nu} [GeV], semi-lep. top decays", 100, 0, 1000);
   neutrino_eta = book<TH1F>("neutrino_eta", "#eta #nu, semi-lep. top decays", 50, -5, 5);
   neutrino_phi = book<TH1F>("neutrino_phi", "#phi #nu, semi-lep. top decays", 50, -4, 4);
   neutrino_px = book<TH1F>("neutrino_px", "p_{x,#nu} [GeV], semi-lep. top decays", 100, 0, 1000);
   neutrino_py = book<TH1F>("neutrino_py", "p_{y,#nu} [GeV], semi-lep. top decays", 100, 0, 1000);
   neutrino_pz = book<TH1F>("neutrino_pz", "p_{z,#nu} [GeV], semi-lep. top decays", 100, 0, 1000);
   
   // DeltaR tops, neutrinos
   DeltaR_bothtops_neutrino = book<TH1F>("DeltaR_bothtops_neutrino", "#DeltaR(tops, #nu)", 50, 0, 5);
   DeltaR_top_neutrino = book<TH1F>("DeltaR_top_neutrino", "#DeltaR(leptonic top, #nu), semi-lep. top decays", 50, 0, 5);
   DeltaR_antitop_neutrino = book<TH1F>("DeltaR_antitop_neutrino", "#DeltaR(leptonic antitop, #nu), semi-lep. top decays", 50, 0, 5);
   
   // DM particles 
   DM_pt = book<TH1F>("DM_pt", "p_{T,DM} [GeV]", 50, 0, 1000);
   DM_eta = book<TH1F>("DM_eta", "#eta DM", 50, -5, 5);
   DM_phi = book<TH1F>("DM_phi", "#phi DM", 50, -4, 4);

   // DeltaR top, DM particles       
   DeltaR_top_DM = book<TH1F>("DeltaR_top_DM", "#DeltaR(top, DM), semi-lep. top decays", 50, 0, 5);
   DeltaR_antitop_DM = book<TH1F>("DeltaR_antitop_DM", "#DeltaR(antitop, DM), semi-lep. top decays", 50, 0, 5);
   DeltaR_neutrino_DM = book<TH1F>("DeltaR_neutrino_DM", "#DeltaR(#nu, DM), semi-lep. top decays", 50, 0, 5);
      
   // invariant masses 
   InvMass_b_lep = book<TH1F>("InvMass_b_lep","invariant mass (b,l)",100,0,300);
   InvMass_neutrino_lep = book<TH1F>("InvMass_neutrino_lep","invariant mass (#nu,l)",60,50,110);
   InvMass_W_closestJet = book<TH1F>("InvMass_W_closestJet","invariant mass (W,clostest jet)",300,0,300);
   
   // leptons
   lepton_pt= book<TH1F>("lepton_pt","lepton p_{T} [GeV]",100,0,1000);
   lepton_eta= book<TH1F>("lepton_eta","lepton #eta",50, -5, 5);
   lepton_phi= book<TH1F>("lepton_phi","lepton #phi",50, -4, 4);
   DeltaR_lepton_genb= book<TH1F>("DeltaR_lepton_genb","#DeltaR(l, b^{gen})", 50, 0, 5); 
   
   // leptons and neutrinos
   DeltaR_lepton_neutrino= book<TH1F>("DeltaR_lepton_neutrino","#DeltaR(l,#nu)", 50, 0, 5);
   Deltaeta_lepton_neutrino= book<TH1F>("Deltaeta_lepton_neutrino","#Delta#eta(l,#nu)", 50, -2.5, 2.5);
   Deltaeta_lepton_neutrino_2= book<TH1F>("Deltaeta_lepton_neutrino_2","#Delta#eta(l,#nu)", 200, 0, 5);
   Deltaphi_lepton_neutrino= book<TH1F>("Deltaphi_lepton_neutrino","#Delta#phi(l,#nu)", 50, -2.5, 2.5);
   Deltaphi_lepton_neutrino_2= book<TH1F>("Deltaphi_lepton_neutrino_2","#Delta#phi(l,#nu)", 200, 0, 5);
   
   DeltaPx_lepton_neutrino = book<TH1F>("DeltaPx_lepton_neutrino","#Deltap_{x}(l,#nu) [GeV]", 50, 0, 5);
   DeltaPy_lepton_neutrino = book<TH1F>("DeltaPy_lepton_neutrino","#Deltap_{y}(l,#nu) [GeV]", 50, 0, 5);
   DeltaPz_lepton_neutrino= book<TH1F>("DeltaPz_lepton_neutrino","#Deltap_{z}(l,#nu) [GeV]", 50, 0, 5);
   DeltaPt_lepton_neutrino= book<TH1F>("DeltaPt_lepton_neutrino","#Deltap_{T}(l,#nu) [GeV]", 50, 0, 5);
   DeltaE_lepton_neutrino= book<TH1F>("DeltaE_lepton_neutrino","#DeltaE(l,#nu)", 100, -1000, 1000);
   Deltax_lepton_neutrino = book<TH1F>("Deltax_lepton_neutrino","#Deltax(l,#nu)", 50, 0, 5);
   Deltay_lepton_neutrino = book<TH1F>("Deltay_lepton_neutrino","#Deltay(l,#nu)", 50, 0, 5);
   Deltaz_lepton_neutrino= book<TH1F>("Deltaz_lepton_neutrino","#Deltaz(l,#nu)", 50, 0, 5);
   DeltaPx_lepton_neutrino_2 = book<TH1F>("DeltaPx_lepton_neutrino_2","#Deltap_{x}(l,#nu) [GeV]", 100, -1000, 1000);
   DeltaPy_lepton_neutrino_2 = book<TH1F>("DeltaPy_lepton_neutrino_2","#Deltap_{y}(l,#nu) [GeV]", 100, -1000, 1000);
   DeltaPz_lepton_neutrino_2= book<TH1F>("DeltaPz_lepton_neutrino_2","#Deltap_{z}(l,#nu) [GeV]", 100, -1000, 1000);
   DeltaPt_lepton_neutrino_2= book<TH1F>("DeltaPt_lepton_neutrino_2","#Deltap_{T}(l,#nu) [GeV]", 100, -1000, 1000);
   Deltax_lepton_neutrino_2 = book<TH1F>("Deltax_lepton_neutrino_2","#Deltax(l,#nu)", 100, -1000, 1000);
   Deltay_lepton_neutrino_2 = book<TH1F>("Deltay_lepton_neutrino_2","#Deltay(l,#nu)",100, -1000, 1000); 
   Deltaz_lepton_neutrino_2= book<TH1F>("Deltaz_lepton_neutrino_2","#Deltaz(l,#nu)", 100, -1000, 1000);


   // likelihood plots, only include if you want to re-do the likelihood
   // h_likelihood = book<TH3F>("likelihood","likelihood",  200, -5, 5, 200, -5, 5, 300,0,300); //eta, phi, mw
   // double binsx2[] = {-5, -3.5, -2, -1.75,  -1.5, -1.25, -1, -0.75, -0.5, -0.25,  0,  0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5};
   // double binsy2[] = {- (TMath::Pi() + (TMath::Pi() -2)),-TMath::Pi(), -2, -1.75,  -1.5, -1.25, -1, -0.75, -0.5, -0.25,  0,  0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, TMath::Pi(),  TMath::Pi() + (TMath::Pi() -2)};
   // double binsz2[] = {0, 20, 40,  60, 65, 70, 75,80, 85, 90, 95, 100, 105, 110, 130, 150, 170, 190, 210,230,  250, 270, 300};
   // h_likelihood_bin2 = book<TH3F>("likelihood_bin2","likelihood_bin2",  20, binsx2, 20, binsy2, 22,binsz2);
   // h_likelihood_bin2_noweights = book<TH3F>("likelihood_bin2_noweights","likelihood_bin2_noweights",  20, binsx2, 20, binsy2, 22,binsz2);
   
   // correlations dPhi, dEta, MW
   h_dphi_deta = book<TH2F>("dphi_deta","#Delta#phi, #Delta#eta", 100, -TMath::Pi(), TMath::Pi(), 100, -5,5);
   h_dW_deta = book<TH2F>("dW_deta","M_{W}, #Delta#eta", 100, 0, 300, 100,-5,5);
   h_dW_dphi = book<TH2F>("dW_dphi","M_{W}, #Delta#phi", 100, 0, 300, 100,-TMath::Pi(),TMath::Pi());

   // met plots
   Resolution_MET = book<TH1F>("Resolution_MET","(met^{rec} - met^{gen}) /( met^{gen})", 50, -5, 5);
   Gen_MET = book<TH1F>("Gen_MET","met^{gen}", 100, 0, 3000);

   //masses
   DM_mass = book<TH1F>("DM_mass","m_{DM}", 1000, 0, 5000);
   DMparticles_mass = book<TH1F>("DMparticles_mass","m_{#chi}", 1000, 0, 5000);
   mediator_mass =   book<TH1F>("mediator_mass","m_{#phi}", 1000, 0, 5000);

   MET_Neutrino_px = book<TH1F>("MET_Neutrino_px","MET_{x}-p_{x}^{nu}", 200, -200, 200);
   MET_Neutrino_py = book<TH1F>("MET_Neutrino_py","MET_{y}-p_{y}^{nu}", 200, -200, 200);
}

void ttDMGenHists::fill(const Event & event){
   
   
   weight = event.weight;
   
   //do not fill histograms if ttbargen information has not been filled
    if(!event.is_valid(h_ttbargen)){
       return;
    }
    const auto & ttbargen = event.get(h_ttbargen);
    
    if(event.isRealData) return;
    assert(event.genparticles);

    //find DM particle
    LorentzVector DM;
    for(unsigned int i=0; i<event.genparticles->size(); ++i) {
       const GenParticle & genp = event.genparticles->at(i);
       if (abs(genp.pdgId())==9100000) 
          {
             DM = genp.v4(); //mediator is stored
             mediator_mass->Fill(DM.M(), weight);
          }
       if (genp.pdgId() == 9100022) { //DM particles are stored
          for(unsigned int j=0; j<event.genparticles->size(); ++j)
             {
                const GenParticle & genp2 = event.genparticles->at(j);
                if (genp2.pdgId() == -9100022) DM = genp.v4() + genp2.v4();
                DMparticles_mass ->Fill(DM.M(), weight);
             }
       }
    }
    DM_mass->Fill(DM.M(), weight);
          
    if (ttbargen.IsSemiLeptonicDecay())
       {
          top_pt->Fill(ttbargen.Top().pt(),weight);
          top_eta->Fill(ttbargen.Top().eta(),weight);
          top_phi->Fill(ttbargen.Top().phi(),weight);
          
          antitop_pt->Fill(ttbargen.Antitop().pt(),weight);
          antitop_eta->Fill(ttbargen.Antitop().eta(),weight);
          antitop_phi->Fill(ttbargen.Antitop().phi(),weight);
          
          bothtops_pt->Fill(ttbargen.Top().pt(),weight);
          bothtops_pt->Fill(ttbargen.Antitop().pt(),weight);
          bothtops_eta->Fill(ttbargen.Top().eta(),weight);
          bothtops_eta->Fill(ttbargen.Antitop().eta(),weight);
          bothtops_phi->Fill(ttbargen.Top().phi(),weight);
          bothtops_phi->Fill(ttbargen.Antitop().phi(),weight);
                             
          neutrino_pt->Fill(ttbargen.Neutrino().pt(),weight);
          neutrino_eta->Fill(ttbargen.Neutrino().eta(),weight);
          neutrino_phi->Fill(ttbargen.Neutrino().phi(),weight);
          neutrino_px->Fill(ttbargen.Neutrino().v4().Px(),weight);
          neutrino_py->Fill(ttbargen.Neutrino().v4().Py(),weight);
          neutrino_pz->Fill(ttbargen.Neutrino().v4().Pz(),weight);
          
          DeltaR_bothtops_neutrino->Fill(uhh2::deltaR(ttbargen.Top(),ttbargen.Neutrino()),weight);
          DeltaR_bothtops_neutrino->Fill(uhh2::deltaR(ttbargen.Antitop(),ttbargen.Neutrino()),weight);
          if (ttbargen.IsTopHadronicDecay()) DeltaR_antitop_neutrino->Fill(uhh2::deltaR(ttbargen.Antitop(),ttbargen.Neutrino()),weight);
          if (ttbargen.IsAntiTopHadronicDecay()) DeltaR_top_neutrino->Fill(uhh2::deltaR(ttbargen.Top(),ttbargen.Neutrino()),weight);
  
          DeltaR_top_DM ->Fill(uhh2::deltaR(ttbargen.Top(),DM),weight);
          DeltaR_antitop_DM ->Fill(uhh2::deltaR(ttbargen.Antitop(),DM),weight);
          DeltaR_neutrino_DM ->Fill(uhh2::deltaR(ttbargen.Neutrino(),DM),weight);
           
          //h_likelihood ->Fill(DeltaEta, DeltaPhi,MW, weight);
          // h_likelihood_bin2->Fill(DeltaEta, DeltaPhi,MW, weight);
          // if (DeltaPhi > 2 && DeltaPhi< TMath::Pi()) h_likelihood_bin2->Fill(DeltaEta, DeltaPhi-2*TMath::Pi() ,MW, weight);
          // if (DeltaPhi < -2 && DeltaPhi > -TMath::Pi()) h_likelihood_bin2->Fill(DeltaEta, DeltaPhi+2*TMath::Pi() ,MW ,weight);
          // h_likelihood_bin2_noweights->Fill(DeltaEta, DeltaPhi,MW, 1.0);
          // if (DeltaPhi > 2 && DeltaPhi< TMath::Pi()) h_likelihood_bin2_noweights->Fill(DeltaEta, DeltaPhi-2*TMath::Pi() ,MW, 1.0);
          // if (DeltaPhi < -2 && DeltaPhi > -TMath::Pi()) h_likelihood_bin2_noweights->Fill(DeltaEta, DeltaPhi+2*TMath::Pi() ,MW ,1.0);

          if (abs(ttbargen.ChargedLepton().pdgId())==11 || abs(ttbargen.ChargedLepton().pdgId())==13)
             {
                lepton_pt-> Fill(ttbargen.ChargedLepton().pt(),weight);
                lepton_eta->Fill(ttbargen.ChargedLepton().eta(), weight);
                lepton_phi->Fill(ttbargen.ChargedLepton().phi(),weight);

                InvMass_b_lep->Fill(inv_mass_safe((ttbargen.ChargedLepton().v4() + ttbargen.BLep().v4()) ),weight);
                InvMass_neutrino_lep->Fill(inv_mass_safe((ttbargen.ChargedLepton().v4() + ttbargen.Neutrino().v4()) ),weight);
                Jet bjet = *nextJet(ttbargen.ChargedLepton(),*event.jets);
                InvMass_W_closestJet ->Fill(inv_mass_safe(ttbargen.ChargedLepton().v4() + ttbargen.Neutrino().v4() + bjet.v4()) ,weight);
                
                DeltaR_lepton_neutrino->Fill(uhh2::deltaR(ttbargen.Neutrino(),ttbargen.ChargedLepton()),weight);
                Deltaeta_lepton_neutrino->Fill(ttbargen.Neutrino().eta()-ttbargen.ChargedLepton().eta(),weight);
                Deltaeta_lepton_neutrino_2->Fill(fabs(ttbargen.Neutrino().eta()-ttbargen.ChargedLepton().eta()),weight);
                Deltaphi_lepton_neutrino->Fill(ttbargen.Neutrino().phi()-ttbargen.ChargedLepton().phi(),weight);
                Deltaphi_lepton_neutrino_2->Fill(fabs(ttbargen.Neutrino().phi()-ttbargen.ChargedLepton().phi()),weight);
                DeltaPx_lepton_neutrino ->Fill(fabs(ttbargen.Neutrino().v4().Px()-ttbargen.ChargedLepton().v4().Px())/ttbargen.Neutrino().v4().Px(),weight);
                DeltaPy_lepton_neutrino ->Fill(fabs(ttbargen.Neutrino().v4().Py()-ttbargen.ChargedLepton().v4().Py())/ttbargen.Neutrino().v4().Py(),weight);
                DeltaPz_lepton_neutrino ->Fill(fabs(ttbargen.Neutrino().v4().Pz()-ttbargen.ChargedLepton().v4().Pz())/ttbargen.Neutrino().v4().Pz(),weight);
                DeltaPt_lepton_neutrino ->Fill(fabs(ttbargen.Neutrino().pt()-ttbargen.ChargedLepton().pt())/ttbargen.Neutrino().pt(),weight);

                DeltaPx_lepton_neutrino_2 ->Fill(ttbargen.Neutrino().v4().Px()-ttbargen.ChargedLepton().v4().Px(),weight);
                DeltaPy_lepton_neutrino_2 ->Fill(ttbargen.Neutrino().v4().Py()-ttbargen.ChargedLepton().v4().Py(),weight);
                DeltaPz_lepton_neutrino_2 ->Fill(ttbargen.Neutrino().v4().Pz()-ttbargen.ChargedLepton().v4().Pz(),weight);
                DeltaPt_lepton_neutrino_2 ->Fill(ttbargen.Neutrino().pt()-ttbargen.ChargedLepton().pt(),weight);
                DeltaE_lepton_neutrino ->Fill(ttbargen.Neutrino().energy()-ttbargen.ChargedLepton().energy(),weight);
                Deltax_lepton_neutrino_2 ->Fill(ttbargen.Neutrino().v4().X()-ttbargen.ChargedLepton().v4().X(),weight);
                Deltay_lepton_neutrino_2 ->Fill(ttbargen.Neutrino().v4().Y()-ttbargen.ChargedLepton().v4().Y(),weight);
                Deltaz_lepton_neutrino_2 ->Fill(ttbargen.Neutrino().v4().Z()-ttbargen.ChargedLepton().v4().Z(),weight);
                
                Deltax_lepton_neutrino ->Fill(fabs(ttbargen.Neutrino().v4().X()-ttbargen.ChargedLepton().v4().X())/ttbargen.Neutrino().v4().X(),weight);
                Deltay_lepton_neutrino ->Fill(fabs(ttbargen.Neutrino().v4().Y()-ttbargen.ChargedLepton().v4().Y())/ttbargen.Neutrino().v4().Y(),weight);
                Deltaz_lepton_neutrino ->Fill(fabs(ttbargen.Neutrino().v4().Z()-ttbargen.ChargedLepton().v4().Z())/ttbargen.Neutrino().v4().Z(),weight);
                MET_Neutrino_px ->Fill(event.met->v4().Px()-ttbargen.Neutrino().v4().Px());
                MET_Neutrino_py ->Fill(event.met->v4().Py()-ttbargen.Neutrino().v4().Py());
             }
          
          LorentzVector lep = ttbargen.ChargedLepton().v4(); //pt,eta,phi,e
          LorentzVector nu= ttbargen.Neutrino().v4();   //pt,eta,phi,e
          Double_t MW = (nu+lep).M();
          Double_t DeltaPhi = nu.Phi()-lep.Phi();
          if(!std::isnan(DeltaPhi))DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi);
          Double_t DeltaEta = nu.Eta()-lep.Eta();
          
          h_dphi_deta ->Fill(DeltaPhi,DeltaEta, weight);
          h_dW_deta ->Fill(MW, DeltaEta, weight);
          h_dW_dphi ->Fill(MW, DeltaPhi, weight);

          Particle lepton = ttbargen.ChargedLepton();
          GenParticle genb = ttbargen.BLep();
          double DeltaR = deltaR(genb, lepton);
          DeltaR_lepton_genb->Fill(DeltaR, weight);

       }
    
    LorentzVector met_gen;
    for(unsigned int i=0; i<event.genparticles->size(); ++i) {
       const GenParticle & genp = event.genparticles->at(i);
       if (abs(genp.pdgId())==12 || abs(genp.pdgId())==14 || abs(genp.pdgId())==16) met_gen = met_gen+genp.v4();
    }     
    met_gen = met_gen+DM;
    Gen_MET -> Fill(met_gen.Pt(),weight);
    Resolution_MET -> Fill((event.met->v4().Pt() - met_gen.Pt())/ met_gen.Pt(),weight);
    
    DM_pt->Fill(DM.pt(),weight);
    DM_phi->Fill(DM.phi(),weight);
    DM_eta->Fill(DM.eta(),weight);
    


    return;
    
}

