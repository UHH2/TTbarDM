#include "UHH2/TTbarDM/include/ttDMReconstructionHists_ScansAndMarginalisation.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TError.h"
using namespace uhh2;


ttDMReconstructionHists_ScansAndMarginalisation::ttDMReconstructionHists_ScansAndMarginalisation(Context & ctx, const std::string & dirname): Hists(ctx, dirname){

   h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

   h_px_py_likelihood_1 = book<TH2F>("px_py_likelihood2_1","px_py_likelihood2_1",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_1 = book<TH2F>("px_pz_likelihood2_1","px_pz_likelihood2_1",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_1 = book<TH2F>("py_pz_likelihood2_1","py_pz_likelihood2_1",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_2 = book<TH2F>("px_py_likelihood2_2","px_py_likelihood2_2",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_2 = book<TH2F>("px_pz_likelihood2_2","px_pz_likelihood2_2",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_2 = book<TH2F>("py_pz_likelihood2_2","py_pz_likelihood2_2",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_3 = book<TH2F>("px_py_likelihood2_3","px_py_likelihood2_3",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_3 = book<TH2F>("px_pz_likelihood2_3","px_pz_likelihood2_3",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_3 = book<TH2F>("py_pz_likelihood2_3","py_pz_likelihood2_3",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_4 = book<TH2F>("px_py_likelihood2_4","px_py_likelihood2_4",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_4 = book<TH2F>("px_pz_likelihood2_4","px_pz_likelihood2_4",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_4 = book<TH2F>("py_pz_likelihood2_4","py_pz_likelihood2_4",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_5 = book<TH2F>("px_py_likelihood2_5","px_py_likelihood2_5",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_5 = book<TH2F>("px_pz_likelihood2_5","px_pz_likelihood2_5",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_5 = book<TH2F>("py_pz_likelihood2_5","py_pz_likelihood2_5",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_6 = book<TH2F>("px_py_likelihood2_6","px_py_likelihood2_6",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_6 = book<TH2F>("px_pz_likelihood2_6","px_pz_likelihood2_6",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_6 = book<TH2F>("py_pz_likelihood2_6","py_pz_likelihood2_6",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_7 = book<TH2F>("px_py_likelihood2_7","px_py_likelihood2_7",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_7 = book<TH2F>("px_pz_likelihood2_7","px_pz_likelihood2_7",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_7 = book<TH2F>("py_pz_likelihood2_7","py_pz_likelihood2_7",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_8 = book<TH2F>("px_py_likelihood2_8","px_py_likelihood2_8",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_8 = book<TH2F>("px_pz_likelihood2_8","px_pz_likelihood2_8",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_8 = book<TH2F>("py_pz_likelihood2_8","py_pz_likelihood2_8",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_9 = book<TH2F>("px_py_likelihood2_9","px_py_likelihood2_9",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_9 = book<TH2F>("px_pz_likelihood2_9","px_pz_likelihood2_9",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_9 = book<TH2F>("py_pz_likelihood2_9","py_pz_likelihood2_9",420,-1050,1050,420,-1050,1050);
   h_px_py_likelihood_10 = book<TH2F>("px_py_likelihood2_10","px_py_likelihood2_10",420,-1050,1050,420,-1050,1050);
   h_px_pz_likelihood_10 = book<TH2F>("px_pz_likelihood2_10","px_pz_likelihood2_10",420,-1050,1050,420,-1050,1050);
   h_py_pz_likelihood_10= book<TH2F>("py_pz_likelihood2_10","py_pz_likelihood2_10",420,-1050,1050,420,-1050,1050);
 
   neutrino_px_1 = book<TH1F>("neutrino_px_1", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_1 = book<TH1F>("neutrino_py_1", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_1 = book<TH1F>("neutrino_pz_1", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_2 = book<TH1F>("neutrino_px_2", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_2 = book<TH1F>("neutrino_py_2", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_2 = book<TH1F>("neutrino_pz_2", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_3 = book<TH1F>("neutrino_px_3", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_3 = book<TH1F>("neutrino_py_3", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_3 = book<TH1F>("neutrino_pz_3", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_4 = book<TH1F>("neutrino_px_4", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_4 = book<TH1F>("neutrino_py_4", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_4 = book<TH1F>("neutrino_pz_4", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_5 = book<TH1F>("neutrino_px_5", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_5 = book<TH1F>("neutrino_py_5", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_5 = book<TH1F>("neutrino_pz_5", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_6 = book<TH1F>("neutrino_px_6", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_6 = book<TH1F>("neutrino_py_6", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_6 = book<TH1F>("neutrino_pz_6", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_7 = book<TH1F>("neutrino_px_7", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_7 = book<TH1F>("neutrino_py_7", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_7 = book<TH1F>("neutrino_pz_7", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_8 = book<TH1F>("neutrino_px_8", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_8 = book<TH1F>("neutrino_py_8", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_8 = book<TH1F>("neutrino_pz_8", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_9 = book<TH1F>("neutrino_px_9", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_9 = book<TH1F>("neutrino_py_9", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_9 = book<TH1F>("neutrino_pz_9", "p_{z} #nu", 100, -1000, 1000);
   neutrino_px_10 = book<TH1F>("neutrino_px_10", "p_{x} #nu", 100, -1000, 1000);
   neutrino_py_10 = book<TH1F>("neutrino_py_10", "p_{y} #nu", 100, -1000, 1000);
   neutrino_pz_10 = book<TH1F>("neutrino_pz_10", "p_{z} #nu", 100, -1000, 1000);

   met_px_1 = book<TH1F>("met_px_1", ";MET [GeV]", 500, -2000, 2000);
   met_py_1 = book<TH1F>("met_py_1", ";MET [GeV]", 500, -2000, 2000);
   met_px_2 = book<TH1F>("met_px_2", ";MET [GeV]", 500, -2000, 2000);
   met_py_2 = book<TH1F>("met_py_2", ";MET [GeV]", 500, -2000, 2000);
   met_px_3 = book<TH1F>("met_px_3", ";MET [GeV]", 500, -2000, 2000);
   met_py_3 = book<TH1F>("met_py_3", ";MET [GeV]", 500, -2000, 2000);
   met_px_4 = book<TH1F>("met_px_4", ";MET [GeV]", 500, -2000, 2000);
   met_py_4 = book<TH1F>("met_py_4", ";MET [GeV]", 500, -2000, 2000);
   met_px_5 = book<TH1F>("met_px_5", ";MET [GeV]", 500, -2000, 2000);
   met_py_5 = book<TH1F>("met_py_5", ";MET [GeV]", 500, -2000, 2000);
   met_px_6 = book<TH1F>("met_px_6", ";MET [GeV]", 500, -2000, 2000);
   met_py_6 = book<TH1F>("met_py_6", ";MET [GeV]", 500, -2000, 2000);
   met_px_7 = book<TH1F>("met_px_7", ";MET [GeV]", 500, -2000, 2000);
   met_py_7 = book<TH1F>("met_py_7", ";MET [GeV]", 500, -2000, 2000);
   met_px_8 = book<TH1F>("met_px_8", ";MET [GeV]", 500, -2000, 2000);
   met_py_8 = book<TH1F>("met_py_8", ";MET [GeV]", 500, -2000, 2000);
   met_px_9 = book<TH1F>("met_px_9", ";MET [GeV]", 500, -2000, 2000);
   met_py_9 = book<TH1F>("met_py_9", ";MET [GeV]", 500, -2000, 2000);
   met_px_10 = book<TH1F>("met_px_10", ";MET [GeV]", 500, -2000, 2000);
   met_py_10 = book<TH1F>("met_py_10", ";MET [GeV]", 500, -2000, 2000);

}

void ttDMReconstructionHists_ScansAndMarginalisation::fill(const Event & event){
   
   weight = event.weight;
   
   if(!event.is_valid(h_ttbargen)){
       return;
    }
   const auto & ttbargen = event.get(h_ttbargen);
   
   if (!ttbargen.IsSemiLeptonicDecay()) return;

   k =k+1;
   if (k==1) {
      neutrino_px_1 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
      neutrino_py_1 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
      neutrino_pz_1 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
      met_px_1 ->Fill(event.met->v4().Px(), weight);
      met_py_1 ->Fill(event.met->v4().Py(), weight);
   }
   if (k==2) {
      neutrino_px_2 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
      neutrino_py_2 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
      neutrino_pz_2 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
      met_px_2 ->Fill(event.met->v4().Px(), weight);
      met_py_2 ->Fill(event.met->v4().Py(), weight);
   }
   if (k==3) {
      neutrino_px_3 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
      neutrino_py_3 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
      neutrino_pz_3 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
      met_px_3 ->Fill(event.met->v4().Px(), weight);
      met_py_3 ->Fill(event.met->v4().Py(), weight);
   }
   if (k==4) {
      neutrino_px_4 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
      neutrino_py_4 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
      neutrino_pz_4 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
      met_px_4 ->Fill(event.met->v4().Px(), weight);
      met_py_4 ->Fill(event.met->v4().Py(), weight);
   }
    if (k==5) {
      neutrino_px_5 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
      neutrino_py_5 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
      neutrino_pz_5 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
      met_px_5 ->Fill(event.met->v4().Px(), weight);
      met_py_5 ->Fill(event.met->v4().Py(), weight);
   }
   //  if (k==6) {
   //    neutrino_px_6 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
   //    neutrino_py_6 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
   //    neutrino_pz_6 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
   //    met_px_6 ->Fill(event.met->v4().Px(), weight);
   //    met_py_6 ->Fill(event.met->v4().Py(), weight);
   // }
   // if (k==7) {
   //    neutrino_px_7 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
   //    neutrino_py_7 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
   //    neutrino_pz_7 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
   //    met_px_7 ->Fill(event.met->v4().Px(), weight);
   //    met_py_7 ->Fill(event.met->v4().Py(), weight);
   // } 
   // if (k==8) {
   //    neutrino_px_8 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
   //    neutrino_py_8 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
   //    neutrino_pz_8 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
   //    met_px_8 ->Fill(event.met->v4().Px(), weight);
   //    met_py_8 ->Fill(event.met->v4().Py(), weight);
   // }
   // if (k==9) {
   //    neutrino_px_9 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
   //    neutrino_py_9 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
   //    neutrino_pz_9 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
   //    met_px_9 ->Fill(event.met->v4().Px(), weight);
   //    met_py_9 ->Fill(event.met->v4().Py(), weight);
   // }
   // if (k==10) {
   //    neutrino_px_10 ->Fill(ttbargen.Neutrino().v4().Px(), weight);
   //    neutrino_py_10 ->Fill(ttbargen.Neutrino().v4().Py(), weight);
   //    neutrino_pz_10 ->Fill(ttbargen.Neutrino().v4().Pz(), weight);
   //    met_px_10 ->Fill(event.met->v4().Px(), weight);
   //    met_py_10 ->Fill(event.met->v4().Py(), weight);
   // }
   
   double Mt0 = 164.716; 
   double sigmaMt = 13.9042; 
   
   Double_t MW;
   Double_t DeltaPhi;
   Double_t DeltaEta;
   Double_t MT;
  
   double px = -1050;
   double py = -1050;
   double pz = ttbargen.Neutrino().v4().Pz();
   double energy;
   double gmax =  0.0191095; 
  
   static TFile *f = new TFile("/nfs/dust/cms/user/mameyer/TTDMDM/selection_MUON_RunII_25ns_v2/Likelihood_VariableBinning_2_Filled.root", "READ");
   static TH3F *h=(TH3F*)f->Get("likelihood3D");
   h->Scale(1./h->Integral(1,20,2,19,1,22,""));
   std::cout<<"Integral:  "<<h->Integral(1,20,2,19,1,22,"")<<std::endl;
  
   Particle lepton;
   assert(event.muons);
   if (event.muons->size()>0) lepton = event.muons->at(0);
   
   LorentzVector nu_scan;
   // Jet bjet2 = *nextJet(lepton,*event.jets);  //here next jet as b-jet
   GenParticle bjet2 = ttbargen.BLep();
   Double_t likelihood;
   
   for (int i=0; i<421;i++) 
      {
         for (int j=0; j<421;j++) 
            {
               energy = TMath::Sqrt(px*px + py*py + pz*pz);
               nu_scan.SetPxPyPzE(px, py, pz, energy); 
               MW = (nu_scan+lepton.v4()).M();
               DeltaPhi = nu_scan.Phi()-lepton.v4().Phi();
               if(!std::isnan(DeltaPhi))DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi);
               DeltaEta = nu_scan.Eta()-lepton.v4().Eta();
               MT = (lepton.v4()+nu_scan+bjet2.v4()).M();
               if (fabs(DeltaEta) < 3.5 && MW > 20 && MW < 270) 
                  {
                     double interpolation = h->Interpolate(DeltaEta,DeltaPhi, MW);
                     if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
                     else likelihood = - TMath::Log(interpolation);
                  }
               else
                  {
                     double val_deta = DEtaParametrization(DeltaEta);
                     double val_dphi = DPhiParametrization(DeltaPhi);                              
                     double val_mw = MWParametrization(MW); 
                     likelihood = gmax*val_deta*val_dphi*val_mw;  
                     if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
                     else likelihood = - TMath::Log(likelihood);
                  }
               likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2));
               //test MET in Chi2
               //likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2)) + (TMath::Power((event.met->v4().Px()-px),2) / TMath::Power(27.3,2)) + (TMath::Power((event.met->v4().Py()-py),2) / TMath::Power(27.4,2));
               
               int binx,biny,binz;
               h_px_py_likelihood_1 ->GetBinXYZ(h_px_py_likelihood_1 ->FindBin(px,py),binx,biny,binz);
               if (k==1) h_px_py_likelihood_1 ->SetBinContent(binx,biny,likelihood);
               if (k==2) h_px_py_likelihood_2 ->SetBinContent(binx,biny,likelihood);
               if (k==3) h_px_py_likelihood_3 ->SetBinContent(binx,biny,likelihood);
               if (k==4) h_px_py_likelihood_4 ->SetBinContent(binx,biny,likelihood);
               if (k==5) h_px_py_likelihood_5 ->SetBinContent(binx,biny,likelihood);
               if (k==6) h_px_py_likelihood_6 ->SetBinContent(binx,biny,likelihood);
               if (k==7) h_px_py_likelihood_7 ->SetBinContent(binx,biny,likelihood);
               if (k==8) h_px_py_likelihood_8 ->SetBinContent(binx,biny,likelihood);
               if (k==9) h_px_py_likelihood_9 ->SetBinContent(binx,biny,likelihood);
               if (k==10) h_px_py_likelihood_10 ->SetBinContent(binx,biny,likelihood);
               py = py+5.0;
            }
         py = -1050;
         px = px+ 5.0;
      }
   
   
   // px = -1050;
   // pz = -1050;
   // py = ttbargen.Neutrino().v4().Py();
   
   // for (int i=0; i<421;i++) 
   //    {
   //       for (int j=0; j<421;j++) 
   //          {
   //             energy = TMath::Sqrt(px*px + py*py + pz*pz);
   //             nu_scan.SetPxPyPzE(px, py, pz, energy); 
   //             MW = (nu_scan+lepton.v4()).M();
   //             DeltaPhi = nu_scan.Phi()-lepton.v4().Phi();
   //             if(!std::isnan(DeltaPhi)) DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi);
   //             DeltaEta = nu_scan.Eta()-lepton.v4().Eta();
   //             MT = (lepton.v4()+nu_scan+bjet2.v4()).M();
   //             if (fabs(DeltaEta) < 3.5 && MW > 20 && MW < 270)
   //                {
   //                   double interpolation = h->Interpolate(DeltaEta,DeltaPhi, MW);
   //                   if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
   //                   else likelihood = - TMath::Log(interpolation);
   //                }
   //             else
   //                {
   //                   // extrapolation in deta
   //                   double val_deta = DEtaParametrization(DeltaEta);
   //                   double val_dphi = DPhiParametrization(DeltaPhi);                             
   //                   double val_mw = MWParametrization(MW); 
   //                   likelihood = gmax*val_deta*val_dphi*val_mw;     
   //                   if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
   //                   else likelihood = - TMath::Log(likelihood);
   //                }
   
   //             likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2));
   //             int binx,biny,binz;
   //             h_px_pz_likelihood_1 ->GetBinXYZ(h_px_pz_likelihood_1 ->FindBin(px,pz),binx,biny,binz);
   //             if (k==1) h_px_pz_likelihood_1 ->SetBinContent(binx,biny,likelihood);
   //             if (k==2) h_px_pz_likelihood_2 ->SetBinContent(binx,biny,likelihood);
   //             if (k==3) h_px_pz_likelihood_3 ->SetBinContent(binx,biny,likelihood);
   //             if (k==4) h_px_pz_likelihood_4 ->SetBinContent(binx,biny,likelihood);
   //             if (k==5) h_px_pz_likelihood_5 ->SetBinContent(binx,biny,likelihood);
   //             if (k==6) h_px_pz_likelihood_6 ->SetBinContent(binx,biny,likelihood);
   //             if (k==7) h_px_pz_likelihood_7 ->SetBinContent(binx,biny,likelihood);
   //             if (k==8) h_px_pz_likelihood_8 ->SetBinContent(binx,biny,likelihood);
   //             if (k==9) h_px_pz_likelihood_9 ->SetBinContent(binx,biny,likelihood);
   //             if (k==10) h_px_pz_likelihood_10 ->SetBinContent(binx,biny,likelihood);
   //             pz = pz+5.0;
   //        }
   //     px = px+ 5.0;
   //     pz = -1050;
   
   //  }
   
   // py = -1050;
   // pz = -1050;
   // px = ttbargen.Neutrino().v4().Px();
   
   // for (int i=0; i<421;i++) 
   //    {
   //       for (int j=0; j<421;j++) 
   //        {
   //           energy = TMath::Sqrt(px*px + py*py + pz*pz);
   //           nu_scan.SetPxPyPzE(px, py, pz, energy); 
   //           MW = (nu_scan+lepton.v4()).M();
   //           DeltaPhi = nu_scan.Phi()-lepton.v4().Phi();
   //           if(!std::isnan(DeltaPhi))DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi);
   //           DeltaEta = nu_scan.Eta()-lepton.v4().Eta();
   //           MT = (lepton.v4()+nu_scan+bjet2.v4()).M();
   //           if (fabs(DeltaEta) < 3.5 && MW > 20 && MW < 270)
   //              {
   //                 double interpolation = h->Interpolate(DeltaEta,DeltaPhi, MW);
   //                 if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
   //                 else likelihood = - TMath::Log(interpolation);
   //              }
   //           else
   //              {
   //                 // extrapolation in deta
   //                 double val_deta = DEtaParametrization(DeltaEta);
   //                 double val_dphi = DPhiParametrization(DeltaPhi);                            
   //                 double val_mw = MWParametrization(MW); 
   //                 likelihood = gmax*val_deta*val_dphi*val_mw;    
   //                 if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
   //                 else likelihood = - TMath::Log(likelihood);
   //              }
   //           likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2));
   //           int binx,biny,binz;
   //           h_py_pz_likelihood_1 ->GetBinXYZ(h_py_pz_likelihood_1 ->FindBin(py,pz),binx,biny,binz);
   //           if (k==1) h_py_pz_likelihood_1 ->SetBinContent(binx,biny,likelihood);
   //           if (k==2) h_py_pz_likelihood_2 ->SetBinContent(binx,biny,likelihood);
   //           if (k==3) h_py_pz_likelihood_3 ->SetBinContent(binx,biny,likelihood);
   //           if (k==4) h_py_pz_likelihood_4 ->SetBinContent(binx,biny,likelihood);
   //           if (k==5) h_py_pz_likelihood_5 ->SetBinContent(binx,biny,likelihood);
   //           if (k==6) h_py_pz_likelihood_6 ->SetBinContent(binx,biny,likelihood);
   //           if (k==7) h_py_pz_likelihood_7 ->SetBinContent(binx,biny,likelihood);
   //           if (k==8) h_py_pz_likelihood_8 ->SetBinContent(binx,biny,likelihood);
   //           if (k==9) h_py_pz_likelihood_9 ->SetBinContent(binx,biny,likelihood);
   //           if (k==10) h_py_pz_likelihood_10 ->SetBinContent(binx,biny,likelihood);
   //           pz = pz+5.0;
   //        }
   //       py = py+ 5.0;
   //       pz = -1050;
   //    }
   
   
   
   // Marginalisation
   
   //          static TFile *f = new TFile("/nfs/dust/cms/user/mameyer/TTDMDM/selection_MUON_RunII_25ns_76X/Marginalisation.root", "RECREATE");
   //          TH2F *hist = new TH2F("hist", "hist",51,-1050,1050,51,-1050,1050);
   //          py = -1050;
   //          px = -1050;
   //          pz =ttbargen.Neutrino().v4().Pz();
   //          for (int i=0; i<52;i++) 
   //             {
   //                for (int j=0; j<52;j++) 
   //                   {
   //                      pz =ttbargen.Neutrino().v4().Pz();
   //                      energy = TMath::Sqrt(px*px + py*py + pz*pz);
   //                      Start.SetPxPyPzE(px, py, pz, energy); 
   //                      TVirtualFitter::SetDefaultFitter("Minuit");
   //                      fit->SetFCN(FCN_Likelihood);
                  
   //                      fit->SetParameter(0,"px_neutrino",Start.Px(),0,-1050,1050); 
   //                      fit->SetParameter(1,"py_neutrino",Start.Py(),0,-1050,1050); 
   //                      fit->SetParameter(2,"pz_neutrino",Start.Pz(),1,-1800,1800);
                  
   //                      double arglist[2]={200,0.01}; //aendern 
   //                      fit->ExecuteCommand("MINUIT", arglist,2);
   //                      energy = TMath::Sqrt(fit->GetParameter(0)*fit->GetParameter(0) + fit->GetParameter(1)*fit->GetParameter(1) + fit->GetParameter(2)*fit->GetParameter(2));
   //                      nu_scan.SetPxPyPzE(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), energy); 
   //                      MW = (nu_scan+lepton.v4()).M();
   //                      DeltaPhi = nu_scan.Phi()-lepton.v4().Phi();
   //                      if(!std::isnan(DeltaPhi))DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi);
   //                      DeltaEta = nu_scan.Eta()-lepton.v4().Eta();
   //                      MT = (lepton.v4()+nu_scan+bjet.v4()).M();
                  
   //                      if (fabs(DeltaEta) < 3.5 && MW > 20 && MW < 270) 
   //                         {
   //                            double interpolation = h->Interpolate(DeltaEta,DeltaPhi, MW);
   //                            if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
   //                            else likelihood = - TMath::Log(interpolation);
   //                         }
   //                      else
   //                         {
   //                            double val_deta = DEtaParametrization(DeltaEta);
   //                            double val_dphi = DPhiParametrization(DeltaPhi);                              
   //                            double val_mw = MWParametrization(MW); 
   //                            likelihood = gmax*val_deta*val_dphi*val_mw;
   //                            if (TMath::Sqrt(px*px + py*py) < 5) likelihood = 500;
    //                            else likelihood = - TMath::Log(likelihood);
   //                         }
   //                      likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2));
   //                      int binx,biny,binz;
   //                      hist ->GetBinXYZ(hist ->FindBin(px,py),binx,biny,binz);
   //                      hist ->SetBinContent(binx,biny,likelihood);
   //                      py = py+ 42.0;
   //                   }
   //                px = px+ 42.0;
   //                py=-1050;
   //             }
   //          hist->Write();
   //          f->Write();
   //          f->Close();
   return;

}


Double_t ttDMReconstructionHists_ScansAndMarginalisation::MWParametrization(Double_t mw)
{
   static TF1* f = new TF1("f", "[0]*TMath::Exp([1]*TMath::Abs(x-[2]))+[3]*TMath::Exp(-0.5*TMath::Power((x-[2])/[4],2))", 50, 110);

   f->SetParameter(0, 1.06420e+03);
   f->SetParameter(1, -1.51973e-01);
   f->SetParameter(2, 8.06140e+01);   
   f->SetParameter(3, 8.57799e+03);
   f->SetParameter(4, 2.34007e+00);

   Double_t y = f->Eval(mw) / 9633.26; // normalize to peak value

   return y;

}

Double_t ttDMReconstructionHists_ScansAndMarginalisation::DPhiParametrization(Double_t dphi)
{
   static TF1* f = new TF1("f", "[0]*TMath::Exp([1]*TMath::Abs(x)) + [3]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)-[2])/[4],2))", -3.1415, 3.1415);

   f->SetParameter(0, 1.40257e+03);
   f->SetParameter(1, -2.21625e+00);
   f->SetParameter(2, 5.18883e-01);   
   f->SetParameter(3, 1.67336e+03);
   f->SetParameter(4, 2.28831e-01);

   Double_t y = f->Eval(dphi) / 2133.87; // normalize to peak value

   return y;

}

Double_t ttDMReconstructionHists_ScansAndMarginalisation::DEtaParametrization(Double_t dphi)
{
   static TF1* f = new TF1("f", "[0]*TMath::Exp([1]*TMath::Abs(x))*TMath::Abs(TMath::ATan(3*x)) + [3]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)-[2])/[4],2))", -2.5, 2.5);

   f->SetParameter(0, 1.85995e+02);
   f->SetParameter(1, -1.65177e+00);
   f->SetParameter(2, 3.25931e-01);   
   f->SetParameter(3, 2.28168e+03);
   f->SetParameter(4, 2.95958e-01);   

   Double_t y = f->Eval(dphi) / 2365.73; // normalize to peak value

   return y;

}
