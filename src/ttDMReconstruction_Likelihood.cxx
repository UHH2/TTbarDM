#include "../include/ttDMReconstruction_Likelihood.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h" 
#include "UHH2/common/include/Utils.h"
#include "TFitter.h" 
#include "TVirtualFitter.h"
#include "TMinuit.h" 
#include "FCN_Likelihood.h"
#include "TH2F.h"
#include "TError.h"


ttDMReconstruction_Likelihood::ttDMReconstruction_Likelihood(uhh2::Context & ctx)
 {
    h_likelihood =ctx.declare_event_output<double>("likelihood");
    h_recneutrino =ctx.declare_event_output<LorentzVector>("rec_neutrino");
    h_bjet =ctx.declare_event_output<Jet>("bjet");
    h_ttbargen =ctx.get_handle<TTbarGen>("ttbargen");
    h_muons=ctx.get_handle<std::vector<Muon>>("h_muons_medium");

    f = new TFile("/nfs/dust/cms/user/mameyer/TTDMDM/selection_MUON_RunII_25ns_v2/Likelihood_VariableBinning_2_Filled.root", "READ");
    h = (TH3F*)f->Get("likelihood3D");
    h->Scale(1./h->Integral(1,20,2,19,1,22,""));
 }

bool ttDMReconstruction_Likelihood::process(uhh2::Event & e)
{

   bool print = false;
   if (print) std::cout<<"Start Minimization"<<std::endl;
   if (print) std::cout<<"integral: "<<h->Integral(1,20,2,19,1,22,"")<<std::endl;

   //assert(e.muons);
   assert(e.electrons);
   std::vector<Muon> muons=e.get(h_muons);

   if (muons.size()>0) lepton = muons.at(0);
   else lepton = e.electrons->at(0);
   indata_likelihood.lep = lepton.v4(); 
   indata_likelihood.met = e.met->v4();
    
   // set up fit
   NFitPar=3;
   TVirtualFitter::SetDefaultFitter("Minuit");
   TFitter *fit=(TFitter*) TVirtualFitter::Fitter(NULL,NFitPar);
   fit->GetMinuit()->SetPrintLevel(-1);
   fit->SetFCN(FCN_Likelihood);
   double arglist[2]={200,0.01}; 
   
   
   if (print) std::cout<<"Perform Fits to find the b-jet and determine pz-coordinate"<<std::endl;
   if (print) std::cout<<"starting point for fit: four momentum of lepton + DeltaEta(leton, neutrino)_Mean/DeltaPhi(lepton, neutrino)_Mean"<<std::endl;
   LorentzVector Start;
   Start= lepton.v4();
   Start.SetEta(lepton.v4().Eta()+0.0006707); 
   Start.SetPhi(TVector2::Phi_mpi_pi(lepton.v4().Phi()-0.002275)); 
   if (print){
      std::cout<<"starting point px: "<<Start.Px()<<std::endl;
      std::cout<<"starting point py: "<<Start.Py()<<std::endl;
      std::cout<<"starting point pz: "<<Start.Pz()<<std::endl;
   }

   LorentzVector nu_scan;
   double MT;
   double Mt0 = 164.716; 
   double sigmaMt = 13.9042; 
   double energy;
   Jet bjet;
   Jet bjet_min;
   LorentzVector bjet_v4;
   double res_min = 1000000;
   double pz_min = 10000;
   
   // loop over jets to find the b-jet, take jet which gives the best top mass constraint
   for (unsigned int i =0; i<e.jets->size();i++)
      {
         bjet = e.jets->at(i);
         indata_likelihood.bjet = bjet.v4();
         bjet_v4 = bjet.v4();
               
         fit->SetParameter(0,"px_neutrino",Start.Px(),10,-1000,1000); 
         fit->SetParameter(1,"py_neutrino",Start.Py(),10,-1050,1050); 
         fit->SetParameter(2,"pz_neutrino",Start.Pz(),10,-1800,1800);
         
         int flag = -9989;
         double gin[NFitPar];
         double res;
         double pars[NFitPar] = {fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2)};
         if (print){
            std::cout << "\n\n\nLikelihood before fit:" << std::endl;
            FCN_Likelihood(NFitPar, gin, res, pars, flag);
         }
         
         fit->ExecuteCommand("MINUIT", arglist, 2);
         
         if(print){
            std::cout << "\n\n\nLikelihood after fit:" << std::endl;
         flag = -9989;
         pars[0] = fit->GetParameter(0);
         pars[1] = fit->GetParameter(1);
         pars[2] = fit->GetParameter(2);
         FCN_Likelihood(NFitPar, gin, res, pars, flag);
         std::cout<<"Results Likelihood= "<<res<<std::endl;
         }
         
         energy = TMath::Sqrt(fit->GetParameter(0)*fit->GetParameter(0) + fit->GetParameter(1)*fit->GetParameter(1)  + fit->GetParameter(2)*fit->GetParameter(2) );
         nu_scan.SetPxPyPzE(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), energy); 
         MT = (indata_likelihood.lep+nu_scan+indata_likelihood.bjet).M();
         double res_T = TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2);
         //double res_T = TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2)+(TMath::Power((indata_likelihood.met.Px()-fit->GetParameter(0)),2) / TMath::Power(27.3,2)) + (TMath::Power((indata_likelihood.met.Py()-fit->GetParameter(1)),2) / TMath::Power(27.4,2));
         if (res_T < res_min)
            {
               res_min =res_T;
               bjet_min =bjet;
               pz_min = fit->GetParameter(2);
            }
       }
    
    indata_likelihood.bjet = bjet_min.v4();
    bjet_v4 = bjet_min.v4();
    
    if (print) {
       std::cout<<"Found b-jet:  "<<bjet_min.v4().pt()<<"   "<<bjet_min.v4().eta()<<"   "<<bjet_min.v4().phi()<<std::endl;
       std::cout<<"After Fit: perform rough scan of likelihood as a function of px and py to find approximate position of the minimum, set pz to value found by the previous fit"<<std::endl;
       std::cout<<"pz was: "<<fit->GetParameter(2)<<std::endl;
    }
    

    //scan with pZ
    double px = -1050;
    double py = -1050;
    double MW;
    double DeltaPhi;
    double DeltaEta;
    double pz = pz_min;
    Double_t likelihood;
    double likelihood_min=100000;
    double px_min=2000;
    double py_min=2000;
    double gmax = 0.0191095; 
  
    for (int i=0; i<52;i++) 
       {
         for (int j=0; j<52;j++) 
            {

               energy = TMath::Sqrt(px*px + py*py + pz*pz);
               nu_scan.SetPxPyPzE(px, py, pz, energy); 
               MW = (nu_scan+lepton.v4()).M();
               DeltaPhi = nu_scan.Phi()-lepton.v4().Phi();
               if(!std::isnan(DeltaPhi))DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi);
               DeltaEta = nu_scan.Eta()-lepton.v4().Eta();
               MT = (lepton.v4()+nu_scan+bjet_v4).M();                  
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
               //likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2)) + (TMath::Power((e.met->v4().Px()-px),2) / TMath::Power(27.3,2)) + (TMath::Power((e.met->v4().Py()-py),2) / TMath::Power(27.4,2));   
               if (likelihood < likelihood_min)
                  {
                     likelihood_min = likelihood;
                     px_min=px;
                     py_min=py;
                  }
               py = py+42.0;
            }
         px = px+ 42.0;
         py = -1050;
      }
    if (print){
       std::cout<<"Scan found minimum at px = "<<px_min<< "and py = "<<py_min<<std::endl;
       std::cout<<"Perform minimization, starting point for px and py from scan, pz from first minimization: "<<std::endl;
       
       std::cout<<"Starting point px = "<<px_min<<std::endl;
       std::cout<<"Starting point py = "<<py_min<<std::endl;
       std::cout<<"Starting point pz = "<<pz<<std::endl;
    }
    
    energy = TMath::Sqrt(px_min*px_min + py_min*py_min + pz*pz);
    Start.SetPxPyPzE(px_min,py_min, pz, energy);
    if (print){
       std::cout<<"Starting point px = "<<px_min<<std::endl;
       std::cout<<"Starting point py = "<<py_min<<std::endl;
       std::cout<<"Starting point pz = "<<pz<<std::endl;
    }
    fit->SetParameter(0,"px_neutrino",Start.Px(),1,-1000,1000); 
    fit->SetParameter(1,"py_neutrino",Start.Py(),1,-1050,1050); 
    fit->SetParameter(2,"pz_neutrino",Start.Pz(),1,-1800,1800);
    fit->ExecuteCommand("MINUIT", arglist,2);
    
    if (print){
       std::cout<<"Result of the fit:  "<<std::endl;
       std::cout<<"               px = "<<fit->GetParameter(0)<<std::endl;
       std::cout<<"               py = "<<fit->GetParameter(1)<<std::endl;
       std::cout<<"               pz = "<<fit->GetParameter(2)<<std::endl;
    }
    
    energy = TMath::Sqrt(fit->GetParameter(0)*fit->GetParameter(0) + fit->GetParameter(1)*fit->GetParameter(1)  + fit->GetParameter(2)*fit->GetParameter(2) );
    nu_scan.SetPxPyPzE(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), energy); 
    MW = (nu_scan+indata_likelihood.lep).M();
    DeltaPhi = nu_scan.Phi()-indata_likelihood.lep.Phi();
    if (!std::isnan(DeltaPhi)) DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi); 
    DeltaEta = nu_scan.Eta()-indata_likelihood.lep.Eta();
    MT = (indata_likelihood.lep+nu_scan+indata_likelihood.bjet).M();
    
    if (fabs(DeltaEta) < 3.5 && MW > 20 && MW < 270 )
       {
          double interpolation = h->Interpolate(DeltaEta,DeltaPhi, MW);
          if (TMath::Sqrt(fit->GetParameter(0)*fit->GetParameter(0) + fit->GetParameter(1)*fit->GetParameter(1)) < 5) likelihood = 500;
          else likelihood = - TMath::Log(interpolation);
       }
    else
       {
          double val_deta = DEtaParametrization(DeltaEta);
          double val_dphi = DPhiParametrization(DeltaPhi);                              
          double val_mw = MWParametrization(MW); 
          likelihood = gmax*val_deta*val_dphi*val_mw; 
          if (TMath::Sqrt(fit->GetParameter(0)*fit->GetParameter(0) + fit->GetParameter(1)*fit->GetParameter(1)) < 5) likelihood = 500;
          else likelihood = - TMath::Log(likelihood);   
       }
    likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2));
    //likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2)) + (TMath::Power((e.met->v4().Px()-fit->GetParameter(0)),2) / TMath::Power(27.3,2)) + (TMath::Power((e.met->v4().Py()-fit->GetParameter(1)),2) / TMath::Power(27.4,2));
    
    if (print){
    std::cout<<"               chi^2 = "<<likelihood <<std::endl;
    std::cout<<"if this was a semi-leptonic event we can compare the result to the generated neutrino: "<<std::endl;
    
    if(!e.is_valid(h_ttbargen)){
       const auto & ttbargen = e.get(h_ttbargen);
       if (ttbargen.IsSemiLeptonicDecay())
          {
             std::cout<<"generated neutrino px: "<<ttbargen.Neutrino().v4().px()<<std::endl;
             std::cout<<"generated neutrino py: "<<ttbargen.Neutrino().v4().py()<<std::endl;
             std::cout<<"generated neutrino pz: "<<ttbargen.Neutrino().v4().pz()<<std::endl;
             
             std::cout<<"pxrec/pxgen: "<< fit->GetParameter(0) / ttbargen.Neutrino().v4().px()<<std::endl;
             std::cout<<"pyrec/pygen: "<< fit->GetParameter(1) / ttbargen.Neutrino().v4().py()<<std::endl;
             std::cout<<"pzrec/pzgen: "<< fit->GetParameter(2) / ttbargen.Neutrino().v4().pz()<<std::endl;
          }
       else std::cout<<"not semi-leptonic"<<std::endl;
    }
    }

    //    std::cout<<"Event done"<<std::endl;

    e.set(h_likelihood, likelihood);
    e.set(h_recneutrino, nu_scan);
    e.set(h_bjet, bjet_min);
   return true;
}
 
 
   

