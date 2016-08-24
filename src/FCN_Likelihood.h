#ifndef FCN_Likelihood__h
#include "TFile.h"
#include "TH3F.h"
#include "TF1.h"
#include "TError.h"
Double_t MWParametrization(Double_t mw);
Double_t DPhiParametrization(Double_t dphi);
Double_t DEtaParametrization(Double_t dphi);

struct fitdata{
   LorentzVector lep;
   LorentzVector bjet;
   LorentzVector met;
};

fitdata indata_likelihood;

void FCN_Likelihood(Int_t &npar, double *gin, double &result, double *par, int iflg){

   static TFile *f = new TFile("/nfs/dust/cms/user/mameyer/TTDMDM/selection_MUON_RunII_25ns_v2/Likelihood_VariableBinning_2_Filled.root", "READ");
   static TH3F *h=(TH3F*)f->Get("likelihood3D");                                                      //Likelihood hist has been divided by bin width and normalized to one. 
   h->Scale(1./h->Integral(1,20,2,19,1,22,""));
   double gmax = 0.0191095;                                                                           //value is needed for correct normalization, has to be updated if likelihood changes
   
   bool print = false;
   if (iflg==-9989) print = true;
  
   double energy = TMath::Sqrt(par[0]*par[0]+par[1]*par[1]+par[2]*par[2]);
   LorentzVector nu;
   nu.SetPxPyPzE(par[0], par[1], par[2], energy);  
   LorentzVector lep = indata_likelihood.lep;  
   LorentzVector bjet = indata_likelihood.bjet;  

   double Mt0 = 164.716; 
   double sigmaMt = 13.9042; 
   Double_t MW = (nu+lep).M();
   Double_t DeltaPhi = nu.Phi()-lep.Phi();
   if (!std::isnan(DeltaPhi)) DeltaPhi= TVector2::Phi_mpi_pi(DeltaPhi);
   Double_t DeltaEta = nu.Eta()-lep.Eta(); 
   Double_t MT = (lep+nu+bjet).M();
   
   if (print){
   std::cout<<"MW: "<<MW<<std::endl;
   std::cout<<"MT: "<<MT<<std::endl;
   std::cout<<"DeltaPhi: "<<DeltaPhi<<std::endl;
   std::cout<<"DeltaEta: "<<DeltaEta<<std::endl;
   }   
   
   Double_t likelihood;
   
   if (fabs(DeltaEta) < 3.5 && MW > 20 && MW < 270 )
      {
         double interpolation = h->Interpolate(DeltaEta,DeltaPhi, MW);
         if (TMath::Sqrt(par[0]*par[0] + par[1]*par[1]) < 5) likelihood = 500;
         else likelihood = - TMath::Log(interpolation);
      }
   else
      {
         double val_deta = DEtaParametrization(DeltaEta);
         double val_dphi = DPhiParametrization(DeltaPhi);                              
         double val_mw = MWParametrization(MW); 
         
         likelihood = gmax*val_deta*val_dphi*val_mw;   
         if (TMath::Sqrt(par[0]*par[0] + par[1]*par[1]) < 5) likelihood = 500;
         else likelihood = - TMath::Log(likelihood);
      }
   if (print){
      std::cout << "L = " << likelihood << std::endl;
      std::cout << "L(mt) = " << (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2)) << std::endl;
   }
   //calculate chi2
   likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2)) ;
   //likelihood = likelihood*2 + (TMath::Power((MT-Mt0),2)/TMath::Power(sigmaMt,2)) + (TMath::Power((indata_likelihood.met.Px()-par[0]),2) / TMath::Power(27.3,2)) + (TMath::Power((indata_likelihood.met.Py()-par[1]),2) / TMath::Power(27.4,2));   
   result =likelihood;
   return;
}

Double_t MWParametrization(Double_t mw)
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

Double_t DPhiParametrization(Double_t dphi)
{
   TF1* f = new TF1("f", "[0]*TMath::Exp([1]*TMath::Abs(x)) + [3]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)-[2])/[4],2))", -3.1415, 3.1415);

   f->SetParameter(0, 1.40257e+03);
   f->SetParameter(1, -2.21625e+00);
   f->SetParameter(2, 5.18883e-01);   
   f->SetParameter(3, 1.67336e+03);
   f->SetParameter(4, 2.28831e-01);

   Double_t y = f->Eval(dphi) / 2133.87; // normalize to peak value

   return y;

}

Double_t DEtaParametrization(Double_t dphi)
{
   TF1* f = new TF1("f", "[0]*TMath::Exp([1]*TMath::Abs(x))*TMath::Abs(TMath::ATan(3*x)) + [3]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)-[2])/[4],2))", -2.5, 2.5);

   f->SetParameter(0, 1.85995e+02);
   f->SetParameter(1, -1.65177e+00);
   f->SetParameter(2, 3.25931e-01);   
   f->SetParameter(3, 2.28168e+03);
   f->SetParameter(4, 2.95958e-01);   

   Double_t y = f->Eval(dphi) / 2365.73; // normalize to peak value

   return y;

}



#endif

