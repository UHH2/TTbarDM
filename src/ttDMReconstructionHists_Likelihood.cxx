#include "UHH2/TTbarDM/include/ttDMReconstructionHists_Likelihood.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/TTbarDM/include/MT2Utility.h"
#include "UHH2/TTbarDM/include/ttDMSemiLeptonicUtils.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TError.h"
using namespace uhh2;
double METUncertainty_px(MET *met);
double METUncertainty_py(MET *met);
double CalcUncertainty(MET *met, bool px, double up, double down);


ttDMReconstructionHists_Likelihood::ttDMReconstructionHists_Likelihood(Context & ctx, const std::string & dirname): Hists(ctx, dirname){

   h_likelihood = ctx.get_handle<double>("likelihood");
   h_recneutrino =ctx.get_handle<LorentzVector>("rec_neutrino");
   h_ttbargen =ctx.get_handle<TTbarGen>("ttbargen");
   h_bjets=ctx.get_handle<Jet>("bjet");
   h_heptopjets_WP3=ctx.get_handle<std::vector<TopJet>>("h_heptopjets_WP3");

   hist_chi2 = book<TH1F>("chi2","#chi^{2}", 100, 0, 100);
   
   hist_pxrec_pxgen = book<TH1F>("pxrec_pxgen","p_{x}^{rec} - p_{x}^{gen}", 100, -500, 500);
   hist_pyrec_pygen = book<TH1F>("pyrec_pygen","p_{y}^{rec} - p_{y}^{gen}", 100, -500, 500);
   hist_pzrec_pzgen = book<TH1F>("pzrec_pzgen","p_{z}^{rec} - p_{z}^{gen}", 100, -500, 500);
   
   hist_DM_MET = book<TH1F>("DM_MET","|#vec{MET} - #vec{p_{T,#nu}^{rec}}| [GeV]", 500, 0, 5000);
   hist_DM_MET_chi2 = book<TH1F>("DM_MET_chi2","DM_MET_chi2", 50, 0, 20000);
   
   hist_unc_px = book<TH1F>("unc_px","unc_px", 30, 0, 300);
   hist_unc_py = book<TH1F>("unc_py","unc_py", 30, 0, 300);
   hist_relunc_px = book<TH1F>("relunc_px","relunc_px", 10, 0, 100);
   hist_relunc_py = book<TH1F>("relunc_py","relunc_py", 10, 0, 100);
   
   //test different variables
   hist_chi2_METpT_components = book<TH1F>("chi2_METpT_components","chi2_METpT_components", 50, 0, 1000);    
   hist_pTmisspT = book<TH1F>("pTmisspT","pTmisspT", 100, -2000, 2000);         
   hist_pTmisspT_fabs = book<TH1F>("pTmisspT_abs","|pTmisspT|", 50, 0, 2000);         //done 
   hist_pTmisspT_quotient = book<TH1F>("pTmisspT_quotient","pTmisspT_quotient", 50, 0, 50);        
   hist_2D_ptmiss_pT=book<TH2F>("hist_2D_ptmiss_pT","hist_2D_ptmiss_pT", 150, 0, 1500, 150,0, 1500);                 //done
   hist_2D_METx_px =book<TH2F>("hist_2D_METx_px","hist_2D_METx_px", 200, -1000, 1000, 200,-1000, 1000);              //done
   hist_2D_METy_py=book<TH2F>("hist_2D_METy_py","hist_2D_METy_py",200, -1000, 1000, 200,-1000, 1000);                //done
    
   hist_Metxpx= book<TH1F>("hist_Metxpx","hist_Metxpx", 100, -2000, 2000); 
   hist_Metypy = book<TH1F>("hist_Metypy","hist_Metypy", 100, -2000, 2000);
   hist_Metxpx_quad = book<TH1F>("hist_Metxpx_quad","hist_Metxpx_quad", 50, 0, 500000); 
   hist_Metypy_quad = book<TH1F>("hist_Metypy_quad","hist_Metypy_quad", 50, 0, 500000);
   hist_Metxpx_Metypy_quad = book<TH1F>("hist_Metxpx_Metypy_quad","hist_Metxpx_Metypy_quad", 50, 0, 1000000); 

   //DM_MET with generated neutrino
   hist_DM_MET_gen = book<TH1F>("DM_MET_gen","|#vec{MET} - #vec{p_{T,#nu}^{gen}}| [GeV]", 500, 0, 5000);
   
   //DM_MET and chi2 only for semi-leptonic decays
   hist_DM_MET_semilept = book<TH1F>("DM_MET_semi-lept","|#vec{MET} - #vec{p_{T,#nu}^{rec}}|, only semi-leptonic decays [GeV]", 500, 0, 5000);
   hist_DM_MET_chi2_semilept = book<TH1F>("DM_MET_chi2_semi-lept","DM_MET_chi2_semi-lept", 250, 0, 1000);
   
   //plot distance of closest jet to lepton and generated b-jet
   hist_DeltaR_genb_nextjet = book<TH1F>("DeltaR_genb_nextjet","#Delta R(gen b-jet, closest jet to lepton)", 100, -5, 5);
   hist_DeltaR_genb_nextjet_lowchi2 = book<TH1F>("DeltaR_genb_nextjet_lowchi2","#Delta R(gen b-jet, closest jet to lepton), low #chi^{2}", 100, -5, 5);
   hist_DeltaR_genb_nextjet_highchi2 = book<TH1F>("DeltaR_genb_nextjet_highchi2","#Delta R(gen b-jet, closest jet to lepton), high #chi^{2}", 100, -5, 5);
   
   //compatiblity MET with generated neutrino
   hist_MET_Neutrino_px = book<TH1F>("MET_Neutrino_px","(MET_{x} - p_{x, #nu}^{gen})/ p_{x, #nu}^{gen}", 50, -5, 5);
   hist_MET_Neutrino_py = book<TH1F>("MET_Neutrino_py","(MET_{y} - p_{y, #nu}^{gen})/ p_{y, #nu}^{gen}", 50, -5, 5);

   //compatiblity MET with generated neutrino + generated DM particle
   hist_MET_Neutrino_DM_px = book<TH1F>("MET_Neutrino_DM_px","(MET_{x} - (p_{x, #nu}^{gen} + p_{x, DM}^{gen}))/ MET_{x}", 50, -5, 5);
   hist_MET_Neutrino_DM_py = book<TH1F>("MET_Neutrino_DM_py","(MET_{y} - (p_{y, #nu}^{gen} + p_{y, DM}^{gen}))/ MET_{y}", 50, -5, 5);
   
   //compatibility MET with reconstructed neutrino
   hist_MET_RecNeutrino_px = book<TH1F>("MET_RecNeutrino_px","(MET_{x} - p_{x, #nu}^{rec})/ MET_{x}", 50, -5, 5);
   hist_MET_RecNeutrino_py = book<TH1F>("MET_RecNeutrino_py","(MET_{y} - p_{y, #nu}^{rec})/ MET_{y}", 50, -5, 5);
   hist_MET_RecNeutrino_px_abs = book<TH1F>("MET_RecNeutrino_px_abs","|(MET_{x} - p_{x, #nu}^{rec})/ MET_{x}|", 50, 0, 5);
   hist_MET_RecNeutrino_py_abs = book<TH1F>("MET_RecNeutrino_py_abs","|(MET_{y} - p_{y, #nu}^{rec})/ MET_{y}|", 50, 0, 5);
   
   //compatiblity MET with reconstructed neutrino + generated DM particle
   hist_MET_RecNeutrino_DM_px = book<TH1F>("MET_RecNeutrino_DM_px","(MET_{x} - (p_{x, #nu}^{rec} + p_{x, DM}^{rec}))/ MET_{x}", 50, -5, 5);
   hist_MET_RecNeutrino_DM_py = book<TH1F>("MET_RecNeutrino_DM_py","(MET_{y} - (p_{y, #nu}^{rec} + p_{y, DM}^{rec}))/ MET_{y}", 50, -5, 5);

   hist_DMMET_GenDMpT = book<TH1F>("hist_DMMET_GenDMpT","|#vec{MET} - #vec{p_{T,#nu}^{rec}}| - |p_{T}^{DM}|", 100, -500, 500);
   hist_pxrec_pxgen_100 = book<TH1F>("pxrec_pxgen_100","p_{x}^{rec}/p_{x}^{gen}, DMMET > 100 GeV", 100, -10, 10);
   hist_pyrec_pygen_100 = book<TH1F>("pyrec_pygen_100","p_{y}^{rec}/p_{y}^{gen}, DMMET > 100 GeV", 100, -10, 10);
   hist_DMMET_GenDMpT_vec = book<TH1F>("hist_DMMET_GenDMpT_vec","|#vec{MET} - #vec{p_{T,#nu}^{rec}} - #vec{p_{T}^{DM}}|", 100, -500, 500);
   hist_DMMET_GenDMpT_px = book<TH1F>("hist_DMMET_GenDMpT_px","|#vec{MET_{x}} - #vec{p_{T,#nu,x}^{rec}}| - |p_{T,x}^{DM}|", 100, -500, 500);
   hist_DMMET_GenDMpT_py = book<TH1F>("hist_DMMET_GenDMpT_py","|#vec{MET_{y}} - #vec{p_{T,#nu,y}^{rec}}| - |p_{T,y}^{DM}|", 100, -500, 500);
   hist_DMMET_GenDMpT_genmet = book<TH1F>("hist_DMMET_GenDMpT_genmet","|#vec{MET^{gen}} - #vec{p_{T,#nu}^{rec}}| - |p_{T}^{DM}|", 100, -500, 500);
   hist_DMMETgen_GenDMpT = book<TH1F>("hist_DMMETgen_GenDMpT","|#vec{MET} - #vec{p_{T,#nu}^{gen}}| - |p_{T}^{DM}|, gen", 100, -500, 500);
   hist_DMMET_GenDMpT_semilept = book<TH1F>("hist_DMMET_GenDMpT_semilept","|#vec{MET} - #vec{p_{T,#nu}^{rec}}| - |p_{T}^{DM}|, semi-lept.", 100, -500, 500);
   hist_DMMETgen_GenDMpT_only = book<TH1F>("hist_DMMET_GenDMpT_only","|#vec{MET^{gen}} - #vec{p_{T,#nu}^{gen}}| - |p_{T}^{DM}|, only gen", 100, -500, 500);

   hist_n_nu = book<TH1F>("hist_n_nu","N_{#nu}", 100, 0, 10);
   hist_neutrino_pT = book<TH1F>("hist_neutrino_pT","p_{T,rec}^{#nu}",100,0,1000);
   hist_neutrino_pT_gen = book<TH1F>("hist_neutrino_pT_gen","p_{T,gen}^{#nu}",100,0,1000);
   hist_neutrino_phi = book<TH1F>("hist_neutrino_phi","#phi^{rec, #nu}",50,-TMath::Pi(),TMath::Pi());
   hist_neutrino_phi_gen = book<TH1F>("hist_neutrino_phi_gen","#phi^{gen, #nu}",50,-TMath::Pi(),TMath::Pi());
   hist_neutrino_gen_pT = book<TH1F>("hist_neutrino_gen_pT","p_{T,rec}^{#nu} - p_{T,gen}^{#nu}",100,-500,500);
   hist_neutrino_gen_phi = book<TH1F>("hist_neutrino_gen_phi","#phi^{rec, #nu}-#phi^{gen, #nu}",50,-TMath::Pi(),TMath::Pi());
   
   //ptrec-ptgen binnend in ptgen
   hist_neutrino_pT_120 = book<TH1F>("hist_neutrino_pT_120","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, p_{T,gen}^{#nu}<120 GeV",100,-500,500);
   hist_neutrino_phi_120 = book<TH1F>("hist_neutrino_phi_120","#phi^{rec, #nu}-#phi^{gen, #nu}, p_{T,gen}^{#nu}<120 GeV",50,-TMath::Pi(),TMath::Pi());
   hist_neutrino_pT_120_220 = book<TH1F>("hist_neutrino_pT_120_220","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, 120 GeV < p_{T,gen}^{#nu} < 220 GeV",100,-500,500);
   hist_neutrino_phi_120_220 = book<TH1F>("hist_neutrino_phi_120_220","#phi^{rec, #nu}-#phi^{gen, #nu}, 120 GeV < p_{T,gen}^{#nu} < 220 GeV",50,-TMath::Pi(),TMath::Pi());
   hist_neutrino_pT_220_320 = book<TH1F>("hist_neutrino_pT_220_320","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, 220 GeV < p_{T,gen}^{#nu} < 320 GeV",100,-500,500);
   hist_neutrino_phi_220_320 = book<TH1F>("hist_neutrino_phi_220_320","#phi^{rec, #nu}-#phi^{gen, #nu}, 220 GeV < p_{T,gen}^{#nu} < 320 GeV",50,-TMath::Pi(),TMath::Pi());
   hist_neutrino_pT_320 = book<TH1F>("hist_neutrino_pT_320","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, p_{T,gen}^{#nu} > 320 GeV",100,-500,500);
   hist_neutrino_phi_320 = book<TH1F>("hist_neutrino_phi_320","#phi^{rec, #nu}-#phi^{gen, #nu}, p_{T,gen}^{#nu} > 320 GeV",50,-TMath::Pi(),TMath::Pi());

   //ptrec-ptgen binnend in MET
   hist_neutrino_pT_met160_320 = book<TH1F>("hist_neutrino_pT_met160_320","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, 160 GeV < met < 320 GeV",100,-500,500);
   hist_neutrino_phi_met160_320= book<TH1F>("hist_neutrino_phi_met160_320","#phi^{rec, #nu}-#phi^{gen, #nu}, 160 GeV < met < 320 GeV",50,-TMath::Pi(),TMath::Pi());
   hist_neutrino_pT_met320_500 = book<TH1F>("hist_neutrino_pT_met320_500","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, 320 GeV < met < 500 GeV",100,-500,500);
   hist_neutrino_phi_met320_500= book<TH1F>("hist_neutrino_phi_met320_500","#phi^{rec, #nu}-#phi^{gen, #nu}, 320 GeV < met < 500 GeV",50,-TMath::Pi(),TMath::Pi());
   hist_neutrino_pT_500= book<TH1F>("hist_neutrino_pT_met500","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, met > 500 GeV",100,-500,500);
   hist_neutrino_phi_500= book<TH1F>("hist_neutrino_phi_met500","#phi^{rec, #nu}-#phi^{gen, #nu}, met > 500 GeV",50,-TMath::Pi(),TMath::Pi());


   //2D plots
   hist_DeltaX_DeltaY = book<TH2F>("hist_DeltaX_DeltaY","hist_DeltaX_DeltaY",100,-500,500, 100,-500,500);
   hist_pX_pY_gen = book<TH2F>("hist_pX_pY_gen","hist_pX_pY_gen",100,-500,500, 100,-500,500);
   hist_pX_pY_rec = book<TH2F>("hist_pX_pY_rec","hist_pX_pY_rec",100,-500,500, 100,-500,500);
 

   //re-weight signal with Mchi =1 and Mmed=20 to ttbar
   hist_neutrino_gen_pT_reweighted=book<TH1F>("hist_neutrino_gen_pT_reweighted","p_{T,rec}^{#nu} - p_{T,gen}^{#nu}, reweighted",100,-500,500);
   hist_neutrino_pT_120_reweighted = book<TH1F>("hist_neutrino_pT_120_reweighted","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, p_{T,gen}^{#nu}<120 GeV, reweighted",100,-500,500);
   hist_neutrino_pT_120_220_reweighted = book<TH1F>("hist_neutrino_pT_120_220_reweighted","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, 120 GeV < p_{T,gen}^{#nu} < 220 GeV reweighted",100,-500,500);
   hist_neutrino_pT_220_320_reweighted = book<TH1F>("hist_neutrino_pT_220_320_reweighted","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, 220 GeV < p_{T,gen}^{#nu} < 320 GeV reweighted",100,-500,500);
   hist_neutrino_pT_320_reweighted = book<TH1F>("hist_neutrino_pT_320_reweighted","p_{T,rec}^{#nu}-p_{T,gen}^{#nu}, p_{T,gen}^{#nu} > 320 GeV reweighted",100,-500,500);

   hist_2d_DMMET_MET =book<TH2F>("hist_2d_DMMET_MET","MET", 250, 0, 5000, 250, 0, 5000);
   hist_2d_DMMET_MT2W=book<TH2F>(" hist_2d_DMMET_MT2W","MT2W", 250, 0, 5000, 50, 0, TMath::Pi());   
   hist_2d_DMMET_deltaphijetmet12=book<TH2F>("hist_2d_DMMET_deltaphijetmet12","#Delta#Phi(MET, jet12)", 250, 0, 5000, 50, 0, TMath::Pi());
   hist_2d_DMMET_deltaphijet123=book<TH2F>("hist_2d_DMMET_deltaphijet123","#Delta#Phi(MET, jet123)", 250, 0, 5000, 50, 0, TMath::Pi());
   hist_2d_DMMET_deltphilep=book<TH2F>("hist_2d_DMMET_deltphilep","#Delta#Phi(MET, lep)", 250, 0, 5000, 50, 0, TMath::Pi());
   hist_2d_DMMET_mtlep=book<TH2F>("hist_2d_DMMET_mtlep","m_{T}^{lep}", 250, 0, 5000, 50, 0, 500);

   hist_pTttbar=book<TH1F>("hist_pTttbar","p_T (t#bar{t})",100,0,1000);
   hist_deltaphi_tj_met=book<TH1F>("hist_deltaphi_tj_met","#Delta#phi(tagged jet, MET)",50,0,TMath::Pi());
   hist_deltaphi_tj_lep=book<TH1F>("hist_deltaphi_tj_lep","#Delta#phi(tagged jet, lep)",50,0,TMath::Pi());
   hist_deltaphi_neutrino_lep=book<TH1F>("hist_deltaphi_neutrino_lep","#Delta#phi(neutrino, lep)",50,0,TMath::Pi());
   hist_deltaphi_tj_neutrino=book<TH1F>("hist_deltaphi_tj_neutrino","#Delta#phi(tagged jet, neutrino)",50,0,TMath::Pi());
   hist_deltaphi_neutrino_met=book<TH1F>("hist_deltaphi_neutrino_met","#Delta#phi(neutrino,MET)",50,0,TMath::Pi());
   hist_deltaphi_thad_tlep=book<TH1F>("hist_deltaphi_thad_tlep","#Delta#phi(thad,tlep)",50,0,TMath::Pi());
   hist_deltaphi_thad_MET=book<TH1F>("hist_deltaphi_thad_MET","#Delta#phi(thad,met)",50,0,TMath::Pi());
}

void ttDMReconstructionHists_Likelihood::fill(const Event & event){
 
   double chi2 = event.get(h_likelihood);
   LorentzVector neutrino = event.get(h_recneutrino);
   std::vector<TopJet> heptopjets_WP3;
 
   if (event.is_valid(h_heptopjets_WP3)) heptopjets_WP3=event.get(h_heptopjets_WP3);
   Jet b_jet;
   if (event.is_valid(h_bjets))b_jet = event.get(h_bjets);
   hist_chi2->Fill(chi2, event.weight);
   hist_chi2_METpT_components->Fill(chi2 + TMath::Power((event.met->v4().Px()-neutrino.Px()),2)/TMath::Power(122.2,2) + TMath::Power((event.met->v4().Py()-neutrino.Py()),2)/TMath::Power(137.4,2), event.weight);
   
   double DM_MET = std::sqrt((event.met->v4().Px()-neutrino.Px())*(event.met->v4().Px()-neutrino.Px())+(event.met->v4().Py()-neutrino.Py())*(event.met->v4().Py()-neutrino.Py()));
   hist_DM_MET->Fill(DM_MET,event.weight);
      
   hist_neutrino_pT->Fill(neutrino.Pt(),event.weight);
   hist_neutrino_phi->Fill(TVector2::Phi_mpi_pi(neutrino.phi()),event.weight);
 
   hist_pTmisspT->Fill(event.met->v4().Pt()-neutrino.Pt(),event.weight);
   hist_pTmisspT_fabs->Fill(fabs(event.met->v4().Pt()-neutrino.Pt()),event.weight);
   hist_pTmisspT_quotient->Fill(event.met->v4().Pt()/neutrino.Pt(),event.weight);
   hist_2D_ptmiss_pT ->Fill(event.met->v4().Pt(),neutrino.Pt(),event.weight);
   hist_2D_METx_px ->Fill(event.met->v4().Px(),neutrino.Px(),event.weight);
   hist_2D_METy_py ->Fill(event.met->v4().Py(),neutrino.Py(),event.weight);

   double unc_px = METUncertainty_px(event.met);
   double unc_py = METUncertainty_py(event.met);
   
   hist_unc_px ->Fill(unc_px, event.weight);
   hist_unc_py ->Fill(unc_py, event.weight);
   hist_relunc_px ->Fill(unc_px*100/event.met->v4().Px(), event.weight);
   hist_relunc_py ->Fill(unc_py*100/event.met->v4().Py(), event.weight);
   

   double chi2_comparisonMetNu = (event.met->v4().Px() - neutrino.Px()) * (event.met->v4().Px() - neutrino.Px()) / (unc_px*unc_px) + (event.met->v4().Py() - neutrino.Py()) * (event.met->v4().Py() - neutrino.Py()) / (unc_py*unc_py);
   hist_DM_MET_chi2->Fill(chi2_comparisonMetNu, event.weight);
      
   hist_MET_RecNeutrino_px->Fill((event.met->v4().Px()-neutrino.Px())/event.met->v4().Px(), event.weight);
   hist_MET_RecNeutrino_py->Fill((event.met->v4().Py()-neutrino.Py())/event.met->v4().Py(), event.weight);
   hist_MET_RecNeutrino_px_abs->Fill((fabs(event.met->v4().Px()-neutrino.Px())/event.met->v4().Px()), event.weight);
   hist_MET_RecNeutrino_py_abs->Fill((fabs(event.met->v4().Py()-neutrino.Py())/event.met->v4().Py()), event.weight);
    
   hist_Metxpx  ->Fill((event.met->v4().Px()-neutrino.Px()),event.weight);
   hist_Metypy  ->Fill((event.met->v4().Py()-neutrino.Py()),event.weight);
   hist_Metxpx_quad ->Fill(TMath::Power((event.met->v4().Px()-neutrino.Px()),2),event.weight);
   hist_Metypy_quad ->Fill(TMath::Power((event.met->v4().Py()-neutrino.Py()),2),event.weight);
   hist_Metxpx_Metypy_quad ->Fill((TMath::Power((event.met->v4().Px()-neutrino.Px()),2) + TMath::Power((event.met->v4().Py()-neutrino.Py()),2)),event.weight);

   hist_2d_DMMET_MET ->Fill(DM_MET, event.met->pt(), event.weight); 
   Muon *lep1 = &event.muons->at(0);
   if (lep1 &&  event.jets->size()> 1) hist_2d_DMMET_MT2W->Fill(DM_MET,CalculateMT2W(event) , event.weight);

   double mindeltaphi = 9999;
   if(event.jets->size()>1) {
      for(size_t i=0; i<2; i++){
         double deltaphi = uhh2::deltaPhi(*event.met, event.jets->at(i));
         if (deltaphi < mindeltaphi) mindeltaphi=deltaphi;
      }
   }
   hist_2d_DMMET_deltaphijetmet12->Fill(DM_MET, mindeltaphi, event.weight); 

   mindeltaphi = 9999;
   if(event.jets->size()>2) {
      for(size_t i=0; i<3; i++){
         double deltaphi = uhh2::deltaPhi(*event.met, event.jets->at(i));
         if (deltaphi < mindeltaphi) mindeltaphi=deltaphi;
      }
   }
   hist_2d_DMMET_deltaphijet123->Fill(DM_MET, mindeltaphi, event.weight); 
   if (event.muons->size()) hist_2d_DMMET_deltphilep->Fill(DM_MET,uhh2::deltaPhi(*event.met, event.muons->at(0)) , event.weight); 
   if (lep1)hist_2d_DMMET_mtlep->Fill(DM_MET,sqrt(2*event.met->pt()*lep1->pt()*(1-cos(uhh2::deltaPhi(*event.met, *lep1)))) , event.weight);

   Muon lepton = event.muons->at(0); 
   double deltaphi_neutrino_lep = uhh2::deltaPhi(neutrino,lepton);
   hist_deltaphi_neutrino_lep->Fill(deltaphi_neutrino_lep, event.weight);
   double deltaphi_neutrino_met = uhh2::deltaPhi(neutrino,*event.met);
   hist_deltaphi_neutrino_met->Fill(deltaphi_neutrino_met, event.weight);

   if(heptopjets_WP3.size()>0)
      {
         TopJet taggedjet=heptopjets_WP3.at(0);
         LorentzVector ttbar = taggedjet.v4()+neutrino+lepton.v4()+b_jet.v4();
         hist_pTttbar->Fill(ttbar.pt(), event.weight);
         double deltaphi_tj_met = uhh2::deltaPhi(taggedjet,*event.met);
         hist_deltaphi_tj_met->Fill(deltaphi_tj_met, event.weight); 
         double deltaphi_tj_lep = uhh2::deltaPhi(taggedjet,lepton);
         hist_deltaphi_tj_lep->Fill(deltaphi_tj_lep, event.weight);
         double deltaphi_tj_neutrino = uhh2::deltaPhi(taggedjet,neutrino);
         hist_deltaphi_tj_neutrino->Fill(deltaphi_tj_neutrino, event.weight);
         double deltaphi_thad_tlep = uhh2::deltaPhi(taggedjet,neutrino+lepton.v4()+b_jet.v4());
         hist_deltaphi_thad_tlep->Fill(deltaphi_thad_tlep, event.weight);
         double deltaphi_thad_MET = uhh2::deltaPhi(taggedjet,*event.met);
         hist_deltaphi_thad_MET->Fill(deltaphi_thad_MET, event.weight);
      }
   
   if (event.isRealData) return;

   TTbarGen ttbargen = event.get(h_ttbargen);
    //find DM particle
   LorentzVector DM;
   for(unsigned int i=0; i<event.genparticles->size(); ++i) {
      const GenParticle & genp = event.genparticles->at(i);
      if (abs(genp.pdgId())==9100000) DM = genp.v4(); //mediator is stored
      if (genp.pdgId() == 9100022) { //DM particles are stored
         for(unsigned int j=0; j<event.genparticles->size(); ++j)
            {
               const GenParticle & genp2 = event.genparticles->at(j);
               if (genp2.pdgId() == -9100022) DM = genp.v4() + genp2.v4();
            }
      }
   }
   
   double x = DM_MET -fabs(DM.Pt());
   hist_DMMET_GenDMpT ->Fill(x, event.weight);
   double DMMET_GenDMpT_x = fabs(event.met->v4().Px()-neutrino.Px())-fabs(DM.Px());
   double DMMET_GenDMpT_y = fabs(event.met->v4().Py()-neutrino.Py())-fabs(DM.Py());
   
   LorentzVector DMMET_GenDMpT_vec = (event.met->v4() - neutrino - DM);


   hist_DMMET_GenDMpT_vec ->Fill( TMath::Sqrt(DMMET_GenDMpT_vec.Px()*DMMET_GenDMpT_vec.Px()+DMMET_GenDMpT_vec.Py()*DMMET_GenDMpT_vec.Py()), event.weight);
   hist_DMMET_GenDMpT_px ->Fill(DMMET_GenDMpT_x,event.weight);
   hist_DMMET_GenDMpT_py->Fill(DMMET_GenDMpT_y,event.weight);

   hist_MET_RecNeutrino_DM_px->Fill((event.met->v4().Px()-(neutrino.Px()+ DM.Px()))/event.met->v4().Px(), event.weight);
   hist_MET_RecNeutrino_DM_py->Fill((event.met->v4().Py()-(neutrino.Py() + DM.Py()))/event.met->v4().Py(), event.weight);
  
   LorentzVector met_gen;
   int n_nu = 0;
   for(unsigned int i=0; i<event.genparticles->size(); ++i) {
      const GenParticle & genp = event.genparticles->at(i);
      if (abs(genp.pdgId())==12 || abs(genp.pdgId())==14 || abs(genp.pdgId())==16) {
         met_gen = met_gen+genp.v4();
         n_nu = n_nu+1;
      }
   }
   met_gen = met_gen + DM;
   DM_MET = std::sqrt((met_gen.Px()-neutrino.Px())*(met_gen.Px()-neutrino.Px())+(met_gen.Py()-neutrino.Py())*(met_gen.Py()-neutrino.Py()));
   x = DM_MET -fabs(DM.Pt());
   hist_DMMET_GenDMpT_genmet -> Fill(x, event.weight);

    
   if(ttbargen.IsSemiLeptonicDecay())
      {
         hist_MET_Neutrino_px->Fill((event.met->v4().Px()-ttbargen.Neutrino().v4().Px())/ttbargen.Neutrino().v4().Px(), event.weight);
         hist_MET_Neutrino_py->Fill((event.met->v4().Py()-ttbargen.Neutrino().v4().Py())/ttbargen.Neutrino().v4().Py(), event.weight);
         
         hist_MET_Neutrino_DM_px->Fill((event.met->v4().Px()-(ttbargen.Neutrino().v4().Px() + DM.Px()))/(ttbargen.Neutrino().v4().Px()+ DM.Px()), event.weight);
         hist_MET_Neutrino_DM_py->Fill((event.met->v4().Py()-(ttbargen.Neutrino().v4().Py() + DM.Py()))/(ttbargen.Neutrino().v4().Py()+ DM.Py()), event.weight);
         
         hist_pxrec_pxgen->Fill(neutrino.Px()-ttbargen.Neutrino().v4().Px(),event.weight);
         hist_pyrec_pygen->Fill(neutrino.Py()-ttbargen.Neutrino().v4().Py(),event.weight);
         hist_pzrec_pzgen->Fill(neutrino.Pz()-ttbargen.Neutrino().v4().Pz(),event.weight);
         
         if (DM_MET>100) hist_pxrec_pxgen_100->Fill(neutrino.Px()/ttbargen.Neutrino().v4().Px(),event.weight);
         if (DM_MET>100) hist_pyrec_pygen_100->Fill(neutrino.Py()/ttbargen.Neutrino().v4().Py(),event.weight);

         hist_DM_MET_chi2_semilept->Fill(chi2_comparisonMetNu, event.weight);
         hist_DM_MET_semilept->Fill(DM_MET,  event.weight);
         hist_DMMET_GenDMpT_semilept ->Fill(x, event.weight);
         
         double DM_MET_gen = std::sqrt((event.met->v4().Px()-ttbargen.Neutrino().v4().Px())*(event.met->v4().Px()-ttbargen.Neutrino().v4().Px())+(event.met->v4().Py()-ttbargen.Neutrino().v4().Py())*(event.met->v4().Py()-ttbargen.Neutrino().v4().Py()));
         hist_DM_MET_gen->Fill(DM_MET_gen,  event.weight);
         
         Particle lepton;
         if (event.muons->size()>0) lepton = event.muons->at(0);
         else lepton = event.electrons->at(0);
         Jet bjet = *nextJet(lepton,*event.jets);
         GenParticle genb = ttbargen.BLep();
         double DeltaR = deltaR(genb, bjet);
         hist_DeltaR_genb_nextjet ->Fill(DeltaR, event.weight);
         if (chi2<15) hist_DeltaR_genb_nextjet_lowchi2 ->Fill(DeltaR, event.weight);
         else hist_DeltaR_genb_nextjet_highchi2->Fill(DeltaR, event.weight);
      
         x = DM_MET_gen - fabs(DM.Pt());
         hist_DMMETgen_GenDMpT ->Fill(x, event.weight);
         
         DM_MET_gen = std::sqrt((met_gen.Px()-ttbargen.Neutrino().v4().Px())*(met_gen.Px()-ttbargen.Neutrino().v4().Px())+(met_gen.Py()-ttbargen.Neutrino().v4().Py())*(met_gen.Py()-ttbargen.Neutrino().v4().Py()));
         x = DM_MET_gen - fabs(DM.Pt());
         hist_DMMETgen_GenDMpT_only ->Fill(x, event.weight);
         hist_n_nu -> Fill(n_nu, event.weight);

         hist_neutrino_gen_pT->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
         hist_neutrino_gen_phi->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         hist_neutrino_pT_gen->Fill(ttbargen.Neutrino().v4().Pt(),event.weight);
         hist_neutrino_phi_gen->Fill(TVector2::Phi_mpi_pi(ttbargen.Neutrino().v4().Phi()),event.weight);

         if (ttbargen.Neutrino().v4().Pt() < 120) {
            hist_neutrino_pT_120->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
            hist_neutrino_phi_120->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         }
         if (ttbargen.Neutrino().v4().Pt() >= 120 && neutrino.Pt() <220) {
            hist_neutrino_pT_120_220->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
            hist_neutrino_phi_120_220->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         }
         if (ttbargen.Neutrino().v4().Pt() >= 220 && neutrino.Pt() <320) {
            hist_neutrino_pT_220_320->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
            hist_neutrino_phi_220_320->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         }
         if (ttbargen.Neutrino().v4().Pt() >= 320) {
            hist_neutrino_pT_320->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
            hist_neutrino_phi_320->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         }

         if (event.met->pt()>160 && event.met->pt()<=320){
            hist_neutrino_pT_met160_320 ->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
            hist_neutrino_phi_met160_320->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         }
         if (event.met->pt()>320 && event.met->pt()<=500){
            hist_neutrino_pT_met320_500  ->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
            hist_neutrino_phi_met320_500->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         }
         if (event.met->pt()>=500){
            hist_neutrino_pT_500->Fill(neutrino.Pt()-ttbargen.Neutrino().v4().Pt(),event.weight);
            hist_neutrino_phi_500->Fill(TVector2::Phi_mpi_pi(neutrino.phi()-ttbargen.Neutrino().v4().Phi()),event.weight);
         }

         hist_DeltaX_DeltaY ->Fill(neutrino.Px()-ttbargen.Neutrino().v4().Px(),neutrino.Py()-ttbargen.Neutrino().v4().Py(), event.weight);
         hist_pX_pY_gen ->Fill(ttbargen.Neutrino().v4().Px(),ttbargen.Neutrino().v4().Py(), event.weight);
         hist_pX_pY_rec ->Fill(neutrino.Px(),neutrino.Py(),event.weight);

         //re-weight pT gen signal to bg, does not work properly yet
         
         // static TFile *f =new TFile("/nfs/dust/cms/user/mameyer/TTDMDM/selection_MUON_RunII_25ns_76X_AfterLikelihood/WithLikelihood/ReweightingSignal.root", "READ"); 
         // TH1 *ratio = (TH1*)f->Get("hist_neutrino_pT_gen");
         // double weight_pT = ratio->GetBinContent(ratio->FindBin(ttbargen.Neutrino().v4().Pt()));
         
         // hist_neutrino_gen_pT_reweighted->Fill((neutrino.Pt()-ttbargen.Neutrino().v4().Pt())*weight_pT,event.weight);
         // if (ttbargen.Neutrino().v4().Pt() < 120)  hist_neutrino_pT_120_reweighted->Fill((neutrino.Pt()-ttbargen.Neutrino().v4().Pt())*weight_pT,event.weight);
         // if (ttbargen.Neutrino().v4().Pt() >= 120 && neutrino.Pt() <220) hist_neutrino_pT_120_220_reweighted->Fill((neutrino.Pt()-ttbargen.Neutrino().v4().Pt())*weight_pT,event.weight);
         // if (ttbargen.Neutrino().v4().Pt() >= 220 && neutrino.Pt() <320) hist_neutrino_pT_220_320_reweighted->Fill((neutrino.Pt()-ttbargen.Neutrino().v4().Pt())*weight_pT,event.weight);
         // if (ttbargen.Neutrino().v4().Pt() >= 320) hist_neutrino_pT_320_reweighted->Fill((neutrino.Pt()-ttbargen.Neutrino().v4().Pt())*weight_pT,event.weight);
      }

   

   return;
}

double METUncertainty_px(MET *met)
{
   double uncertainty_up_px[6] = {met->shiftedPx_JetEnUp(), met->shiftedPx_JetResUp(),  met->shiftedPx_UnclusteredEnUp(), met->shiftedPx_ElectronEnUp(),  met->shiftedPx_TauEnUp(), met->shiftedPx_MuonEnUp()}; 
   double uncertainty_down_px[6] = {met->shiftedPx_JetEnDown(), met->shiftedPx_JetResDown(),  met->shiftedPx_UnclusteredEnDown(), met->shiftedPx_ElectronEnDown(),  met->shiftedPx_TauEnDown(), met->shiftedPx_MuonEnDown()}; 

   double total_unc =0;
   for (int i = 0;i<6;i++){   
      double unc = CalcUncertainty(met, true, uncertainty_up_px[i], uncertainty_down_px[i]);
      total_unc = total_unc+ unc*unc; 
   }
   total_unc = sqrt(total_unc);
   return total_unc;
}   

double METUncertainty_py(MET *met)
{
   double uncertainty_up_py[6] = {met->shiftedPy_JetEnUp(), met->shiftedPy_JetResUp(),  met->shiftedPy_UnclusteredEnUp(), met->shiftedPy_ElectronEnUp(),  met->shiftedPy_TauEnUp(), met->shiftedPy_MuonEnUp()}; 
   double uncertainty_down_py[6] = {met->shiftedPy_JetEnDown(), met->shiftedPy_JetResDown(),  met->shiftedPy_UnclusteredEnDown(), met->shiftedPy_ElectronEnDown(),  met->shiftedPy_TauEnDown(), met->shiftedPy_MuonEnDown()};  

   double total_unc =0;
   for (int i = 0;i<6;i++){   
      double unc = CalcUncertainty(met, false, uncertainty_up_py[i], uncertainty_down_py[i]);
      total_unc = total_unc+ unc*unc; 
   }
   total_unc = sqrt(total_unc);
   return total_unc;
}

double CalcUncertainty(MET *met, bool px, double up, double down)
{
   double unc =0;
   double met_unc_up;
   double met_unc_down;
   if (px){
      met_unc_up = fabs(up - met->v4().Px());
      met_unc_down = fabs(met->v4().Px() -down);
   }
   else {
      met_unc_up = fabs(up - met->v4().Py());
      met_unc_down = fabs(met->v4().Py()-down);
   }
   if (met_unc_up >= met_unc_down ) unc=met_unc_up;
   else unc = met_unc_down;
   return unc;
}
