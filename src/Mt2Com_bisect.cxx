//#include "EventReco/interface/Mt2Com_bisect.h"
#include "UHH2/TTbarDM/include/Mt2Com_bisect.h"

Mt2Com_bisect::Mt2Com_bisect(){}
Mt2Com_bisect::~Mt2Com_bisect(){}

//###########################################
//    Global Functions
//###########################################

// double Mt2Com_bisect::calculateMT2w(const std::vector<IPHCTree::NTJet>& jets,
//                                     const std::vector<IPHCTree::NTJet>& bjets,
//                                     const std::vector<IPHCTree::NTMuon>& muon,
//                                     const std::vector<IPHCTree::NTElectron>& electron,
//                                     TVector2 vecMet,
//                                     string mt2type)
// {

//     TLorentzVector tmpLepton;

//          if (muon.size()     == 0) tmpLepton = electron[0].p4;
//     else if (electron.size() == 0) tmpLepton = muon[0].p4;


//     return Mt2Com_bisect::calculateMT2w(jets, bjets, tmpLepton, vecMet, mt2type);
// }


double Mt2Com_bisect::calculateMT2w(const std::vector<TLorentzVector>& jets,
				    const std::vector<TLorentzVector>& bjets,
				    const TLorentzVector leptonInput,
				    TVector2 vecMet,
				    string mt2type)
 {

// double Mt2Com_bisect::calculateMT2w(const std::vector<IPHCTree::NTJet>& jets,
//                                     const std::vector<IPHCTree::NTJet>& bjets,
//                                     const TLorentzVector leptonInput,
//                                     TVector2 vecMet,
//                                     string mt2type)

// {

    //Variables
    metval = vecMet.Mod();
    metphi = vecMet.Phi();
    lep = leptonInput;

    // I am asumming that jets is sorted by Pt
    //assert ( jets.size() == btag.size() );
    // require at least 2 jets
    if ( jets.size()<2 ) return 99999.;

    n_btag = (int) bjets.size();

    // We do different things depending on the number of b-tagged jets
    // arXiv:1203.4813 recipe

    nMax=3;
    ///nMax=2; //try 2

    int n_jets = int(jets.size());

//     if(entry_ == 58469) {
//       cout<<"check all jets/bjet " << jets.size()<<" "<<bjets.size()<<endl;
//       for(int j=0; j<int(jets.size()); j++){
// 	cout<<" jets " << jets[j].X()<<" "<< jets[j].Y()<<" "<<jets[j].Z()<<" "<<jets[j].E()<<endl;
//       }
//       for(int j=0; j<int(bjets.size()); j++){
//         cout<<" jets " << bjets[j].X()<<" "<< bjets[j].Y()<<" "<<bjets[j].Z()<<" "<<bjets[j].E()<<endl;
//       }

//     }


    // -------
    // 0 bTag
    // -------

    if (n_btag == 0)
    {
      // If no b-jets select the minimum of the mt2w from all combinations with
      // the three leading jets
      float min_mt2w = 9999;

      for (int i=0; i<nMax && i < n_jets; i++)  //fixed by Yang (it<n_jets)
	for (int j=0; j<nMax && j< n_jets; j++)
	  {
          if (i == j) continue;
          float c_mt2w = Mt2Com_bisect::mt2wWrapper(lep, jets[i], jets[j], metval, metphi, mt2type);
          if (c_mt2w < min_mt2w)
          min_mt2w = c_mt2w;
      }
      return min_mt2w;
    }
    // -------
    // 1 bTag
    // -------
    else if (n_btag == 1 )
    {
      // if only one b-jet choose the three non-b leading jets and choose the smaller
      float min_mt2w = 9999;
      int it(0);
      int ctr(0);

      while (ctr<nMax && it < n_jets ) //fixed by Yang (it<n_jets)
	{
	if (bjets[0]!=jets[it])
	  {
	    ///if(entry_ == 58469) cout<<"it " << it <<endl;

	    float c_mt2w = Mt2Com_bisect::mt2wWrapper(lep, bjets[0], jets[it], metval, metphi, mt2type);
	    if (c_mt2w < min_mt2w) min_mt2w = c_mt2w;
	    ctr ++;
	    it ++;
	  }
	else it ++;
      }
      it=0;
      ctr=0;
      while (ctr<nMax && it < n_jets )//fixed by Yang (it<n_jets)
      {
            if (bjets[0]!=jets[it])
            {
	      ///if(entry_ == 58469) cout<<"it2 " << it <<endl;

	      float c_mt2w = Mt2Com_bisect::mt2wWrapper(lep, jets[it], bjets[0], metval, metphi, mt2type);
                if (c_mt2w < min_mt2w) min_mt2w = c_mt2w;
                ctr ++;
                it ++;
            }
            else it ++;
      }
      return min_mt2w;
    }
    // -------
    // 2 bTag
    // -------
    else if (n_btag >= 2)
    {
      // if 3 or more b-jets the paper says ignore b-tag and do like 0-bjets
      // but we are going to make the combinations with the b-jets
      float min_mt2w = 9999;
      for (int i=0; i<n_btag; i++)
        for (int j=0; j<n_btag; j++)
        {
          if (i == j) continue;
          float c_mt2w = Mt2Com_bisect::mt2wWrapper(lep, bjets[i], bjets[j], metval, metphi, mt2type);
	  //	  cout<<" c_mt2w " << c_mt2w <<" "<<metphi<<endl;

          if (c_mt2w < min_mt2w)
          min_mt2w = c_mt2w;
        }
      return min_mt2w;
    }

    return -1.;
}




// This funcion is a wrapper for mt2w_bisect etc that takes LorentzVectors instead of doubles
double Mt2Com_bisect::mt2wWrapper(TLorentzVector& lep, const TLorentzVector& jet_o,
                                                                    const TLorentzVector& jet_b, float& metval, float& metphi,
                                                                    string& mt2type){

    // same for all MT2x variables
    metx = metval * cos( metphi );
    mety = metval * sin( metphi );


//     if(entry_ == 58469) {
//       cout<<" metx/y " << metx <<" "<<mety <<endl;
//       cout<<" mt2wWrapper  jet_b " <<  jet_b.X()<<" "<<jet_b.Y()<<" "<<jet_b.Z()<<" "<<jet_b.E()<<endl;
//     }


    pl[0]= lep.E(); pl[1]= lep.Px(); pl[2]= lep.Py(); pl[3]= lep.Pz();
    pb1[1] = jet_o.Px();  pb1[2] = jet_o.Py();   pb1[3] = jet_o.Pz();
    pb2[1] = jet_b.Px();  pb2[2] = jet_b.Py();   pb2[3] = jet_b.Pz();
    pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

    // specifics for each variable
    if (mt2type == "MT2b")
    {
      pmiss_lep[0] = 0.;
      pmiss_lep[1] = pmiss[1]+pl[1]; pmiss_lep[2] = pmiss[2]+pl[2];

      pb1[0] = jet_o.M();
      pb2[0] = jet_b.M();

      //mt2_bisect::mt2 mt2_event;
      mt2::set_momenta( pb1, pb2, pmiss_lep );
      mt2::set_mn( 80.385 );   // Invisible particle mass == W mass
      return mt2::get_mt2();
    }

    else
    {
      pb1[0] = jet_o.E();
      pb2[0] = jet_b.E();

      if (mt2type == "MT2bl")
      {

	//	if(entry_ == 58469) cout<<" check "<< entry_ <<endl;

        mt2bl::set_momenta(pl, pb1, pb2, pmiss);
        return mt2bl::get_mt2bl();
      }
      else if (mt2type == "MT2w")
      {
        mt2w::set_momenta(pl, pb1, pb2, pmiss);
        return mt2w::get_mt2w();
      }
    }

    return -1.;
}

