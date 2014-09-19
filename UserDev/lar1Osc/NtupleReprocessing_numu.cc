//
//   Author: Corey Adams, Yale (corey.adams@yale.edu), in the past
//   Stolen and broken by Joseph Zennamo, UChicago (jzennamo@uchicago.edu), Aug 2014
//
//
 
#define NtupleReprocessing_numu_cxx
#include "NtupleReprocessing_numu.hh"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>

#define print(a) ( std::cout << a << std::endl)

namespace lar1{

  void NtupleReprocessing_numu::Loop( std::string signal,
                                 Int_t iDet,
                                 Long64_t max_entry,
                                 bool verbose,
                                 double scale ){

    // Signal options are:
    // numu - muon neutrino disappearance 
    // 
    // Will want to add NC background in the future

    //---------------------------------------------
    // some things to be configured by the user:

    // Check to make sure the signal is one of the above:
    if (signal != "numu"){
      std::cout << "Error: incorrect signal \"" << signal << "\"  used." << std::endl;
      std::cout << "This Only does numu's for the time being..." << std::endl;
      return;
    }


    double muonIDeff     = 0.8;
    Double_t prot_thresh = 0.02;

    Double_t detect_dist = 0;   // 10000=ND, 47000=MicroBooNE, 70000=FD
    if (iDet == 0) detect_dist = 10000.0;
    else if (iDet == 1) detect_dist = 47000.0;
    else if (iDet == 2) detect_dist = 70000.0;
    else if (iDet == 3) detect_dist = 70000.0; // 3 is MB using FD monte carlo
    else if (iDet == 4) detect_dist = 70000.0; // 4 is IC using FD monte carlo
    else if (iDet == 6) detect_dist = 60000.0; // 4 is IC @ 600m using FD monte carlo
    else if (iDet == 7) detect_dist = 80000.0; // 4 is IC @ 600m using FD monte carlo
    else if (iDet == 8) detect_dist = 15000.0; // alternative near detector locations
    else if (iDet == 9) detect_dist = 17500.0; // alternative near detector locations
    else if (iDet == 10) detect_dist = 20000.0; //alternative near detector locations

    //---------------------------------------------

    //   In a ROOT session, you can do:
    //      Root > .L NtupleReprocessing_numu.C+ (+ to compile - this is necessary)
    //      Root > NtupleReprocessing_numu t
    //      Root > t.GetEntry(12); // Fill t data members with entry number 12
    //      Root > t.Show();       // Show values of entry 12
    //      Root > t.Show(16);     // Read and show values of entry 16
    //      Root > t.Loop();       // Loop on all entries
    //
       
    /**
       In this section, I describe how this object creates ntuples that could be used in sensitivity analysis.
       
       Each event is examined and given a weight. We then build a variety of histograms which are saved and 
       passed to the sensitivity code. There are two main histograms:
       
       The null hypothesis histograms and the prediction histogram.

       The prediction takes into account part of the oscillation probability, completing part of the Raster 
       scan. 

       I will expand here about how we select events.

    **/

    if (fChain == 0) return;

    //----------------------------------
    // Create the output file
    //----------------------------------
    std::cout << "input file: " << InFile().Remove(InFile().Length()-5) << std::endl;
    // TString sample = "ccinc";
    // if( ccqe_only ) sample = "ccqe";
    TString outfile = InFile().Remove(InFile().Length()-5) + "_processed_";
    if (iDet == 3) outfile += "MB_";
    if (iDet == 4) outfile += "IC_";
    if (iDet == 6) outfile += "IC_600_";
    if (iDet == 7) outfile += "IC_800_";
    if (scale != 1.0) {
      outfile += "scale_";
      outfile += scale;
      outfile += "_";
    }
    if (specialNameText != "") outfile = outfile + specialNameText + "_";
    outfile = outfile + signal + ".root";

    std::cout << "outfile file: " << outfile << std::endl;
    TFile *f = new TFile(outfile, "RECREATE");
    if (f->IsZombie()){
      std::cout << "Error, couldn't create output file." << std::endl;
      return;
    }  

    double efficiency = 0.0;
    double fluxweight = 1.0;

    // If necessary, scale the far det:
    if (scale != 1.0) utils.ScaleFarDet(scale/100.0); 

    //Get the POT weight.  Need to get iflux
    b_iflux->GetEntry(0); //iflux should now be filled
    double potweight;
    potweight = utils.GetPOTNorm( iflux, iDet );


    const static int npoints = 500; //Number of grid points in Raster scan of parameter space
    Double_t dm2min = 0.01;                       //eV**2
    Double_t dm2max = 100.;                       //eV**2


    double emin = 0.0, emax = 3.0;
    double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};

    TH1D *numuCC = new TH1D("NumuCC","NumuCC;Reconstructed Neutrino Energy [GeV];Events",19, bins);
    TH1D *numuCC_truth = new TH1D("NumuCC_True","NumuCC;True Neutrino Energy [GeV];Events",19, bins);
    TH1D *numuCC_cont = new TH1D("NumuCC_Contained","NumuCC;Reconstructed Neutrino Energy [GeV];Events",19, bins);

    TH1D *numuCC_Osc[npoints+1];
    TH1D *numuCC_Osc_truth[npoints+1];
    TH1D *numuCC_Osc_cont[npoints+1];


    for(int i = 0; i <= npoints; i++){
      numuCC_Osc[i] = (TH1D*)numuCC->Clone();    
      numuCC_Osc_truth[i] = (TH1D*)numuCC_truth->Clone();
      numuCC_Osc_cont[i] = (TH1D*)numuCC_cont->Clone();

    }


    double mu_L = 0;
    double DMpoint;
    double mu_E = 0;
    double mu_Theta = 0;
    double smeared_mu_E = 0;

    bool isFid, isActive;

    TVector3 MuonExitPos(0,0,0), MuonInitialMom(0,0,0);

    double EnuReco = 0.;

    double nu_LE_wgt = 0;

    bool contained = true;
    int iTop = 0, iChan = 0;

    //====================================================
    // Loop over entries in incoming ntuple
    //====================================================
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    if( max_entry == -1 ) max_entry = nentries;

    for( Long64_t jentry=0; jentry<nentries && jentry<max_entry; jentry++ ){
    // for( Long64_t jentry=430000; jentry<nentries && jentry<max_entry; jentry++ ){

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb;
      iTop = 0;
      iChan = 0;

      if (verbose || jentry % 50000 == 0) std::cout << "========== Event " << ientry << " of " << nentries << " ==========" << std::endl;
      
      // Clear ntuple variables
      efficiency = 0.0;
      wgt = 0.0;

      // Calculate the FS lepton's 4-vector and return neutrino path length
      nuleng = CalcLepton( detect_dist );

      // Get flux weight from FluxRW utilities
      fluxweight = 1.0;

      // Lump in the Pot weight to flux weight so it isn't lost:
      fluxweight *= potweight;
        
      // Is the interaction in the fiducial volume?
      TVector3 vtx(Vx, Vy, Vz);
      isFid = utils.IsFiducial( iDet, vtx );
      isActive = utils.IsActive( iDet, vtx );

      if( !isActive ) continue;

      Double_t photon_energy = utils.TotalPhotonEnergy( iDet, p1PhotonConversionPos, p1PhotonConversionMom, 
                    p2PhotonConversionPos, p2PhotonConversionMom );

      nu_LE_wgt = (nuleng/100.)/(enugen*1000.);

      //------------------------------------------
      // Analyze CC events
      //------------------------------------------

      // electron CC events
      if( isFid && isCC && abs(inno) == 14){

      //----------------------------------------------
      // Track the muon till it exits the detector:
      //----------------------------------------------
      
      unsigned int muonPosIndex = 0;
      // Make sure default values are 0:
      MuonExitPos *= 0.0;

      while ( contained && muonPosIndex < MuonPos->size() ){
	TVector3 pos(MuonPos->at(muonPosIndex).X(), MuonPos->at(muonPosIndex).Y(), MuonPos->at(muonPosIndex).Z());
	contained = utils.IsActive(iDet,pos,5); //fiducial cut at 5cm
	muonPosIndex++;
      }
      
      // Take the data from the previous point.  It's now exiting the detector.
      // Or has stopped!
      MuonExitPos = TVector3(MuonPos->at(muonPosIndex-1).X(), MuonPos->at(muonPosIndex-1).Y(), MuonPos->at(muonPosIndex-1).Z());
      
      //Unsmeared Lepton Energy
      mu_E = MuonMom->at(0).E();
   
      //internal muon length
      mu_L = 0;
      mu_L = sqrt(pow((vtx.X()-MuonExitPos.X()),2)+pow((vtx.Y()-MuonExitPos.Y()),2)+pow((vtx.Z()-MuonExitPos.Z()),2));

      smeared_mu_E = utils.MuE_Res(contained, mu_E, mu_L);

      if(smeared_mu_E < 0) smeared_mu_E = 0;

      if(smeared_mu_E > 0){
	EnuReco = utils.NuEnergyCaloSmeared(  GeniePDG,
					      GenieE, false, false,
					      prot_thresh, false,
					      smeared_mu_E,
					      mu_E,
					      "muon")+photon_energy;
      }
      else{EnuReco = 0;}

     }// accepted neutrinos 
      
      
      numuCC->Fill(EnuReco,wgt);
      numuCC_truth->Fill(enugen, wgt);
      if(contained){numuCC_cont->Fill(EnuReco,wgt);}

      //Partial oscillation computation
      for(int dm2 = 0; dm2 <= npoints; dm2++){	  
	DMpoint = pow(10., (TMath::Log10(dm2min)+ (dm2 * (1./npoints))*TMath::Log10(dm2max/dm2min)) );
	numuCC_Osc[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));
	numuCC_Osc_truth[dm2]->Fill(enugen,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));
	if(contained){numuCC_Osc_cont[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));}
      }

    }// loop over all events 

    
    numuCC->Write();
    numuCC_truth->Write();
    numuCC_cont->Write();

    for(int dm2 = 0; dm2 <= npoints; dm2++){
      std::string dmpoint = std::to_string(dm2);
      std::string name = "Osc_";
      std::string true_name = "Osc_True_";
      std::string cont_name = "Osc_Cont_";

      name = name+dmpoint;
      true_name = true_name+dmpoint;
      cont_name = cont_name+dmpoint;

      numuCC_Osc[dm2]->Write(name.c_str());
      numuCC_Osc_truth[dm2]->Write(true_name.c_str());
      numuCC_Osc_cont[dm2]->Write(cont_name.c_str());

    }


    f->Write();
    f->Close();


  }   


}// namespace lar1

