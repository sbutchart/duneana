 ////////////////////////////////////////////////////////////////////////
// Class:       SupernovaAna
// Module Type: analyzer
// File:        SupernovaAna_module.cc
//
// Generated at Mon Jul 11 21:36:48 2016 by Michael Baird using the old
// copy and paste...
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <string>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dunecore/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dunecore/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"

class SupernovaAna2;

class SupernovaAna2 : public art::EDAnalyzer {

public:

  explicit SupernovaAna2(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  SupernovaAna2(SupernovaAna2 const &) = delete;
  SupernovaAna2(SupernovaAna2 &&) = delete;
  SupernovaAna2 & operator = (SupernovaAna2 const &) = delete;
  SupernovaAna2 & operator = (SupernovaAna2 &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;

  void endJob() override;

private:
 
  std::string fTruthLabel; 
  std::string fBgLabel;

  TCanvas *MarleyEC;
  TCanvas *MarleyNuEC;
  TCanvas *BgEC;
  TCanvas *BgTC;
  TCanvas *BgTrackXC;
  TCanvas *BgTrackYC;
  TCanvas *BgTrackZC;


  TH1F *MarleyE;
  TH1F *MarleyNuE;
  TH1F *BgE;
  TH1F *BgT;

  TCanvas *Bg_PosXC;
  TH2F	  *Bg_PosX;

  TCanvas *Bg_PosZC;
  TH2F	  *Bg_PosZ;

  TGraph *BgTrackX;
  
  std::vector<float> BgEs;
  std::vector<float>Bg_PosXs;
  std::vector<float>Bg_PosZs;
  std::vector<float>Bg_PosYs;
  std::vector<float>BgTs;
  std::vector<float> Bg_TrackTs;
  std::vector<float> Bg_TrackXs;
  std::vector<float> Bg_TrackYs;
  std::vector<float> Bg_TrackZs;

};



//......................................................
SupernovaAna2::SupernovaAna2(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void SupernovaAna2::reconfigure(fhicl::ParameterSet const & p)
{

  fTruthLabel = p.get<std::string> ("TruthLabel");
  fBgLabel   = p.get<std::string> ("BackgroundLabel");

}



//......................................................
void SupernovaAna2::beginJob()
{

  MarleyEC = new TCanvas ("MarleyEC", "", 1024, 768);
  MarleyNuEC = new TCanvas ("MarleyNuEC", "", 1024, 768);
  BgEC = new TCanvas ("BgEC", "", 1024, 768);
  BgTC = new TCanvas("BgTC", "", 1024, 768);
  BgTrackXC = new TCanvas("BgTrackXC", "", 1024, 768);

  MarleyE = new TH1F ("MarleyE", "", 100, 0, 0.05);
  MarleyNuE = new TH1F ("MarleyNuE", "", 100, 0, 0.07);  

  Bg_PosXC = new TCanvas("Bg_PosXC", "", 1024, 768);
  Bg_PosZC = new TCanvas("Bg_PosZC", "", 1024, 768);

}

//......................................................
void SupernovaAna2::endJob()
{

  //Make plot of supernova electron energies
  MarleyEC->cd();
  MarleyE->SetTitle("Marley events: e energy");
  MarleyE->GetXaxis()->SetTitle("Energy [GeV]");
  MarleyE->GetYaxis()->SetTitle("count");
  MarleyE->Draw();
  MarleyEC->SaveAs("MarleyE.png");

  //Make plot of supernova nu energies
  MarleyNuEC->cd();
  MarleyNuE->SetTitle("Marley events: Nu energy");
  MarleyNuE->GetXaxis()->SetTitle("Energy [GeV]");
  MarleyNuE->GetYaxis()->SetTitle("count");
  MarleyNuE->Draw();
  MarleyNuEC->SaveAs("MarleyNuE.png");

  //Fill histogram with background energies (all particles)
  BgE = new TH1F ("BgE", "", 150, *min_element( BgEs.begin(), BgEs.end() )-0.0001, *max_element( BgEs.begin(), BgEs.end() )+0.0001);
  for(int i = 0; i < (int)BgEs.size(); i++)
  {
    BgE->Fill( BgEs[i] );
  }

  //Make plot of selected background energies (all particles)
  BgEC->cd();
  std::stringstream fEtitle;
  fEtitle << "Background events: " << fBgLabel << " energy";
  BgE->SetTitle(fEtitle.str().c_str());
  BgE->GetXaxis()->SetTitle("Energy [GeV]");
  BgE->GetYaxis()->SetTitle("count");
  BgE->GetYaxis()->SetRange(0, 20000000);  
  BgE->GetXaxis()->SetMaxDigits(2);
  BgE->Draw();
  std::stringstream fEname;
  fEname << "BgE_" << fBgLabel << ".png"; 
  BgEC->SaveAs(fEname.str().c_str());


  //Make and fill histogram with background counts (all particles)
  BgT = new TH1F("BgT", "", 150, *min_element( BgTs.begin(), BgTs.end() )-0.0001, *max_element( BgTs.begin(), BgTs.end() )+0.0001);
  for(int i = 0; i < (int)BgTs.size(); i++)
  {
    BgT->Fill( BgTs[i] );
  }

  BgTC->cd();
  std::stringstream fTtitle;
  fTtitle << "Background events: " << fBgLabel << " count rate";
  BgT->SetTitle(fTtitle.str().c_str());
  BgT->GetXaxis()->SetTitle("Time");
  BgT->GetYaxis()->SetTitle("count");
  BgT->GetXaxis()->SetMaxDigits(2);
  BgT->Draw();
  std::stringstream fTname;
  fTname << "BgT_" << fBgLabel << ".png";
  BgTC->SaveAs(fTname.str().c_str());


 

  //TGraph *BgTrackX = new TGraph (n, &Bg_TrackTs[0], &Bg_TrackXs[0]);
  //TGraph *BgTrackY = new TGraph (n, &Bg_TrackTs[0], &Bg_TrackYs[0]);
  //TGraph *BgTrackZ = new TGraph (n, &Bg_TrackTs[0], &Bg_TrackZs[0]);

  BgTrackX = new TGraph();
  for(int i = 0; i < (int)Bg_TrackXs.size(); i++)
  {
    BgTrackX->SetPoint( i, Bg_TrackTs[i], Bg_TrackXs[i] );
    std::cout << i << ", " << Bg_TrackTs[i] << ", " << Bg_TrackXs[i] << std::endl;
  }

  BgTrackXC->cd();
  std::stringstream fTXtitle;
  fTXtitle << "Background event: " << fBgLabel << " X track";
  BgTrackX->SetTitle(fTXtitle.str().c_str());
  BgTrackX->GetXaxis()->SetTitle("Time");
  BgTrackX->GetYaxis()->SetTitle("X position");
  BgTrackX->GetXaxis()->SetMaxDigits(2);
  BgTrackX->Draw("L");
  std::stringstream fTXname;
  fTXname << "BgTrackX_" << fBgLabel << ".png";
  BgTrackXC->SaveAs(fTXname.str().c_str());


  //Bg_PosX = new TH2F ("Bg_PosX", "", 100, *min_element( Bg_PosXs.begin(), Bg_PosXs.end() )-10, *max_element( Bg_PosXs.begin(), Bg_PosXs.end() )+10, 100, *min_element( Bg_PosYs.begin(), Bg_PosYs.end() )-10, *max_element( Bg_PosYs.begin(), Bg_PosYs.end() )+10);
  Bg_PosX = new TH2F ("Bg_PosX", "", 100, -360, 360, 100, -700, 700);
  for(int i = 0; i < (int)Bg_PosYs.size(); i++)
  {
    Bg_PosX->Fill( Bg_PosXs[i], Bg_PosYs[i] );
  }

  //Bg_PosZ = new TH2F ("Bg_PosZ", "", 100, *min_element( Bg_PosZs.begin(), Bg_PosZs.end() )-10, *max_element( Bg_PosZs.begin(), Bg_PosZs.end() )+10, 100, *min_element( Bg_PosYs.begin(), Bg_PosYs.end() )-10, *max_element( Bg_PosYs.begin(), Bg_PosYs.end() )+10);
  Bg_PosZ = new TH2F ("Bg_PosZ", "", 100, -100, 1500, 100, -700, 700);
  for(int i = 0; i < (int)Bg_PosYs.size(); i++)
  {
    Bg_PosZ->Fill( Bg_PosZs[i], Bg_PosYs[i] );
  }

  Bg_PosXC->cd();
  Bg_PosX->SetTitle("Bg Position XY");
  Bg_PosX->GetXaxis()->SetTitle("StartX");
  Bg_PosX->GetYaxis()->SetTitle("StartY");
  Bg_PosX->Draw("COLZP");
  Bg_PosXC->SaveAs("Bg_Position_histogram_XY.png");
  Bg_PosXC->Close();

  Bg_PosZC->cd();
  Bg_PosZ->SetTitle("Bg Position ZY");
  Bg_PosZ->GetXaxis()->SetTitle("StartZ");
  Bg_PosZ->GetYaxis()->SetTitle("StartY");
  Bg_PosZ->Draw("COLZP");
  Bg_PosZC->SaveAs("Bg_Position_histogram_ZY.png");
  Bg_PosZC->Close();


}

//......................................................
void SupernovaAna2::analyze(art::Event const & e)
{

  auto MarlTrue = e.getHandle< std::vector< simb::MCTruth > >(fTruthLabel);
  if (MarlTrue) {
    //std::cout << "MARLEY event!" << std::endl;
    for(size_t i = 0; i < MarlTrue->size(); i++)
    {
      //std::cout << "lepton E: " <<  MarlTrue->at(i).GetNeutrino().Lepton().E() << std::endl;
      //std::cout << "Nu E: " <<  MarlTrue->at(i).GetNeutrino().Nu().E() << std::endl;
      MarleyE->Fill( MarlTrue->at(i).GetNeutrino().Lepton().E() - MarlTrue->at(i).GetNeutrino().Lepton().Mass() );
      MarleyNuE->Fill( MarlTrue->at(i).GetNeutrino().Nu().E() );
    }
  }

  auto BgTrue = e.getHandle< std::vector<simb::MCTruth> >(fBgLabel);
  if (BgTrue) {
    //std::cout << "BG event!" << std::endl;  
    for(size_t i = 0; i < BgTrue->size(); i++)
    {
      for(int j = 0; j < BgTrue->at(i).NParticles(); j++)
      {       
        //std::cout << "E: " <<  BgTrue->at(i).GetParticle(j).E() << std::endl;
        BgEs.push_back( BgTrue->at(i).GetParticle(j).E() - BgTrue->at(i).GetParticle(j).Mass() );        
        BgTs.push_back( BgTrue->at(i).GetParticle(j).T() );
        
        if (i==0) {
          Bg_TrackTs.push_back( BgTrue->at(i).GetParticle(j).T() );
          Bg_TrackXs.push_back( BgTrue->at(i).GetParticle(j).Vx());
          Bg_TrackYs.push_back( BgTrue->at(i).GetParticle(j).Vy());
          Bg_TrackZs.push_back( BgTrue->at(i).GetParticle(j).Vz());
        }

	if (j==0) {
	  Bg_PosXs.push_back( BgTrue->at(i).GetParticle(j).Vx());
	  Bg_PosZs.push_back( BgTrue->at(i).GetParticle(j).Vz());

	  Bg_PosYs.push_back( BgTrue->at(i).GetParticle(j).Vy());

	}
      }
    }
  }

  std::cout << " " << std::endl;
}
DEFINE_ART_MODULE(SupernovaAna2)
