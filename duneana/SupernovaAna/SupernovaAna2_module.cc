///////////////////////////////////////////////////////////////////////
// Class:       SupernovaAna
// Module Type: analyzer
// File:        SupernovaAna_module.cc
//
// Generated at Mon Jul 11 21:36:48 2016 by Michael Baird using the old
// copy and paste...
////////////////////////////////////////////////////////////////////////

// C++ includes

// ROOT includes
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTimeStamp.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/EventPrincipal.h"
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

#include "art/Framework/Core/EDProducer.h"

// DUNETPC specific includes
#include "dunecore/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dunecore/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SupernovaTruth.h"
#include "lardataobj/RecoBase/Hit.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/OutputModule.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "cetlib/column_width.h"
#include "cetlib/lpad.h"
#include "cetlib/rpad.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

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

  bool GenOnly;

  bool firstEv;
  uint64_t firstEvTime; 
  std::vector< std::string > labels; 
  
  //what information we want to read in from the events
  struct data {
    std::string                 label;
    std::vector <double>        E;
    std::vector <std::uint64_t> Time;
    std::vector <double>        SNTime;
    std::vector <double>        XPosn;
    std::vector <double>        YPosn;
    std::vector <double>        ZPosn;
  };
  std::vector < data > all_data;

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
  GenOnly = p.get<bool>("GenOnly",0);
}

//......................................................
void SupernovaAna2::beginJob()
{

  firstEv = true;
  firstEvTime = 0;

}

//......................................................
void SupernovaAna2::endJob()
{

  // Summary
  std::cout << "Final counts: " << std::endl;
  for (auto & data_label : all_data) {
    std::cout << data_label.label << " " << data_label.E.size() << std::endl;
  }

  art::ServiceHandle<art::TFileService> tfs;
  TTree* fEnergies;
  TTree* fXPosns;
  fEnergies = tfs->make<TTree>("Energies", "Energy plots");
  fXPosns = tfs->make<TTree>("XY_positions", "XY positions");

  // Make plots here
  for (auto & data_label : all_data) {

    // Energy
    if (data_label.E.size() > 0) {
      // prep plot
      TCanvas *EnergyC;
      TH1F *Energy;
      std::stringstream name;
      name << data_label.label << "_energy";
      EnergyC = new TCanvas (name.str().c_str(), "", 1024, 768);
      const auto [min, max] = std::minmax_element(begin(data_label.E), end(data_label.E));
      const auto scale = 0.1*(*max - *min);
      Energy = new TH1F (name.str().c_str(), "", 100, *min-scale, *max+scale); 
      // fill data
      for (auto & energy : data_label.E) {
        Energy->Fill( energy );
      }
      // plot and save
      EnergyC->cd();
      std::stringstream title;
      title << "Energy: " << data_label.label;
      Energy->SetTitle(title.str().c_str());
      Energy->GetXaxis()->SetTitle("Energy [GeV]");
      Energy->GetYaxis()->SetTitle("count");
      Energy->Draw();
      fEnergies->Branch(name.str().c_str(), &Energy, "Energy/I");
      Energy->Write();
      name << ".root";
      EnergyC->SaveAs(name.str().c_str());
      delete EnergyC; delete Energy;
    }
  
  
   // PosX
    if (data_label.XPosn.size() > 0) {
      TCanvas *PosXC;
      TH2F *PosX;
      std::stringstream name;
      name << data_label.label << "_posX";
      PosXC = new TCanvas (name.str().c_str(), "", 1024, 768);
      PosX = new TH2F (name.str().c_str(), "", 100, -400, 400, 100, -700, 700);

      // fill data
      for (int i = 0; i < (int)data_label.YPosn.size(); i++) {
      PosX->Fill(data_label.XPosn[i], data_label.YPosn[i]);
      }

      // plot and save
      PosXC->cd();
      std::stringstream title;
      title << "Position XY: " << data_label.label;
      PosX->SetTitle(title.str().c_str());
      PosX->GetXaxis()->SetTitle("X position");
      PosX->GetYaxis()->SetTitle("Y position");
      PosX->SetOption("COLZP");
      fXPosns->Branch(name.str().c_str(), &PosX, "PosX/I");
      PosX->Write();
      name << ".png";
      PosXC->SaveAs(name.str().c_str());
      delete PosXC; delete PosX;
    }


   //PosZ
    if (data_label.ZPosn.size() > 0) {
      TCanvas *PosZC;
      TH2F *PosZ;
      std::stringstream name;
      name << data_label.label << "_posZ";
      PosZC = new TCanvas (name.str().c_str(), "", 1024, 768);
      PosZ = new TH2F (name.str().c_str(), "", 100, -400, 1500, 100, -700, 700);

    //fill data
    for (int i = 0; i < (int)data_label.YPosn.size(); i++) {
    PosZ->Fill(data_label.ZPosn[i], data_label.YPosn[i]);
    }

    //plot and save
    PosZC->cd();
    std::stringstream title; 
    title << "Position ZY: " << data_label.label;
    PosZ->SetTitle(title.str().c_str());
    PosZ->GetXaxis()->SetTitle("Z position");
    PosZ->GetYaxis()->SetTitle("Y position");
    PosZ->SetOption("COLZP");
    //fSupernovaAnaTree->Branch("PosZ", &PosZ, "PosZ/I");
    PosZ->Write();
    name << ".png";
    PosZC->SaveAs(name.str().c_str());
    delete PosZC; delete PosZ;
  }



   // Time
    if (data_label.Time.size() > 0) {
      // prep plot
      TCanvas *TimeC;
      TH1F *Time; 
      std::stringstream name;
      name << data_label.label << "_time";
      TimeC = new TCanvas (name.str().c_str(), "", 1024, 768); 
      //for (auto & time : data_label.Time) { 
      const auto [min, max] = std::minmax_element(begin(data_label.Time), end(data_label.Time));
      const auto scale = 0.1*(*max - *min);
      
      Time = new TH1F (name.str().c_str(), "", 100, *min-scale, *max+scale);
      //Time = new TH1F (name.str().c_str(), "", 100, -10, 1500000000000);

      // fill data
      for (auto & time : data_label.Time) {
        Time->Fill( time );
      }
      // plot and save
      TimeC->cd();
      std::stringstream title;
      title << "Time: " << data_label.label;
      Time->SetTitle(title.str().c_str());
      Time->GetXaxis()->SetTitle("Time [ns]");
      Time->GetYaxis()->SetTitle("count");
      Time->Write();
      name << ".png";
      TimeC->SaveAs(name.str().c_str());
      delete TimeC; delete Time;
    }

    // SN Truth stuff
    if (data_label.label == "marley") {
      if (data_label.SNTime.size() > 0) {
        TCanvas *SNTimeC;
        TH1F *SNTime; 
        std::stringstream name;
        name << data_label.label << "_sntime";
        SNTimeC = new TCanvas (name.str().c_str(), "", 1024, 768); 
        const auto [min, max] = std::minmax_element(begin(data_label.SNTime), end(data_label.SNTime));
        const auto scale = 0.1*(*max - *min);
        SNTime = new TH1F (name.str().c_str(), "", 100, *min-scale, *max+scale);
        // fill data
        for (auto & sntime : data_label.SNTime) {
          SNTime->Fill( sntime );
        }
        // plot and save
        SNTimeC->cd();
        std::stringstream title;
        title << "SN Time: " << data_label.label;
        SNTime->SetTitle(title.str().c_str());
        SNTime->GetXaxis()->SetTitle("Time [s]");
        SNTime->GetYaxis()->SetTitle("count");
        SNTime->Draw();
        name << ".png";
        SNTimeC->SaveAs(name.str().c_str());
        delete SNTimeC; delete SNTime;
      }
    }
  } // label

} // end of run

//......................................................
void SupernovaAna2::analyze(art::Event const & e)
{

  // For first event, get labels, store them
  if (firstEv) {
 
    auto mcHandles = e.getMany<std::vector<simb::MCTruth>>();
    for (auto const& mcHandle : mcHandles) {
      const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
      labels.push_back( sModuleLabel );
    }

    // set up data    
    for (auto const& label : labels) {
      data temp;
      temp.label = label;
      all_data.push_back( temp ); 
    }
  
    firstEv = false;
    firstEvTime = e.time().value();

    std::cout << "Available labels: " << std::endl;
    for (auto & label : labels) {       
      std::cout << label << std::endl;
    }

    std::cout << "First Ev time: " << firstEvTime << std::endl;

  }

  if (GenOnly) {

    // Now for each event, get the type, and do work
    for (auto & label : labels) {
      if (label != "marley") {
        auto eventLabel = e.getHandle< std::vector< simb::MCTruth > >( label );
        if (eventLabel) {
          for(size_t i = 0; i < eventLabel->size(); i++) {
            for(int j = 0; j < eventLabel->at(i).NParticles(); j++) {
              for (auto & data_label : all_data) {
                if (data_label.label == label) {
                  data_label.E.push_back( eventLabel->at(i).GetParticle(j).E() - eventLabel->at(i).GetParticle(j).Mass() );

                  //exclude anomalous event times
                  if (( e.time().value() + static_cast<uint64_t>(eventLabel->at(i).GetParticle(j).T()) - firstEvTime ) < 1500000000000 ) {
	            data_label.Time.push_back( e.time().value() + static_cast<uint64_t>(eventLabel->at(i).GetParticle(j).T()) - firstEvTime );
                  }  
                 
                  data_label.XPosn.push_back ( eventLabel->at(i).GetParticle(j).Vx());
                  data_label.YPosn.push_back ( eventLabel->at(i).GetParticle(j).Vy());
                  data_label.ZPosn.push_back ( eventLabel->at(i).GetParticle(j).Vz());
  
                }
              }
            }
          }
        }
      } else { // marley
        auto eventLabel = e.getHandle< std::vector< simb::MCTruth > >( label );
        if (eventLabel) {
          for  (size_t i = 0; i < eventLabel->size(); i++) {
            for (auto & data_label : all_data) {
              if (data_label.label == label) {
                data_label.E.push_back( eventLabel->at(i).GetNeutrino().Nu().E() - eventLabel->at(i).GetNeutrino().Nu().Mass() );
                data_label.Time.push_back( e.time().value() - firstEvTime );
                
                // SN Truth
                art::FindManyP<sim::SupernovaTruth> SNTruth(eventLabel, e, label);
                for (size_t j = 0; j < SNTruth.at(i).size(); j++) {
                  const sim::SupernovaTruth ThisTr = (*SNTruth.at(i).at(j));
                  data_label.SNTime.push_back( ThisTr.SupernovaTime );
                }
              }
            }      
          }
        }
      }
    }

  } else { // gen only

    // Now for each event, get the type, and do work
    for (auto & label : labels) {
      if (label != "marley") {
        auto eventLabel = e.getHandle< std::vector< simb::MCTruth > >( label );
        art::FindManyP<simb::MCParticle> Assn(eventLabel,e,"largeant");
        if (eventLabel) {
          for ( size_t L1=0; L1 < eventLabel->size(); ++L1 ) {
            for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
              for (auto & data_label : all_data) {
                if (data_label.label == label) {
                  const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
                  data_label.E.push_back( ThisPar.E() - ThisPar.Mass() );
                  data_label.Time.push_back( e.time().value() + static_cast<uint64_t>(ThisPar.T()) - firstEvTime );
                  data_label.XPosn.push_back ( ThisPar.Vx());
                  data_label.YPosn.push_back ( ThisPar.Vy());
                  data_label.ZPosn.push_back ( ThisPar.Vz());
                }
              }
            }
          }
        }
      } else { // marley
        auto eventLabel = e.getHandle< std::vector< simb::MCTruth > >( label );
        art::FindManyP<simb::MCParticle> Assn(eventLabel,e,"largeant");
        if (eventLabel) {
          for ( size_t L1=0; L1 < eventLabel->size(); ++L1 ) {
            for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
              for (auto & data_label : all_data) {
                if (data_label.label == label) {
                  const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
                  data_label.E.push_back( ThisPar.E() - ThisPar.Mass() );
                  data_label.Time.push_back( e.time().value() - firstEvTime );
                }
            } 
              // SN Truth
              for (auto & data_label : all_data) {
                if (data_label.label == label) {
                  art::FindManyP<sim::SupernovaTruth> SNTruth(eventLabel, e, label);
                  for (size_t j = 0; j < SNTruth.at(L1).size(); j++) {
                    const sim::SupernovaTruth ThisTr = (*SNTruth.at(L1).at(j));
                    data_label.SNTime.push_back( ThisTr.SupernovaTime );
                  }
                }
              }
            }     
          }
        }
      }
    }
  }

}
DEFINE_ART_MODULE(SupernovaAna2)
