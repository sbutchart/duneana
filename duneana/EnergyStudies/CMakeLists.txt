
cet_build_plugin(dunezsanalysis art::module
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RawData
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except
  ROOT::Core ROOT::Hist ROOT::Tree
  BASENAME_ONLY
  )

cet_build_plugin(dEdx art::module
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RawData
  lardataobj::RecoBase
  lardata::Utilities
  larsim::Simulation 
  nug4::ParticleNavigation 
  lardataobj::Simulation
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except
  ROOT::Core ROOT::Hist ROOT::Tree
  BASENAME_ONLY
  )

cet_build_plugin(EMEnergyCalib art::module
  lardataalg::DetectorInfo
  lardataobj::RawData
  lardataobj::RecoBase
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larsim::Simulation 
  nug4::ParticleNavigation 
  lardataobj::Simulation
  lardata::Utilities
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except
  ROOT::Core ROOT::Hist ROOT::Tree
  BASENAME_ONLY
  )

cet_build_plugin(EMPi0Energy art::module
  lardataobj::RawData
  lardataobj::RecoBase
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larsim::Simulation nug4::ParticleNavigation lardataobj_Simulation
  lardata::Utilities
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except
  ROOT::Core ROOT::Hist ROOT::Tree
  BASENAME_ONLY
  )

install_headers()
install_fhicl()
install_source()
install_scripts()
