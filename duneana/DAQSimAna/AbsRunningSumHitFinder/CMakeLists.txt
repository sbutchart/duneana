make_dune_tool_builder(dunecore::AbsRunningSumTPFinderTool)

art_make( BASENAME_ONLY MODULE_LIBRARIES
  dunepdlegacy_Services_ChannelMap_PdspChannelMapService_service
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  larcorealg::Geometry
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  dunecore::DuneObj
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larsim::MCCheater_PhotonBackTrackerService_service
  larsim::Simulation nug4::ParticleNavigation lardataobj_Simulation
  lardata::ArtDataHelper
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
  EXCLUDE read_packed.cxx AbsRunningSumTPFinderPass1_tool.cc
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

cet_build_plugin(AbsRunningSumTPFinderPass1 art::tool
  fhiclcpp cetlib::cetlib cetlib_except

)

install_fhicl()
install_headers()
install_source()
