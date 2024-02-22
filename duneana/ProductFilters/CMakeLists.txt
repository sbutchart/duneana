
art_make( BASENAME_ONLY
         MODULE_LIBRARIES
         lardataobj::RecoBase
         lardata::ArtDataHelper
         larcorealg::Geometry
         larcore::Geometry_Geometry_service
         larsim::MCCheater_BackTrackerService_service
         larsim::MCCheater_ParticleInventoryService_service
         lardata::Utilities
	 nurandom::RandomUtils_NuRandomService_service
	 nusimdata::SimulationBase
 	 art::Framework_Core
	 art::Framework_Principal
	 art::Framework_Services_Registry
	 art_root_io::tfile_support
	 ROOT::Core
	 art_root_io::TFileService_service
	 art::Persistency_Common canvas
	 art::Utilities 
	 canvas
	 messagefacility::MF_MessageLogger		   
	 fhiclcpp::fhiclcpp
	 cetlib::cetlib 
	 cetlib_except
         ROOT::Core ROOT::Hist ROOT::Tree
         ROOT::Geom
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
