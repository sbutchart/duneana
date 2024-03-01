art_make(BASENAME_ONLY
	 MODULE_LIBRARIES
	 larreco::RecoAlg
	 lardataobj::RecoBase
         lardataobj::RawData
         lardata::ArtDataHelper
	 larcore::Geometry_Geometry_service
         larcorealg::Geometry
	 lardataalg::DetectorInfo
         nusimdata::SimulationBase
	 larsim::MCCheater_BackTrackerService_service
         larsim::MCCheater_ParticleInventoryService_service
	 dunereco::HitFinderDUNE
         art::Framework_Core
	 art::Framework_Principal
	 art::Framework_Services_Registry
         art_root_io::tfile_support
         ROOT::Core
         art_root_io::TFileService_service
 	 art::Persistency_Common
	 art::Persistency_Provenance 
	 art::Utilities canvas::canvas
	 messagefacility::MF_MessageLogger	 
	 cetlib::cetlib cetlib_except::cetlib_except
	 ROOT::Core ROOT::Hist ROOT::Tree
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_scripts()
