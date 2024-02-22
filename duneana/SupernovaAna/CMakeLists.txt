art_make(BASENAME_ONLY
         MODULE_LIBRARIES
         dunecore::DAQTriggerSim_TriggerDataProducts
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
         art::Persistency_Common
         art::Persistency_Provenance
         art::Utilities 
	 canvas
         art::Persistency_Common
         art::Persistency_Provenance
         messagefacility::MF_MessageLogger
         cetlib::cetlib 
	 cetlib_except
         ROOT::Core ROOT::Hist ROOT::Tree
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_fhicl()
install_headers()
install_source()
