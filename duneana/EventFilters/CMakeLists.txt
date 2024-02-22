
cet_build_plugin(DuneEventFilter art::module
			dunecore::ArtSupport
			lardataobj::RawData
                        lardataobj::RecoBase
	                larcorealg::Geometry
	                larcore::Geometry_Geometry_service
                        lardata::Utilities    
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
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
			ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
