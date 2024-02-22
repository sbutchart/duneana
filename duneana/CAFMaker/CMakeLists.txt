if( DEFINED ENV{TENSORFLOW_DIR} )
# nusystematics
##cet_find_library( SYSTTOOLS NAMES systematicstools_interface systematicstools_interpreters systematicstools_utility PATHS ENV SYSTEMATICSTOOLS_LIB NO_DEFAULT_PATH )

cet_build_plugin(CAFMaker   art::module
                        duneanaobj_StandardRecord
                        duneanaobj_StandardRecordFlat
                        nugen::NuReweight
                        nugen::NuReweight_art
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
			${GENIE_LIB_LIST}
                        systematicstools::interface
			systematicstools::interpreters
			systematicstools::utility
                        dunereco::CVN_func
              BASENAME_ONLY
)
endif()

install_headers()
install_fhicl()
install_source()
install_scripts()
