include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make(BASENAME_ONLY
   LIB_LIBRARIES
    ${ROOT_CORE}
    ${ROOT_EVE}
    cetlib_except::cetlib_except
    cetlib::cetlib 
    larcorealg::Geometry
    lardataob::_RawData
    art::Persistency_Common canvas::canvas
   SERVICE_LIBRARIES
    duneana::EventDisplay3D
    ${ART_ROOT_IO_ROOTINPUT_SOURCE}
   MODULE_LIBRARIES
    nusimdata::SimulationBase
    larreco::RecoAlg
    duneana::EventDisplay3D
    EvtDisplayService_service
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
