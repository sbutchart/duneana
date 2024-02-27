# ======================================================================
#  duneana main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.20 FATAL_ERROR) ### MIGRATE-CHANGED (migrate-3.22.02) - VERSION 3.19 -> 3.20
# cmake_policy(VERSION 3.19) ### MIGRATE-ACTION-COMPATIBILITY (migrate-3.22.02) - Uncomment to preserve compatibility with older CMake versions

#project(duneana VERSION 09.48.00.d00 LANGUAGES CXX)
find_package(cetmodules REQUIRED) ### MIGRATE-ADDED (migrate-3.22.02) - find_package(cetmodules) must precede project()
project(duneana LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.83.01d00)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags() ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - add args: REPORT_THRESHOLD VERBOSE

find_package( art REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art ...) -> find_package(art REQUIRED ...)
find_package( art_root_io REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art_root_io ...) -> find_package(art_root_io REQUIRED ...)
find_package( canvas_root_io REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(canvas_root_io ...) -> find_package(canvas_root_io REQUIRED ...)
find_package( Boost REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package(ROOT REQUIRED) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package( dunecore REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunecore ...) -> find_package(dunecore REQUIRED ...)
find_package( dunereco REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunereco ...) -> find_package(dunereco REQUIRED ...)
find_package( duneopdet REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(duneopdet ...) -> find_package(duneopdet REQUIRED ...)
find_package( duneanaobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(duneanaobj ...) -> find_package(duneanaobj REQUIRED ...)
find_package( larcoreobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcoreobj ...) -> find_package(larcoreobj REQUIRED ...)
find_package( larreco REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larreco ...) -> find_package(larreco REQUIRED ...)
find_package( larsim REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larsim ...) -> find_package(larsim REQUIRED ...)
find_package( larevt REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larevt ...) -> find_package(larevt REQUIRED ...)
find_package( larpandora REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larpandora ...) -> find_package(larpandora REQUIRED ...)
#find_package( srproxy REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(srproxy ...) -> find_package(srproxy REQUIRED ...)
find_package( systematicstools REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(systematicstools ...) -> find_package(systematicstools REQUIRED ...)
find_package(GENIE REQUIRED EXPORT)
#find_ups_product( dunepdlegacy )
find_package( Geant4 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(duneana)

cet_cmake_config() ### MIGRATE-ADDED (migrate-3.22.02) - required cmd
