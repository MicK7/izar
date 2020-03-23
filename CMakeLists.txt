cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(PluginName Izar)
project(${PluginName})

include(GNUInstallDirs)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(BUILD_SHARED_LIBS ON)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
endif()

# define a warning option
option(ENABLE_SUPPORT_WARNING "Check this to disable the warnings at the filters loading" ON)
mark_as_advanced(FORCE ENABLE_SUPPORT_WARNING)
if(${ENABLE_SUPPORT_WARNING})
	add_definitions(-DPRINT_SUPPORT_WARNING)
endif()

find_package(yaml-cpp REQUIRED)
find_package(Boost REQUIRED)

include_directories(${YAML_CPP_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
set(LIBS ${LIBS} ${YAML_CPP_LIBRARIES} ${Boost_LIBRARIES})

set(CMAKE_CXX_STANDARD 11)

# Scan the plugin file in order to set up internal data structures for building
# plugins.
paraview_plugin_scan(
  # The `paraview.plugin` file describing the plugin.
  PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/paraview.plugin"
  # A result variable for the (enabled) plugins found during the scan.
  PROVIDES_PLUGINS plugins
  # Enable plugins during this scan by default.
  ENABLE_BY_DEFAULT ON)

paraview_plugin_build(
  PLUGINS ${plugins})
