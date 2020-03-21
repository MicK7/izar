paraview_add_plugin(${PluginName}
  REQUIRED_ON_SERVER
  VERSION "0.2"
  MODULES Izar::IzarModule
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/module/vtk.module")
