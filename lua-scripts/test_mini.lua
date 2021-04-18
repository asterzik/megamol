mmCreateView("GraphEntry_1","View3DGL","::View3DGL_1") 

mmCreateModule("PDBLoader","::PDBLoader_1") 
mmCreateModule("MoleculeSESRenderer","::MoleculeSESRenderer_1") 

mmCreateCall("MolecularDataCall","::MoleculeSESRenderer_1::getData","::PDBLoader_1::dataout")
mmCreateCall("CallRender3DGL","::View3DGL_1::rendering","::MoleculeSESRenderer_1::rendering")

mmSetParamValue("::PDBLoader_1::pdbFilename",[=[C:/Users/asterzik/bin/megamol/data/6wx4-DESRES-Trajectory_sarscov2-11441075-no-water-zinc-glueCA/sarscov2-11441075-no-water-zinc-glueCA/protein_static.pdb]=])
mmSetParamValue("::PDBLoader_1::xtcFilename",[=[]=])
mmSetParamValue("::PDBLoader_1::capFilename",[=[]=])
mmSetParamValue("::PDBLoader_1::maxFrames",[=[500]=])
mmSetParamValue("::PDBLoader_1::strideFlag",[=[true]=])
mmSetParamValue("::PDBLoader_1::solventResidues",[=[]=])
mmSetParamValue("::PDBLoader_1::calcBBoxPerFrame",[=[false]=])
mmSetParamValue("::PDBLoader_1::calculateBonds",[=[true]=])
mmSetParamValue("::PDBLoader_1::recomputeSTRIDEeachFrame",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::probeRadius",[=[1.400000]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::coloringMode0",[=[Element]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::coloringMode1",[=[Molecule]=])
mmSetParamValue("::MoleculeSESRenderer_1::displayedProperty",[=[Contour]=])
mmSetParamValue("::MoleculeSESRenderer_1::curvatureMode",[=[Prantl2Radial]=])
mmSetParamValue("::MoleculeSESRenderer_1::contourMode",[=[ShadingAndCurvature]=])
mmSetParamValue("::MoleculeSESRenderer_1::blur parameter",[=[DepthSensitive]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::colorWeighting",[=[0.500000]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::minGradColor",[=[#146496]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::midGradColor",[=[#f0f0f0]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::maxGradColor",[=[#ae3b32]=])
mmSetParamValue("::MoleculeSESRenderer_1::molIdxList",[=[0]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::colorTableFilename",[=[colors.txt]=])
mmSetParamValue("::MoleculeSESRenderer_1::SCRadius",[=[10]=])
mmSetParamValue("::MoleculeSESRenderer_1::SCNeighbourThreshold",[=[0.600000]=])
mmSetParamValue("::MoleculeSESRenderer_1::SCDiffThreshold",[=[0.050000]=])
mmSetParamValue("::MoleculeSESRenderer_1::SCMedianFilter",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::SCCircularNeighborhood",[=[true]=])
mmSetParamValue("::MoleculeSESRenderer_1::testcase",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::cylinder",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::cut off point for contours",[=[10.000000]=])
mmSetParamValue("::MoleculeSESRenderer_1::depthDiff",[=[0.100000]=])
mmSetParamValue("::MoleculeSESRenderer_1::# normal blurring iterations",[=[10.000000]=])
mmSetParamValue("::MoleculeSESRenderer_1::# position blurring iterations",[=[10.000000]=])
mmSetParamValue("::MoleculeSESRenderer_1::# curvature blurring iterations",[=[10.000000]=])
mmSetParamValue("::MoleculeSESRenderer_1::white background",[=[true]=])
mmSetParamValue("::MoleculeSESRenderer_1::overlay",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::curvatureDiff",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::drawSES",[=[true]=])
mmSetParamValue("::MoleculeSESRenderer_1::drawSAS",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::offscreenRendering",[=[true]=])
mmSetParamValue("::View3DGL_1::camstore::settings",[=[]=])
mmSetParamValue("::View3DGL_1::camstore::overrideSettings",[=[false]=])
mmSetParamValue("::View3DGL_1::camstore::autoSaveSettings",[=[false]=])
mmSetParamValue("::View3DGL_1::camstore::autoLoadSettings",[=[true]=])
mmSetParamValue("::View3DGL_1::resetViewOnBBoxChange",[=[false]=])
mmSetParamValue("::View3DGL_1::anim::play",[=[false]=])
mmSetParamValue("::View3DGL_1::anim::speed",[=[4.000000]=])
mmSetParamValue("::View3DGL_1::anim::time",[=[78.112000]=])
mmSetParamValue("::View3DGL_1::backCol",[=[#000020]=])
mmSetParamValue("::View3DGL_1::showLookAt",[=[false]=])
mmSetParamValue("::View3DGL_1::viewKey::MoveStep",[=[0.500000]=])
mmSetParamValue("::View3DGL_1::viewKey::RunFactor",[=[2.000000]=])
mmSetParamValue("::View3DGL_1::viewKey::AngleStep",[=[90.000000]=])
mmSetParamValue("::View3DGL_1::viewKey::FixToWorldUp",[=[true]=])
mmSetParamValue("::View3DGL_1::viewKey::MouseSensitivity",[=[3.000000]=])
mmSetParamValue("::View3DGL_1::viewKey::RotPoint",[=[Look-At]=])
mmSetParamValue("::View3DGL_1::hookOnChange",[=[false]=])
mmSetParamValue("::View3DGL_1::cam::position",[=[-88.73172;-36.7335663;33.7506866]=])
mmSetParamValue("::View3DGL_1::cam::orientation",[=[0.568206131;0.0919057652;-0.801451862;0.162421152]=])
mmSetParamValue("::View3DGL_1::cam::projectiontype",[=[Perspective]=])
mmSetParamValue("::View3DGL_1::cam::convergenceplane",[=[0.000000]=])
mmSetParamValue("::View3DGL_1::cam::centeroffset",[=[0;0]=])
mmSetParamValue("::View3DGL_1::cam::halfaperturedegrees",[=[15.000000]=])
mmSetParamValue("::View3DGL_1::cam::halfdisparity",[=[0.025000]=])
mmSetParamValue("::View3DGL_1::cam::ovr::up",[=[0;0;0]=])
mmSetParamValue("::View3DGL_1::cam::ovr::lookat",[=[0;0;0]=])
mmSetParamValue("::View3DGL_1::view::defaultView",[=[Front]=])
mmSetParamValue("::View3DGL_1::view::defaultOrientation",[=[Top]=])
mmSetParamValue("::View3DGL_1::view::cubeOrientation",[=[0.568206131;0.0919057652;-0.801451862;0.162421152]=])
mmSetParamValue("::View3DGL_1::view::showViewCube",[=[false]=])

mmSetGUIVisible(true)
mmSetGUIScale(1.500000)
mmSetGUIState([=[{"ConfiguratorState":{"module_list_sidebar_width":250.0,"show_module_list_sidebar":true},"GUIState":{"font_file_name":"","font_size":19,"imgui_settings":"[Window][DockSpaceViewport_11111111]\nPos=0,26\nSize=574,335\nCollapsed=0\n\n[Window][Debug##Default]\nPos=60,60\nSize=400,400\nCollapsed=0\n\n[Window][Log Console     F9]\nPos=1400,1052\nSize=1920,264\nCollapsed=0\n\n[Window][Parameters     F10]\nPos=-15,333\nSize=574,335\nCollapsed=0\n\n[Window][Save Project (.lua)]\nPos=9,-372\nSize=600,750\nCollapsed=0\n\n[Window][Select File]\nPos=661,165\nSize=600,750\nCollapsed=0\n\n[Window][Select Filename for Screenshot (.png)]\nPos=1629,658\nSize=600,750\nCollapsed=0\n\n[Docking][Data]\nDockSpace ID=0x8B93E3BD Window=0xA787BDB4 Pos=0,26 Size=574,335 CentralNode=1\n\n","menu_visible":true,"style":2},"GraphStates":{"Project":{"Modules":{"::MoleculeSESRenderer_1":{"graph_position":[3.4028234663852886e+38,3.4028234663852886e+38]},"::PDBLoader_1":{"graph_position":[3.4028234663852886e+38,3.4028234663852886e+38]},"::View3DGL_1":{"graph_position":[3.4028234663852886e+38,3.4028234663852886e+38]}},"canvas_scrolling":[0.0,0.0],"canvas_zooming":1.0,"param_extended_mode":false,"parameter_sidebar_width":300.0,"params_readonly":false,"params_visible":true,"project_name":"Project_1","show_call_label":true,"show_call_slots_label":false,"show_grid":false,"show_module_label":true,"show_parameter_sidebar":true,"show_slot_label":false}},"ParameterStates":{"::MoleculeSESRenderer_1::# curvature blurring iterations":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::# normal blurring iterations":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::# position blurring iterations":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::SCCircularNeighborhood":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::SCDiffThreshold":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::SCMedianFilter":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::SCNeighbourThreshold":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::SCRadius":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::ViewType":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::blur parameter":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::color::colorTableFilename":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::color::colorWeighting":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::color::coloringMode0":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::color::coloringMode1":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::color::maxGradColor":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::color::midGradColor":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::color::minGradColor":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::contourMode":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::curvatureDiff":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::curvatureMode":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::cut off point for contours":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::cylinder":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::depthDiff":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::displayedProperty":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::drawSAS":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::drawSES":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::extendContours":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::molIdxList":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::offscreenRendering":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::orthographic projection":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::overlay":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::probeRadius":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::smoothCurvature":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::smoothNormals":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::smoothPositions":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::smoothTimesteps":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::testcase":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::MoleculeSESRenderer_1::white background":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::calcBBoxPerFrame":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::calculateBonds":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::capFilename":{"gui_presentation_mode":16,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::maxFrames":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::pdbFilename":{"gui_presentation_mode":16,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::recomputeSTRIDEeachFrame":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::solventResidues":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::strideFlag":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::PDBLoader_1::xtcFilename":{"gui_presentation_mode":16,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::ParameterGroup::anim":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::ParameterGroup::view":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::anim::SpeedDown":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::anim::SpeedUp":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::anim::play":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::anim::speed":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::anim::time":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::anim::togglePlay":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::backCol":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::centeroffset":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::convergenceplane":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::halfaperturedegrees":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::halfdisparity":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::orientation":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::ovr::lookat":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::ovr::override":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::ovr::up":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::position":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::cam::projectiontype":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::camstore::autoLoadSettings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::camstore::autoSaveSettings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::camstore::overrideSettings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::camstore::restorecam":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::camstore::settings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::camstore::storecam":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::enableMouseSelection":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::hookOnChange":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::resetViewOnBBoxChange":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::showLookAt":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::view::cubeOrientation":{"gui_presentation_mode":2,"gui_read-only":true,"gui_visibility":false},"::View3DGL_1::view::defaultOrientation":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::view::defaultView":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::view::resetView":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::view::showViewCube":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::viewKey::AngleStep":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::viewKey::FixToWorldUp":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::viewKey::MouseSensitivity":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::viewKey::MoveStep":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::viewKey::RotPoint":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::View3DGL_1::viewKey::RunFactor":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true}},"WindowConfigurations":{"Configurator":{"fpsms_max_value_count":20,"fpsms_mode":0,"fpsms_refresh_rate":2.0,"fpsms_show_options":false,"log_force_open":true,"log_level":4294967295,"param_extended_mode":false,"param_module_filter":0,"param_modules_list":[],"param_show_hotkeys":false,"tfe_active_param":"","tfe_view_minimized":false,"tfe_view_vertical":false,"win_callback":6,"win_collapsed":false,"win_flags":1032,"win_hotkey":[300,0],"win_position":[0.0,0.0],"win_reset_position":[0.0,0.0],"win_reset_size":[1920.0,1080.0],"win_show":false,"win_size":[1920.0,1080.0]},"Log Console":{"fpsms_max_value_count":20,"fpsms_mode":0,"fpsms_refresh_rate":2.0,"fpsms_show_options":false,"log_force_open":true,"log_level":4294967295,"param_extended_mode":false,"param_module_filter":0,"param_modules_list":[],"param_show_hotkeys":false,"tfe_active_param":"","tfe_view_minimized":false,"tfe_view_vertical":false,"win_callback":7,"win_collapsed":false,"win_flags":3072,"win_hotkey":[298,0],"win_position":[1400.0,1052.0],"win_reset_position":[0.0,904.0],"win_reset_size":[1920.0,176.0],"win_show":false,"win_size":[1280.0,176.0]},"Parameters":{"fpsms_max_value_count":20,"fpsms_mode":0,"fpsms_refresh_rate":2.0,"fpsms_show_options":false,"log_force_open":true,"log_level":4294967295,"param_extended_mode":false,"param_module_filter":0,"param_modules_list":[],"param_show_hotkeys":false,"tfe_active_param":"","tfe_view_minimized":false,"tfe_view_vertical":false,"win_callback":1,"win_collapsed":false,"win_flags":8,"win_hotkey":[299,0],"win_position":[-15.0,333.0],"win_reset_position":[0.0,0.0],"win_reset_size":[400.0,500.0],"win_show":true,"win_size":[382.6666564941406,223.3333282470703]},"Performance Metrics":{"fpsms_max_value_count":20,"fpsms_mode":0,"fpsms_refresh_rate":2.0,"fpsms_show_options":false,"log_force_open":true,"log_level":4294967295,"param_extended_mode":false,"param_module_filter":0,"param_modules_list":[],"param_show_hotkeys":false,"tfe_active_param":"","tfe_view_minimized":false,"tfe_view_vertical":false,"win_callback":3,"win_collapsed":false,"win_flags":65,"win_hotkey":[296,0],"win_position":[960.0,0.0],"win_reset_position":[960.0,0.0],"win_reset_size":[0.0,0.0],"win_show":false,"win_size":[0.0,0.0]},"Transfer Function Editor":{"fpsms_max_value_count":20,"fpsms_mode":0,"fpsms_refresh_rate":2.0,"fpsms_show_options":false,"log_force_open":true,"log_level":4294967295,"param_extended_mode":false,"param_module_filter":0,"param_modules_list":[],"param_show_hotkeys":false,"tfe_active_param":"","tfe_view_minimized":false,"tfe_view_vertical":false,"win_callback":5,"win_collapsed":false,"win_flags":64,"win_hotkey":[297,0],"win_position":[400.0,0.0],"win_reset_position":[400.0,0.0],"win_reset_size":[0.0,0.0],"win_show":false,"win_size":[0.0,0.0]}}}]=])
