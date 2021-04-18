mmCreateView("GraphEntry_1","View3DGL","::View3DGL_1") 

mmCreateModule("PDBLoader","::PDBLoader_1") 
mmCreateModule("MoleculeSESRenderer","::MoleculeSESRenderer_1") 

mmCreateCall("MolecularDataCall","::MoleculeSESRenderer_1::getData","::PDBLoader_1::dataout")
mmCreateCall("CallRender3DGL","::View3DGL_1::rendering","::MoleculeSESRenderer_1::rendering")

mmSetParamValue("::PDBLoader_1::pdbFilename",[=[/home/anna/bin/megamol/data/SARSCov2/DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA.pdb]=])
mmSetParamValue("::PDBLoader_1::xtcFilename",[=[/home/anna/bin/megamol/data/SARSCov2/DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-no-water-no-ion-glueCA/sarscov2-10880334-0000.xtc]=])
mmSetParamValue("::PDBLoader_1::capFilename",[=[]=])
mmSetParamValue("::PDBLoader_1::maxFrames",[=[500]=])
mmSetParamValue("::PDBLoader_1::strideFlag",[=[true]=])
mmSetParamValue("::PDBLoader_1::solventResidues",[=[]=])
mmSetParamValue("::PDBLoader_1::calcBBoxPerFrame",[=[false]=])
mmSetParamValue("::PDBLoader_1::calculateBonds",[=[true]=])
mmSetParamValue("::PDBLoader_1::recomputeSTRIDEeachFrame",[=[false]=])
mmSetParamValue("::MoleculeSESRenderer_1::probeRadius",[=[1.400000]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::coloringMode0",[=[Chain]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::coloringMode1",[=[Element]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::colorWeighting",[=[0.500000]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::minGradColor",[=[#146496]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::midGradColor",[=[#f0f0f0]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::maxGradColor",[=[#ae3b32]=])
mmSetParamValue("::MoleculeSESRenderer_1::molIdxList",[=[0]=])
mmSetParamValue("::MoleculeSESRenderer_1::color::colorTableFilename",[=[colors.txt]=])
mmSetParamValue("::MoleculeSESRenderer_1::renderingMode",[=[GPU Ray Casting]=])
mmSetParamValue("::MoleculeSESRenderer_1::postProcessingMode",[=[None]=])
mmSetParamValue("::MoleculeSESRenderer_1::silhouetteColor",[=[255255255]=])
mmSetParamValue("::MoleculeSESRenderer_1::SSAOsigma",[=[5.000000]=])
mmSetParamValue("::MoleculeSESRenderer_1::SSAOlambda",[=[10.000000]=])
mmSetParamValue("::MoleculeSESRenderer_1::fogStart",[=[0.500000]=])
mmSetParamValue("::MoleculeSESRenderer_1::drawRS",[=[false]=])
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
mmSetParamValue("::View3DGL_1::anim::time",[=[0.040009]=])
mmSetParamValue("::View3DGL_1::backCol",[=[#000020]=])
mmSetParamValue("::View3DGL_1::showLookAt",[=[false]=])
mmSetParamValue("::View3DGL_1::viewKey::MoveStep",[=[0.500000]=])
mmSetParamValue("::View3DGL_1::viewKey::RunFactor",[=[2.000000]=])
mmSetParamValue("::View3DGL_1::viewKey::AngleStep",[=[90.000000]=])
mmSetParamValue("::View3DGL_1::viewKey::FixToWorldUp",[=[true]=])
mmSetParamValue("::View3DGL_1::viewKey::MouseSensitivity",[=[3.000000]=])
mmSetParamValue("::View3DGL_1::viewKey::RotPoint",[=[Look-At]=])
mmSetParamValue("::View3DGL_1::hookOnChange",[=[false]=])
mmSetParamValue("::View3DGL_1::cam::position",[=[-150.688477;87.9237289;54.8933563]=])
mmSetParamValue("::View3DGL_1::cam::orientation",[=[-0.357302755;-0.473673195;0.0943411738;0.7994169]=])
mmSetParamValue("::View3DGL_1::cam::projectiontype",[=[Perspective]=])
mmSetParamValue("::View3DGL_1::cam::convergenceplane",[=[0.000000]=])
mmSetParamValue("::View3DGL_1::cam::centeroffset",[=[0;0]=])
mmSetParamValue("::View3DGL_1::cam::halfaperturedegrees",[=[15.000000]=])
mmSetParamValue("::View3DGL_1::cam::halfdisparity",[=[0.025000]=])
mmSetParamValue("::View3DGL_1::cam::ovr::up",[=[0;0;0]=])
mmSetParamValue("::View3DGL_1::cam::ovr::lookat",[=[0;0;0]=])
mmSetParamValue("::View3DGL_1::view::defaultView",[=[Front]=])
mmSetParamValue("::View3DGL_1::view::defaultOrientation",[=[Top]=])
mmSetParamValue("::View3DGL_1::view::cubeOrientation",[=[-0.357302755;-0.473673195;0.0943411738;0.7994169]=])
mmSetParamValue("::View3DGL_1::view::showViewCube",[=[false]=])

