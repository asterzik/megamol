/*
 * MoleculeSESRenderer.cpp
 *
 * Copyright (C) 2009-2021 by Universitaet Stuttgart (VIS). Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define _USE_MATH_DEFINES 1

#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include "../include/protein/magic_enum.hpp"
#include "Color.h"
#include "MoleculeSESRenderer.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/StringParam.h"
#include "mmcore/utility/ColourParser.h"
#include "mmcore/utility/sys/ASCIIFileBuffer.h"
#include "protein/RMSF.h"
#include "vislib/OutOfRangeException.h"
#include "vislib/StringConverter.h"
#include "vislib/StringTokeniser.h"
#include "vislib/Trace.h"
#include "vislib/assert.h"
#include "vislib/graphics/gl/AbstractOpenGLShader.h"
#include "vislib/graphics/gl/IncludeAllGL.h"
#include "vislib/graphics/gl/ShaderSource.h"
#include "vislib/sys/File.h"


using namespace megamol;
using namespace megamol::core;
using namespace megamol::protein;
using namespace megamol::protein_calls;
using namespace megamol::core::utility::log;
MoleculeSESRenderer::MoleculeSESRenderer(void)
        : Renderer3DModuleGL()
        , molDataCallerSlot("getData", "Connects the protein SES rendering with protein data storage")
        , bsDataCallerSlot("getBindingSites", "Connects the molecule rendering with binding site data storage")
        , coloringModeParam0("color::coloringMode0", "The first coloring mode.")
        , coloringModeParam1("color::coloringMode1", "The second coloring mode.")
        , cmWeightParam("color::colorWeighting", "The weighting of the two coloring modes.")
        , minGradColorParam("color::minGradColor", "The color for the minimum value for gradient coloring")
        , midGradColorParam("color::midGradColor", "The color for the middle value for gradient coloring")
        , maxGradColorParam("color::maxGradColor", "The color for the maximum value for gradient coloring")
        , displayedPropertyParam("displayedProperty", "Choose the property to be displayed")
        , curvatureModeParam("curvatureMode", "curvature mode.")
        , contourModeParam("contourMode", "Specify which contour generation method to use")
        , drawSESParam("drawSES", "Draw the SES: ")
        , drawSASParam("drawSAS", "Draw the SAS: ")
        , molIdxListParam("molIdxList", "The list of molecule indices for RS computation:")
        , colorTableFileParam("color::colorTableFilename", "The filename of the color table.")
        , offscreenRenderingParam("offscreenRendering", "Toggle offscreen rendering.")
        , probeRadiusSlot("probeRadius", "The probe radius for the surface computation")
        , smoothNormalsParam("smoothNormals", "Determines whether the normals will be smoothed")
        , smoothPositionsParam("smoothPositions", "Determines whether the positions will be smoothed")
        , smoothCurvatureParam("smoothCurvature", "Determines whether the curvature will be smoothed")
        , numBlurParam("# normal blurring iterations", " How many iterations of blurring for normals?")
        , numPosBlurParam("# position blurring iterations", " How many iterations of blurring for positions?")
        , numCurvBlurParam("# curvature blurring iterations", " How many iterations of blurring for curvature?")
        // , pyramidWeightsParam(
        //       "pyramidWeights", "The factor for the weights in the pull phase of the pull-push algorithm")
        // , pyramidLayersParam("pyramidLayers", "Number of layers in the pull-push normalPyramid")
        // , pyramidGammaParam("pyramidGamma",
        //       "The higher the exponent gamma, the more non-linear the interpolation between points becomes")
        , SCRadiusParam("SCRadius", "Radius to consider around one pixel for SC")
        , SCNeighbourThresholdParam("SCNeighbourThreshold",
              "How many darker pixels are allowed to be in the surrounding to still be rendered")
        , SCDiffThresholdParam(
              "SCDiffThreshold", "How much intensity difference needs to be there, for the pixel to be rendered")
        , SCMedianFilterParam("SCMedianFilter", "Use median filter for suggestive contours?")
        , SCCircularNeighborhoodParam(
              "SCCircularNeighborhood", "Use circular neighborhood for suggestive contours? Alternative is quadratic.")
        , ViewTypeParam("ViewType", "zDir: viewDir = (0,0,1), minusPos: viewdir = normalize(viewPos - fragPos) "
                                    "(viewPos = vec3(0)), interpolate: interpolate between both.")
        , OrthoProjParam("orthographic projection",
              "Was orthographic projection used for the data generation? Other possibility perspective.")
        , TestCaseParam("testcase", "Use the test case? Otherwise real data.")
        , cutOffParam("cut off point for contours",
              "Curvature Contours: How big can the dot product between normal and "
              "viewdir get, such that the point is still considered a contour?")
        , blurParam("blur parameter", " Which blur shader to use?")
        , depthDiffParam("depthDiff", " How big is the z-Position difference of two pixels allowed to be for blurring "
                                      "with the depth sensitive blur shader?")
        , whiteBackgroundParam("white background", "White Background instead of the standard megamol blue.")
        , extendContoursParam("extendContours", "extendContours")
        , dilation2RadiusParam("dilation2Radius", "dilation2Radius")
        , dilation1RadiusParam("dilation1Radius", "dilation1Radius")
        , erosionRadiusParam("erosionRadius", "erosionRadius")
        , smoothTimestepsParam("smoothTimesteps", "smoothTimesteps")
        , computeSesPerMolecule(false) {
#pragma region // Set parameters

    this->molDataCallerSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->molDataCallerSlot);
    this->bsDataCallerSlot.SetCompatibleCall<BindingSiteCallDescription>();
    this->MakeSlotAvailable(&this->bsDataCallerSlot);

    // set epsilon value for float-comparison
    this->epsilon = vislib::math::FLOAT_EPSILON;
    // set probe radius
    this->probeRadius = 1.4f;

    this->probeRadiusSlot.SetParameter(new param::FloatParam(1.4f, 0.1f));
    this->MakeSlotAvailable(&this->probeRadiusSlot);

    // coloring modes
    this->currentColoringMode0 = Color::CHAIN;
    this->currentColoringMode1 = Color::ELEMENT;
    param::EnumParam* cm0 = new param::EnumParam(int(this->currentColoringMode0));
    param::EnumParam* cm1 = new param::EnumParam(int(this->currentColoringMode1));
    MolecularDataCall* mol = new MolecularDataCall();
    BindingSiteCall* bs = new BindingSiteCall();
    unsigned int cCnt;
    Color::ColoringMode cMode;
    for (cCnt = 0; cCnt < Color::GetNumOfColoringModes(mol, bs); ++cCnt) {
        cMode = Color::GetModeByIndex(mol, bs, cCnt);
        cm0->SetTypePair(cMode, Color::GetName(cMode).c_str());
        cm1->SetTypePair(cMode, Color::GetName(cMode).c_str());
    }
    delete mol;
    delete bs;
    this->coloringModeParam0 << cm0;
    this->coloringModeParam1 << cm1;
    this->MakeSlotAvailable(&this->coloringModeParam0);
    this->MakeSlotAvailable(&this->coloringModeParam1);

    // displayedProperty
    this->currentDisplayedProperty = Contour;
    param::EnumParam* dmp = new param::EnumParam(int(this->currentDisplayedProperty));
    constexpr auto& property_entries = magic_enum::enum_entries<displayedProperty>();
    for (int i = 0; i < magic_enum::enum_count<displayedProperty>(); ++i) {
        dmp->SetTypePair((int) property_entries[i].first, std::string(property_entries[i].second).c_str());
    }
    this->displayedPropertyParam << dmp;
    this->MakeSlotAvailable(&this->displayedPropertyParam);


    // curvature modes
    this->currentCurvatureMode = EvansCurvature;
    param::EnumParam* cmp = new param::EnumParam(int(this->currentCurvatureMode));
    constexpr auto& curvature_entries = magic_enum::enum_entries<curvatureMode>();
    for (int i = 0; i < magic_enum::enum_count<curvatureMode>(); ++i) {
        cmp->SetTypePair((int) curvature_entries[i].first, std::string(curvature_entries[i].second).c_str());
    }
    this->curvatureModeParam << cmp;
    this->MakeSlotAvailable(&this->curvatureModeParam);

    // contour modes
    this->currentContourMode = Shading;
    param::EnumParam* cContourp = new param::EnumParam(int(this->currentContourMode));
    constexpr auto& contour_entries = magic_enum::enum_entries<contourMode>();
    for (int i = 0; i < magic_enum::enum_count<contourMode>(); ++i) {
        cContourp->SetTypePair((int) contour_entries[i].first, std::string(contour_entries[i].second).c_str());
    }
    this->contourModeParam << cContourp;
    this->MakeSlotAvailable(&this->contourModeParam);

    // blur modes
    this->currentBlurMode = Gaussian;
    param::EnumParam* cbp = new param::EnumParam(int(this->currentBlurMode));
    constexpr auto& blur_entries = magic_enum::enum_entries<blurMode>();
    for (int i = 0; i < magic_enum::enum_count<blurMode>(); ++i) {
        cbp->SetTypePair((int) blur_entries[i].first, std::string(blur_entries[i].second).c_str());
    }
    this->blurParam << cbp;
    this->MakeSlotAvailable(&this->blurParam);

    // blur modes
    this->currentViewType = minusPos;
    param::EnumParam* cvt = new param::EnumParam(int(this->currentViewType));
    constexpr auto& view_entries = magic_enum::enum_entries<viewType>();
    for (int i = 0; i < magic_enum::enum_count<viewType>(); ++i) {
        cvt->SetTypePair((int) view_entries[i].first, std::string(view_entries[i].second).c_str());
    }
    this->ViewTypeParam << cvt;
    this->MakeSlotAvailable(&this->ViewTypeParam);

    // Color weighting parameter
    this->cmWeightParam.SetParameter(new param::FloatParam(0.5f, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->cmWeightParam);

    // the color for the minimum value (gradient coloring
    this->minGradColorParam.SetParameter(new param::StringParam("#146496"));
    this->MakeSlotAvailable(&this->minGradColorParam);

    // the color for the middle value (gradient coloring
    this->midGradColorParam.SetParameter(new param::StringParam("#f0f0f0"));
    this->MakeSlotAvailable(&this->midGradColorParam);

    // the color for the maximum value (gradient coloring
    this->maxGradColorParam.SetParameter(new param::StringParam("#ae3b32"));
    this->MakeSlotAvailable(&this->maxGradColorParam);

    // ----- draw SES param -----
    this->drawSES = true;
    param::BoolParam* sespm = new param::BoolParam(this->drawSES);
    this->drawSESParam << sespm;

    // ----- draw SAS param -----
    this->drawSAS = false;
    param::BoolParam* saspm = new param::BoolParam(this->drawSAS);
    this->drawSASParam << saspm;

    // ----- ofsfcreen rendering param -----
    this->offscreenRendering = true;
    param::BoolParam* orpm = new param::BoolParam(this->offscreenRendering);
    this->offscreenRenderingParam << orpm;

    // ----- molecular indices list param -----
    this->molIdxList.Add("0");
    this->molIdxListParam.SetParameter(new param::StringParam("0"));
    this->MakeSlotAvailable(&this->molIdxListParam);

    // fill color table with default values and set the filename param
    vislib::StringA filename("colors.txt");
    Color::ReadColorTableFromFile(filename, this->colorLookupTable);
    this->colorTableFileParam.SetParameter(new param::StringParam(A2T(filename)));
    this->MakeSlotAvailable(&this->colorTableFileParam);

#pragma endregion Set parameters

    this->smoothNormals = true;
    this->smoothNormalsParam.SetParameter(new param::BoolParam(this->smoothNormals));
    this->MakeSlotAvailable(&this->smoothNormalsParam);

    this->smoothCurvature = true;
    this->smoothCurvatureParam.SetParameter(new param::BoolParam(this->smoothCurvature));
    this->MakeSlotAvailable(&this->smoothCurvatureParam);

    // this->pyramidWeight = 0.5f;
    // this->pyramidWeightsParam.SetParameter(new param::FloatParam(this->pyramidWeight));
    // this->MakeSlotAvailable(&this->pyramidWeightsParam);

    // normalPyramid.create("fragNormal", this->width, this->height, this->GetCoreInstance(), "pullpush::pullNormal",
    //     "pullpush::pushNormal");
    const int mipmapNumber = (int) glm::log2(glm::max<float>(this->width, this->height)) + 1;
    // this->pyramidLayers = 3;
    // // TODO: somehow this breaks if i do not hardcode the second parameter
    // this->pyramidLayersParam.SetParameter(new param::IntParam(this->pyramidLayers, 1, 11));
    // this->MakeSlotAvailable(&this->pyramidLayersParam);

    // this->pyramidGamma = 1.0f;
    // this->pyramidGammaParam.SetParameter(new param::FloatParam(this->pyramidGamma, 1.0));
    // this->MakeSlotAvailable(&this->pyramidGammaParam);

    this->smoothPositions = true;
    this->smoothPositionsParam.SetParameter(new param::BoolParam(this->smoothPositions));
    this->MakeSlotAvailable(&this->smoothPositionsParam);

    // // Suggestive Contours parameters

    this->SCRadius = 3;
    this->SCRadiusParam.SetParameter(new param::IntParam(this->SCRadius, 1, 10));
    this->MakeSlotAvailable(&this->SCRadiusParam);

    this->SCNeighbourThreshold = 0.2f; // percentage s in original  SC paper //  in original paper 0.2
    this->SCNeighbourThresholdParam.SetParameter(new param::FloatParam(this->SCNeighbourThreshold, 0.0, 1.0));
    this->MakeSlotAvailable(&this->SCNeighbourThresholdParam);

    this->SCDiffThreshold = 0.2f; // threshold d in original SC paper // in original paper 0.25
    this->SCDiffThresholdParam.SetParameter(new param::FloatParam(this->SCDiffThreshold, 0.0, 1.0));
    this->MakeSlotAvailable(&this->SCDiffThresholdParam);

    this->SCMedianFilter = false;
    this->SCMedianFilterParam.SetParameter(new param::BoolParam(this->SCMedianFilter));
    this->MakeSlotAvailable(&this->SCMedianFilterParam);

    this->SCCircularNeighborhood = true;
    this->SCCircularNeighborhoodParam.SetParameter(new param::BoolParam(this->SCCircularNeighborhood));
    this->MakeSlotAvailable(&this->SCCircularNeighborhoodParam);

    this->orthoproj = false;
    this->OrthoProjParam.SetParameter(new param::BoolParam(this->orthoproj));
    this->MakeSlotAvailable(&this->OrthoProjParam);

    this->testcase = false;
    this->TestCaseParam.SetParameter(new param::BoolParam(this->testcase));
    this->MakeSlotAvailable(&this->TestCaseParam);

    this->cutOff = 0.15f;
    this->cutOffParam.SetParameter(new param::FloatParam(this->cutOff));
    this->MakeSlotAvailable(&this->cutOffParam);

    this->depthDiff = 0.3f;
    this->depthDiffParam.SetParameter(new param::FloatParam(this->depthDiff));
    this->MakeSlotAvailable(&this->depthDiffParam);

    this->numBlur = 1;
    this->numBlurParam.SetParameter(new param::FloatParam(this->numBlur, 1.0));
    this->MakeSlotAvailable(&this->numBlurParam);

    this->numPosBlur = 10;
    this->numPosBlurParam.SetParameter(new param::FloatParam(this->numPosBlur, 1.0));
    this->MakeSlotAvailable(&this->numPosBlurParam);

    this->numCurvBlur = 1;
    this->numCurvBlurParam.SetParameter(new param::FloatParam(this->numCurvBlur, 1.0));
    this->MakeSlotAvailable(&this->numCurvBlurParam);

    this->whiteBackground = false;
    this->whiteBackgroundParam.SetParameter(new param::BoolParam(this->whiteBackground));
    this->MakeSlotAvailable(&this->whiteBackgroundParam);

    this->extendContoursBool = false;
    this->extendContoursParam.SetParameter(new param::BoolParam(this->extendContoursBool));
    this->MakeSlotAvailable(&this->extendContoursParam);

    this->dilation2Radius = 4;
    this->dilation2RadiusParam.SetParameter(new param::IntParam(this->dilation2Radius));
    this->MakeSlotAvailable(&this->dilation2RadiusParam);

    this->dilation1Radius = 3;
    this->dilation1RadiusParam.SetParameter(new param::IntParam(this->dilation1Radius, 0));
    this->MakeSlotAvailable(&this->dilation1RadiusParam);

    this->erosionRadius = 2;
    this->erosionRadiusParam.SetParameter(new param::IntParam(this->erosionRadius, 0));
    this->MakeSlotAvailable(&this->erosionRadiusParam);

    this->smoothTimestepsBool = true;
    this->smoothTimestepsParam.SetParameter(new param::BoolParam(this->smoothTimestepsBool));
    this->MakeSlotAvailable(&this->smoothTimestepsParam);

    // fill rainbow color table
    Color::MakeRainbowColorTable(100, this->rainbowColors);

    this->contourFBO = 0;
    this->extendContourFBO[0] = 0;
    this->extendContourFBO[1] = 0;
    this->timestepsFBO[0] = 0;
    this->timestepsFBO[1] = 0;
    this->timestepsFBO[2] = 0;
    this->contourDepthRBO = 0;
    this->curvatureFBO = 0;
    this->positionFBO[0] = 0;
    this->positionFBO[1] = 0;
    this->normalFBO[0] = 0;
    this->normalFBO[1] = 0;
    this->smoothCurvFBO[0] = 0;
    this->smoothCurvFBO[1] = 0;

    // width and height of the screen
    this->width = 0;
    this->height = 0;

    // clear singularity texture
    singularityTexture.clear();
    // set singTexData to 0
    this->singTexData = 0;

    this->preComputationDone = false;

#pragma region // export parameters
    this->MakeSlotAvailable(&this->drawSESParam);
    this->MakeSlotAvailable(&this->drawSASParam);
    this->MakeSlotAvailable(&this->offscreenRenderingParam);
#pragma endregion // export parameters
}

MoleculeSESRenderer::~MoleculeSESRenderer(void) {
    if (contourFBO) {
        glDeleteFramebuffers(1, &contourFBO);
        glDeleteRenderbuffers(1, &contourDepthRBO);
        glDeleteTextures(1, &normalTexture);
        glDeleteTextures(1, &positionTexture);
        glDeleteTextures(1, &objPositionTexture);
    }
    if (extendContourFBO) {
        glDeleteFramebuffers(2, extendContourFBO);
        glDeleteTextures(2, contourTexture);
    }
    if (curvatureFBO) {
        glDeleteFramebuffers(1, &curvatureFBO);
        glDeleteTextures(1, &curvatureTexture);
    }
    if (positionFBO) {
        glDeleteFramebuffers(2, positionFBO);
        glDeleteTextures(2, smoothPositionTexture);
    }
    if (normalFBO) {
        glDeleteFramebuffers(2, normalFBO);
        glDeleteTextures(2, smoothNormalTexture);
    }
    if (smoothCurvFBO) {
        glDeleteFramebuffers(2, smoothCurvFBO);
        glDeleteTextures(2, smoothCurvatureTexture);
    }
    if (timestepsFBO) {
        glDeleteFramebuffers(3, timestepsFBO);
        glDeleteTextures(3, timestepsTexture);
    }
    // delete singularity texture
    for (unsigned int i = 0; i < singularityTexture.size(); ++i)
        glDeleteTextures(1, &singularityTexture[i]);
    // release
    this->sphereShader.Release();
    this->sphericalTriangleShader.Release();
    this->torusShader.Release();
    this->lightShader.Release();
    this->Release();
}

void MoleculeSESRenderer::release(void) {}

bool MoleculeSESRenderer::loadShader(
    vislib::graphics::gl::GLSLShader& Shader, std::string vertex, std::string fragment) {
    using namespace vislib::graphics::gl;

    ShaderSource compSrc;
    ShaderSource vertSrc;
    ShaderSource geomSrc;
    ShaderSource fragSrc;

    CoreInstance* ci = this->GetCoreInstance();

    std::string msg;
    if (!ci->ShaderSourceFactory().MakeShaderSource(vertex.c_str(), vertSrc)) {
        std::ostringstream stream;
        stream << this->ClassName() << ": Unable to load vertex shader source: " << vertex;
        msg = stream.str();
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, msg.c_str());
        return false;
    }
    if (!ci->ShaderSourceFactory().MakeShaderSource(fragment.c_str(), fragSrc)) {
        std::ostringstream fragStream;
        fragStream << this->ClassName() << ": Unable to load fragment shader source: " << fragment;
        msg = fragStream.str();
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, msg.c_str());
        return false;
    }

    try {
        if (!Shader.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        std::ostringstream exStream;
        exStream << this->ClassName() << ": Unable to create shader programm from " << vertex << " and " << fragment
                 << ": " << e.GetMsgA();
        msg = exStream.str();
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, msg.c_str());
        return false;
    }
    return true;
}

bool MoleculeSESRenderer::create(void) {
    if (!ogl_IsVersionGEQ(2, 0) || !areExtsAvailable("GL_EXT_framebuffer_object GL_ARB_texture_float"))
        return false;

    if (!vislib::graphics::gl::GLSLShader::InitialiseExtensions())
        return false;

    cur_timestep = 0;
    // glEnable( GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);

    float spec[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);

    // original MoleculeSESRenderer shaders
    if (!this->loadShader(this->sphereShader, "protein::ses::sphereVertex", "protein::ses::sphereFragment"))
        return false;
    if (!this->loadShader(this->sphereShaderOR, "protein::ses::sphereVertex", "protein::ses::sphereFragmentOR"))
        return false;
    if (!this->loadShader(this->torusShader, "protein::ses::torusVertex", "protein::ses::torusFragment"))
        return false;
    if (!this->loadShader(this->torusShaderOR, "protein::ses::torusVertex", "protein::ses::torusFragmentOR"))
        return false;
    if (!this->loadShader(this->sphericalTriangleShader, "protein::ses::sphericaltriangleVertex",
            "protein::ses::sphericaltriangleFragment"))
        return false;
    if (!this->loadShader(this->sphericalTriangleShaderOR, "protein::ses::sphericaltriangleVertex",
            "protein::ses::sphericaltriangleFragmentOR"))
        return false;

    // shaders for contour drawing
    if (!this->loadShader(this->SC_Shader, "contours::vertex", "contours::contours::SC"))
        return false;
    if (!this->loadShader(this->SC_Curvature_Shader, "contours::vertex", "contours::contours::SC_Curvature"))
        return false;
    if (!this->loadShader(this->C_Shader, "contours::vertex", "contours::contours::C"))
        return false;
    if (!this->loadShader(this->C_Curvature_Shader, "contours::vertex", "contours::contours::C_Curvature"))
        return false;
    // if (!this->loadShader(this->SCfromCurvatureShader, "contours::vertex", "contours::curvature::fragment"))
    //     return false;
    if (!this->loadShader(this->curvatureShader, "contours::vertex", "contours::curvature::evans"))
        return false;
    if (!this->loadShader(this->normalCurvatureShader, "contours::vertex", "contours::curvature::normal"))
        return false;
    if (!this->loadShader(this->nathanReedCurvatureShader, "contours::vertex", "contours::curvature::nathanReed"))
        return false;
    if (!this->loadShader(this->meanCurvatureShader, "contours::vertex", "contours::curvature::mean"))
        return false;
    if (!this->loadShader(this->prantlMeanShader, "contours::vertex", "contours::curvature::prantlmean"))
        return false;
    if (!this->loadShader(this->prantl2MeanShader, "contours::vertex", "contours::curvature::prantl2mean"))
        return false;
    if (!this->loadShader(this->prantlGaussianShader, "contours::vertex", "contours::curvature::prantlgaussian"))
        return false;
    if (!this->loadShader(this->prantl2GaussianShader, "contours::vertex", "contours::curvature::prantl2gaussian"))
        return false;
    if (!this->loadShader(this->prantlRadialShader, "contours::vertex", "contours::curvature::prantlradial"))
        return false;
    if (!this->loadShader(this->prantl2RadialShader, "contours::vertex", "contours::curvature::prantl2radial"))
        return false;
    if (!this->loadShader(this->passThroughShader, "contours::vertex", "contours::passThrough"))
        return false;
    if (!this->loadShader(this->normalizePositionsShader, "contours::vertex", "contours::normalizePositions"))
        return false;
    if (!this->loadShader(this->gaussianBlurShader, "contours::vertex", "contours::postprocessing::blur"))
        return false;
    if (!this->loadShader(
            this->peronaMalikBlurShader, "contours::vertex", "contours::postprocessing::perona-malik-blur"))
        return false;
    if (!this->loadShader(this->depthBlurShader, "contours::vertex", "contours::postprocessing::depthBlur"))
        return false;
    if (!this->loadShader(this->depthGaussBlurShader, "contours::vertex", "contours::postprocessing::gauss-depth-blur"))
        return false;
    if (!this->loadShader(
            this->perspectiveTestCaseShader, "testCase_perspective::vertex", "testCase_perspective::fragment"))
        return false;
    if (!this->loadShader(this->dilationShader, "contours::vertex", "distance::dilation"))
        return false;
    if (!this->loadShader(this->erosionShader, "contours::vertex", "distance::erosion"))
        return false;
    if (!this->loadShader(this->medianShader, "contours::vertex", "distance::median"))
        return false;
    if (!this->loadShader(this->smoothTimestepsShader, "contours::vertex", "timesteps::fragment"))
        return false;

    // // Load test case shader
    using namespace vislib::graphics::gl;

    ShaderSource vertSrc;
    ShaderSource geomSrc;
    ShaderSource fragSrc;

    CoreInstance* ci = this->GetCoreInstance();
    std::string vertex = "testCase::vertex";
    std::string geometry = "testCase::geometry";
    std::string fragment = "testCase::fragment";

    std::string msg;
    if (!ci->ShaderSourceFactory().MakeShaderSource(vertex.c_str(), vertSrc)) {
        std::ostringstream stream;
        stream << this->ClassName() << ": Unable to load vertex shader source: " << vertex;
        msg = stream.str();
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, msg.c_str());
        return false;
    }
    if (!ci->ShaderSourceFactory().MakeShaderSource(geometry.c_str(), geomSrc)) {
        std::ostringstream stream;
        stream << this->ClassName() << ": Unable to load geometry shader source: " << geometry;
        msg = stream.str();
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, msg.c_str());
        return false;
    }
    if (!ci->ShaderSourceFactory().MakeShaderSource(fragment.c_str(), fragSrc)) {
        std::ostringstream fragStream;
        fragStream << this->ClassName() << ": Unable to load fragment shader source: " << fragment;
        msg = fragStream.str();
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, msg.c_str());
        return false;
    }

    try {
        if (!testCaseShader.Compile(
                vertSrc.Code(), vertSrc.Count(), geomSrc.Code(), geomSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
        testCaseShader.Link();
    } catch (vislib::Exception e) {
        std::ostringstream exStream;
        exStream << this->ClassName() << ": Unable to create shader programm from " << vertex << ", " << geometry
                 << " and " << fragment << ": " << e.GetMsgA();
        msg = exStream.str();
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, msg.c_str());
        return false;
    }
    return true;
}

bool MoleculeSESRenderer::GetExtents(view::CallRender3DGL& call) {

    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    if (mol == NULL)
        return false;
    if (!(*mol)(1))
        return false;

    call.AccessBoundingBoxes() = mol->AccessBoundingBoxes();
    call.SetTimeFramesCount(mol->FrameCount());
    return true;
}

bool MoleculeSESRenderer::Render(view::CallRender3DGL& call) {

#pragma region // Set up camera variables
    // temporary variables
    unsigned int cntRS = 0;

    // get camera information
    this->cameraInfo = call.GetCamera();
    cam_type::snapshot_type snapshot;
    cam_type::matrix_type viewTemp, projTemp;
    cameraInfo.calc_matrices(snapshot, viewTemp, projTemp, thecam::snapshot_content::all);
    auto resolution = cameraInfo.resolution_gate();
    glm::mat4 view = viewTemp;
    glm::mat4 proj = projTemp;

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glLoadMatrixf(glm::value_ptr(proj));

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glLoadMatrixf(glm::value_ptr(view));

#pragma endregion // Set up camera variables

    float callTime = call.Time();

#pragma region // Get Protein an Binding Site data
    // get pointer to CallProteinData
    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    // if something went wrong --> return
    if (!mol)
        return false;

    // execute the call
    mol->SetFrameID(static_cast<int>(callTime));
    if (!(*mol)(MolecularDataCall::CallForGetData))
        return false;

    // get pointer to BindingSiteCall
    BindingSiteCall* bs = this->bsDataCallerSlot.CallAs<BindingSiteCall>();
    if (bs) {
        (*bs)(BindingSiteCall::CallForGetData);
    }

#pragma endregion // Get Protein an Binding Site data

    // check parameter
    this->UpdateParameters(mol, bs);

    // ==================== Precomputations ====================


    this->probeRadius = this->probeRadiusSlot.Param<param::FloatParam>()->Value();

    // init the reduced surfaces
    if (this->reducedSurface.empty()) {
        time_t t = clock();
        // create the reduced surface
        unsigned int chainIds;
        if (!this->computeSesPerMolecule) {
            this->reducedSurface.push_back(new ReducedSurface(mol, this->probeRadius));
            this->reducedSurface.back()->ComputeReducedSurface();
        } else {
            // if no molecule indices are given, compute the SES for all molecules
            if (this->molIdxList.IsEmpty()) {
                for (chainIds = 0; chainIds < mol->MoleculeCount(); ++chainIds) {
                    this->reducedSurface.push_back(new ReducedSurface(chainIds, mol, this->probeRadius));
                    this->reducedSurface.back()->ComputeReducedSurface();
                }
            } else {
                // else compute the SES for all selected molecules
                for (chainIds = 0; chainIds < this->molIdxList.Count(); ++chainIds) {
                    this->reducedSurface.push_back(
                        new ReducedSurface(atoi(this->molIdxList[chainIds]), mol, this->probeRadius));
                    this->reducedSurface.back()->ComputeReducedSurface();
                }
            }
        }
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_INFO,
            "%s: RS computed in: %f s\n", this->ClassName(), (double(clock() - t) / double(CLOCKS_PER_SEC)));
    }
    // update the data / the RS
    for (cntRS = 0; cntRS < this->reducedSurface.size(); ++cntRS) {
        if (this->reducedSurface[cntRS]->UpdateData(1.0f, 5.0f)) {
            this->ComputeRaycastingArrays(cntRS);
        }
    }

    if (!this->preComputationDone) {
        // compute the color table
        Color::MakeColorTable(mol, this->currentColoringMode0, this->currentColoringMode1,
            this->cmWeightParam.Param<param::FloatParam>()->Value(),        // weight for the first cm
            1.0f - this->cmWeightParam.Param<param::FloatParam>()->Value(), // weight for the second cm
            this->atomColorTable, this->colorLookupTable, this->rainbowColors,
            this->minGradColorParam.Param<param::StringParam>()->Value(),
            this->midGradColorParam.Param<param::StringParam>()->Value(),
            this->maxGradColorParam.Param<param::StringParam>()->Value(), true, bs);
        this->ComputeRaycastingArrays();
        // set the precomputation of the data as done
        this->preComputationDone = true;
    }

    bool virtualViewportChanged = false;
    if (static_cast<unsigned int>(resolution.width()) != this->width ||
        static_cast<unsigned int>(resolution.height()) != this->height) {
        this->width = static_cast<unsigned int>(resolution.width());
        this->height = static_cast<unsigned int>(resolution.height());
        virtualViewportChanged = true;
    }

    if (virtualViewportChanged) {
        // TODO: All this stuff should only be created if necessary
        normalPyramid.create("fragNormal", this->width, this->height, this->GetCoreInstance(), "pullpush::pullNormal",
            "pullpush::pushNormal");
        positionPyramid.create("fragPosition", this->width, this->height, this->GetCoreInstance(),
            "pullpush::pullPosition", "pullpush::pushPosition");
        depthPyramid.create(
            "fragMaxDepth", this->width, this->height, this->GetCoreInstance(), "pullpush::pullMaxDepth");
        heightPyramid.create("fragMaxY", this->width, this->height, this->GetCoreInstance(), "pullpush::pullMaxY");
        widthPyramid.create("fragMaxX", this->width, this->height, this->GetCoreInstance(), "pullpush::pullMaxX");
        SCpyramid.create(
            "outData", this->width, this->height, this->GetCoreInstance(), "pullpush::pullSC", "pullpush::pushSC");
        curvaturePyramid.create("fragCurvature", this->width, this->height, this->GetCoreInstance(),
            "pullpush::pullCurvature", "pullpush::pushCurvature");
        this->CreateQuadBuffers();
        this->CreateFBO();
    }

    // Scale & Translate
    glPushMatrix();

    // ==================== Start actual rendering ====================

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);

    // get clear color (i.e. background color
    float* clearColor = new float[4];
    glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColor);
    clear_color = glm::vec4(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);
    // delete pointers
    delete[] clearColor;


    // start rendering to frame buffer object
    if (offscreenRendering) {
        glBindFramebuffer(GL_FRAMEBUFFER, this->contourFBO);
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

    // if (!protein::computeRMSF(mol)) {
    //     std::cout << "Computation of RMSF failed!" << std::endl;
    //     return false;
    // }
    float* bfactors = mol->AtomBFactors();

    if (testcase && orthoproj)
        this->RenderTestCase();
    else
        this->RenderSESGpuRaycasting(mol);

    if (offscreenRendering) {

        if (this->smoothPositions) {
            this->SmoothPositions(*blurShaderMap[currentBlurMode]);
        }
        if (this->smoothNormals) {
            this->SmoothNormals(*blurShaderMap[currentBlurMode]);
        }

        if (this->currentDisplayedProperty == Position) {
            this->displayPositions();
        } else if (this->currentDisplayedProperty == NormalizedPosition) {
            this->displayNormalizedPositions();
        } else if (this->currentDisplayedProperty == Normal) {
            this->displayNormals();
        } else if (this->currentDisplayedProperty == Curvature) {
            renderCurvature(*curvatureShaderMap[currentCurvatureMode]);
        } else {
            if (this->currentContourMode == Suggestive || this->currentContourMode == SuggestiveAndCurvature)
                this->SuggestiveContours(*contourShaderMap[currentContourMode]);
            else
                this->Contours(*contourShaderMap[currentContourMode]);
        }
        // stop rendering to frame buffer object
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    glPopMatrix();

    // unlock the current frame
    mol->Unlock();

#pragma region // do some matrix stuff
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
#pragma endregion // do some matrix stuff

    glActiveTexture(GL_TEXTURE0);

    return true;
}

void MoleculeSESRenderer::UpdateParameters(const MolecularDataCall* mol, const BindingSiteCall* bs) {

    // variables
    bool recomputeColors = false;
    // ==================== check parameters ====================
    if (this->coloringModeParam0.IsDirty() || this->coloringModeParam1.IsDirty() || this->cmWeightParam.IsDirty()) {
        this->currentColoringMode0 =
            static_cast<Color::ColoringMode>(this->coloringModeParam0.Param<param::EnumParam>()->Value());
        this->currentColoringMode1 =
            static_cast<Color::ColoringMode>(this->coloringModeParam1.Param<param::EnumParam>()->Value());

        Color::MakeColorTable(mol, this->currentColoringMode0, this->currentColoringMode1,
            this->cmWeightParam.Param<param::FloatParam>()->Value(),        // weight for the first cm
            1.0f - this->cmWeightParam.Param<param::FloatParam>()->Value(), // weight for the second cm
            this->atomColorTable, this->colorLookupTable, this->rainbowColors,
            this->minGradColorParam.Param<param::StringParam>()->Value(),
            this->midGradColorParam.Param<param::StringParam>()->Value(),
            this->maxGradColorParam.Param<param::StringParam>()->Value(), true, bs);

        this->preComputationDone = false;
        this->coloringModeParam0.ResetDirty();
        this->coloringModeParam1.ResetDirty();
        this->cmWeightParam.ResetDirty();
    }
    if (this->curvatureModeParam.IsDirty()) {
        this->currentCurvatureMode =
            static_cast<curvatureMode>(this->curvatureModeParam.Param<param::EnumParam>()->Value());
        this->curvatureModeParam.ResetDirty();
    }
    if (this->contourModeParam.IsDirty()) {
        this->currentContourMode = static_cast<contourMode>(this->contourModeParam.Param<param::EnumParam>()->Value());
        this->contourModeParam.ResetDirty();
    }
    if (this->blurParam.IsDirty()) {
        this->currentBlurMode = static_cast<blurMode>(this->blurParam.Param<param::EnumParam>()->Value());
        this->blurParam.ResetDirty();
    }
    if (this->displayedPropertyParam.IsDirty()) {
        this->currentDisplayedProperty =
            static_cast<displayedProperty>(this->displayedPropertyParam.Param<param::EnumParam>()->Value());
        this->displayedPropertyParam.ResetDirty();
    }
    if (this->drawSESParam.IsDirty()) {
        this->drawSES = this->drawSESParam.Param<param::BoolParam>()->Value();
        this->drawSESParam.ResetDirty();
    }
    if (this->drawSASParam.IsDirty()) {
        this->drawSAS = this->drawSASParam.Param<param::BoolParam>()->Value();
        this->drawSASParam.ResetDirty();
        this->preComputationDone = false;
    }
    if (this->offscreenRenderingParam.IsDirty()) {
        this->offscreenRendering = this->offscreenRenderingParam.Param<param::BoolParam>()->Value();
        this->offscreenRenderingParam.ResetDirty();
    }
    if (this->molIdxListParam.IsDirty()) {
        vislib::StringA tmpStr(this->molIdxListParam.Param<param::StringParam>()->Value());
        this->molIdxList = vislib::StringTokeniser<vislib::CharTraitsA>::Split(tmpStr, ';', true);
        this->molIdxListParam.ResetDirty();
    }
    // color table param
    if (this->colorTableFileParam.IsDirty()) {
        Color::ReadColorTableFromFile(
            this->colorTableFileParam.Param<param::StringParam>()->Value(), this->colorLookupTable);
        this->colorTableFileParam.ResetDirty();
        recomputeColors = true;
    }
    if (this->probeRadiusSlot.IsDirty()) {
        this->probeRadius = this->probeRadiusSlot.Param<param::FloatParam>()->Value();
        this->reducedSurface.clear();
        this->preComputationDone = false;
        this->probeRadiusSlot.ResetDirty();
    }
    if (this->smoothNormalsParam.IsDirty()) {
        this->smoothNormals = this->smoothNormalsParam.Param<param::BoolParam>()->Value();
        this->smoothNormalsParam.ResetDirty();
    }
    if (this->smoothCurvatureParam.IsDirty()) {
        this->smoothCurvature = this->smoothCurvatureParam.Param<param::BoolParam>()->Value();
        this->smoothCurvatureParam.ResetDirty();
    }
    if (this->smoothPositionsParam.IsDirty()) {
        this->smoothPositions = this->smoothPositionsParam.Param<param::BoolParam>()->Value();
        this->smoothPositionsParam.ResetDirty();
    }
    // if (this->pyramidWeightsParam.IsDirty()) {
    //     this->pyramidWeight = this->pyramidWeightsParam.Param<param::FloatParam>()->Value();
    //     this->pyramidWeightsParam.ResetDirty();
    // }
    // if (this->pyramidLayersParam.IsDirty()) {
    //     this->pyramidLayers = this->pyramidLayersParam.Param<param::IntParam>()->Value();
    //     this->pyramidLayersParam.ResetDirty();
    // }
    // if (this->pyramidGammaParam.IsDirty()) {
    //     this->pyramidGamma = this->pyramidGammaParam.Param<param::FloatParam>()->Value();
    //     this->pyramidGammaParam.ResetDirty();
    // }
    if (this->SCRadiusParam.IsDirty()) {
        this->SCRadius = this->SCRadiusParam.Param<param::IntParam>()->Value();
        this->SCRadiusParam.ResetDirty();
    }
    if (this->SCNeighbourThresholdParam.IsDirty()) {
        this->SCNeighbourThreshold = this->SCNeighbourThresholdParam.Param<param::FloatParam>()->Value();
        this->SCNeighbourThresholdParam.ResetDirty();
    }
    if (this->SCDiffThresholdParam.IsDirty()) {

        this->SCDiffThreshold = this->SCDiffThresholdParam.Param<param::FloatParam>()->Value();
        this->SCDiffThresholdParam.ResetDirty();
    }
    if (this->SCMedianFilterParam.IsDirty()) {
        this->SCMedianFilter = this->SCMedianFilterParam.Param<param::BoolParam>()->Value();
        this->SCMedianFilterParam.ResetDirty();
    }
    if (this->SCCircularNeighborhoodParam.IsDirty()) {
        this->SCCircularNeighborhood = this->SCCircularNeighborhoodParam.Param<param::BoolParam>()->Value();
        this->SCCircularNeighborhoodParam.ResetDirty();
    }
    if (this->ViewTypeParam.IsDirty()) {
        this->currentViewType = static_cast<viewType>(this->ViewTypeParam.Param<param::EnumParam>()->Value());
        this->ViewTypeParam.ResetDirty();
    }
    if (this->OrthoProjParam.IsDirty()) {
        this->orthoproj = this->OrthoProjParam.Param<param::BoolParam>()->Value();
        this->OrthoProjParam.ResetDirty();
    }
    if (this->TestCaseParam.IsDirty()) {
        this->testcase = this->TestCaseParam.Param<param::BoolParam>()->Value();
        this->TestCaseParam.ResetDirty();
    }
    if (this->cutOffParam.IsDirty()) {
        this->cutOff = this->cutOffParam.Param<param::FloatParam>()->Value();
        this->cutOffParam.ResetDirty();
    }
    if (this->numBlurParam.IsDirty()) {
        this->numBlur = this->numBlurParam.Param<param::FloatParam>()->Value();
        this->numBlurParam.ResetDirty();
    }
    if (this->numPosBlurParam.IsDirty()) {
        this->numPosBlur = this->numPosBlurParam.Param<param::FloatParam>()->Value();
        this->numPosBlurParam.ResetDirty();
    }
    if (this->numCurvBlurParam.IsDirty()) {
        this->numCurvBlur = this->numCurvBlurParam.Param<param::FloatParam>()->Value();
        this->numCurvBlurParam.ResetDirty();
    }
    if (this->depthDiffParam.IsDirty()) {
        this->depthDiff = this->depthDiffParam.Param<param::FloatParam>()->Value();
        this->depthDiffParam.ResetDirty();
    }
    if (this->whiteBackgroundParam.IsDirty()) {
        this->whiteBackground = this->whiteBackgroundParam.Param<param::BoolParam>()->Value();
        this->whiteBackgroundParam.ResetDirty();
    }
    if (this->extendContoursParam.IsDirty()) {
        this->extendContoursBool = this->extendContoursParam.Param<param::BoolParam>()->Value();
        this->extendContoursParam.ResetDirty();
    }
    if (this->dilation2RadiusParam.IsDirty()) {
        this->dilation2Radius = this->dilation2RadiusParam.Param<param::IntParam>()->Value();
        this->dilation2RadiusParam.ResetDirty();
    }
    if (this->dilation1RadiusParam.IsDirty()) {
        this->dilation1Radius = this->dilation1RadiusParam.Param<param::IntParam>()->Value();
        this->dilation1RadiusParam.ResetDirty();
    }
    if (this->erosionRadiusParam.IsDirty()) {
        this->erosionRadius = this->erosionRadiusParam.Param<param::IntParam>()->Value();
        this->erosionRadiusParam.ResetDirty();
    }
    if (this->smoothTimestepsParam.IsDirty()) {
        this->smoothTimestepsBool = this->smoothTimestepsParam.Param<param::BoolParam>()->Value();
        this->smoothTimestepsParam.ResetDirty();
    }
    if (recomputeColors) {
        this->preComputationDone = false;
    }
}
void MoleculeSESRenderer::calculateTextureBBX() {
    /*
     * Find max/min depth using pull phase of pull-push algorithm
     */
    depthPyramid.pullShaderProgram.Enable();
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, positionTexture);
    glUniform1i(depthPyramid.pullShaderProgram.ParameterLocation("inputTex_fragPosition"), 1);

    depthPyramid.clear();
    depthPyramid.pull();
    /*
     * Find max/min height using pull phase of pull-push algorithm
     */
    heightPyramid.pullShaderProgram.Enable();
    glUniform1i(heightPyramid.pullShaderProgram.ParameterLocation("inputTex_fragPosition"), 1);

    heightPyramid.clear();
    heightPyramid.pull();

    /*
     * Find max/min width using pull phase of pull-push algorithm
     */
    widthPyramid.pullShaderProgram.Enable();
    glUniform1i(widthPyramid.pullShaderProgram.ParameterLocation("inputTex_fragPosition"), 1);

    widthPyramid.clear();
    widthPyramid.pull();

    this->bbx_levelMax = widthPyramid.getMipmapNumber();
}

void MoleculeSESRenderer::SmoothNormals(vislib::graphics::gl::GLSLShader& Shader) {
    glDisable(GL_DEPTH_TEST);
    Shader.Enable();
    glUniform1i(Shader.ParameterLocation("positionTexture"), 1);
    glUniform1i(Shader.ParameterLocation("screenTexture"), 0);
    glUniform1f(Shader.ParameterLocation("depthdiff"), this->depthDiff);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, positionTexture);
    glActiveTexture(GL_TEXTURE0);

    normal_horizontal = true;
    bool first_iteration = true;
    for (unsigned int i = 0; i < this->numBlur; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, normalFBO[normal_horizontal]);
        glUniform1i(Shader.ParameterLocation("horizontal"), normal_horizontal);
        glBindTexture(GL_TEXTURE_2D, first_iteration ? normalTexture : smoothNormalTexture[!normal_horizontal]);
        glBindVertexArray(quadVAO);
        glClear(GL_COLOR_BUFFER_BIT);

        glDrawArrays(GL_TRIANGLES, 0, 6);
        normal_horizontal = !normal_horizontal;
        if (first_iteration)
            first_iteration = false;
    }
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    Shader.Disable();
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    // /*
    //  * Execute Pull-Push algorithm for smoothing
    //  */
    // glDisable(GL_DEPTH_TEST);
    // this->calculateTextureBBX();
    // normalPyramid.pullShaderProgram.Enable();
    // glActiveTexture(GL_TEXTURE1);
    // glBindTexture(GL_TEXTURE_2D, this->normalTexture);
    // glUniform1i(normalPyramid.pullShaderProgram.ParameterLocation("inputTex_fragNormal"), 1);
    // glActiveTexture(GL_TEXTURE2);
    // glBindTexture(GL_TEXTURE_2D, this->positionTexture);
    // glUniform1i(normalPyramid.pullShaderProgram.ParameterLocation("inputTex_fragPosition"), 2);
    // glActiveTexture(GL_TEXTURE3);
    // glBindTexture(GL_TEXTURE_2D, depthPyramid.get("fragMaxDepth"));
    // glUniform1i(normalPyramid.pullShaderProgram.ParameterLocation("maxDepth_texture"), 3);
    // glUniform1f(normalPyramid.pullShaderProgram.ParameterLocation("weightFactor"), this->pyramidWeight);
    // normalPyramid.pushShaderProgram.Enable();
    // glUniform1f(normalPyramid.pushShaderProgram.ParameterLocation("gamma"), this->pyramidGamma);

    // normalPyramid.clear();
    // normalPyramid.pull_until(this->pyramidLayers);
    // normalPyramid.push_from(this->pyramidLayers);
    // glEnable(GL_DEPTH_TEST);
}
void MoleculeSESRenderer::SmoothCurvature(vislib::graphics::gl::GLSLShader& Shader) {
    glDisable(GL_DEPTH_TEST);
    Shader.Enable();
    glUniform1i(Shader.ParameterLocation("positionTexture"), 1);
    glUniform1i(Shader.ParameterLocation("screenTexture"), 0);
    glUniform1f(Shader.ParameterLocation("depthdiff"), this->depthDiff);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, positionTexture);
    glActiveTexture(GL_TEXTURE0);

    curv_horizontal = true;
    bool first_iteration = true;
    for (unsigned int i = 0; i < this->numCurvBlur; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, smoothCurvFBO[curv_horizontal]);
        glUniform1i(Shader.ParameterLocation("horizontal"), curv_horizontal);
        glBindTexture(GL_TEXTURE_2D, first_iteration ? curvatureTexture : smoothCurvatureTexture[!curv_horizontal]);
        glBindVertexArray(quadVAO);
        glClear(GL_COLOR_BUFFER_BIT);

        glDrawArrays(GL_TRIANGLES, 0, 6);
        curv_horizontal = !curv_horizontal;
        if (first_iteration)
            first_iteration = false;
    }
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    Shader.Disable();
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void MoleculeSESRenderer::SmoothPositions(vislib::graphics::gl::GLSLShader& Shader) {
    /*
     * Execute Pull-Push algorithm for smoothing
     */
    glDisable(GL_DEPTH_TEST);
    Shader.Enable();
    glActiveTexture(GL_TEXTURE1);
    glUniform1i(Shader.ParameterLocation("positionTexture"), 1);
    glBindTexture(GL_TEXTURE_2D, positionTexture);
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(Shader.ParameterLocation("screenTexture"), 0);

    pos_horizontal = true;
    bool first_iteration = true;
    for (unsigned int i = 0; i < this->numPosBlur; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, positionFBO[pos_horizontal]);
        glUniform1i(Shader.ParameterLocation("horizontal"), pos_horizontal);
        glBindTexture(GL_TEXTURE_2D, first_iteration ? positionTexture : smoothPositionTexture[!pos_horizontal]);
        glBindVertexArray(quadVAO);
        glClear(GL_COLOR_BUFFER_BIT);

        glDrawArrays(GL_TRIANGLES, 0, 6);
        pos_horizontal = !pos_horizontal;
        if (first_iteration)
            first_iteration = false;
    }
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    Shader.Disable();
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    // this->calculateTextureBBX();
    // positionPyramid.pullShaderProgram.Enable();
    // glActiveTexture(GL_TEXTURE2);
    // glBindTexture(GL_TEXTURE_2D, positionTexture);
    // glUniform1i(positionPyramid.pullShaderProgram.ParameterLocation("inputTex_fragPosition"), 2);
    // glActiveTexture(GL_TEXTURE3);
    // glBindTexture(GL_TEXTURE_2D, depthPyramid.get("fragMaxDepth"));
    // glUniform1i(positionPyramid.pullShaderProgram.ParameterLocation("maxDepth_texture"), 3);
    // glUniform1f(positionPyramid.pullShaderProgram.ParameterLocation("weightFactor"), this->pyramidWeight);
    // positionPyramid.pushShaderProgram.Enable();
    // glUniform1f(positionPyramid.pushShaderProgram.ParameterLocation("gamma"), this->pyramidGamma);

    // positionPyramid.clear();
    // positionPyramid.pull_until(this->pyramidLayers);
    // positionPyramid.push_from(this->pyramidLayers);
    // glEnable(GL_DEPTH_TEST);
}

void MoleculeSESRenderer::SuggestiveContours(vislib::graphics::gl::GLSLShader& Shader) {


    calculateCurvature(*curvatureShaderMap[this->currentCurvatureMode]);
    glDisable(GL_DEPTH_TEST);
    Shader.Enable();
    glUniform1i(Shader.ParameterLocation("radius"), this->SCRadius);
    glUniform1f(Shader.ParameterLocation("neighbourThreshold"), this->SCNeighbourThreshold);
    glUniform1f(Shader.ParameterLocation("intensityDiffThreshold"), this->SCDiffThreshold);
    glUniform1i(Shader.ParameterLocation("medianFilter"), this->SCMedianFilter);
    glUniform1i(Shader.ParameterLocation("circularNeighborhood"), this->SCCircularNeighborhood);
    glUniform1i(Shader.ParameterLocation("viewType"), this->currentViewType);
    glUniform1f(Shader.ParameterLocation("cutOff"), this->cutOff);
    glUniform1i(Shader.ParameterLocation("whiteBackground"), this->whiteBackground);
    glActiveTexture(GL_TEXTURE1);
    if (smoothNormals) {
        glBindTexture(GL_TEXTURE_2D, smoothNormalTexture[!pos_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, normalTexture);
    }
    glUniform1i(Shader.ParameterLocation("normalTexture"), 1);
    glActiveTexture(GL_TEXTURE2);
    if (smoothPositions) {
        glBindTexture(GL_TEXTURE_2D, smoothPositionTexture[!pos_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, positionTexture);
    }
    glUniform1i(Shader.ParameterLocation("positionTexture"), 2);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, curvatureTexture);
    glUniform1i(Shader.ParameterLocation("curvatureTexture"), 3);
    glGetError();
    glBindFramebuffer(GL_FRAMEBUFFER, 1);
    if (this->whiteBackground) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
    } else {
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    Shader.Disable();
}
void MoleculeSESRenderer::Contours(vislib::graphics::gl::GLSLShader& Shader) {

    calculateTextureBBX();
    calculateCurvature(*curvatureShaderMap[this->currentCurvatureMode]);
    if (smoothCurvature) {
        SmoothCurvature(*blurShaderMap[this->currentBlurMode]);
    }
    glDisable(GL_DEPTH_TEST);
    // auto shader = *contourShaderMap[this->currentContourMode];
    Shader.Enable();
    glActiveTexture(GL_TEXTURE1);
    if (smoothNormals) {
        glBindTexture(GL_TEXTURE_2D, smoothNormalTexture[!normal_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, normalTexture);
    }
    // glBindTexture(GL_TEXTURE_2D, *this->cur_normalTexture);
    glUniform1i(Shader.ParameterLocation("normalTexture"), 1);
    glActiveTexture(GL_TEXTURE2);
    if (smoothPositions) {
        glBindTexture(GL_TEXTURE_2D, smoothPositionTexture[!pos_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, positionTexture);
    }
    // glBindTexture(GL_TEXTURE_2D, *this->cur_positionTexture);
    glUniform1i(Shader.ParameterLocation("positionTexture"), 2);
    glActiveTexture(GL_TEXTURE3);
    if (smoothCurvature) {
        glBindTexture(GL_TEXTURE_2D, smoothCurvatureTexture[!curv_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, curvatureTexture);
    }
    glUniform1i(Shader.ParameterLocation("curvatureTexture"), 3);
    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_2D, depthPyramid.get("fragMaxDepth"));
    glUniform1i(Shader.ParameterLocation("depthTexture"), 4);
    glUniform1f(Shader.ParameterLocation("cutOff"), cutOff);
    glUniform1i(Shader.ParameterLocation("viewType"), this->currentViewType);
    glUniform1i(Shader.ParameterLocation("orthoproj"), this->orthoproj);
    float near_plane = this->cameraInfo.near_clipping_plane();
    glUniform1f(Shader.ParameterLocation("near_plane"), near_plane);
    glUniform1i(Shader.ParameterLocation("level_max"), this->bbx_levelMax);
    glUniform1i(Shader.ParameterLocation("whiteBackground"), this->whiteBackground);
    if (!extendContoursBool && !smoothTimestepsBool) {
        glBindFramebuffer(GL_FRAMEBUFFER, 1);
        if (whiteBackground) {
            glClearColor(1.0, 1.0, 1.0, 0.0);
        } else {
            glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        }
    } else if (extendContoursBool) {
        glBindFramebuffer(GL_FRAMEBUFFER, extendContourFBO[0]);
        glClearColor(0.0, 0.0, 0.0, 0.0);
    } else {
        auto test = cur_timestep % 3;
        glBindFramebuffer(GL_FRAMEBUFFER, timestepsFBO[cur_timestep % 3]);
        glClearColor(0.0, 0.0, 0.0, 0.0);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    Shader.Disable();
    if (extendContoursBool) {
        extendContours();
    } else if (smoothTimestepsBool) {
        smoothTimesteps();
    }
}
void MoleculeSESRenderer::extendContours() {

    glDisable(GL_DEPTH_TEST);

    // dilation
    dilationShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, contourTexture[0]);
    glUniform1i(dilationShader.ParameterLocation("contourTexture"), 0);
    glUniform1i(dilationShader.ParameterLocation("whiteBackground"), this->whiteBackground);
    glUniform1i(dilationShader.ParameterLocation("radius"), this->dilation1Radius);
    glBindFramebuffer(GL_FRAMEBUFFER, extendContourFBO[1]);
    glClearColor(0, 0, 0, 0);
    // glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    dilationShader.Disable();

    // erosion
    erosionShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, contourTexture[1]);
    glUniform1i(erosionShader.ParameterLocation("contourTexture"), 0);
    glUniform1i(erosionShader.ParameterLocation("whiteBackground"), this->whiteBackground);
    glUniform1i(erosionShader.ParameterLocation("radius"), this->erosionRadius);
    glBindFramebuffer(GL_FRAMEBUFFER, extendContourFBO[0]);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    erosionShader.Disable();

    // erosion
    erosionShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, contourTexture[0]);
    glUniform1i(erosionShader.ParameterLocation("contourTexture"), 0);
    glUniform1i(erosionShader.ParameterLocation("whiteBackground"), this->whiteBackground);
    glUniform1i(erosionShader.ParameterLocation("radius"), this->erosionRadius);
    glBindFramebuffer(GL_FRAMEBUFFER, extendContourFBO[1]);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    erosionShader.Disable();

    // dilation
    dilationShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, contourTexture[1]);
    glUniform1i(dilationShader.ParameterLocation("contourTexture"), 0);
    glUniform1i(dilationShader.ParameterLocation("whiteBackground"), this->whiteBackground);
    glUniform1i(dilationShader.ParameterLocation("radius"), this->dilation2Radius);
    if (smoothTimestepsBool) {
        glBindFramebuffer(GL_FRAMEBUFFER, timestepsFBO[cur_timestep % 3]);
        glClearColor(0.0, 0.0, 0.0, 0.0);
    } else {
        glBindFramebuffer(GL_FRAMEBUFFER, 1);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    dilationShader.Disable();
    if (smoothTimestepsBool) {
        smoothTimesteps();
    }
}
void MoleculeSESRenderer::smoothTimesteps() {
    if (cur_timestep < 2) {
        cur_timestep++;
        return;
    }
    smoothTimestepsShader.Enable();

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, timestepsTexture[(cur_timestep - 2) % 3]);
    glUniform1i(smoothTimestepsShader.ParameterLocation("first"), 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, timestepsTexture[(cur_timestep - 1) % 3]);
    glUniform1i(smoothTimestepsShader.ParameterLocation("second"), 1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, timestepsTexture[cur_timestep % 3]);
    glUniform1i(smoothTimestepsShader.ParameterLocation("third"), 2);

    // glUniform1i(smoothTimestepsShader.ParameterLocation("whiteBackground"), this->whiteBackground);
    glBindFramebuffer(GL_FRAMEBUFFER, 1);
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
    smoothTimestepsShader.Disable();
    cur_timestep++;

    passThroughShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(passThroughShader.ParameterLocation("screenTexture"), 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
    passThroughShader.Disable();
}
void MoleculeSESRenderer::displayPositions() {

    glDisable(GL_DEPTH_TEST);

    passThroughShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    if (smoothPositions) {
        glBindTexture(GL_TEXTURE_2D, smoothPositionTexture[!pos_horizontal]); // positionPyramid.get("fragPosition"));
    } else {
        glBindTexture(GL_TEXTURE_2D, positionTexture);
    }
    // glBindTexture(GL_TEXTURE_2D, *this->cur_positionTexture);
    glUniform1i(passThroughShader.ParameterLocation("screenTexture"), 0);
    glGetError();


    glBindFramebuffer(GL_FRAMEBUFFER, 1);
    if (this->whiteBackground) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
    } else {
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
}
void MoleculeSESRenderer::displayNormalizedPositions() {

    this->calculateTextureBBX();
    glDisable(GL_DEPTH_TEST);
    glBindFramebuffer(GL_FRAMEBUFFER, 1);
    if (this->whiteBackground) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
    } else {
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    normalizePositionsShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    if (smoothPositions) {
        glBindTexture(GL_TEXTURE_2D, smoothPositionTexture[!pos_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, positionTexture);
    }
    // glBindTexture(GL_TEXTURE_2D, *cur_positionTexture);
    glUniform1i(normalizePositionsShader.ParameterLocation("positionTexture"), 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, widthPyramid.get("fragMaxX"));
    glUniform1i(normalizePositionsShader.ParameterLocation("widthTexture"), 1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, heightPyramid.get("fragMaxY"));
    glUniform1i(normalizePositionsShader.ParameterLocation("heightTexture"), 2);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, depthPyramid.get("fragMaxDepth"));
    glUniform1i(normalizePositionsShader.ParameterLocation("depthTexture"), 3);
    glUniform1i(normalizePositionsShader.ParameterLocation("level_max"), this->bbx_levelMax);
    glGetError();

    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
    normalizePositionsShader.Disable();
}
void MoleculeSESRenderer::displayNormals() {

    glDisable(GL_DEPTH_TEST);

    passThroughShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    if (smoothNormals) {
        glBindTexture(GL_TEXTURE_2D, smoothNormalTexture[!normal_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, normalTexture);
    }
    // glBindTexture(GL_TEXTURE_2D, *this->cur_normalTexture);
    glUniform1i(passThroughShader.ParameterLocation("screenTexture"), 0);
    glGetError();

    glBindFramebuffer(GL_FRAMEBUFFER, 1);
    if (this->whiteBackground) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
    } else {
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
}
void MoleculeSESRenderer::calculateCurvature(vislib::graphics::gl::GLSLShader& Shader) {

    glDisable(GL_DEPTH_TEST);

    Shader.Enable();
    glActiveTexture(GL_TEXTURE0);
    if (smoothPositions) {
        glBindTexture(GL_TEXTURE_2D, smoothPositionTexture[!pos_horizontal]); // positionPyramid.get("fragPosition"));
    } else {
        glBindTexture(GL_TEXTURE_2D, positionTexture);
    }
    glUniform1i(Shader.ParameterLocation("tex_fragPosition"), 0);
    glActiveTexture(GL_TEXTURE1);
    if (smoothNormals) {
        glBindTexture(GL_TEXTURE_2D, smoothNormalTexture[!normal_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, normalTexture);
    }
    // glBindTexture(GL_TEXTURE_2D, *this->cur_normalTexture);
    glUniform1i(Shader.ParameterLocation("tex_fragNormal"), 1);
    glGetError();

    glBindFramebuffer(GL_FRAMEBUFFER, curvatureFBO);
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
}

void MoleculeSESRenderer::renderCurvature(vislib::graphics::gl::GLSLShader& Shader) {

    this->calculateCurvature(Shader);
    if (smoothCurvature) {
        SmoothCurvature(*blurShaderMap[this->currentBlurMode]);
    }

    glDisable(GL_DEPTH_TEST);
    this->passThroughShader.Enable();
    glActiveTexture(GL_TEXTURE0);
    if (smoothCurvature) {
        glBindTexture(GL_TEXTURE_2D, smoothCurvatureTexture[!curv_horizontal]);
    } else {
        glBindTexture(GL_TEXTURE_2D, curvatureTexture);
    }
    glUniform1i(passThroughShader.ParameterLocation("screenTexture"), 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 1);
    if (this->whiteBackground) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
    } else {
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glEnable(GL_DEPTH_TEST);
    glBindVertexArray(0);
}

void MoleculeSESRenderer::CreateFBO() {
    if (contourFBO) {
        glDeleteFramebuffers(1, &contourFBO);
        glDeleteRenderbuffers(1, &contourDepthRBO);
        glDeleteTextures(1, &normalTexture);
        glDeleteTextures(1, &positionTexture);
        glDeleteTextures(1, &objPositionTexture);
    }
    if (extendContourFBO) {
        glDeleteFramebuffers(2, extendContourFBO);
        glDeleteTextures(2, contourTexture);
    }
    if (timestepsFBO) {
        glDeleteFramebuffers(3, extendContourFBO);
        glDeleteTextures(3, contourTexture);
    }
    if (curvatureFBO) {
        glDeleteFramebuffers(1, &curvatureFBO);
        glDeleteTextures(1, &curvatureTexture);
    }
    if (positionFBO) {
        glDeleteFramebuffers(2, positionFBO);
        glDeleteTextures(2, smoothPositionTexture);
    }
    if (normalFBO) {
        glDeleteFramebuffers(2, normalFBO);
        glDeleteTextures(2, smoothNormalTexture);
    }
    if (smoothCurvFBO) {
        glDeleteFramebuffers(2, smoothCurvFBO);
        glDeleteTextures(2, smoothCurvatureTexture);
    }
    glGenFramebuffers(1, &contourFBO);
    glGenFramebuffers(2, extendContourFBO);
    glGenTextures(2, contourTexture);
    glGenFramebuffers(3, timestepsFBO);
    glGenTextures(3, timestepsTexture);
    glGenFramebuffers(1, &curvatureFBO);
    glGenFramebuffers(2, positionFBO);
    glGenTextures(2, smoothPositionTexture);
    glGenFramebuffers(2, normalFBO);
    glGenTextures(2, smoothNormalTexture);
    glGenFramebuffers(2, smoothCurvFBO);
    glGenTextures(2, smoothCurvatureTexture);
    glGenTextures(1, &normalTexture);
    glGenTextures(1, &curvatureTexture);
    glGenTextures(1, &positionTexture);
    glGenTextures(1, &objPositionTexture);
    glGenRenderbuffers(1, &contourDepthRBO);

    // contour FBO
    glBindFramebuffer(GL_FRAMEBUFFER, this->contourFBO);

    // texture for normals
    glBindTexture(GL_TEXTURE_2D, this->normalTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, normalTexture, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    // texture for positions
    glBindTexture(GL_TEXTURE_2D, this->positionTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, positionTexture, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    // texture for object positions
    glBindTexture(GL_TEXTURE_2D, this->objPositionTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, objPositionTexture, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    // Depth RBO
    glBindRenderbuffer(GL_RENDERBUFFER, contourDepthRBO);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, this->width, this->height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, contourDepthRBO);

    GLuint attachments[3] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2};
    glDrawBuffers(3, attachments);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%: Unable to complete contourFBO", this->ClassName());
    }

    // extend Contour FBO
    for (unsigned int i = 0; i < 2; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, this->extendContourFBO[i]);

        // texture for contours
        glBindTexture(GL_TEXTURE_2D, this->contourTexture[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, contourTexture[i], 0);
        glBindTexture(GL_TEXTURE_2D, 0);

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%: Unable to complete extendContourFBO", this->ClassName());
        }
    }
    // timesteps FBO
    for (unsigned int i = 0; i < 3; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, this->timestepsFBO[i]);

        // texture for contours
        glBindTexture(GL_TEXTURE_2D, this->timestepsTexture[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, timestepsTexture[i], 0);
        glBindTexture(GL_TEXTURE_2D, 0);

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%: Unable to complete timestepsFBO", this->ClassName());
            std::cout << i << std::endl;
        }
    }


    glBindFramebuffer(GL_FRAMEBUFFER, this->curvatureFBO);

    // texture for curvature
    glBindTexture(GL_TEXTURE_2D, this->curvatureTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, curvatureTexture, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%: Unable to complete curvatureFBO", this->ClassName());
    }


    for (unsigned int i = 0; i < 2; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, positionFBO[i]);
        glBindTexture(GL_TEXTURE_2D, smoothPositionTexture[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, smoothPositionTexture[i], 0);
    }
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%: Unable to complete positionFBO", this->ClassName());
    }
    for (unsigned int i = 0; i < 2; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, normalFBO[i]);
        glBindTexture(GL_TEXTURE_2D, smoothNormalTexture[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, smoothNormalTexture[i], 0);
    }
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%: Unable to complete normalFBO", this->ClassName());
    }

    for (unsigned int i = 0; i < 2; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, smoothCurvFBO[i]);
        glBindTexture(GL_TEXTURE_2D, smoothCurvatureTexture[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, smoothCurvatureTexture[i], 0);
    }
    glBindTexture(GL_TEXTURE_2D, 0);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%: Unable to complete curvatureFBO", this->ClassName());
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void MoleculeSESRenderer::CreateQuadBuffers() {

    glDeleteVertexArrays(1, &quadVAO);
    glDeleteBuffers(1, &quadVBO);
    // Create VAO and VBO for screen filling quad (contour generation)
    float quadVertices[] = {
        // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
        // positions   // texCoords
        -1.0f, 1.0f, 0.0f, 1.0f, -1.0f, -1.0f, 0.0f, 0.0f, 1.0f, -1.0f, 1.0f, 0.0f,

        -1.0f, 1.0f, 0.0f, 1.0f, 1.0f, -1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &quadVBO);
    glBindVertexArray(quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*) 0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*) (2 * sizeof(float)));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    std::cout << "quadBuffer" << std::endl;
}

void MoleculeSESRenderer::RenderSESGpuRaycasting(const MolecularDataCall* mol) {

    auto resolution = cameraInfo.resolution_gate();

    bool virtualViewportChanged = false;
    if (static_cast<unsigned int>(resolution.width()) != this->width ||
        static_cast<unsigned int>(resolution.height()) != this->height) {
        this->width = static_cast<unsigned int>(resolution.width());
        this->height = static_cast<unsigned int>(resolution.height());
        virtualViewportChanged = true;
    }

#pragma region // set viewport
    glm::vec4 viewportStuff;
    viewportStuff[0] = 0.0f;
    viewportStuff[1] = 0.0f;
    viewportStuff[2] = static_cast<float>(resolution.width());
    viewportStuff[3] = static_cast<float>(resolution.height());
    if (viewportStuff[2] < 1.0f)
        viewportStuff[2] = 1.0f;
    if (viewportStuff[3] < 1.0f)
        viewportStuff[3] = 1.0f;
    viewportStuff[2] = 2.0f / viewportStuff[2];
    viewportStuff[3] = 2.0f / viewportStuff[3];
#pragma endregion // set viewport

#pragma region // set camerastuff
    glm::vec4 camdir = cameraInfo.view_vector();
    glm::vec4 right = cameraInfo.right_vector();
    glm::vec4 up = cameraInfo.up_vector();
    float nearplane = cameraInfo.near_clipping_plane();
    float farplane = cameraInfo.far_clipping_plane();
    auto campos = cameraInfo.position();
#pragma endregion // set camerastuff

    unsigned int cntRS;

    for (cntRS = 0; cntRS < this->reducedSurface.size(); ++cntRS) {
        //////////////////////////////////
        // ray cast the tori on the GPU //
        //////////////////////////////////
        GLuint attribInParams;
        GLuint attribQuatC;
        GLuint attribInSphere;
        GLuint attribInColors;
        GLuint attribInCuttingPlane;

        if (this->drawSES) {
            if (!testcase) {
#pragma region // torus shader
                if (offscreenRendering) {
                    this->torusShaderOR.Enable();

#pragma region // set shader variables and attributes
                    glUniform4fvARB(
                        this->torusShaderOR.ParameterLocation("viewAttr"), 1, glm::value_ptr(viewportStuff));
                    glUniform3fvARB(this->torusShaderOR.ParameterLocation("camIn"), 1, glm::value_ptr(camdir));
                    glUniform3fvARB(this->torusShaderOR.ParameterLocation("camRight"), 1, glm::value_ptr(right));
                    glUniform3fvARB(this->torusShaderOR.ParameterLocation("camUp"), 1, glm::value_ptr(up));
                    glUniform3fARB(this->torusShaderOR.ParameterLocation("zValues"), 0.0, nearplane, farplane);

                    attribInParams = glGetAttribLocationARB(this->torusShaderOR, "inParams");
                    attribQuatC = glGetAttribLocationARB(this->torusShaderOR, "quatC");
                    attribInSphere = glGetAttribLocationARB(this->torusShaderOR, "inSphere");
                    attribInColors = glGetAttribLocationARB(this->torusShaderOR, "inColors");
                    attribInCuttingPlane = glGetAttribLocationARB(this->torusShaderOR, "inCuttingPlane");
#pragma endregion // set shader variables and attributes
                } else {
                    this->torusShader.Enable();
                    // set shader variables


                    glUniform4fvARB(this->torusShader.ParameterLocation("viewAttr"), 1, glm::value_ptr(viewportStuff));
                    glUniform3fvARB(this->torusShader.ParameterLocation("camIn"), 1, glm::value_ptr(camdir));
                    glUniform3fvARB(this->torusShader.ParameterLocation("camRight"), 1, glm::value_ptr(right));
                    glUniform3fvARB(this->torusShader.ParameterLocation("camUp"), 1, glm::value_ptr(up));
                    glUniform3fARB(this->torusShader.ParameterLocation("zValues"), 0.0, nearplane, farplane);
                    // get attribute locations
                    attribInParams = glGetAttribLocationARB(this->torusShader, "inParams");
                    attribQuatC = glGetAttribLocationARB(this->torusShader, "quatC");
                    attribInSphere = glGetAttribLocationARB(this->torusShader, "inSphere");
                    attribInColors = glGetAttribLocationARB(this->torusShader, "inColors");
                    attribInCuttingPlane = glGetAttribLocationARB(this->torusShader, "inCuttingPlane");
                }

                // set color to orange
                glColor3f(1.0f, 0.75f, 0.0f);
                glEnableClientState(GL_VERTEX_ARRAY);

#pragma region // enable attributes and set pointer
               // enable vertex attribute arrays for the attribute locations
                glEnableVertexAttribArrayARB(attribInParams);
                glEnableVertexAttribArrayARB(attribQuatC);
                glEnableVertexAttribArrayARB(attribInSphere);
                glEnableVertexAttribArrayARB(attribInColors);
                glEnableVertexAttribArrayARB(attribInCuttingPlane);
                // set vertex and attribute pointers and draw them
                glVertexAttribPointerARB(
                    attribInParams, 3, GL_FLOAT, 0, 0, this->torusInParamArray[cntRS].PeekElements());
                glVertexAttribPointerARB(attribQuatC, 4, GL_FLOAT, 0, 0, this->torusQuatCArray[cntRS].PeekElements());
                glVertexAttribPointerARB(
                    attribInSphere, 4, GL_FLOAT, 0, 0, this->torusInSphereArray[cntRS].PeekElements());
                glVertexAttribPointerARB(attribInColors, 4, GL_FLOAT, 0, 0, this->torusColors[cntRS].PeekElements());
                glVertexAttribPointerARB(
                    attribInCuttingPlane, 3, GL_FLOAT, 0, 0, this->torusInCuttingPlaneArray[cntRS].PeekElements());
#pragma endregion // enable attirbutes and set pointers

                glVertexPointer(3, GL_FLOAT, 0, this->torusVertexArray[cntRS].PeekElements());
                glDrawArrays(GL_POINTS, 0, ((unsigned int) this->torusVertexArray[cntRS].Count()) / 3);

#pragma region // disable vertex attribute arrays for the attribute locations
                glDisableVertexAttribArrayARB(attribInParams);
                glDisableVertexAttribArrayARB(attribQuatC);
                glDisableVertexAttribArrayARB(attribInSphere);
                glDisableVertexAttribArrayARB(attribInColors);
                glDisableVertexAttribArrayARB(attribInCuttingPlane);
                glDisableClientState(GL_VERTEX_ARRAY);
#pragma endregion // disable vertex attribute arrays for the attribute locations

                if (offscreenRendering) {
                    this->torusShaderOR.Disable();
                } else {
                    this->torusShader.Disable();
                }

#pragma endregion // torus shader

#pragma region // spherical triangles
               /////////////////////////////////////////////////
               // ray cast the spherical triangles on the GPU //
               /////////////////////////////////////////////////
                GLuint attribVec1;
                GLuint attribVec2;
                GLuint attribVec3;
                GLuint attribTexCoord1;
                GLuint attribTexCoord2;
                GLuint attribTexCoord3;
                GLuint attribColors;

                // bind texture
                glBindTexture(GL_TEXTURE_2D, singularityTexture[cntRS]);
                // enable spherical triangle shader
                if (offscreenRendering) {
                    this->sphericalTriangleShaderOR.Enable();
#pragma region // set shader variables and get attribute locations
               // set shader variables
                    glUniform4fvARB(this->sphericalTriangleShaderOR.ParameterLocation("viewAttr"), 1,
                        glm::value_ptr(viewportStuff));
                    glUniform3fvARB(
                        this->sphericalTriangleShaderOR.ParameterLocation("camIn"), 1, glm::value_ptr(camdir));
                    glUniform3fvARB(
                        this->sphericalTriangleShaderOR.ParameterLocation("camRight"), 1, glm::value_ptr(right));
                    glUniform3fvARB(this->sphericalTriangleShaderOR.ParameterLocation("camUp"), 1, glm::value_ptr(up));
                    glUniform3fARB(
                        this->sphericalTriangleShaderOR.ParameterLocation("zValues"), 0.0, nearplane, farplane);
                    glUniform2fARB(this->sphericalTriangleShaderOR.ParameterLocation("texOffset"),
                        1.0f / (float) this->singTexWidth[cntRS], 1.0f / (float) this->singTexHeight[cntRS]);
                    // get attribute locations
                    attribVec1 = glGetAttribLocationARB(this->sphericalTriangleShaderOR, "attribVec1");
                    attribVec2 = glGetAttribLocationARB(this->sphericalTriangleShaderOR, "attribVec2");
                    attribVec3 = glGetAttribLocationARB(this->sphericalTriangleShaderOR, "attribVec3");
                    attribTexCoord1 = glGetAttribLocationARB(this->sphericalTriangleShaderOR, "attribTexCoord1");
                    attribTexCoord2 = glGetAttribLocationARB(this->sphericalTriangleShaderOR, "attribTexCoord2");
                    attribTexCoord3 = glGetAttribLocationARB(this->sphericalTriangleShaderOR, "attribTexCoord3");
                    attribColors = glGetAttribLocationARB(this->sphericalTriangleShaderOR, "attribColors");
#pragma endregion // set shader variables and get attribute locations
                } else {
                    this->sphericalTriangleShader.Enable();
                    // set shader variables
                    glUniform4fvARB(
                        this->sphericalTriangleShader.ParameterLocation("viewAttr"), 1, glm::value_ptr(viewportStuff));
                    glUniform3fvARB(
                        this->sphericalTriangleShader.ParameterLocation("camIn"), 1, glm::value_ptr(camdir));
                    glUniform3fvARB(
                        this->sphericalTriangleShader.ParameterLocation("camRight"), 1, glm::value_ptr(right));
                    glUniform3fvARB(this->sphericalTriangleShader.ParameterLocation("camUp"), 1, glm::value_ptr(up));
                    glUniform3fARB(
                        this->sphericalTriangleShader.ParameterLocation("zValues"), 0.0, nearplane, farplane);
                    glUniform2fARB(this->sphericalTriangleShader.ParameterLocation("texOffset"),
                        1.0f / (float) this->singTexWidth[cntRS], 1.0f / (float) this->singTexHeight[cntRS]);
                    // get attribute locations
                    attribVec1 = glGetAttribLocationARB(this->sphericalTriangleShader, "attribVec1");
                    attribVec2 = glGetAttribLocationARB(this->sphericalTriangleShader, "attribVec2");
                    attribVec3 = glGetAttribLocationARB(this->sphericalTriangleShader, "attribVec3");
                    attribTexCoord1 = glGetAttribLocationARB(this->sphericalTriangleShader, "attribTexCoord1");
                    attribTexCoord2 = glGetAttribLocationARB(this->sphericalTriangleShader, "attribTexCoord2");
                    attribTexCoord3 = glGetAttribLocationARB(this->sphericalTriangleShader, "attribTexCoord3");
                    attribColors = glGetAttribLocationARB(this->sphericalTriangleShader, "attribColors");
                }

                // set color to turquoise
                glColor3f(0.0f, 0.75f, 1.0f);
                glEnableClientState(GL_VERTEX_ARRAY);
#pragma region // vertex attributes
               // enable vertex attribute arrays for the attribute locations
                glEnableVertexAttribArrayARB(attribVec1);
                glEnableVertexAttribArrayARB(attribVec2);
                glEnableVertexAttribArrayARB(attribVec3);
                glEnableVertexAttribArrayARB(attribTexCoord1);
                glEnableVertexAttribArrayARB(attribTexCoord2);
                glEnableVertexAttribArrayARB(attribTexCoord3);
                glEnableVertexAttribArrayARB(attribColors);
                // set vertex and attribute pointers and draw them
                glVertexAttribPointerARB(attribVec1, 4, GL_FLOAT, 0, 0, this->sphericTriaVec1[cntRS].PeekElements());
                glVertexAttribPointerARB(attribVec2, 4, GL_FLOAT, 0, 0, this->sphericTriaVec2[cntRS].PeekElements());
                glVertexAttribPointerARB(attribVec3, 4, GL_FLOAT, 0, 0, this->sphericTriaVec3[cntRS].PeekElements());
                glVertexAttribPointerARB(
                    attribTexCoord1, 3, GL_FLOAT, 0, 0, this->sphericTriaTexCoord1[cntRS].PeekElements());
                glVertexAttribPointerARB(
                    attribTexCoord2, 3, GL_FLOAT, 0, 0, this->sphericTriaTexCoord2[cntRS].PeekElements());
                glVertexAttribPointerARB(
                    attribTexCoord3, 3, GL_FLOAT, 0, 0, this->sphericTriaTexCoord3[cntRS].PeekElements());
                glVertexAttribPointerARB(
                    attribColors, 3, GL_FLOAT, 0, 0, this->sphericTriaColors[cntRS].PeekElements());
#pragma endregion // vertex attributes
                glVertexPointer(4, GL_FLOAT, 0, this->sphericTriaVertexArray[cntRS].PeekElements());
                glDrawArrays(GL_POINTS, 0, ((unsigned int) this->sphericTriaVertexArray[cntRS].Count()) / 4);
#pragma region // disable vertex attributes
               // disable vertex attribute arrays for the attribute locations
                glDisableVertexAttribArrayARB(attribVec1);
                glDisableVertexAttribArrayARB(attribVec2);
                glDisableVertexAttribArrayARB(attribVec3);
                glDisableVertexAttribArrayARB(attribTexCoord1);
                glDisableVertexAttribArrayARB(attribTexCoord2);
                glDisableVertexAttribArrayARB(attribTexCoord3);
                glDisableVertexAttribArrayARB(attribColors);
                glDisableClientState(GL_VERTEX_ARRAY);
#pragma endregion // disable vertex attributes


                // disable spherical triangle shader
                if (offscreenRendering) {
                    this->sphericalTriangleShaderOR.Disable();
                } else {
                    this->sphericalTriangleShader.Disable();
                }
                // unbind texture
                glBindTexture(GL_TEXTURE_2D, 0);
#pragma endregion // spherical triangles
            }

#pragma region // sphere shader
            /////////////////////////////////////
            // ray cast the spheres on the GPU //
            /////////////////////////////////////
            // enable sphere shader
            if (offscreenRendering) {
                this->sphereShaderOR.Enable();
#pragma region // set shader variables
                glUniform4fvARB(this->sphereShaderOR.ParameterLocation("viewAttr"), 1, glm::value_ptr(viewportStuff));
                glUniform3fvARB(this->sphereShaderOR.ParameterLocation("camIn"), 1, glm::value_ptr(camdir));
                glUniform3fvARB(this->sphereShaderOR.ParameterLocation("camRight"), 1, glm::value_ptr(right));
                glUniform3fvARB(this->sphereShaderOR.ParameterLocation("camUp"), 1, glm::value_ptr(up));
                glUniform3fARB(this->sphereShaderOR.ParameterLocation("zValues"), 0.0, nearplane, farplane);
#pragma endregion // set shader variables
            } else {
                this->sphereShader.Enable();
                // set shader variables


                glUniform4fvARB(this->sphereShader.ParameterLocation("viewAttr"), 1, glm::value_ptr(viewportStuff));
                glUniform3fvARB(this->sphereShader.ParameterLocation("camIn"), 1, glm::value_ptr(camdir));
                glUniform3fvARB(this->sphereShader.ParameterLocation("camRight"), 1, glm::value_ptr(right));
                glUniform3fvARB(this->sphereShader.ParameterLocation("camUp"), 1, glm::value_ptr(up));
                glUniform3fARB(this->sphereShader.ParameterLocation("zValues"), 0.0, nearplane, farplane);
            }
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            // set vertex and color pointers and draw them

            if (testcase) {
                float colors[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
                float vertices[12] = {10.328, 20.121, -4.012, 2, 8.233, 20.042, -0.768, 0.5, 8.13, 20.952, -7.128, 0.5};
                glColorPointer(3, GL_FLOAT, 0, colors);
                glVertexPointer(4, GL_FLOAT, 0, vertices);
                glDrawArrays(GL_POINTS, 0, 3);
            } else {
                glColorPointer(3, GL_FLOAT, 0, this->sphereColors[cntRS].PeekElements());
                glVertexPointer(4, GL_FLOAT, 0, this->sphereVertexArray[cntRS].PeekElements());
                glDrawArrays(GL_POINTS, 0, ((unsigned int) this->sphereVertexArray[cntRS].Count()) / 4);
            }
            // disable sphere shader
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_VERTEX_ARRAY);


            if (offscreenRendering) {
                this->sphereShaderOR.Disable();
            } else {
                this->sphereShader.Disable();
            }
#pragma endregion // sphere shader
            // Offscreen Rendering
        }
    }
}

/*
 * Compute the vertex and attribute arrays for the raycasting shaders
 * (spheres, spherical triangles & tori)
 */
void MoleculeSESRenderer::ComputeRaycastingArrays() {
    // time_t t = clock();

    unsigned int cntRS;
    unsigned int i;

    // resize lists of vertex, attribute and color arrays
    this->sphericTriaVertexArray.resize(this->reducedSurface.size());
    this->sphericTriaVec1.resize(this->reducedSurface.size());
    this->sphericTriaVec2.resize(this->reducedSurface.size());
    this->sphericTriaVec3.resize(this->reducedSurface.size());
    this->sphericTriaTexCoord1.resize(this->reducedSurface.size());
    this->sphericTriaTexCoord2.resize(this->reducedSurface.size());
    this->sphericTriaTexCoord3.resize(this->reducedSurface.size());
    this->sphericTriaColors.resize(this->reducedSurface.size());
    this->torusVertexArray.resize(this->reducedSurface.size());
    this->torusInParamArray.resize(this->reducedSurface.size());
    this->torusQuatCArray.resize(this->reducedSurface.size());
    this->torusInSphereArray.resize(this->reducedSurface.size());
    this->torusColors.resize(this->reducedSurface.size());
    this->torusInCuttingPlaneArray.resize(this->reducedSurface.size());
    this->sphereVertexArray.resize(this->reducedSurface.size());
    this->sphereColors.resize(this->reducedSurface.size());


    // compute singulatity textures
    this->CreateSingularityTextures();

    for (cntRS = 0; cntRS < this->reducedSurface.size(); ++cntRS) {
        ///////////////////////////////////////////////////////////////////////
        // compute arrays for ray casting the spherical triangles on the GPU //
        ///////////////////////////////////////////////////////////////////////
        vislib::math::Vector<float, 3> tmpVec;
        vislib::math::Vector<float, 3> tmpDualProbe(1.0f, 1.0f, 1.0f);
        float dualProbeRad = 0.0f;

        this->sphericTriaVertexArray[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 4);
        this->sphericTriaVec1[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 4);
        this->sphericTriaVec2[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 4);
        this->sphericTriaVec3[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 4);
        this->sphericTriaTexCoord1[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 3);
        this->sphericTriaTexCoord2[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 3);
        this->sphericTriaTexCoord3[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 3);
        this->sphericTriaColors[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSFaceCount() * 3);

        this->probeRadius = this->probeRadiusSlot.Param<param::FloatParam>()->Value();

        // loop over all RS-faces
        for (i = 0; i < this->reducedSurface[cntRS]->GetRSFaceCount(); ++i) {
            // if the face has a dual face --> store the probe of this face
            if (this->reducedSurface[cntRS]->GetRSFace(i)->GetDualFace() != NULL) {
                tmpDualProbe = this->reducedSurface[cntRS]->GetRSFace(i)->GetDualFace()->GetProbeCenter();
                dualProbeRad = this->probeRadius;
            }
            // first RS-vertex
            tmpVec = this->reducedSurface[cntRS]->GetRSFace(i)->GetVertex1()->GetPosition() -
                     this->reducedSurface[cntRS]->GetRSFace(i)->GetProbeCenter();
            this->sphericTriaVec1[cntRS][i * 4 + 0] = tmpVec.GetX();
            this->sphericTriaVec1[cntRS][i * 4 + 1] = tmpVec.GetY();
            this->sphericTriaVec1[cntRS][i * 4 + 2] = tmpVec.GetZ();
            this->sphericTriaVec1[cntRS][i * 4 + 3] = 1.0f;
            // second RS-vertex
            tmpVec = this->reducedSurface[cntRS]->GetRSFace(i)->GetVertex2()->GetPosition() -
                     this->reducedSurface[cntRS]->GetRSFace(i)->GetProbeCenter();
            this->sphericTriaVec2[cntRS][i * 4 + 0] = tmpVec.GetX();
            this->sphericTriaVec2[cntRS][i * 4 + 1] = tmpVec.GetY();
            this->sphericTriaVec2[cntRS][i * 4 + 2] = tmpVec.GetZ();
            this->sphericTriaVec2[cntRS][i * 4 + 3] = 1.0f;
            // third RS-vertex
            tmpVec = this->reducedSurface[cntRS]->GetRSFace(i)->GetVertex3()->GetPosition() -
                     this->reducedSurface[cntRS]->GetRSFace(i)->GetProbeCenter();
            this->sphericTriaVec3[cntRS][i * 4 + 0] = tmpVec.GetX();
            this->sphericTriaVec3[cntRS][i * 4 + 1] = tmpVec.GetY();
            this->sphericTriaVec3[cntRS][i * 4 + 2] = tmpVec.GetZ();
            this->sphericTriaVec3[cntRS][i * 4 + 3] = dualProbeRad * dualProbeRad;
            // store number of cutting probes and texture coordinates for each edge
            this->sphericTriaTexCoord1[cntRS][i * 3 + 0] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge1()->cuttingProbes.size();
            this->sphericTriaTexCoord1[cntRS][i * 3 + 1] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge1()->GetTexCoordX();
            this->sphericTriaTexCoord1[cntRS][i * 3 + 2] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge1()->GetTexCoordY();
            this->sphericTriaTexCoord2[cntRS][i * 3 + 0] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge2()->cuttingProbes.size();
            this->sphericTriaTexCoord2[cntRS][i * 3 + 1] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge2()->GetTexCoordX();
            this->sphericTriaTexCoord2[cntRS][i * 3 + 2] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge2()->GetTexCoordY();
            this->sphericTriaTexCoord3[cntRS][i * 3 + 0] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge3()->cuttingProbes.size();
            this->sphericTriaTexCoord3[cntRS][i * 3 + 1] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge3()->GetTexCoordX();
            this->sphericTriaTexCoord3[cntRS][i * 3 + 2] =
                (float) this->reducedSurface[cntRS]->GetRSFace(i)->GetEdge3()->GetTexCoordY();
            // colors
            this->sphericTriaColors[cntRS][i * 3 + 0] = CodeColor(
                &this->atomColorTable[this->reducedSurface[cntRS]->GetRSFace(i)->GetVertex1()->GetIndex() * 3]);
            this->sphericTriaColors[cntRS][i * 3 + 1] = CodeColor(
                &this->atomColorTable[this->reducedSurface[cntRS]->GetRSFace(i)->GetVertex2()->GetIndex() * 3]);
            this->sphericTriaColors[cntRS][i * 3 + 2] = CodeColor(
                &this->atomColorTable[this->reducedSurface[cntRS]->GetRSFace(i)->GetVertex3()->GetIndex() * 3]);
            // sphere center
            this->sphericTriaVertexArray[cntRS][i * 4 + 0] =
                this->reducedSurface[cntRS]->GetRSFace(i)->GetProbeCenter().GetX();
            this->sphericTriaVertexArray[cntRS][i * 4 + 1] =
                this->reducedSurface[cntRS]->GetRSFace(i)->GetProbeCenter().GetY();
            this->sphericTriaVertexArray[cntRS][i * 4 + 2] =
                this->reducedSurface[cntRS]->GetRSFace(i)->GetProbeCenter().GetZ();
            this->sphericTriaVertexArray[cntRS][i * 4 + 3] = this->GetProbeRadius();
        }

        ////////////////////////////////////////////////////////
        // compute arrays for ray casting the tori on the GPU //
        ////////////////////////////////////////////////////////
        vislib::math::Quaternion<float> quatC;
        vislib::math::Vector<float, 3> zAxis, torusAxis, rotAxis, P, X1, X2, C, planeNormal;
        zAxis.Set(0.0f, 0.0f, 1.0f);
        float distance, d;
        vislib::math::Vector<float, 3> tmpDir1, tmpDir2, tmpDir3, cutPlaneNorm;

        this->torusVertexArray[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSEdgeCount() * 3);
        this->torusInParamArray[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSEdgeCount() * 3);
        this->torusQuatCArray[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSEdgeCount() * 4);
        this->torusInSphereArray[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSEdgeCount() * 4);
        this->torusColors[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSEdgeCount() * 4);
        this->torusInCuttingPlaneArray[cntRS].SetCount(this->reducedSurface[cntRS]->GetRSEdgeCount() * 3);

        // loop over all RS-edges
        for (i = 0; i < this->reducedSurface[cntRS]->GetRSEdgeCount(); ++i) {
            // get the rotation axis of the torus
            torusAxis = this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetPosition() -
                        this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter();
            torusAxis.Normalise();
            // get the axis for rotating the torus rotations axis on the z-axis
            rotAxis = torusAxis.Cross(zAxis);
            rotAxis.Normalise();
            // compute quaternion
            quatC.Set(torusAxis.Angle(zAxis), rotAxis);
            // compute the tangential point X2 of the spheres
            P = this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter() +
                rotAxis * this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusRadius();

            X1 = P - this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetPosition();
            X1.Normalise();
            X1 *= this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetRadius();
            X2 = P - this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetPosition();
            X2.Normalise();
            X2 *= this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetRadius();
            d = (X1 + this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetPosition() -
                 this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter())
                    .Dot(torusAxis);

            C = this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetPosition() -
                this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetPosition();
            C = ((P - this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetPosition()).Length() /
                    ((P - this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetPosition()).Length() +
                        (P - this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetPosition()).Length())) *
                C;
            distance = (X2 - C).Length();
            C = (C + this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetPosition()) -
                this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter();

            // compute normal of the cutting plane
            tmpDir1 = this->reducedSurface[cntRS]->GetRSEdge(i)->GetFace1()->GetProbeCenter();
            tmpDir2 = this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetPosition() - tmpDir1;
            tmpDir2.Normalise();
            tmpDir2 *= this->probeRadius;
            tmpDir2 = tmpDir2 + tmpDir1 - this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter();
            tmpDir3 = this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetPosition() - tmpDir1;
            tmpDir3.Normalise();
            tmpDir3 *= this->probeRadius;
            tmpDir3 = tmpDir3 + tmpDir1 - this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter();
            // tmpDir2 and tmpDir3 now store the position of the intersection points for face 1
            cutPlaneNorm = tmpDir2 - tmpDir3;
            // cutPlaneNorm now stores the vector between the two intersection points for face 1
            tmpDir1 = this->reducedSurface[cntRS]->GetRSEdge(i)->GetFace2()->GetProbeCenter();
            tmpDir2 = this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetPosition() - tmpDir1;
            tmpDir2.Normalise();
            tmpDir2 *= this->probeRadius;
            tmpDir2 = tmpDir2 + tmpDir1 - this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter();
            // tmpDir2 now stores the position of the intersection point 1 for face 2
            tmpDir2 = tmpDir2 - tmpDir3;
            // tmpDir2 and tmpDir3 now span the plane containing the four intersection points
            cutPlaneNorm = cutPlaneNorm.Cross(tmpDir2);
            cutPlaneNorm = torusAxis.Cross(cutPlaneNorm);
            cutPlaneNorm.Normalise();

            // attributes
            this->torusInParamArray[cntRS][i * 3 + 0] = this->probeRadius;
            this->torusInParamArray[cntRS][i * 3 + 1] = this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusRadius();
            this->torusInParamArray[cntRS][i * 3 + 2] = this->reducedSurface[cntRS]->GetRSEdge(i)->GetRotationAngle();
            this->torusQuatCArray[cntRS][i * 4 + 0] = quatC.GetX();
            this->torusQuatCArray[cntRS][i * 4 + 1] = quatC.GetY();
            this->torusQuatCArray[cntRS][i * 4 + 2] = quatC.GetZ();
            this->torusQuatCArray[cntRS][i * 4 + 3] = quatC.GetW();
            this->torusInSphereArray[cntRS][i * 4 + 0] = C.GetX();
            this->torusInSphereArray[cntRS][i * 4 + 1] = C.GetY();
            this->torusInSphereArray[cntRS][i * 4 + 2] = C.GetZ();
            this->torusInSphereArray[cntRS][i * 4 + 3] = distance;
            // colors
            this->torusColors[cntRS][i * 4 + 0] = CodeColor(
                &this->atomColorTable[this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex1()->GetIndex() * 3]);
            this->torusColors[cntRS][i * 4 + 1] = CodeColor(
                &this->atomColorTable[this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetIndex() * 3]);
            this->torusColors[cntRS][i * 4 + 2] = d;
            // this->torusColors[cntRS][i*4+3] = ( X2 - X1).Length();
            this->torusColors[cntRS][i * 4 + 3] =
                (X2 + this->reducedSurface[cntRS]->GetRSEdge(i)->GetVertex2()->GetPosition() -
                    this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter())
                    .Dot(torusAxis) -
                d;
            // cutting plane
            this->torusInCuttingPlaneArray[cntRS][i * 3 + 0] = cutPlaneNorm.GetX();
            this->torusInCuttingPlaneArray[cntRS][i * 3 + 1] = cutPlaneNorm.GetY();
            this->torusInCuttingPlaneArray[cntRS][i * 3 + 2] = cutPlaneNorm.GetZ();
            // torus center
            this->torusVertexArray[cntRS][i * 3 + 0] =
                this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter().GetX();
            this->torusVertexArray[cntRS][i * 3 + 1] =
                this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter().GetY();
            this->torusVertexArray[cntRS][i * 3 + 2] =
                this->reducedSurface[cntRS]->GetRSEdge(i)->GetTorusCenter().GetZ();
        }

        ///////////////////////////////////////////////////////////
        // compute arrays for ray casting the spheres on the GPU //
        ///////////////////////////////////////////////////////////
        /*
        this->sphereVertexArray[cntRS].SetCount( this->reducedSurface[cntRS]->GetRSVertexCount() * 4);
        this->sphereColors[cntRS].SetCount( this->reducedSurface[cntRS]->GetRSVertexCount() * 3);
        */
        this->sphereVertexArray[cntRS].AssertCapacity(this->reducedSurface[cntRS]->GetRSVertexCount() * 4);
        this->sphereVertexArray[cntRS].Clear();
        this->sphereColors[cntRS].AssertCapacity(this->reducedSurface[cntRS]->GetRSVertexCount() * 3);
        this->sphereColors[cntRS].Clear();

        // loop over all RS-vertices (i.e. all protein atoms)
        for (i = 0; i < this->reducedSurface[cntRS]->GetRSVertexCount(); ++i) {
            // add only surface atoms (i.e. with not buried RS-vertices)
            if (this->reducedSurface[cntRS]->GetRSVertex(i)->IsBuried())
                continue;
            // set vertex color
            this->sphereColors[cntRS].Append(
                this->atomColorTable[this->reducedSurface[cntRS]->GetRSVertex(i)->GetIndex() * 3 + 0]);
            this->sphereColors[cntRS].Append(
                this->atomColorTable[this->reducedSurface[cntRS]->GetRSVertex(i)->GetIndex() * 3 + 1]);
            this->sphereColors[cntRS].Append(
                this->atomColorTable[this->reducedSurface[cntRS]->GetRSVertex(i)->GetIndex() * 3 + 2]);
            // set vertex position
            this->sphereVertexArray[cntRS].Append(this->reducedSurface[cntRS]->GetRSVertex(i)->GetPosition().GetX());
            this->sphereVertexArray[cntRS].Append(this->reducedSurface[cntRS]->GetRSVertex(i)->GetPosition().GetY());
            this->sphereVertexArray[cntRS].Append(this->reducedSurface[cntRS]->GetRSVertex(i)->GetPosition().GetZ());
            if (this->drawSAS) {
                this->sphereVertexArray[cntRS].Append(
                    this->reducedSurface[cntRS]->GetRSVertex(i)->GetRadius() + this->probeRadius);
            } else {
                this->sphereVertexArray[cntRS].Append(this->reducedSurface[cntRS]->GetRSVertex(i)->GetRadius());
            }
        }
    }
    // print the time of the computation
    // std::cout << "computation of arrays for GPU ray casting finished:" << ( double( clock() - t) / double(
    // CLOCKS_PER_SEC) ) << std::endl;
}


/*
 * Compute the vertex and attribute arrays for the raycasting shaders
 * (spheres, spherical triangles & tori)
 */
void MoleculeSESRenderer::ComputeRaycastingArrays(unsigned int idxRS) {
    // do nothing if the given index is out of bounds
    if (idxRS > this->reducedSurface.size())
        return;

    this->probeRadius = this->probeRadiusSlot.Param<param::FloatParam>()->Value();

    // check if all arrays have the correct size
    if (this->sphericTriaVertexArray.size() != this->reducedSurface.size() ||
        this->sphericTriaVec1.size() != this->reducedSurface.size() ||
        this->sphericTriaVec2.size() != this->reducedSurface.size() ||
        this->sphericTriaVec3.size() != this->reducedSurface.size() ||
        this->sphericTriaTexCoord1.size() != this->reducedSurface.size() ||
        this->sphericTriaTexCoord2.size() != this->reducedSurface.size() ||
        this->sphericTriaTexCoord3.size() != this->reducedSurface.size() ||
        this->sphericTriaColors.size() != this->reducedSurface.size() ||
        this->torusVertexArray.size() != this->reducedSurface.size() ||
        this->torusInParamArray.size() != this->reducedSurface.size() ||
        this->torusQuatCArray.size() != this->reducedSurface.size() ||
        this->torusInSphereArray.size() != this->reducedSurface.size() ||
        this->torusColors.size() != this->reducedSurface.size() ||
        this->torusInCuttingPlaneArray.size() != this->reducedSurface.size() ||
        this->sphereVertexArray.size() != this->reducedSurface.size() ||
        this->sphereColors.size() != this->reducedSurface.size()) {
        // recompute everything if one of the arrays has the wrong size
        // ComputeRaycastingArrays();
        this->preComputationDone = false;
        return;
    }

    unsigned int i;

    // compute singulatity textures
    this->CreateSingularityTexture(idxRS);

    ///////////////////////////////////////////////////////////////////////
    // compute arrays for ray casting the spherical triangles on the GPU //
    ///////////////////////////////////////////////////////////////////////
    vislib::math::Vector<float, 3> tmpVec;
    vislib::math::Vector<float, 3> tmpDualProbe(1.0f, 1.0f, 1.0f);
    float dualProbeRad = 0.0f;

    this->sphericTriaVertexArray[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 4);
    this->sphericTriaVec1[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 4);
    this->sphericTriaVec2[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 4);
    this->sphericTriaVec3[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 4);
    this->sphericTriaTexCoord1[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 3);
    this->sphericTriaTexCoord2[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 3);
    this->sphericTriaTexCoord3[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 3);
    this->sphericTriaColors[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSFaceCount() * 3);

    // loop over all RS-faces
    for (i = 0; i < this->reducedSurface[idxRS]->GetRSFaceCount(); ++i) {
        // if the face has a dual face --> store the probe of this face
        if (this->reducedSurface[idxRS]->GetRSFace(i)->GetDualFace() != NULL) {
            tmpDualProbe = this->reducedSurface[idxRS]->GetRSFace(i)->GetDualFace()->GetProbeCenter();
            dualProbeRad = this->probeRadius;
        }
        // first RS-vertex
        tmpVec = this->reducedSurface[idxRS]->GetRSFace(i)->GetVertex1()->GetPosition() -
                 this->reducedSurface[idxRS]->GetRSFace(i)->GetProbeCenter();
        this->sphericTriaVec1[idxRS][i * 4 + 0] = tmpVec.GetX();
        this->sphericTriaVec1[idxRS][i * 4 + 1] = tmpVec.GetY();
        this->sphericTriaVec1[idxRS][i * 4 + 2] = tmpVec.GetZ();
        this->sphericTriaVec1[idxRS][i * 4 + 3] = 1.0f;
        // second RS-vertex
        tmpVec = this->reducedSurface[idxRS]->GetRSFace(i)->GetVertex2()->GetPosition() -
                 this->reducedSurface[idxRS]->GetRSFace(i)->GetProbeCenter();
        this->sphericTriaVec2[idxRS][i * 4 + 0] = tmpVec.GetX();
        this->sphericTriaVec2[idxRS][i * 4 + 1] = tmpVec.GetY();
        this->sphericTriaVec2[idxRS][i * 4 + 2] = tmpVec.GetZ();
        this->sphericTriaVec2[idxRS][i * 4 + 3] = 1.0f;
        // third RS-vertex
        tmpVec = this->reducedSurface[idxRS]->GetRSFace(i)->GetVertex3()->GetPosition() -
                 this->reducedSurface[idxRS]->GetRSFace(i)->GetProbeCenter();
        this->sphericTriaVec3[idxRS][i * 4 + 0] = tmpVec.GetX();
        this->sphericTriaVec3[idxRS][i * 4 + 1] = tmpVec.GetY();
        this->sphericTriaVec3[idxRS][i * 4 + 2] = tmpVec.GetZ();
        this->sphericTriaVec3[idxRS][i * 4 + 3] = dualProbeRad * dualProbeRad;
        // store number of cutting probes and texture coordinates for each edge
        this->sphericTriaTexCoord1[idxRS][i * 3 + 0] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge1()->cuttingProbes.size();
        this->sphericTriaTexCoord1[idxRS][i * 3 + 1] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge1()->GetTexCoordX();
        this->sphericTriaTexCoord1[idxRS][i * 3 + 2] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge1()->GetTexCoordY();
        this->sphericTriaTexCoord2[idxRS][i * 3 + 0] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge2()->cuttingProbes.size();
        this->sphericTriaTexCoord2[idxRS][i * 3 + 1] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge2()->GetTexCoordX();
        this->sphericTriaTexCoord2[idxRS][i * 3 + 2] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge2()->GetTexCoordY();
        this->sphericTriaTexCoord3[idxRS][i * 3 + 0] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge3()->cuttingProbes.size();
        this->sphericTriaTexCoord3[idxRS][i * 3 + 1] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge3()->GetTexCoordX();
        this->sphericTriaTexCoord3[idxRS][i * 3 + 2] =
            (float) this->reducedSurface[idxRS]->GetRSFace(i)->GetEdge3()->GetTexCoordY();
        // colors
        this->sphericTriaColors[idxRS][i * 3 + 0] =
            CodeColor(&this->atomColorTable[this->reducedSurface[idxRS]->GetRSFace(i)->GetVertex1()->GetIndex() * 3]);
        this->sphericTriaColors[idxRS][i * 3 + 1] =
            CodeColor(&this->atomColorTable[this->reducedSurface[idxRS]->GetRSFace(i)->GetVertex2()->GetIndex() * 3]);
        this->sphericTriaColors[idxRS][i * 3 + 2] =
            CodeColor(&this->atomColorTable[this->reducedSurface[idxRS]->GetRSFace(i)->GetVertex3()->GetIndex() * 3]);
        // sphere center
        this->sphericTriaVertexArray[idxRS][i * 4 + 0] =
            this->reducedSurface[idxRS]->GetRSFace(i)->GetProbeCenter().GetX();
        this->sphericTriaVertexArray[idxRS][i * 4 + 1] =
            this->reducedSurface[idxRS]->GetRSFace(i)->GetProbeCenter().GetY();
        this->sphericTriaVertexArray[idxRS][i * 4 + 2] =
            this->reducedSurface[idxRS]->GetRSFace(i)->GetProbeCenter().GetZ();
        this->sphericTriaVertexArray[idxRS][i * 4 + 3] = this->GetProbeRadius();
    }

    ////////////////////////////////////////////////////////
    // compute arrays for ray casting the tori on the GPU //
    ////////////////////////////////////////////////////////
    vislib::math::Quaternion<float> quatC;
    vislib::math::Vector<float, 3> zAxis, torusAxis, rotAxis, P, X1, X2, C, planeNormal;
    zAxis.Set(0.0f, 0.0f, 1.0f);
    float distance, d;
    vislib::math::Vector<float, 3> tmpDir1, tmpDir2, tmpDir3, cutPlaneNorm;

    this->torusVertexArray[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSEdgeCount() * 3);
    this->torusInParamArray[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSEdgeCount() * 3);
    this->torusQuatCArray[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSEdgeCount() * 4);
    this->torusInSphereArray[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSEdgeCount() * 4);
    this->torusColors[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSEdgeCount() * 4);
    this->torusInCuttingPlaneArray[idxRS].SetCount(this->reducedSurface[idxRS]->GetRSEdgeCount() * 3);

    // loop over all RS-edges
    for (i = 0; i < this->reducedSurface[idxRS]->GetRSEdgeCount(); ++i) {
        // get the rotation axis of the torus
        torusAxis = this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetPosition() -
                    this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter();
        torusAxis.Normalise();
        // get the axis for rotating the torus rotations axis on the z-axis
        rotAxis = torusAxis.Cross(zAxis);
        rotAxis.Normalise();
        // compute quaternion
        quatC.Set(torusAxis.Angle(zAxis), rotAxis);
        // compute the tangential point X2 of the spheres
        P = this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter() +
            rotAxis * this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusRadius();

        X1 = P - this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetPosition();
        X1.Normalise();
        X1 *= this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetRadius();
        X2 = P - this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetPosition();
        X2.Normalise();
        X2 *= this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetRadius();
        d = (X1 + this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetPosition() -
             this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter())
                .Dot(torusAxis);

        C = this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetPosition() -
            this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetPosition();
        C = ((P - this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetPosition()).Length() /
                ((P - this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetPosition()).Length() +
                    (P - this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetPosition()).Length())) *
            C;
        distance = (X2 - C).Length();
        C = (C + this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetPosition()) -
            this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter();

        // compute normal of the cutting plane
        tmpDir1 = this->reducedSurface[idxRS]->GetRSEdge(i)->GetFace1()->GetProbeCenter();
        tmpDir2 = this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetPosition() - tmpDir1;
        tmpDir2.Normalise();
        tmpDir2 *= this->probeRadius;
        tmpDir2 = tmpDir2 + tmpDir1 - this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter();
        tmpDir3 = this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetPosition() - tmpDir1;
        tmpDir3.Normalise();
        tmpDir3 *= this->probeRadius;
        tmpDir3 = tmpDir3 + tmpDir1 - this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter();
        // tmpDir2 and tmpDir3 now store the position of the intersection points for face 1
        cutPlaneNorm = tmpDir2 - tmpDir3;
        // cutPlaneNorm now stores the vector between the two intersection points for face 1
        tmpDir1 = this->reducedSurface[idxRS]->GetRSEdge(i)->GetFace2()->GetProbeCenter();
        tmpDir2 = this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetPosition() - tmpDir1;
        tmpDir2.Normalise();
        tmpDir2 *= this->probeRadius;
        tmpDir2 = tmpDir2 + tmpDir1 - this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter();
        // tmpDir2 now stores the position of the intersection point 1 for face 2
        tmpDir2 = tmpDir2 - tmpDir3;
        // tmpDir2 and tmpDir3 now span the plane containing the four intersection points
        cutPlaneNorm = cutPlaneNorm.Cross(tmpDir2);
        cutPlaneNorm = torusAxis.Cross(cutPlaneNorm);
        cutPlaneNorm.Normalise();

        // attributes
        this->torusInParamArray[idxRS][i * 3 + 0] = this->probeRadius;
        this->torusInParamArray[idxRS][i * 3 + 1] = this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusRadius();
        this->torusInParamArray[idxRS][i * 3 + 2] = this->reducedSurface[idxRS]->GetRSEdge(i)->GetRotationAngle();
        this->torusQuatCArray[idxRS][i * 4 + 0] = quatC.GetX();
        this->torusQuatCArray[idxRS][i * 4 + 1] = quatC.GetY();
        this->torusQuatCArray[idxRS][i * 4 + 2] = quatC.GetZ();
        this->torusQuatCArray[idxRS][i * 4 + 3] = quatC.GetW();
        this->torusInSphereArray[idxRS][i * 4 + 0] = C.GetX();
        this->torusInSphereArray[idxRS][i * 4 + 1] = C.GetY();
        this->torusInSphereArray[idxRS][i * 4 + 2] = C.GetZ();
        this->torusInSphereArray[idxRS][i * 4 + 3] = distance;
        // colors
        this->torusColors[idxRS][i * 4 + 0] =
            CodeColor(&this->atomColorTable[this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex1()->GetIndex() * 3]);
        this->torusColors[idxRS][i * 4 + 1] =
            CodeColor(&this->atomColorTable[this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetIndex() * 3]);
        this->torusColors[idxRS][i * 4 + 2] = d;
        this->torusColors[idxRS][i * 4 + 3] =
            (X2 + this->reducedSurface[idxRS]->GetRSEdge(i)->GetVertex2()->GetPosition() -
                this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter())
                .Dot(torusAxis) -
            d;
        // cutting plane
        this->torusInCuttingPlaneArray[idxRS][i * 3 + 0] = cutPlaneNorm.GetX();
        this->torusInCuttingPlaneArray[idxRS][i * 3 + 1] = cutPlaneNorm.GetY();
        this->torusInCuttingPlaneArray[idxRS][i * 3 + 2] = cutPlaneNorm.GetZ();
        // torus center
        this->torusVertexArray[idxRS][i * 3 + 0] = this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter().GetX();
        this->torusVertexArray[idxRS][i * 3 + 1] = this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter().GetY();
        this->torusVertexArray[idxRS][i * 3 + 2] = this->reducedSurface[idxRS]->GetRSEdge(i)->GetTorusCenter().GetZ();
    }

    ///////////////////////////////////////////////////////////
    // compute arrays for ray casting the spheres on the GPU //
    ///////////////////////////////////////////////////////////
    this->sphereVertexArray[idxRS].AssertCapacity(this->reducedSurface[idxRS]->GetRSVertexCount() * 4);
    this->sphereVertexArray[idxRS].Clear();
    this->sphereColors[idxRS].AssertCapacity(this->reducedSurface[idxRS]->GetRSVertexCount() * 3);
    this->sphereColors[idxRS].Clear();

    // loop over all RS-vertices (i.e. all protein atoms)
    for (i = 0; i < this->reducedSurface[idxRS]->GetRSVertexCount(); ++i) {
        // add only surface atoms (i.e. with not buried RS-vertices)
        if (this->reducedSurface[idxRS]->GetRSVertex(i)->IsBuried())
            continue;
        // set vertex color
        this->sphereColors[idxRS].Append(
            this->atomColorTable[this->reducedSurface[idxRS]->GetRSVertex(i)->GetIndex() * 3 + 0]);
        this->sphereColors[idxRS].Append(
            this->atomColorTable[this->reducedSurface[idxRS]->GetRSVertex(i)->GetIndex() * 3 + 1]);
        this->sphereColors[idxRS].Append(
            this->atomColorTable[this->reducedSurface[idxRS]->GetRSVertex(i)->GetIndex() * 3 + 2]);
        // set vertex position
        this->sphereVertexArray[idxRS].Append(this->reducedSurface[idxRS]->GetRSVertex(i)->GetPosition().GetX());
        this->sphereVertexArray[idxRS].Append(this->reducedSurface[idxRS]->GetRSVertex(i)->GetPosition().GetY());
        this->sphereVertexArray[idxRS].Append(this->reducedSurface[idxRS]->GetRSVertex(i)->GetPosition().GetZ());
        if (this->drawSAS) {
            this->sphereVertexArray[idxRS].Append(
                this->reducedSurface[idxRS]->GetRSVertex(i)->GetRadius() + this->probeRadius);
        } else {
            this->sphereVertexArray[idxRS].Append(this->reducedSurface[idxRS]->GetRSVertex(i)->GetRadius());
        }
    }
}


/*
 * code a rgb-color into one float
 */
float MoleculeSESRenderer::CodeColor(const float* col) const {
    return float((int) (col[0] * 255.0f) * 1000000 // red
                 + (int) (col[1] * 255.0f) * 1000  // green
                 + (int) (col[2] * 255.0f));       // blue
}


/*
 * decode a coded color to the original rgb-color
 */
vislib::math::Vector<float, 3> MoleculeSESRenderer::DecodeColor(int codedColor) const {
    int col = codedColor;
    vislib::math::Vector<float, 3> color;
    float red, green;
    if (col >= 1000000)
        red = floor((float) col / 1000000.0f);
    else
        red = 0.0;
    col = col - int(red * 1000000.0f);
    if (col > 1000)
        green = floor((float) col / 1000.0f);
    else
        green = 0.0;
    col = col - int(green * 1000.0f);
    // color.Set( red / 255.0f, green / 255.0f, float(col) / 255.0f);
    color.Set(std::min(1.0f, std::max(0.0f, red / 255.0f)), std::min(1.0f, std::max(0.0f, green / 255.0f)),
        std::min(1.0f, std::max(0.0f, col / 255.0f)));
    return color;
}

void MoleculeSESRenderer::CreateSingularityTextures() {
    // time_t t = clock();
    unsigned int cnt1, cnt2, cntRS;

    // delete old singularity textures
    for (cnt1 = 0; cnt1 < this->singularityTexture.size(); ++cnt1) {
        glDeleteTextures(1, &singularityTexture[cnt1]);
    }
    // check if the singularity texture has the right size
    if (this->reducedSurface.size() != this->singularityTexture.size()) {
        // store old singularity texture size
        unsigned int singTexSizeOld = (unsigned int) this->singularityTexture.size();
        // resize singularity texture to fit the number of reduced surfaces
        this->singularityTexture.resize(this->reducedSurface.size());
        // generate a new texture for each new singularity texture
        for (cnt1 = singTexSizeOld; cnt1 < singularityTexture.size(); ++cnt1) {
            glGenTextures(1, &singularityTexture[cnt1]);
        }
    }
    // resize singularity texture dimension arrays
    this->singTexWidth.resize(this->reducedSurface.size());
    this->singTexHeight.resize(this->reducedSurface.size());

    // get maximum texture size
    GLint texSize;
    glGetIntegerv(GL_MAX_TEXTURE_SIZE, &texSize);

    unsigned int numProbes = 16;

    for (cntRS = 0; cntRS < this->reducedSurface.size(); ++cntRS) {
        // set width and height of texture
        if ((unsigned int) texSize < this->reducedSurface[cntRS]->GetCutRSEdgesCount()) {
            this->singTexHeight[cntRS] = texSize;
            this->singTexWidth[cntRS] =
                numProbes * (int) ceil(double(this->reducedSurface[cntRS]->GetCutRSEdgesCount()) / (double) texSize);
        } else {
            this->singTexHeight[cntRS] = this->reducedSurface[cntRS]->GetCutRSEdgesCount();
            this->singTexWidth[cntRS] = numProbes;
        }
        // generate float-array for texture with the appropriate dimension
        if (this->singTexData)
            delete[] this->singTexData;
        this->singTexData = new float[this->singTexWidth[cntRS] * this->singTexHeight[cntRS] * 3];
        // write probes to singularity texture
        unsigned int coordX = 0;
        unsigned int coordY = 0;
        unsigned int counter = 0;
        for (cnt1 = 0; cnt1 < this->reducedSurface[cntRS]->GetRSEdgeCount(); ++cnt1) {
            if (this->reducedSurface[cntRS]->GetRSEdge(cnt1)->cuttingProbes.empty()) {
                this->reducedSurface[cntRS]->GetRSEdge(cnt1)->SetTexCoord(0, 0);
            } else {
                // set texture coordinates
                this->reducedSurface[cntRS]->GetRSEdge(cnt1)->SetTexCoord(coordX, coordY);
                // compute texture coordinates for next entry
                coordY++;
                if (coordY == this->singTexHeight[cntRS]) {
                    coordY = 0;
                    coordX = coordX + numProbes;
                }
                // write probes to texture
                for (cnt2 = 0; cnt2 < numProbes; ++cnt2) {
                    if (cnt2 < this->reducedSurface[cntRS]->GetRSEdge(cnt1)->cuttingProbes.size()) {
                        singTexData[counter] =
                            this->reducedSurface[cntRS]->GetRSEdge(cnt1)->cuttingProbes[cnt2]->GetProbeCenter().GetX();
                        counter++;
                        singTexData[counter] =
                            this->reducedSurface[cntRS]->GetRSEdge(cnt1)->cuttingProbes[cnt2]->GetProbeCenter().GetY();
                        counter++;
                        singTexData[counter] =
                            this->reducedSurface[cntRS]->GetRSEdge(cnt1)->cuttingProbes[cnt2]->GetProbeCenter().GetZ();
                        counter++;
                    } else {
                        singTexData[counter] = 0.0f;
                        counter++;
                        singTexData[counter] = 0.0f;
                        counter++;
                        singTexData[counter] = 0.0f;
                        counter++;
                    }
                }
            }
        }
        // texture generation
        glBindTexture(GL_TEXTURE_2D, singularityTexture[cntRS]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F_ARB, this->singTexWidth[cntRS], this->singTexHeight[cntRS], 0, GL_RGB,
            GL_FLOAT, this->singTexData);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    // std::cout << "Create texture: " << ( double( clock() - t) / double( CLOCKS_PER_SEC) ) << std::endl;
}


void MoleculeSESRenderer::CreateSingularityTexture(unsigned int idxRS) {
    // do nothing if the index is out of bounds
    if (idxRS > this->reducedSurface.size())
        return;

    // check if all arrays have the appropriate size
    if (this->singularityTexture.size() != this->reducedSurface.size() ||
        this->singTexWidth.size() != this->reducedSurface.size() ||
        this->singTexHeight.size() != this->reducedSurface.size()) {
        // create all singularity textures
        CreateSingularityTextures();
        return;
    }

    unsigned int cnt1, cnt2;

    // delete old singularity texture
    glDeleteTextures(1, &singularityTexture[idxRS]);

    // get maximum texture size
    GLint texSize;
    glGetIntegerv(GL_MAX_TEXTURE_SIZE, &texSize);

    unsigned int numProbes = 16;

    // set width and height of texture
    if ((unsigned int) texSize < this->reducedSurface[idxRS]->GetCutRSEdgesCount()) {
        this->singTexHeight[idxRS] = texSize;
        this->singTexWidth[idxRS] =
            numProbes * (int) ceil(double(this->reducedSurface[idxRS]->GetCutRSEdgesCount()) / (double) texSize);
    } else {
        this->singTexHeight[idxRS] = this->reducedSurface[idxRS]->GetCutRSEdgesCount();
        this->singTexWidth[idxRS] = numProbes;
    }
    // generate float-array for texture with the appropriate dimension
    if (this->singTexData)
        delete[] this->singTexData;
    this->singTexData = new float[this->singTexWidth[idxRS] * this->singTexHeight[idxRS] * 3];
    // write probes to singularity texture
    unsigned int coordX = 0;
    unsigned int coordY = 0;
    unsigned int counter = 0;
    for (cnt1 = 0; cnt1 < this->reducedSurface[idxRS]->GetRSEdgeCount(); ++cnt1) {
        if (this->reducedSurface[idxRS]->GetRSEdge(cnt1)->cuttingProbes.empty()) {
            this->reducedSurface[idxRS]->GetRSEdge(cnt1)->SetTexCoord(0, 0);
        } else {
            // set texture coordinates
            this->reducedSurface[idxRS]->GetRSEdge(cnt1)->SetTexCoord(coordX, coordY);
            // compute texture coordinates for next entry
            coordY++;
            if (coordY == this->singTexHeight[idxRS]) {
                coordY = 0;
                coordX = coordX + numProbes;
            }
            // write probes to texture
            for (cnt2 = 0; cnt2 < numProbes; ++cnt2) {
                if (cnt2 < this->reducedSurface[idxRS]->GetRSEdge(cnt1)->cuttingProbes.size()) {
                    singTexData[counter] =
                        this->reducedSurface[idxRS]->GetRSEdge(cnt1)->cuttingProbes[cnt2]->GetProbeCenter().GetX();
                    counter++;
                    singTexData[counter] =
                        this->reducedSurface[idxRS]->GetRSEdge(cnt1)->cuttingProbes[cnt2]->GetProbeCenter().GetY();
                    counter++;
                    singTexData[counter] =
                        this->reducedSurface[idxRS]->GetRSEdge(cnt1)->cuttingProbes[cnt2]->GetProbeCenter().GetZ();
                    counter++;
                } else {
                    singTexData[counter] = 0.0f;
                    counter++;
                    singTexData[counter] = 0.0f;
                    counter++;
                    singTexData[counter] = 0.0f;
                    counter++;
                }
            }
        }
    }
    // texture generation
    glBindTexture(GL_TEXTURE_2D, singularityTexture[idxRS]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F_ARB, this->singTexWidth[idxRS], this->singTexHeight[idxRS], 0, GL_RGB,
        GL_FLOAT, this->singTexData);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void MoleculeSESRenderer::deinitialise(void) {
    if (contourFBO) {
        glDeleteFramebuffers(1, &contourFBO);
        glDeleteRenderbuffers(1, &contourDepthRBO);
        glDeleteTextures(1, &normalTexture);
        glDeleteTextures(1, &positionTexture);
        glDeleteTextures(1, &objPositionTexture);
    }
    if (curvatureFBO) {
        glDeleteFramebuffers(1, &curvatureFBO);
        glDeleteTextures(1, &curvatureTexture);
    }
    if (positionFBO) {
        glDeleteFramebuffers(2, positionFBO);
        glDeleteTextures(2, smoothPositionTexture);
    }
    if (normalFBO) {
        glDeleteFramebuffers(2, normalFBO);
        glDeleteTextures(2, smoothNormalTexture);
    }
    if (smoothCurvFBO) {
        glDeleteFramebuffers(2, smoothCurvFBO);
        glDeleteTextures(2, smoothCurvatureTexture);
    }
    // delete singularity texture
    for (unsigned int i = 0; i < singularityTexture.size(); ++i)
        glDeleteTextures(1, &singularityTexture[i]);
    // release shaders
    this->sphereShader.Release();
    this->sphericalTriangleShader.Release();
    this->torusShader.Release();
    this->lightShader.Release();
}


/*
 * returns the color of the atom 'idx' for the current coloring mode
 */
vislib::math::Vector<float, 3> MoleculeSESRenderer::GetProteinAtomColor(unsigned int idx) {
    if (idx < this->atomColorTable.Count() / 3)
        // return this->atomColorTable[idx];
        return vislib::math::Vector<float, 3>(
            this->atomColorTable[idx * 3 + 0], this->atomColorTable[idx * 3 + 1], this->atomColorTable[idx * 3 + 0]);
    else
        return vislib::math::Vector<float, 3>(0.5f, 0.5f, 0.5f);
}

/*
 * Render two spheres, one big one small for testing purposes
 */
void MoleculeSESRenderer::RenderTestCase() {
    unsigned int VBO, VAO;
    glGenBuffers(1, &VBO);
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 3 * 4 * sizeof(float), NULL, GL_STATIC_DRAW);
    // fill buffer
    float pos[9] = {-0.2, 0, 1, 0.7, 0, 1, 0.2, 0, -2};
    float col[3] = {0.2, 0.7, 0.2};
    float* positions = &pos[0];
    float* colors = &col[0];
    int posSize = 9 * sizeof(float);
    glBufferSubData(GL_ARRAY_BUFFER, 0, posSize, positions);
    glBufferSubData(GL_ARRAY_BUFFER, posSize, 3 * sizeof(float), colors);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*) (intptr_t)(posSize));
    glBindVertexArray(0);
    // scale coordinates such that they lie in cube with edges [-1,1].
    glm::vec3 viewPos = glm::vec3(0.0f, 0.0f, -3.0f);
    glm::mat4 view = glm::mat4(1.0f);
    view = glm::translate(view, glm::vec3(viewPos.x, viewPos.y, viewPos.z));
    // float radiusSphere = 0.2; // radius for the spheres representing the atoms
    testCaseShader.Enable();
    // ------
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glm::mat4 proj;
    proj = glm::ortho(-1.0f, 3.0f, -1.0f, 1.0f, 1.0f, 10.0f);

    glUniformMatrix4fv(testCaseShader.ParameterLocation("view"), 1, false, &view[0][0]);
    glUniformMatrix4fv(testCaseShader.ParameterLocation("proj"), 1, false, &proj[0][0]);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBindVertexArray(VAO);
    glDrawArrays(GL_POINTS, 0, 3);
}
/*
 * Render two spheres, one big one small for testing purposes
 */
void MoleculeSESRenderer::RenderPerspectiveTestCase() {
    unsigned int VBO, VAO;
    glGenBuffers(1, &VBO);
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    // fill buffer
    float positions[] = {
        1.0f, 1.0f, 0.0f,   // top right
        1.0f, -1.0f, 0.0f,  // bottom right
        -1.0f, 1.0f, 0.0f,  // top left
                            // second triangle
        1.0f, -1.0f, 0.0f,  // bottom right
        -1.0f, -1.0f, 0.0f, // bottom left
        -1.0f, 1.0f, 0.0f   // top left
    };
    glBufferData(GL_ARRAY_BUFFER, sizeof(positions), positions, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);

    // scale coordinates such that they lie in cube with edges [-1,1].
    glm::vec3 viewPos = glm::vec3(0.0f, 0.0f, -3.0f);
    glm::mat4 view = glm::mat4(1.0f);
    view = glm::translate(view, glm::vec3(viewPos.x, viewPos.y, viewPos.z));
    // float radiusSphere = 0.2; // radius for the spheres representing the atoms
    perspectiveTestCaseShader.Enable();
    // ------
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glm::mat4 proj;
    proj = glm::perspective(glm::radians(45.0f), (float) this->width / (float) this->height, 0.1f, 10.0f);
    glUniformMatrix4fv(perspectiveTestCaseShader.ParameterLocation("view"), 1, false, &view[0][0]);
    glUniformMatrix4fv(perspectiveTestCaseShader.ParameterLocation("proj"), 1, false, &proj[0][0]);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBindVertexArray(VAO);

    // first sphere
    float radius = 0.2;
    glm::mat4 model = glm::mat4(1.0);
    model = glm::translate(model, glm::vec3(-0.2, 0.0, 1.0));
    model = glm::scale(model, glm::vec3(radius));
    glUniformMatrix4fv(perspectiveTestCaseShader.ParameterLocation("model"), 1, false, &model[0][0]);
    glUniform1f(perspectiveTestCaseShader.ParameterLocation("radius"), radius);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    // second sphere
    radius = 0.7;
    model = glm::mat4(1.0);
    model = glm::translate(model, glm::vec3(0.7, 0.0, 1.0));
    model = glm::scale(model, glm::vec3(radius));
    glUniformMatrix4fv(perspectiveTestCaseShader.ParameterLocation("model"), 1, false, &model[0][0]);
    glUniform1f(perspectiveTestCaseShader.ParameterLocation("radius"), radius);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    // third sphere
    radius = 0.2;
    model = glm::mat4(1.0);
    model = glm::translate(model, glm::vec3(0.2, 0.0, -2.0));
    model = glm::scale(model, glm::vec3(radius));
    glUniformMatrix4fv(perspectiveTestCaseShader.ParameterLocation("model"), 1, false, &model[0][0]);
    glUniform1f(perspectiveTestCaseShader.ParameterLocation("radius"), radius);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}
